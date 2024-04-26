import os
import glob
import subprocess
import gzip
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_histogram(data, column_name, title, xlabel, ylabel, ax, log_scale=False):
    n, bins, patches = ax.hist(data[column_name], bins=100, edgecolor='black', log=log_scale, color='white', lw=1.2, zorder=2)
    median_value = data[column_name].median()
    rounded_median = round(median_value) if column_name == "sequence_length_template" else round(median_value, 2)
    ax.set_title(f"{title}", loc='left')  # Align title to the left
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if log_scale:
        ax.set_yscale('log')
    ax.grid(axis='y', linestyle='-', alpha=0.5, color='lightgrey', zorder=1)  # Add light-grey gridlines

    # Plot median as a vertical line
    ax.axvline(x=median_value, color='red', linestyle='dashed', linewidth=2, label=f"Median: {rounded_median}")
    ax.legend()

def generate_sequencing_summary_pdf(directory):
    pdf_path = os.path.join(directory, "sequencing_summary.pdf")

    # Check if the sequencing_summary.pdf file already exists
    if os.path.exists(pdf_path):
        print("The sequencing_summary.pdf file has already been created.")
        return pdf_path

    # Find the file with the name "sequencing_summary_*.txt"
    file_pattern = os.path.join(directory, "sequencing_summary_*.txt")
    files = glob.glob(file_pattern)

    if not files:
        print("The sequencing_summary.txt file is not available yet.")
        return pdf_path

    # Take the first matching file (you can modify this logic if needed)
    file_path = files[0]

    # Read the data from the file
    data = pd.read_csv(file_path, sep='\t')

    # Create histograms and save in one PDF file
    with PdfPages(pdf_path) as pdf:
        # Plot for all data
        fig, axs = plt.subplots(figsize=(8, 6))
        plot_histogram(data, "mean_qscore_template", "Mean Q score (all reads)",
                       "Mean Q score", "Frequency", axs)

        pdf.savefig(fig)
        plt.close(fig)

        # Filter based on the "passes_filtering" column
        filtered_data = data[data["passes_filtering"] == True]

        # Plot for filtered data
        fig, axs = plt.subplots(figsize=(8, 6))  # Create a new figure
        plot_histogram(filtered_data, "mean_qscore_template", "Mean Q score (quality-filtered reads)",
                       "Mean Q score", "Frequency", axs)

        pdf.savefig(fig)
        plt.close(fig)

        # Plot for sequence length (all reads)
        fig, axs = plt.subplots(figsize=(8, 6))  # Create a new figure
        plot_histogram(data, "sequence_length_template", "Sequence length (all reads)",
                       "Sequence length", "Frequency", axs, log_scale=True)

        pdf.savefig(fig)
        plt.close(fig)

        # Plot for sequence length (quality-filtered reads)
        fig, axs = plt.subplots(figsize=(8, 6))  # Create a new figure
        plot_histogram(filtered_data, "sequence_length_template", "Sequence length (quality-filtered reads)",
                       "Sequence length", "Frequency", axs, log_scale=True)

        pdf.savefig(fig)
        plt.close(fig)

    print(f"Plots saved in {pdf_path}")
    return pdf_path

def run_kraken(directory, kraken_db, analyze_bacterial):
    # Path to the output .fastq file
    output_fastq = os.path.join(directory, "concatenated.fastq")
    
    # Collect fastq.gz files in the directory
    fastq_files = [file for file in os.listdir(directory) if file.endswith('.fastq.gz')]
    
    # Check if there are fastq.gz files in the directory
    if fastq_files:
        # Open the output file in binary mode and concatenate the content of fastq.gz files
        with open(output_fastq, 'wb') as output_handle:
            for file in fastq_files:
                with open(os.path.join(directory, file), 'rb') as input_handle:
                    output_handle.write(input_handle.read())
        
        # Run Kraken2 on the concatenated .fastq file
        output_kraken = os.path.join(directory, f"{os.path.basename(directory)}.kreport.txt")
        kraken_command = f"kraken2 --db {kraken_db} --output {output_fastq.replace('.fastq', '.kraken')} --report {output_kraken} {output_fastq}"
        subprocess.run(kraken_command, shell=True)
        
        # Add a column with title to the kreport.txt file
        title = os.path.basename(directory)
        with open(output_kraken, 'r') as f:
            lines = f.readlines()
        lines = [f"{title}\t{line}" for line in lines]  # Prepend the title column
        with open(output_kraken, 'w') as f:
            f.writelines(lines)
        
        # Cleanup: Remove the concatenated .fastq file after Kraken2 analysis
        os.remove(output_fastq)
    else:
        print(f"No fastq.gz files found in {directory}")

def main():
    # Use the current working directory as the analysis directory
    analysis_directory = os.getcwd()

    # Ask if bacterial reads analysis is desired
    analyze_bacterial = input("Do you wish to analyze bacterial reads? (yes/y or no/n): ").lower().strip()

    # Generate sequencing_summary.pdf
    pdf_path = generate_sequencing_summary_pdf(analysis_directory)

    # Set Kraken2 database path based on user input
    if analyze_bacterial in ['yes', 'y']:
        kraken_db_path = "/data/kraken_databases/k2_pluspf_08gb_20231009/"
        output_file_name = "virus_bacteria.kraken.txt"
    else:
        kraken_db_path = "/data/kraken_databases/k2_human-viral_20240111/"
        output_file_name = "virus.kraken.txt"

    # List to store individual kreport.txt DataFrames
    kreport_dfs = []

    # Navigate through each sub-folder in "fastq_pass" and run Kraken2 analysis
    fastq_pass_directory = os.path.join(analysis_directory, "fastq_pass")
    for subdir in os.listdir(fastq_pass_directory):
        subdir_path = os.path.join(fastq_pass_directory, subdir)
        if os.path.isdir(subdir_path) and any(os.scandir(subdir_path)):
            run_kraken(subdir_path, kraken_db_path, analyze_bacterial)
            output_kraken = os.path.join(subdir_path, f"{os.path.basename(subdir_path)}.kreport.txt")
            kreport_df = pd.read_csv(output_kraken, sep='\t', header=None)
            kreport_dfs.append(kreport_df)
        else:
            # For real-time analysis when the folder already exists but no file yet
            print(f"Skipping empty subdirectory: {subdir}")

    # Combine all kreport.txt DataFrames into one DataFrame
    combined_kreport = pd.concat(kreport_dfs)

    # Duplicate and modify the DataFrame for "all" rows
    combined_kreport_all = combined_kreport.copy()
    combined_kreport_all[0] = "all"

    # Concatenate the original DataFrame and the modified DataFrame
    combined_kreport_final = pd.concat([combined_kreport, combined_kreport_all])

    # Save the combined DataFrame with the dynamically generated output file name
    combined_output_path = os.path.join(analysis_directory, output_file_name)
    combined_kreport_final.to_csv(combined_output_path, sep='\t', header=False, index=False)

    print(f"Kraken analysis results saved in {combined_output_path}")

if __name__ == "__main__":
    main()
