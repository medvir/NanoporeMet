import os
import subprocess
import glob
import matplotlib.pyplot as plt

def create_output_directory(directory, reference_name):
    output_directory = os.path.join(directory, reference_name)
    os.makedirs(output_directory, exist_ok=True)
    return output_directory

def find_reference_fasta(directory):
    os.chdir(directory)
    reference_files = [file for file in os.listdir() if file.endswith('.fasta') and not file.startswith('._')]
    
    if not reference_files:
        raise FileNotFoundError("No reference fasta file found in the specified directory.")
    
    if len(reference_files) > 1:
        raise ValueError("Multiple reference fasta files found. Please ensure only one reference file is present.")
    
    return reference_files[0]

def concatenate_fastq(directory):
    os.chdir(directory)
    
    concatenated_file = 'concatenated.fastq'

    if os.path.exists(concatenated_file):
        print("Using previously generated concatenated.fastq file.")
        return
    
    all_fastq_files = glob.glob('**/*.fastq.gz', recursive=True)

    with open(concatenated_file, 'wb') as outfile:
        for fastq_file in all_fastq_files:
            with open(fastq_file, 'rb') as infile:
                outfile.write(infile.read())

def run_minimap2(reference_file, concatenated_file, output_directory, reference_name):
    sam_output = os.path.join(output_directory, f'{reference_name}.sam')
    subprocess.run(['minimap2', '-ax', 'map-ont', reference_file, concatenated_file], stdout=open(sam_output, 'w'))

def run_samtools_sort(output_directory, reference_name):
    bam_output = os.path.join(output_directory, f'{reference_name}.bam')
    subprocess.run(['samtools', 'sort', os.path.join(output_directory, f'{reference_name}.sam'), '-o', bam_output])

def run_samtools_depth(output_directory, reference_name):
    subprocess.run(['samtools', 'depth', '-a', os.path.join(output_directory, f'{reference_name}.bam')],
                   stdout=open(os.path.join(output_directory, f'{reference_name}.coverage'), 'w'))

def plot_coverage(output_directory, reference_name):
    data = {'x': [], 'y': []}
    with open(os.path.join(output_directory, f'{reference_name}.coverage'), 'r') as coverage_file:
        for line in coverage_file:
            columns = line.strip().split('\t')
            data['x'].append(int(columns[1]))
            data['y'].append(int(columns[2]))

    horizontal_coverage = sum(1 for val in data['y'] if val > 0) / len(data['y'])
    vertical_coverage = sum(data['y']) / len(data['y'])

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(data['x'], data['y'], color='black', linewidth=1.5)
    ax.set_yscale('log')
    ax.set_xlabel('Genome position')
    ax.set_ylabel('Coverage (reads)')
    ax.set_title(f"{columns[0]} ({len(data['x'])} nt)", fontweight='bold', loc='left', y=1.08)
    ax.grid(which='major', linestyle='-', color='lightgrey')
    ax.minorticks_on()
    ax.grid(which='minor', linestyle='-', color='lightgrey', alpha=0.5)

    subtitle = f"Horizontal coverage: {horizontal_coverage:.2%} | Mean vertical coverage: {vertical_coverage:.2f}X"
    fig.text(0.125, 0.92, subtitle, ha='left', va='center', fontweight='normal', fontsize=10)

    ax.tick_params(axis='both', which='both', length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    pdf_path = os.path.join(output_directory, f'{reference_name}.pdf')
    plt.savefig(pdf_path, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.show()

    print(f"Plot saved as {pdf_path}")

if __name__ == '__main__':
    
    fastq_pass_directory = os.getcwd()
    concatenate_fastq(fastq_pass_directory)

    reference_directory = input('Enter the path to the reference sequence directory: ')
    reference_file = find_reference_fasta(reference_directory)
    
    reference_name = os.path.splitext(reference_file)[0]

    output_directory = create_output_directory(fastq_pass_directory, reference_name)

    os.chdir(fastq_pass_directory)
    run_minimap2(os.path.join(reference_directory, reference_file), 'concatenated.fastq', output_directory, reference_name)

    run_samtools_sort(output_directory, reference_name)

    run_samtools_depth(output_directory, reference_name)

    plot_coverage(output_directory, reference_name)
