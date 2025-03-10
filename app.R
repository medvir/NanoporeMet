library(shiny)
library(tidyverse)

ui <- fluidPage(
  titlePanel("NanoporeMet"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload virus.kraken.txt or virus_bacteria.kraken.txt"),
      selectInput("barcode", "Select barcode", choices = NULL),
      selectInput("category", "Select domain", choices = c("", "Virus", "Bacteria")),
      selectInput("taxonomy_level", "Select taxonomy level", choices = c("", "Species", "Genus")),
      checkboxInput("remove_retrovirus", "Hide phages and endogenous retroviruses"),
      checkboxInput("hide_blocklisted_viruses", "Hide blocklisted viruses")
    ),
    mainPanel(
      textOutput("text"),
      plotOutput("barplot"),
      tableOutput("dataTable")
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$file, {
    req(input$file$datapath)
    data <- read.table(input$file$datapath, sep = "\t", header = FALSE)
    updateSelectInput(session, "barcode", choices = sort(unique(data$V1)))
  })
  
  output$text <- renderText({
    req(input$file$datapath, input$barcode)
    
    data <- read.table(input$file$datapath, sep = "\t", header = FALSE)
    
    # Divide by 2 because dataset is duplicated to create "all" barcode
    total_analyzed_reads <- sum(data$V3[data$V7 %in% c("unclassified", "root")])/2
    
    selected_rows <- data[data$V1 == input$barcode & data$V7 %in% c("root", "unclassified"), ]
    total_reads <- sum(selected_rows$V3)
    
    return(paste0("Analyzed reads total: ", total_analyzed_reads, ", ", input$barcode, ": ", total_reads))
  })
  
  output$barplot <- renderPlot({
    req(input$barcode)
    
    data <- read.table(input$file$datapath, sep = "\t", header = FALSE)
    data <- data[data$V1 == input$barcode, ]
    data$V7 <- trimws(data$V7)
    data$V7 <- recode(data$V7, "Homo sapiens" = "Human",
                      "Bacteria" = "Bacterial",
                      "Fungi" = "Fungal",
                      "Viruses" = "Viral",
                      "unclassified" = "Unclassified")
    filtered_data <- subset(data, data$V7 %in% c("Human", "Bacterial", "Fungal", "Viral", "Unclassified"))
    filtered_data$V7 <- factor(filtered_data$V7, levels = c("Human", "Bacterial", "Fungal", "Viral", "Unclassified"))
    
    # Sum up the reads for each category
    aggregated_data <- filtered_data %>%
      group_by(V7) %>%
      summarize(Reads = sum(V3))
    
    ggplot(aggregated_data, aes(x = V7, y = Reads, fill = V7)) +
      geom_bar(position = "dodge", stat = "identity") +
      labs(x = "", y = "reads", fill = "") +
      scale_y_log10(breaks = 10^(-10:10), minor_breaks = rep(1:9, 21) * (10^rep(-10:10, each = 9))) +
      theme_minimal() +
      theme(axis.text.x = element_blank()) +
      theme(plot.title = element_text(face = "bold", size = 14)) +
      scale_fill_manual("", values = c("Human" = "#143642",
                                       "Bacterial" = "#0F8B8D",
                                       "Fungal" = "#EC9A29",
                                       "Viral" = "#A8201A",
                                       "Unclassified" = "#A69CAC"),
                        limits = force)
  })
  
  output$dataTable <- renderTable({
    req(input$barcode)
    
    data <- read.table(input$file$datapath, sep = "\t", header = FALSE)
    data <- data[data$V1 == input$barcode, ]
    names(data) <- paste0("X", 1:ncol(data))
    
    data <- data %>%
      mutate(X8 = case_when(
        grepl("Viruses", X7) ~ "Virus",
        grepl("Bacteria", X7) ~ "Bacteria",
        X7 == "unclassified" ~ "unclassified",
        X7 == "root" ~ "root",
        grepl("cellular organisms", X7) ~ "Human",
        grepl("Fungi", X7) ~ "Fungi",
        TRUE ~ NA_character_
      ))
    
    # Go through every row from top to bottom
    for (i in 2:nrow(data)) {
      # When row X8 in row n has a value, leave it
      if (!is.na(data$X8[i])) {
        next
      }
      
      # When row n+1 is empty, take value from row n
      if (is.na(data$X8[i]) & !is.na(data$X8[i - 1])) {
        data$X8[i] <- data$X8[i - 1]
      }
    }
    
    # Check and merge rows with the same values in columns 6 and 7 among "all" rows
    if ("all" %in% data$X1) {
      taxonomy_filter <- ifelse(input$taxonomy_level == "Species", "S", "G")
      
      data <- data %>%
        filter(X1 == "all" & X5 == taxonomy_filter) %>%
        group_by(X6, X7) %>%
        mutate(X3 = sum(X3)) %>%
        distinct(X6, X7, .keep_all = TRUE) %>%
        ungroup() %>%
        bind_rows(data %>% filter(X1 != "all" | X5 != taxonomy_filter))
    }
    
    if (input$remove_retrovirus) {
      data <- data %>%
        filter(!(grepl("phage|escherichia|streptococcus|staphylococcus|bacillus|actinomyces|ostreococcus|myoviridae|clostridium|shigella|haemophilus|endogenous retrovirus|endogenous", X7, ignore.case = TRUE) & X8 == "Virus"))
    }
    
    if (input$hide_blocklisted_viruses) {
      data <- data %>%
        filter(!(grepl("Spodoptera litura nucleopolyhedrovirus II|baculovirus|nudivirus|Badaguanvirus IME347|Dolusvirus pacificense|Cyvirus cyprinidallo3|Mollivirus sibericum|
                       |pandoravirus|Thaumetopoea pityocampa iflavirus 1|Caulimovirus tessellodahliae|Feravirus neuropterus|Leptopilina boulardi filamentous virus|
                       |Aresaunavirus RsoM1USA|Spbetavirus SPbeta|Eganvirus SW9|Groundnut rosette assistor virus|Erannis ankeraria nucleopolyhedrovirus|Kostyavirus porky|
                       |Lacusarxvirus lacusarx|Daphnia iridovirus 1|Anayavirus chris|Kinglevirus lutadaptatum|Fowl aviadenovirus C|Megavirus chilense|Twortvirus twort|
                       |Peduovirus P22H4|Spodoptera exigua multiple nucleopolyhedrovirus|Ugandan cassava brown streak virus|Mardivirus columbidalpha1|Cosavirus D|
                       |Noni mosaic virus|Tomato leaf curl Yemen betasatellite|Alphaportoglobovirus SPV2|Yellowstone lake phycodnavirus 3|Cymopoleiavirus swam2|
                       |Ostreavirus ostreidmalaco1|Moumouvirus australiense|Canhaevirus hiberniae|Plateaulakevirus pv4L372XY|Franklinbayvirus fv9A|Palaemonvirus pssm7|
                       |Gemykroznavirus hydro1|Acinetobacter virus Acj61|BeAn 58058 virus|Cyvirus cyprinidallo1|Betatectivirus AP50|Pahexavirus PHL171M01|Haetaevirus PBC2|
                       |Hubei picorna-like virus 45|Mason-Pfizer monkey virus|Pinnievirus moorethemaryer|Fadolivirus algeromassiliense|Phaeocystis globosa virus|
                       |Eceepunavirus EcP1|Orpheovirus IHUMI-LCC2|Theiavirus salishense|Gibbon ape leukemia virus|Emiliania huxleyi virus 86|Laroyevirus laroye|
                       |Elephant endotheliotropic herpesvirus 4|Ictavirus acipenseridallo2|Firehammervirus CP21|Caeruleovirus Bcp1|White spot syndrome virus|
                       |Cohcovirus hiberniae|Sinsheimervirus phiX174|Sida golden mottle virus|Muldoonvirus muldoon|Anaposvirus socalone|Cyvirus cyprinidallo2|
                       |Yellowseavirus thirtyeight|Dishui lake phycodnavirus 1|Anaposvirus socalone|Haifavirus tim68|Eponavirus epona|Muromegalovirus muridbeta1|
                       |Snuvirus SNUABM7|Baltimorevirus DFL12|Percavirus equidgamma5|Acanthamoeba castellanii medusavirus|Hanrivervirus slyngel|Donellivirus gee|
                       |", X7, ignore.case = TRUE) & X8 == "Virus"))
    }
    
    category_selected <- input$category
    taxonomy_level <- ifelse(input$taxonomy_level == "Species", "S", "G")
    
    filtered_data <- data %>%
      filter(X8 == category_selected & X5 == taxonomy_level) %>%
      select("Taxonomy" = X7, "NCBI taxonomy ID" = X6, "Reads" = X3)
    
    names(filtered_data)[1] <- ifelse(input$taxonomy_level == "Species", "Species", "Genus")
    
    filtered_data$RPM <- as.integer(round(filtered_data$Reads * 1000000 / sum(data$X3[data$X7 %in% c("root", "unclassified")])))
    
    filtered_data <- arrange(filtered_data, desc(Reads))
    
    filtered_data <- filtered_data %>%
      mutate_at(vars(1), ~if_else(grepl("T1|Tuna|MS2|Levi|Emes|zinderi", .), sprintf('<span style="font-weight: bold">%s</span>', .), .))
    
    return(filtered_data)
  }, sanitize.text.function = function(x) x)
}

shinyApp(ui = ui, server = server)

