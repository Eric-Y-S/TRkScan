library(shiny)
library(ggplot2)
library(igraph)
library(visNetwork)
library(stringr)

setwd("D:/MyFile/git_repo/TRkScan/testData/")

ui <- fluidPage(
  
  ###################################################
  # CSS setting
  ###################################################
  tags$style(HTML("
    .title {
      display: flex;
      justify-content: center;  /* horizontally center */
      align-items: center;      /* vertically center */
      font-size: 30px;
      font-weight: bold;
      height: 80px;
      /* background: lightgrey; */
      margin-bottom: 20px;
    }
    .title p {
      height: auto;
      color: RoyalBlue;
      margin-bottom: 0;
    }
    .readme {
      display: block;
      color: black;
      font-size: 14px;
      background: WhiteSmoke;
      margin-bottom: 20px;
      padding-left: 15px;
      padding-top: 5px;
      padding-bottom: 5px;
    }
    .readme p {
      color: black;
      margin-bottom: 0;
    }
    .module-name {
      display: flex;
      justify-content: center;  /* horizontally center */
      align-items: center;      /* vertically center */
      font-size: 24px;
      height: 40px;
      background: WhiteSmoke;
      margin-top: 30px;
      margin-bottom: 10px;
    }
  ")),
  
  ###################################################
  # title and notes
  ###################################################
  div(class = "title", p("TRkScan")),
  div(class = "readme", 
        p("readme readme"),
        p("readme readme")
      ),
  
  ###################################################
  # global parameters
  ###################################################
  div(class = "module-name", "Data and Global Parameters"),
  # data to read and global parameters
  fluidRow(     
    column(6, 
          fileInput("upload", NULL, buttonLabel = "Upload...", multiple = TRUE),
    ),
    column(6, 
           # uiOutput("region_to_show_ui"),
           # sliderInput("region to show", "Range", value = c(10, 20), min = 0, max = 100),
           selectInput( 
             "select", 
             "Color Palette", 
             list("Choice 1A" = "1A", "Choice 1B" = "1B") 
           )
    )
  ),
  
  # uploaded data
  fluidRow(     
    column(12, 
           tableOutput("files"),
           tableOutput("test")
    )
  ),
  

  ###################################################
  # MODULE 1: motif count and cluster
  ###################################################
  div(class = "module-name", "Motif Count and Cluster"),
  
  fluidRow(     # motif count and cluster
    column(4, 
      uiOutput("num_of_motif_ui"),
      checkboxInput("merge_rc", "merge_rc", FALSE)
    ),
    column(8,
      visNetworkOutput("network")
    )
  ),
  
  ###################################################
  # MODULE 2: motif annotation
  ###################################################
  div(class = "module-name", "Motif Annotation Visualization"),
  
  fluidRow(    
    column(4, 
       selectInput( 
         "x-axis", 
         "X-axis unit", 
         list("base pair" = "base-pair", "repetition number" = "rep-num") 
       )
    ),
    column(8, 
       plotOutput("motif_vis")
    )
  ),
  
  ###################################################
  # MODULE 3: xxxxxxxxxxxxxxx
  ###################################################
  div(class = "module-name", "xxxxx"),
  
  fluidRow(     
    column(4, 
           
    ),
    column(8, 
           
    )
  ),
)

server <- function(input, output, session) {

  ###################################################
  # get data 
  ###################################################
  # all data information
  uploaded_data <- reactive({
    req(input$upload)
    uploaded_data <- input$upload
  })

  output$files <- renderTable({
    req(uploaded_data())
    uploaded_data()
  })
  
  output$test <- renderTable({
    req(dist_rc())
    dist_rc()
    
  })

  # concise.tsv
  concise <- reactive({
    req(uploaded_data())
    filename <- uploaded_data()[grepl(".concise.tsv",uploaded_data()$name),][1,"datapath"]
    read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  
  # annotation.tsv
  annotation <- reactive({
    req(uploaded_data())
    filename <- uploaded_data()[grepl(".annotation.tsv",uploaded_data()$name),][1,"datapath"]
    read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  
  # motif.tsv
  motif <- reactive({
    req(uploaded_data())
    filename <- uploaded_data()[grepl(".motif.tsv",uploaded_data()$name),][1,"datapath"]
    read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })

  # motif2id
  motif2id <- reactive({
    req(motif())

  })
  
  # distance (not rc)
  dist <- reactive({
    req(motif())

    tmp = motif()
    motif_list = tmp$motif
    # tmp = tmp[motif_list]
    result <- data.frame(motif1 = character(0), motif2 = character(0), distance = integer(0))
    for(i in 1:(length(motif_list) - 1)){
      for(j in (i+1):length(motif_list)){
        # print(i,j)
        motif1 <- motif_list[i]
        motif2 <- motif_list[j]
        distance1 <- as.integer(str_extract(tmp[tmp$motif == motif1, motif2], "^\\d+"))
        result <- rbind(result, data.frame(motif1 = motif1, motif2 = motif2, distance = distance1))
      }
    }
    result = result[order(result$distance), ]
    rownames(result) <- 1:nrow(result)
    result
  })
  
  # distance (rc)
  dist_rc <- reactive({
    req(motif())
    
    tmp = motif()
    motif_list = tmp$motif
    # tmp = tmp[motif_list]
    result <- data.frame(motif1 = character(0), motif2 = character(0), distance = integer(0))
    for(i in 1:(length(motif_list) - 1)){
      for(j in (i+1):length(motif_list)){
        # print(i,j)
        motif1 <- motif_list[i]
        motif2 <- motif_list[j]
        distance1 <- as.integer(str_extract(tmp[tmp$motif == motif1, motif2], "(?<=,)(\\d+)"))
        result <- rbind(result, data.frame(motif1 = motif1, motif2 = motif2, distance = distance1))
      }
    }
    result = result[order(result$distance), ]
    rownames(result) <- 1:nrow(result)
    result
  })
  
  #########################
  output$file_preview <- renderText({
    tmp = uploaded_data()
    tmp = tmp[grepl(".concise.tsv",tmp$name),]
    filename = tmp[1,"name"]
    
    concise_tsv <- reactive({
      read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    })
    
    paste("You are reading:",filename)
    
    # df = read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # df
  })
  
  
  ###################################################
  # global parameters
  ###################################################
  # output$region_to_show_ui <- renderUI({
  #   req(ncol(concise) > 0)  # ensure file has data
  #   
  #   column_data <- motif[[1]]
  #   
  #   sliderInput("slider", "选择值", 
  #               min = min(column_data, na.rm = TRUE), 
  #               max = max(column_data, na.rm = TRUE), 
  #               value = c(min(column_data, na.rm = TRUE), max(column_data, na.rm = TRUE)))
  # })
  
  ###################################################
  # MODULE 1: motif count and cluster
  ###################################################
  output$num_of_motif_ui <- renderUI({
    req(motif())
    max_num <- length(unique(motif()$motif))

    sliderInput("num_of_motifs", "Number of motifs",
                min = 1,
                max = max_num,
                value = max_num,
                step = 1)
  })
  
  motif_cluster <- reactive({
    req(motif())
    req(input$num_of_motifs)
    
  })
  
  
  
  output$network <- renderVisNetwork({
    req(motif())
    motif_list <- unique(motif()$motif)
    
    
    
    nodes <- data.frame(id = 1:length(motif_list), label = motif()$motif, value = motif()$rep_num)
    edges <- data.frame(from = integer(0), to = integer(0), value = integer(0))
    
    ########!!!!!!!!!!!!!!!!!!
    distance_df <- dist()
    distance_df <- distance_df[order(distance_df$distance), ]

    # create a empty graph
    g <- make_empty_graph(directed = FALSE)
    g <- add_vertices(g, length(motif_list), name = motif_list)

    # add edge in cycle，until number of cluster = X
    num_components <- components(g)$no
    if(num_components > input$num_of_motifs) {
      for(i in 1:nrow(distance_df)) {
        motif1 <- distance_df$motif1[i]
        motif2 <- distance_df$motif2[i]
        from_id <- which(motif_list == motif1)
        to_id <- which(motif_list == motif2)
        
        edges <- rbind(edges, data.frame(from = from_id, to = to_id, value = distance_df$distance[i]))
        g <- add_edges(g, c(from_id, to_id))
        num_components <- components(g)$no
        if(num_components <= input$num_of_motifs) {
          break
        }
      }
    }
    
    
    visNetwork(nodes, edges, directed = FALSE)
  })
  
  ###################################################
  # MODULE 2: motif annotation
  ###################################################
  
  
  # output$motif_vis <- renderPlot({
  #   seq(concise())
  #   seq_list = concise()$seq
  #   
  #   ggplot() + 
  #     geom_tile() +
  #     geom_text(aes(label = seq_list), data = 1:length(seq_list)) +
  #     theme_classic()
  # })
  
  
  ###################################################
  # MODULE 3: xxxxxxxxxxxxxxx
  ###################################################
  
 
  
  
}


shinyApp(ui, server)