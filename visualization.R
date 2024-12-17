library(shiny)
library(stringr)
library(dplyr)
library(future)
library(promises)
# plot package
library(ggplot2)
library(igraph)
library(visNetwork)

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
             "color", 
             "Color Palette", 
             list("contrast color palette" = "contrast", "rainbow" = "rainbow") 
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

  contrast <- c("#D00000","#3185FC","#FFBA08","#5D2E8C","#8FE388","#FF9B85","#8c5c2b","#696663") # 7 colors
  rainbow <- c()
  ###################################################
  # get data 
  ###################################################
  # all data information
  uploaded_data <- reactive({
    req(input$upload)
    input$upload
  })

  output$files <- renderTable({
    req(uploaded_data())
    uploaded_data()
  })
  
  output$test <- renderTable({
    req(dist())
    head(dist())
    ### head(dist_rc(), n = 3)
  })

  # concise.tsv
  concise <- reactive({
    req(uploaded_data())
    filename <- uploaded_data()[grepl(".concise.tsv",uploaded_data()$name),][1,"datapath"]
    read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # print("concise done")
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
    # print("concise done")
  })

  # dist.tsv
  dist <- reactive({
    req(uploaded_data())
    filename <- uploaded_data()[grepl(".dist.tsv",uploaded_data()$name),][1,"datapath"]
    read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # print("concise done")
  })
  
  motif2id <- reactive({
    req(motif())
    tmp = setNames(motif()$id, motif()$motif)
    print(tmp)
    tmp
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
  output$region_to_show_ui <- renderUI({
    req(ncol(concise) > 0)  # ensure file has data

    column_data <- motif[[1]]

    sliderInput("slider", "选择值",
                min = min(column_data, na.rm = TRUE),
                max = max(column_data, na.rm = TRUE),
                value = c(min(column_data, na.rm = TRUE), max(column_data, na.rm = TRUE)))
  })
  
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
  
  graph_data <- reactive({
    req(motif())

    motif_list <- unique(motif()$motif)

    # 创建空图并添加顶点
    g <- make_empty_graph(directed = FALSE)
    g <- add_vertices(g, length(motif_list), name = motif_list)

    # 计算距离矩阵
    if(input$merge_rc){
      distance_df <- rbind(dist(), dist_rc())
    } else {
      distance_df <- dist()
    }
    distance_df <- distance_df[order(distance_df$distance), ]

    tmp_edges <- data.frame(from = integer(0), to = integer(0), value = integer(0))

    # 添加边并计算连通分量
    num_components <- components(g)$no
    if(num_components > input$num_of_motifs) {
      for(i in 1:nrow(distance_df)) {
        motif1 <- distance_df$motif1[i]
        motif2 <- distance_df$motif2[i]
        from_id <- which(motif_list == motif1)
        to_id <- which(motif_list == motif2)

        tmp_edges <- rbind(tmp_edges, data.frame(from = from_id, to = to_id,
                                                 value = min(length(motif1), length(motif2)) - distance_df$distance[i],
                                                 label = as.character(distance_df$distance[i]),
                                                 dashes = distance_df$rc[i]))
        g <- add_edges(g, c(from_id, to_id))

        num_components <- components(g)$no
        if(num_components <= input$num_of_motifs) {
          break
        }
      }
    }


    # 返回图和边的数据
    list(graph = g, edges = tmp_edges, motif_list = motif_list, distance_df = distance_df)
  })

  nodes <- reactive({
    req(graph_data())  # 使用共享的 graph_data

    if(input$color == 'contrast'){
      color_list = contrast
    }else{
      color_list = rainbow
    }

    g <- graph_data()$graph
    motif_list <- graph_data()$motif_list
    motif_rep_num <- motif()$rep_num

    # 获取连通分量和计算每个子图的权重总和
    comp <- components(g)
    print("#######")
    ### print(comp)
    component_sizes <- sapply(unique(comp$membership), function(x) sum(motif_rep_num[comp$membership == x]) )
    tmp <- data.frame( total_rep_num = component_sizes, old_cluster_id = unique(comp$membership))
    tmp = tmp[order(tmp$total_rep_num, decreasing = TRUE),]
    tmp$new_cluster_id = 1:nrow(tmp)
    print(tmp)

    old_cluster_id <- data.frame(old_cluster_id  = comp$membership)
    code_color_id <- old_cluster_id %>%
      left_join(tmp, by = "old_cluster_id") %>%
      mutate(color_id = ifelse(new_cluster_id < length(color_list), new_cluster_id, length(color_list))) %>%
      pull(color_id)
    print(code_color_id)

    # 生成排序后的节点信息
    data.frame(id = 1:length(motif_list),
               label = motif_list,
               value = motif_rep_num,
               color = color_list[code_color_id])

  })

  edges <- reactive({
    req(graph_data())  # 使用共享的 graph_data
    graph_data()$edges
  })

  output$network <- renderVisNetwork({
    req(motif())
    motif_list <- unique(motif()$motif)


    visNetwork(nodes(), edges(), directed = FALSE)
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