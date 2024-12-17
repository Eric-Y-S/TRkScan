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
       uiOutput("region_to_show_ui"),
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
  edge_scaling_factor <- 1
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
    tmp = read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    tmp$id <- as.integer(tmp$id)
    tmp
  })

  # dist.tsv
  dist <- reactive({
    req(uploaded_data())
    filename <- uploaded_data()[grepl(".dist.tsv",uploaded_data()$name),][1,"datapath"]
    tmp = read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    tmp$ref <- as.integer(tmp$ref)
    tmp$query <- as.integer(tmp$query)
    transfer = c('True' = TRUE, 'False' = FALSE)
    tmp$is_rc = transfer[tmp$is_rc]
    tmp
  })
  
  id2motif <- reactive({
    req(motif())
    setNames(motif()$motif ,motif()$id)
  })
  
  motif2id <- reactive({
    req(motif())
    setNames(motif()$id, motif()$motif)
  })
  
  ###################################################
  # global parameters
  ###################################################

  
  ###################################################
  # MODULE 1: motif count and cluster
  ###################################################
  output$num_of_motif_ui <- renderUI({
    req(motif())
    max_num <- length(motif()$motif)

    sliderInput("num_of_motifs", "Number of motifs",
                min = 1,
                max = max_num,
                value = max_num,
                step = 1)
  })
  
  graph_data <- reactive({
    req(motif())
    req(dist())

    # create empty graph and add nodes
    g <- make_empty_graph(directed = FALSE)
    g <- add_vertices(g, nrow(motif()), name = motif()$id)
    ### print(typeof(motif()$id))
    # get distance data
    distance_df <- dist()
    if(! input$merge_rc){
      distance_df <- distance_df[! distance_df$is_rc,]
    }

    # add edges and calculate connected components
    tmp_edges <- data.frame(from = integer(0), to = integer(0), value = integer(0))
    num_components <- components(g)$no
    if(num_components > input$num_of_motifs) {
      for(i in 1:nrow(distance_df)) {
        id1 <- distance_df$ref[i]
        id2 <- distance_df$query[i]

        tmp_edges <- rbind(tmp_edges, data.frame(from = id1, to = id2,
                                                 value = pmax(edge_scaling_factor / distance_df$dist[i], 0.1),
                                                 label = as.character(distance_df$dist[i]),
                                                 dashes = distance_df$is_rc[i]))
        ### print(c(id1,id2))
        g <- add_edges(g, c(id1 + 1, id2 + 1))

        num_components <- components(g)$no
        if(num_components <= input$num_of_motifs) {
          break
        }
      }
    }


    # return data
    list(graph = g, edges = tmp_edges, distance_df = distance_df)
  })

  nodes <- reactive({
    req(graph_data())  # use shared graph_data

    if(input$color == 'contrast'){
      color_list = contrast
    }else{
      color_list = rainbow
    }

    g <- graph_data()$graph
    motif_rep_num <- motif()$rep_num

    # get connected components and calculate the sum weight of each sub graph
    comp <- components(g)
    component_sizes <- sapply(unique(comp$membership), function(x) sum(motif_rep_num[comp$membership == x]) )
    tmp <- data.frame( total_rep_num = component_sizes, old_cluster_id = unique(comp$membership))
    tmp = tmp[order(tmp$total_rep_num, decreasing = TRUE),]
    tmp$new_cluster_id = 1:nrow(tmp)

    old_cluster_id <- data.frame(old_cluster_id  = comp$membership)
    code_color_id <- old_cluster_id %>%
      left_join(tmp, by = "old_cluster_id") %>%
      mutate(color_id = ifelse(new_cluster_id < length(color_list), new_cluster_id, length(color_list))) %>%
      pull(color_id)
    ### print(code_color_id)

    # 生成排序后的节点信息
    data.frame(id = motif()$id,
               label = motif()$motif,
               value = motif()$rep_num,
               color = color_list[code_color_id])

  })

  edges <- reactive({
    req(graph_data())  # 使用共享的 graph_data
    graph_data()$edges
  })

  output$network <- renderVisNetwork({
    req(nodes())
    req(edges())
    visNetwork(nodes(), edges(), directed = FALSE)
  })
  
  ###################################################
  # MODULE 2: motif annotation
  ###################################################
  output$region_to_show_ui <- renderUI({
    req(concise())  # ensure file has data
    
    length <- concise()[1,'length']
    
    sliderInput("rgn2show", "region to show",
                min = 0,
                max = length,
                value = c(0, length),
                step = 1)
  })
  
  split_cigar <- function(seq, length, start, end, motif, cigar) {
    segment_data <- data.frame()
    
    new_start = start
    cur = new_start
    cigar_str = ''
    num_str = ''
    # print(length(cigar))
    for (i in 1:nchar(cigar)){
      symbol = substring(cigar, i, i)
      # print(symbol)
      if (symbol %in% c('=','X','I') ){
        ### print("#############33")
        cur = cur + as.integer(num_str)
        cigar_str =  paste0(cigar_str, num_str, symbol)
        num_str = ''
      } else if (symbol == 'D'){
        cigar_str =  paste0(cigar_str, num_str, symbol)
        num_str = ''
      } else if (symbol == 'N'){
        segment_data <- rbind(segment_data, data.frame(
          seq = seq,
          length = length,
          start = new_start,
          end = cur,
          motif = motif,  
          CIGAR = cigar_str
        ))
        new_start = cur + as.integer(num_str)
        cur = new_start
        cigar_str = ''
        num_str = ''
      } else if (symbol == '/'){
        cigar_str = paste0(cigar_str, '/')
      } else {
        num_str = paste0(num_str, symbol)
      }
    }
    if (length(cigar_str) > 0){
      segment_data <- rbind(segment_data, data.frame(
        seq = seq,
        length = length,
        start = new_start,
        end = cur,
        motif = motif,  
        CIGAR = cigar_str
      ))
    }
    return(segment_data)
  }
  
  base_pair_concise_annotation <- reactive({
    seq(concise())
    seq(nodes())
    
    df = concise()
    # split with N character
    new_df <- data.frame()
    for (i in 1:nrow(df)) {
      seq <- df$seq[i]
      length <- df$length[i]
      start <- df$start[i]
      end <- df$end[i]
      motif <- df$motif[i]
      cigar <- df$CIGAR[i]
      # print(typeof(cigar))
      # print(grepl("N", cigar))
      if (grepl("N", cigar)) {
        # print(typeof(start))
        # print(typeof(end))
        # print(typeof(cigar))
        split_result <- split_cigar(seq, length, start, end, motif, cigar)
        # print(split_result)
        new_df <- rbind(new_df, split_result)
      } else {
        new_df <- rbind(new_df, df[i, c('seq','length','start','end','motif','CIGAR')])
      }
    }
    
    colormap = setNames(nodes()$color, nodes()$label)
    new_df$color = colormap[new_df$motif]
    
    df = new_df
    new_df = data.frame()
    ### print(df)
    for (i in 1:(nrow(df)-1)){
      if (i == 1){
        start = df[1,'start']
        end = df[1,'end']
      }
      if (df[i,'color'] == df[i+1, 'color'] && df[i,'end'] == df[i+1,'start']){
        end = df[i+1,'end']
      } else {
        new_df <- rbind(new_df, data.frame(
          seq = df[i,'seq'],
          length = df[i,'seq'],
          start = start,
          end = end,
          color = df[i,'color']
        ))
        start = df[i+1,'start']
        end = df[i+1,'end']
      }
    }
    new_df <- rbind(new_df, data.frame(
      seq = df[i,'seq'],
      length = df[i,'seq'],
      start = start,
      end = end,
      color = df[i,'color']
    ))
    
    
    # give y axis
    seq_list = unique(new_df$seq)
    seq2ymin = setNames(-1*(1:length(seq_list)),seq_list)
    new_df$y_min = seq2ymin[new_df$seq]
    new_df$y_max = new_df$y_min + 0.8
    
    new_df
  })
  
  output$motif_vis <- renderPlot({
    seq(base_pair_concise_annotation())
    
    seq_list = unique(base_pair_concise_annotation()$seq)
    seq_pos = unique(base_pair_concise_annotation()$y_min) + 0.4
    color_list = unique(base_pair_concise_annotation()$color)

    if(input$"x-axis" == 'base-pair'){
      data4plot = base_pair_concise_annotation()
      # print(names(colormap))
      # print(unique(data4plot$motif) %in% names(colormap) ) 
      # print(colormap[unique(data4plot$motif)]) 
      # print(tail(data4plot))
      ggplot() +
        #geom_rect(data = data4plot, aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, fill = color)) +  
        geom_tile(data = data4plot, aes(x = (start + end) / 2, y = (y_min + y_max) / 2, width = end - start, height = y_max - y_min, fill = color)) +
        geom_text(aes(label = seq_list), x = rep(-1,length(seq_list)), y = seq_pos) +
        theme_minimal() + 
        ### labs(title = "Tandem Repeat Units", x = "Position", y = "") +
        scale_fill_manual(values = setNames(color_list, color_list)) +
        coord_cartesian(xlim = input$rgn2show, expand = TRUE) +
        theme(axis.text.y = element_blank(),
              legend.position = "none")  # 隐藏 y 轴标签
    } else {
      
    }

    
  })
  
  
  ###################################################
  # MODULE 3: xxxxxxxxxxxxxxx
  ###################################################
  
 
  
  
}


shinyApp(ui, server)