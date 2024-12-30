library(shiny)
library(stringr)
library(tidyr)
library(dplyr)
library(future)
library(promises)
# plot package
library(ggplot2)
library(igraph)
library(visNetwork)
library(htmlwidgets)

setwd("D:/MyFile/git_repo/TRkScan/")
source("./plotScripts/triangular_heatmap.R")

options(shiny.maxRequestSize=30*1024^2)

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
      checkboxInput("merge_rc", "merge_rc", FALSE),
      downloadButton("dl_cluster_df", "save cluster result")
    ),
    column(8,
      visNetworkOutput("network"),
      downloadButton("dl_network", "save")
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
       ),
       checkboxInput("align", "align", FALSE)
    ),
    column(8, 
       plotOutput("motif_vis"),
       downloadButton("dl_anno", "save")
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
        plotOutput("id_heatmap")
    )
  ),
)

server <- function(input, output, session) {

  contrast <- c("lightgrey","#3185FC","#D00000","#FFBA08","#5D2E8C","#8FE388","#FF9B85","#8c5c2b","#696663") # 7 colors
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
  
  output$dl_cluster_df <- downloadHandler(
    filename = function() {
      paste("dataframe-", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      write.table(nodes(), file, sep = '\t', row.names = FALSE)
    }
  )
  
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

    print('graph done!')

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

    print(tmp)
    
    old_cluster_id <- data.frame(old_cluster_id  = comp$membership)
    code_color_id <- old_cluster_id %>%
      left_join(tmp, by = "old_cluster_id") %>%
      mutate(color_id = ifelse(new_cluster_id < length(color_list), new_cluster_id, length(color_list))) %>%
      pull(color_id)
    
    print(code_color_id)
    # decide direction
    dir_list = c('+')
    dist = dist()
    motif_list = motif()$id

    if (length(motif_list) > 1){
      for (i in 2:length(motif_list)){
        dir = ''
        for (j in 1:(i-1)) {
          if (code_color_id[i] != code_color_id[j]){
            next
          }
          ### print(dist[dist$ref == motif_list[i] & dist$query == motif_list[j], ])
          if (dist[dist$ref == motif_list[i] & dist$query == motif_list[j],][1, 'is_rc']){
            dir = ifelse(dir_list[j] == '+', '-', '+')
            break
          } else {
            dir = ifelse(dir_list[j] == '+', '+', '-')
            break
          }
        }
      }
      dir_list = c(dir_list, c(ifelse(dir != '', dir, '+')))

    }
    
    print('nodes done!')
    
    data.frame(id = motif()$id,
               label = motif()$label,
               motif = motif()$motif,
               value = motif()$rep_num,
               color = color_list[code_color_id],
               dir = dir_list)
  })

  edges <- reactive({
    req(graph_data())  # use shared graph_data
    print('edges done!')
    graph_data()$edges
  })

  networkPlot <- reactive({
    req(nodes())
    req(edges())
    visNetwork(nodes(), edges(), directed = FALSE) %>%
      visNodes(font = list(size = 20)) %>%
      visEdges(font = list(size = 15)) %>%
      visOptions(highlightNearest = TRUE)
  })
  
  output$network <- renderVisNetwork({
    req(networkPlot())
    networkPlot()
  })
  
  output$dl_network <- downloadHandler(
    filename = function() {
      paste("network-", Sys.Date(),".html", sep = "")
    },
    content = function(file) {
      saveWidget(networkPlot(),
                 file,
                 selfcontained = TRUE)  # 保存为自包含的 HTML 文件
    }
  )
  
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
  
  
 # solve split motif
  split_cigar <- function(row, name) {
    result <- data.frame()  # seq, bp_x, height, motif_id, name
    
    seq = row[1, 'seq']
    start = row[1, 'start']
    cigar = row[1, 'CIGAR']
    motif_id = motif2id()[row[1, 'motif']]
    cur = start
    length = 0
    sub_name = 1
    
    result <- rbind(result, data.frame(
      seq = seq,
      bp_x = cur,
      height = length,
      motif_id = motif_id,
      name = paste(name, as.character(sub_name), sep = '-')
    ))
    

    num_str = ''
    for (i in 1:nchar(cigar)){
      symbol = substring(cigar, i, i)
      if (symbol %in% c('=','X','I') ){
        cur = cur + as.integer(num_str)
        length = length + as.integer(num_str)
        num_str = ''
      } else if (symbol == 'D'){
        num_str = ''
      } else if (symbol == 'N'){
        result <- rbind(result, data.frame(
          seq = seq,
          bp_x = cur,
          height = length,
          motif_id = motif_id,
          name = paste(name, as.character(sub_name), sep = '-')
        ))
        sub_name = sub_name + 1
        cur = cur + as.integer(num_str)
        num_str = ''
        if (i != nchar(cigar)){
          result <- rbind(result, data.frame(
            seq = seq,
            bp_x = cur,
            height = length,
            motif_id = motif_id,
            name = paste(name, as.character(sub_name), sep = '-')
          ))
        }

      } else {  # 0-9
        num_str = paste0(num_str, symbol)
      }
    }
    if (symbol != 'N'){
      result <- rbind(result, data.frame(
        seq = seq,
        bp_x = cur,
        height = length,
        motif_id = motif_id,
        name = paste(name, as.character(sub_name), sep = '-')
      ))
    }

    result$height = (length - result$height) / length

    return(result)
  }
  
  # get primary data4plot
  annoData4plot <- reactive({
    seq(annotation())
    df = annotation()

    df_noedit = df[! grepl('N', df$CIGAR), c('seq', 'start', 'end', 'motif')]

    new_df = data.frame(
      seq = rep(df_noedit$seq, 2),
      bp_x = c(df_noedit$start, df_noedit$end),
      height = c(rep(1, nrow(df_noedit)), rep(0, nrow(df_noedit))),
      motif_id = motif2id()[rep(df_noedit$motif, 2)],
      name = rep(as.character(1:nrow(df_noedit)), 2)
    )

    cur_name = nrow(new_df) + 1
    # make primary data for annotation plot
    df_edit = df[grepl('N', df$CIGAR), ]
    row.names(df_edit) = 1:nrow(df_edit)
    for (i in 1:nrow(df_edit)) {
      tmp = split_cigar(df_edit[i, ], cur_name)
      new_df = rbind(new_df, tmp)
      cur_name = cur_name + 1
    }
    rownames(new_df) = 1:nrow(new_df)
    
    new_df
  })
  
  # get final data4plot
  annoData4plot_final <- reactive({
    seq(annoData4plot)
    seq(nodes())
    
    df = annoData4plot()
    
    # add color
    colormap = setNames(nodes()$color, as.character(nodes()$id))
    df$color = colormap[as.character(df$motif_id)]
    # add direction
    print(nodes()$dir)
    dirmap = setNames(nodes()$dir, as.character(nodes()$id))
    df$dir = dirmap[as.character(df$motif_id)]
    # give y axis
    seq_list = unique(df$seq)
    seq2y = setNames(-1*(1:length(seq_list)), seq_list)


    # make final data based on direction
    df <- df %>%
      uncount(2) %>%
      mutate(height = ifelse(dir == '+', height, 1 - height)) %>%
      mutate(y = ifelse(row_number() %% 2 == 1, seq2y[seq] - height * 0.8, seq2y[seq] + height * 0.8))
    
   
    # add 
    # # combine same line
    # df = new_df
    # new_df = data.frame()
    # for (i in 1:(nrow(df)-1)){
    #   if (i == 1){
    #     start = df[1,'start']
    #     end = df[1,'end']
    #   }
    #   if (df[i,'color'] == df[i+1, 'color'] && df[i,'end'] == df[i+1,'start']){
    #     end = df[i+1,'end']
    #   } else {
    #     new_df <- rbind(new_df, data.frame(
    #       seq = df[i,'seq'],
    #       length = df[i,'seq'],
    #       start = start,
    #       end = end,
    #       color = df[i,'color']
    #     ))
    #     start = df[i+1,'start']
    #     end = df[i+1,'end']
    #   }
    # }
    # new_df <- rbind(new_df, data.frame(
    #   seq = df[i,'seq'],
    #   length = df[i,'seq'],
    #   start = start,
    #   end = end,
    #   color = df[i,'color']
    # ))
    df <- df[!duplicated(df), ]
    print('final annotation data4plot done!')
    print(df)
    df
  })
  
  annoPlot <- reactive({
    seq(annoData4plot_final())
    
    seq_list = unique(annoData4plot_final()$seq)
    seq_pos = -1 * (1:length(seq_list))
    color_list = unique(annoData4plot_final()$color)
    
    if(input$"x-axis" == 'base-pair'){
      data4plot = annoData4plot_final()
      ggplot() +
        geom_polygon(data = data4plot, aes(x = bp_x, y = y, group = name, fill = color)) +
        #geom_rect(data = data4plot, aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, fill = color)) +  
        #geom_tile(data = data4plot, aes(x = (start + end) / 2, y = seq, width = end - start, height = y_max - y_min, fill = color)) +
        geom_text(aes(label = seq_list), x = rep(-1,length(seq_list)), y = seq_pos, size = 10, hjust = 0) +
        theme_minimal() + 
        ### labs(title = "Tandem Repeat Units", x = "Position", y = "") +
        scale_fill_manual(values = setNames(color_list, color_list)) +
        coord_cartesian(xlim = input$rgn2show, expand = FALSE) +
        xlab('base pair (bp)') +
        ylab('sequences') +
        theme(axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14, color = 'black'),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title = element_text(size = 20),
              legend.position = "none")  # 隐藏 y 轴标签
      
    } else {
      
    }
  }) 
  
  output$motif_vis <- renderPlot({
    req(annoPlot())
    annoPlot()
  })
  
  output$dl_anno <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = annoPlot(), device = "pdf", width = 10, height = 6)
    }
  )
  
  ###################################################
  # MODULE 3: xxxxxxxxxxxxxxx
  ###################################################
  output$id_heatmap <- renderPlot({
    
  })
 
  
  
}


shinyApp(ui, server)