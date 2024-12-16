library(ggplot2)
library(igraph)
library(visNetwork)
library(stringr)
library(dplyr)

setwd('D:/MyFile/git_repo/TRkScan')

contrast <- c("#D00000","#3185FC","#FFBA08","#5D2E8C","#8FE388","#FF9B85","#8c5c2b","#696663") # 7 colors


filename <- "testData/bo_chr2a_14110848-14463880.concise.tsv"
concise = read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# annotation.tsv

filename <- "testData/bo_chr2a_14110848-14463880.annotation.tsv"
annotation = read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# motif.tsv

filename <- "testData/bo_chr2a_14110848-14463880.motif.tsv"
motif = read.table(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


tmp = motif
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
result$rc = FALSE
dist = result


# distance (rc)
tmp = motif
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
result$rc = TRUE
dist_rc = result





motif_list <- unique(motif$motif)

# 创建空图并添加顶点
g <- make_empty_graph(directed = FALSE)
g <- add_vertices(g, length(motif_list), name = motif_list)

# 计算距离矩阵

distance_df <- rbind(dist, dist_rc)

distance_df <- distance_df[order(distance_df$distance), ]

tmp_edges <- data.frame(from = integer(0), to = integer(0), value = integer(0))

# 添加边并计算连通分量
num_of_motifs = 100
num_components <- components(g)$no
if(num_components > num_of_motifs) {
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
    if(num_components <= num_of_motifs) {
      break
    }
  }
}


# 返回图和边的数据
graph_data = list(graph = g, edges = tmp_edges, motif_list = motif_list, distance_df = distance_df)




color_list = contrast


g <- graph_data$graph
motif_list <- graph_data$motif_list
motif_rep_num <- motif$rep_num

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
nodes = data.frame(id = 1:length(motif_list),
           label = motif_list,
           value = motif_rep_num,
           color = color_list[code_color_id])

edges <- graph_data$edges


visNetwork(nodes, edges, directed = FALSE)

concise$color = 'grey'
for(i in 1:nrow(concise)){
  motif = concise[i,'motif']
  concise[i,'color'] = nodes[nodes$label == motif, 'color']
}

concise$x_mid = (concise$start + concise$end)/2
concise$x_width = concise$end - concise$start

ggplot(data = concise, aes(xmin = start, xmax = end, ymin = 1, ymax = 3, fill = color)) +
  geom_rect() +
  scale_fill_manual(values = c("#D00000" = "#D00000","#3185FC" = "#3185FC","#FFBA08" = "#FFBA08","#5D2E8C" = "#5D2E8C","#8FE388" = "#8FE388",
                               "#FF9B85" = "#FF9B85","#8c5c2b" = "#8c5c2b","#696663" = "#696663")) + 
  theme_void()


