#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


.diamond <- function(
    row
) {
  
  side_length <- as.numeric(row[".window"])
  x <- as.numeric(row[".w"])
  y <- as.numeric(row[".z"])
  
  base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) * sqrt(2) / 2
  trans <- (base * side_length) + c(x, y)
  df <- as.data.frame(t(trans))
  colnames(df) <- c(".w", ".z")
  df$.color <- row[".color"]
  df$.group <- as.numeric(row[".group"])
  
  return(df)
}


triangular_heatmap <- function(
    df,
    reference_start_col = "reference_start",
    reference_end_col = "reference_end",
    query_start_col = "query_start",
    query_end_col = "query_end",
    identity_col = "identity",
    color_col = "auto",
    coordinate_boundary = c(NA, NA),
    title = NULL
) {
  
  #'
  #' Generate a triangular heatmap based on ggplot2
  #' @param color_col if `color_col` is `auto`, the color gradient will be determined by the `identity_col`.
  #' @return a ggplot2 object
  #'
  
  sdf <- copy(df)

  if (color_col == "auto") {
    sdf <- sdf %>%
      mutate(
        .color = case_when(
          !!sym(identity_col) >= 70 & !!sym(identity_col) < 80 ~ "#873f82",
          !!sym(identity_col) >= 80 & !!sym(identity_col) < 90 ~ "#3f6698",
          !!sym(identity_col) >= 90 & !!sym(identity_col) < 91 ~ "#0098c9",
          !!sym(identity_col) >= 91 & !!sym(identity_col) < 92 ~ "#36b3a6",
          !!sym(identity_col) >= 92 & !!sym(identity_col) < 93 ~ "#90cd9e",
          !!sym(identity_col) >= 93 & !!sym(identity_col) < 94 ~ "#b6dd8e",
          !!sym(identity_col) >= 94 & !!sym(identity_col) < 95 ~ "#fddf86",
          !!sym(identity_col) >= 95 & !!sym(identity_col) < 96 ~ "#f8ac68",
          !!sym(identity_col) >= 96 & !!sym(identity_col) < 97 ~ "#f17746",
          !!sym(identity_col) >= 97 & !!sym(identity_col) < 98 ~ "#d45b59",
          !!sym(identity_col) >= 98 & !!sym(identity_col) < 99 ~ "#b23363",
          !!sym(identity_col) >= 99 & !!sym(identity_col) <= 100 ~ "#b90064",
          TRUE ~ NA_character_
        )
      ) %>%
      na.omit()
  } else {
    sdf <- sdf %>% mutate(.color = !!sym(color_col))
  }
  
  if (is.na(coordinate_boundary[1])) {
    coordinate_boundary[1] <- min(sdf[[query_start_col]], sdf[[reference_start_col]])
  }
  if (is.na(coordinate_boundary[2])) {
    coordinate_boundary[2] <- max(sdf[[query_end_col]], sdf[[reference_end_col]])
  }
  
  window <- max(sdf[[query_end_col]] - sdf[[query_start_col]])
  sdf$.first_pos <- sdf[[query_start_col]] / window
  sdf$.second_pos <- sdf[[reference_start_col]] / window
  
  sdf$.w <- sdf$.first_pos + sdf$.second_pos
  sdf$.z <- - sdf$.first_pos + sdf$.second_pos
  triangle_scale <- max(sdf[[query_start_col]]) / max(sdf$.w)
  sdf$.window <- max(sdf[[query_end_col]] - sdf[[query_start_col]]) / triangle_scale
  sdf$.group <- seq(nrow(sdf))
  
  df.d <- rbindlist(apply(sdf, 1, .diamond))
  
  plot <- ggplot(df.d) +
    geom_polygon(mapping = aes(x = .w * triangle_scale, y = .z * window, group = .group, fill = .color)) +
    scale_x_continuous(limits = c(NA, NA), expand = c(0, 0), name = NULL, breaks = coordinate_boundary, labels = scales::comma) +
    scale_y_continuous(limits = c(NA, NA), expand = c(0, 0), name = NULL, breaks = coordinate_boundary, labels = scales::comma) +
    coord_cartesian(xlim = coordinate_boundary) +
    scale_fill_identity() +
    ggtitle(title) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(plot)
  
}
