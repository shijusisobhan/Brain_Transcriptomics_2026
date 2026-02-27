rm(list = ls())

# This is the code for plotting behavior (sleep) vs time of day

library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Please set the path where all you data are stored
setwd("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/Brain_Transcriptomics_2026/Data")

file_path  <- "Sleep_data_all.xlsx"
sheet_names <- excel_sheets(file_path)

max_sleep  <- 30
time_hours <- seq(1, 96)

# LD cycle generator
make_ld_cycle <- function(order = c("Light","Dark")) {
  data.frame(
    xmin = c(0,24,48,72),
    xmax = c(24,48,72,96),
    ymin = -5,
    ymax = 0,
    phase = rep(order, length.out = 4)
  )
}

# Highlight region generator

make_highlight <- function(xmin = 0, xmax = 0) {
  data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = 0,
    ymax = 100
  )
}


# Read + compute sleep %

compute_sleep <- function(sheet_name) {

df <- read_excel(file_path, sheet = sheet_name)

df <- df[, -1]
df <- df[, 1:96]

avg_sleep <- colMeans(df, na.rm = TRUE)

data.frame(
  Time = time_hours,
  Sleep = avg_sleep / max_sleep * 100
)
}


# SINGLE plotting function

make_sleep_plot <- function(sheet_index,
                            ld_order,
                            highlight_range) {
  
  sheet_name <- sheet_names[sheet_index]
  
  sleep_df <- compute_sleep(sheet_name)
  
  ld_cycle_df <- make_ld_cycle(ld_order)
  highlight_df <- make_highlight(highlight_range[1],
                                 highlight_range[2])
  
  ggplot() +
    geom_rect(data = ld_cycle_df,
              aes(xmin=xmin,xmax=xmax,
                  ymin=ymin,ymax=ymax,fill=phase),
              color="black") +
    scale_fill_manual(values=c("Light"="white","Dark"="black")) +
    
    geom_rect(data = highlight_df,
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill="grey20", alpha=.3) +
    
    geom_line(data=sleep_df,
              aes(Time,Sleep),
              size=1,color="blue") +
    
    coord_cartesian(ylim=c(-5,100)) +
    scale_y_continuous(breaks=seq(0,100,20)) +
    
    theme_minimal(base_size=14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position="none"
    ) +
    labs(title = sheet_name)
}


generate_plot_list <- function(indices,
ld_order,
highlight_range) {
  
  lapply(indices, make_sleep_plot,
         ld_order = ld_order,
         highlight_range = highlight_range)
}


plot_list1 <- generate_plot_list(
c(1,4,7,10),
ld_order = c("Dark","Light"),
highlight_range = c(0,0)
)

plot_list2 <- generate_plot_list(
  c(16),
  ld_order = c("Light","Dark"),
  highlight_range = c(42,48)
)

plot_list3 <- generate_plot_list(
  c(13,19),
  ld_order = c("Light","Dark"),
  highlight_range = c(0,0)
)


plot_list11 <- generate_plot_list(
  c(2,5,8,11),
  c("Dark","Light"),
  c(66,72)
)

plot_list22 <- generate_plot_list(
  c(17),
  c("Light","Dark"),
  c(36,48)
)

plot_list33 <- generate_plot_list(
  c(14),
  c("Light","Dark"),
  c(48,72)
)


plot_list111 <- generate_plot_list(
  c(3,6,9,12),
  c("Dark","Light"),
  c(66,72)
)

plot_list222 <- generate_plot_list(
  c(18),
  c("Light","Dark"),
  c(24,48)
)

plot_list333 <- generate_plot_list(
  c(15),
  c("Light","Dark"),
  c(48,72)
)


Plot_arrangement <- c(
  plot_list3[2],
  plot_list2[1], plot_list22[1], plot_list222[1],
  plot_list3[1], plot_list33[1], plot_list333[1],
  plot_list1[1], plot_list11[1], plot_list111[1],
  plot_list1[2], plot_list11[2], plot_list111[2],
  plot_list1[3], plot_list11[3], plot_list111[3],
  plot_list1[4], plot_list11[4], plot_list111[4]
)

layout <- rbind(
  c(1, NA, NA),
  matrix(2:19, ncol = 3, byrow = TRUE)
)

grid.arrange(grobs = Plot_arrangement,
             layout_matrix = layout)