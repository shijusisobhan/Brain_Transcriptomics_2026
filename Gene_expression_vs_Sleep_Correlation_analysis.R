
#----------------------------------------------------------------------------
# correltion between gene expression and sleep histoy & Sleep rebound
# This code includes PCA plot, corelation analysis, 
# heat map of correlation, Clustering of correlated samples

#-------------------------------------------------------------------------
rm(list=ls())

library(dplyr)
library(readxl)
library(tibble)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
#library(gplots)
#library(dendextend)

# Please set the path where all you data are stored
# Download the data from git hub and stored in you lcal machine and change the following path

setwd("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/Brain_Transcriptomics_2026/Data")
#-----------------------------------------------------------------------------------------------------------------------------

df<-read.csv('TPM_all_samples.csv') # gene expression data

# Data pre filtering : exclude genes whose TPM<1 in more than 50% samples

#Set threshold
threshold <- 1
max_allowed <- floor(0.5 * 72)  # 50% of 72 samples = 36 → floor to 9

# Subset only the numeric sample columns
tpm_data <- df[-1]

# Count number of samples with TPM < 1 for each gene
low_tpm_counts <- rowSums(tpm_data < threshold)

# Exclude genes with less than 1 TPM in more than 9 samples
filtered_df <- df[low_tpm_counts <= max_allowed, ]
filtered_df<-na.omit(filtered_df)
filtered_df <- filtered_df[!duplicated(filtered_df$ext_gene), ]


rownames(filtered_df) <- filtered_df[, 1]   # set first column as rownames
data <- filtered_df[, -1]            # remove the first column

#------------------------------------------------------------------------------------

# set sample name to a group for better visualization in latter
sample_name<-colnames(data)
# Create the data frame
sample_condition <- data.frame(
  sample = sample_name,
  condition = sub("_[0-9]+$", "", sample_name),
  stringsAsFactors = FALSE
)


condition<-sample_condition$condition

sample_condition<-sample_condition %>% 
  mutate(Exp=ifelse(grepl("SD|29", sample), "Activated", "CTL"))



library(dplyr)

sample_condition <- sample_condition %>%
  mutate(
    
    group = case_when(
      
      # 1. starts with R and ends with TrpA1_29 → R_TrpA1_29
      grepl("TrpA1_29", condition) ~ "GAL4/TrpA1 29",
      
      # 2. contains +29 → Gal4s29
      grepl("\\_29", condition) ~ "GAL4/+29",
      
      # 3. ZT_0 → ZT_0
      grepl("^ZT_0$", condition) ~ "ZT0",
      
      # 4. ZT_anynumber except 0 → ZT
      grepl("^ZT_[1-9][0-9]*$", condition) ~ "ZT4-20",
      
      # 5. contains MechSD → MechSD
      grepl("MechSD", condition) ~ "MechSD",
      
      # 6. contains TrpA1_21 → gal4s_TrpA1_21
      grepl("TrpA1_21", condition) ~ "GAL4/TrpA1 21",
      
      TRUE ~ NA_character_   # default if no rule matches
    )
  )


#-------------------------------------------------------------------------------------------------------
                        ############## PCA plot #######################
#-------------------------------------------------------------------------------------------------------
data_t <- t(data)  # Transpose data so samples are rows

# Perform PCA
pca_result <- prcomp(data_t, scale. = TRUE)
# Create a DataFrame for PCA results
pca_df <- as.data.frame(pca_result$x)
pca_df$condition <- factor(sample_condition$condition)  # Adjust for your genotypes
pca_df$sample <- factor(sample_condition$sample)  # Adjust for your genotypes
pca_df$Exp<- factor(sample_condition$Exp)  # Adjust for your genotypes
pca_df$group<- factor(sample_condition$group)  # Adjust for your genotypes

library(dplyr)
library(ggplot2)

# Define colors

group_colors <- c(
  "GAL4/TrpA1 29" = "darkgreen",
  "GAL4/+29"      = "brown1",
  "ZT0"           = "blue",
  "ZT4-20"        = "magenta",
  "MechSD"        = "cyan",
  "GAL4/TrpA1 21" = "goldenrod4"
)



# PCA plot with labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = condition)) +
  geom_point(size = 4) +                          # Points
  scale_color_manual(values = group_colors) +
  geom_text(vjust = -1, size = 4) +               # Text labels above points
  # geom_label(vjust = -0.5, size = 4) +          # Use this instead if you want boxed labels
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),   # axis label font size
        axis.text  = element_text(size = 18)) +   # tick label font size) + # Removes all grid lines
  labs(
    x = "PC1",
    y = "PC2",
    color = "Group"
  )


ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 18)
  ) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Group"
  )

#---------------------------------------------------------------------------------------------------------------
   # Correlation between gene expression & Sleep history
#---------------------------------------------------------------------------------------------------------------

 # Step-1 - Find the average gene expression across replicates

TPM_log2<-log2(data+1)

# Mean_mat=matrix(0,14145,20)
Mean_mat=matrix(0,nrow(TPM_log2),24)
for (i in 1:24) {
  
  Mean_mat[,i]<-rowMeans(TPM_log2[,(((i-1)*3)+1):(((i-1)*3)+3)], na.rm = F)
  
}

Mean_mat<-as.data.frame(Mean_mat)
rownames(Mean_mat)<-rownames(TPM_log2)
colnames(Mean_mat) <- sample_condition$condition[seq(1, nrow(sample_condition), by = 3)]



#------------- step #2 -Get sleep history-------------


# Path to your Excel file for sleep history which derived from Sleep_data_all.xlsx
file_path <- "Sleep_before_sampling.xlsx"
# Get sheet names
sheets <- excel_sheets(file_path)


# To store results for all sheets
all_results <- list()

for (sheet in sheets) {
  # Read data (first row is header)
  data_SH <- read_excel(file_path, sheet = sheet)
  
  # Remove first column
  data_trimmed <- data_SH[ , -1]
  
  # Ensure it's numeric (remove factors or characters if any)
  data_trimmed <- as.data.frame(lapply(data_trimmed, as.numeric))
  
  # Get number of columns
  n_cols <- ncol(data_trimmed)
  
  # Initialize vector to store means
  cum_means <- numeric()
  
  # Loop for last 2, 4, ..., 24 columns (12 steps)
  for (i in seq(2, 24, by = 2)) {
    # Get the last i columns
    selected_cols <- data_trimmed[ , (n_cols - i + 1):n_cols]
    
    # Row-wise sum
    row_sums <- rowSums(selected_cols, na.rm = TRUE)
    
    # Average across all rows
    cum_means <- c(cum_means, mean(row_sums, na.rm = TRUE))
  }
  
  # Store result for this sheet
  all_results[[sheet]] <- cum_means
}

# Convert results to data frame
result_sleep_history <- do.call(rbind, all_results)
colnames(result_sleep_history) <- paste0(seq(1,12,by=1), "_hr")
rownames(result_sleep_history) <- sheets



# step 4- prepare sleep history data and gene expression data for correlation analysis

sleep_data<-result_sleep_history
row_names<-rownames(Mean_mat)
column_names<-colnames(sleep_data)
gene_data<-t(Mean_mat)
# Check if all rownames are the same
setequal(rownames(gene_data), rownames(sleep_data))
# Row names in df1 but not in df2
missing_in_gene_data <- setdiff(rownames(sleep_data), rownames(gene_data))
missing_in_gene_data

identical(rownames(gene_data), rownames(sleep_data)) # Returns TRUE → if both have exactly the same row names in the same order

sleep_data <- sleep_data[rownames(gene_data), , drop = FALSE]

#- Step 5 pearson correlation analysis -----------------------------------------------------

# Calculate correlations and p-values for each gene and sleep interval
correlation_matrix <- matrix(nrow = ncol(gene_data), ncol = ncol(sleep_data))
p_value_matrix <- matrix(nrow = ncol(gene_data), ncol = ncol(sleep_data))

for (i in 1:ncol(gene_data)) {
  for (j in 1:ncol(sleep_data)) {
    cor_test <- cor.test(gene_data[,i], sleep_data[,j])
    correlation_matrix[i, j] <- cor_test$estimate
    p_value_matrix[i, j] <- cor_test$p.value
  }
}

# Function to calculate adjusted p-values for each row
adjusted_p_values <- apply(p_value_matrix, 1, function(row) {
  p.adjust(row, method = "BH")
})

adjusted_p_values_t<-round(t(adjusted_p_values), digits = 3)



results_df <- data.frame(ext_gene = rep(row_names, each=12),
                         Sleep = rep(column_names, times=ncol(gene_data)),
                         Cor_val =unlist(lapply(1:nrow(correlation_matrix), function(i) correlation_matrix[i,])),
                         t_P = unlist(lapply(1:nrow(p_value_matrix), function(i) p_value_matrix[i,])),
                         BH_Q=unlist(lapply(1:nrow(adjusted_p_values_t), function(i) adjusted_p_values_t[i,]))
)




Correlaton_results<-results_df[which(results_df$BH_Q<0.05),] # select only significant correlations

positive_cor_results<-Correlaton_results[which(Correlaton_results$Cor_val>0),c("ext_gene", "Sleep")]

Negetive_cor_results<-Correlaton_results[which(Correlaton_results$Cor_val<0),c("ext_gene", "Sleep")]


#
sig_cor_genes<-unique(Correlaton_results$ext_gene)
sig_pos_cor<-unique(positive_cor_results$ext_gene)
sig_neg_cor<-unique(Negetive_cor_results$ext_gene)


rownames(correlation_matrix)<-row_names
sig_correlation_matrix<-correlation_matrix[which(rownames(correlation_matrix) %in% sig_cor_genes),]
sig_correlation_matrix_pos<-correlation_matrix[which(rownames(correlation_matrix) %in% sig_pos_cor),]
sig_correlation_matrix_neg<-correlation_matrix[which(rownames(correlation_matrix) %in% sig_neg_cor),]

#-------------------------------------------------------------------------------------------------------

gene_matrix<-t(gene_data)


sig_gene_matrix<-gene_matrix[which(rownames(gene_matrix) %in% sig_cor_genes),]
sig_gene_matrix_pos<-gene_matrix[which(rownames(gene_matrix) %in% sig_pos_cor),]
sig_gene_matrix_neg<-gene_matrix[which(rownames(gene_matrix) %in% sig_neg_cor),]

#--------------------------------------------------------------------------------------------------------------------------------



x_labels <- as.character(1:12)

#Positive correlation plot
sig_correlation_matrix_df<-as.data.frame(sig_correlation_matrix_pos)
sig_correlation_matrix_df$MaxCol <- unname(apply(sig_correlation_matrix_df, 1, which.max))
sig_correlation_matrix_df <- sig_correlation_matrix_df[order(sig_correlation_matrix_df$MaxCol), ]
sig_correlation_matrix_df<-as.matrix(sig_correlation_matrix_df[-13])

heatmap.2(sig_correlation_matrix_df,Rowv=FALSE,Colv=FALSE,col=rev(brewer.pal(9,"RdBu")), trace='none',
          dendrogram = "none",
          labCol = x_labels,
          key = TRUE,
          density.info = "none")

#Negative correlation plot
sig_correlation_matrix_df<-as.data.frame(sig_correlation_matrix_neg)
sig_correlation_matrix_df$MaxCol <- unname(apply(sig_correlation_matrix_df, 1, which.min))
sig_correlation_matrix_df <- sig_correlation_matrix_df[order(sig_correlation_matrix_df$MaxCol), ]
sig_correlation_matrix_df<-as.matrix(sig_correlation_matrix_df[-13])

heatmap.2(sig_correlation_matrix_df,Rowv=FALSE,Colv=FALSE,col=rev(brewer.pal(9,"RdBu")), trace='none',
          dendrogram = "none",
          labCol = x_labels,
          key = TRUE,
          density.info = "none")

#----------------------------------------------------------------------------------------

# Sort in decreasing order of Cor_val and pick top 5 unique ext_gene
corlated_values <- results_df %>%
  arrange(desc(abs(Cor_val))) %>%
  distinct(ext_gene, .keep_all = TRUE) 

mt_genes<- c('CG9377', 'AstA-R2') 

mt_corlated_values<-corlated_values[which(corlated_values$ext_gene %in% mt_genes),]



plots <- list()

# Loop over each gene in top5_high
for (i in 1:nrow(mt_corlated_values)) {
  gene <- mt_corlated_values$ext_gene[i]
  sleep_interval <- mt_corlated_values$Sleep[i]
  
  r_val <- round(mt_corlated_values$Cor_val[i], 3)
  q_val <- formatC(mt_corlated_values$BH_Q[i],3)  # scientific notation
  
  
  # Extract sleep values for the given interval
  sleep_vec <- sleep_data[,sleep_interval]
  
  # Extract gene expression values
  expr_vec <- as.numeric(gene_data[,gene])
  
  # Combine into a data frame for plotting
  plot_df <- tibble(
    Name = rownames(gene_data),
    Sleep = sleep_vec,
    Expression = expr_vec
  )
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = Sleep, y = Expression, color = Name)) +
    geom_point(size = 3) +
    #geom_text(aes(label = Name), vjust = -1, size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.7) +
    ggtitle(paste0(gene),
            subtitle = paste0("r = ", r_val, ", q = ", q_val)) +
    theme_minimal() +
    theme(
      legend.position = "none",   # keep legend visible
      panel.grid = element_blank()
    ) +
    xlab("sleep history (min)") +
    ylab("Gene Expression (log2)")
  
  plots[[i]] <- p
}


# Combine all plots into one panel with 3 rows
final_plot <- wrap_plots(plots, nrow = 2)
print(final_plot)

#-----------------------------------------------------------------------------------------------------------------
                               #Cluster Heat Map Sleep History Genes
#------------------------------------------------------------------------------------------------------------------




gene_matrix<-TPM_log2

sig_gene_matrix<-gene_matrix[which(rownames(gene_matrix) %in% sig_cor_genes),]
sig_gene_matrix_pos<-gene_matrix[which(rownames(gene_matrix) %in% sig_pos_cor),]
sig_gene_matrix_neg<-gene_matrix[which(rownames(gene_matrix) %in% sig_neg_cor),]


# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")





#--------------------------------------------------------


all_genes <- rownames(sig_gene_matrix_z)


gene_sign <- ifelse(all_genes %in% sig_pos_cor, "pos",
                    ifelse(all_genes %in% sig_neg_cor, "neg", NA))


gene_clusters <- ifelse(gene_sign == "pos", 1,
                        ifelse(gene_sign == "neg", 2, NA))

row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z))), method = "complete")
row_dend <- as.dendrogram(row_clust)

row_dend_colored <- color_branches(
  row_dend,
  clusters = gene_clusters,
  col = c("red", "blue"),  # 1 = pos, 2 = neg
  height = max(height(row_dend)) * 0   # color only lowest 10%
)


# Remove .number suffix from column names
colnames(sig_gene_matrix_z) <- sub("\\.\\d+$", "", colnames(sig_gene_matrix_z))



# assign colors based on group
col_side_colors <- group_colors[sample_condition$group]

# FIXED (use sample instead of condition)
col_side_colors <- col_side_colors[
  match(colnames(sig_gene_matrix_z), sample_condition$sample)
]



clean_labels <- sub("\\.\\d+$", "", colnames(sig_gene_matrix_z))

P1 <- heatmap.2(
  sig_gene_matrix_z,
  Rowv = row_dend_colored,
  Colv = as.dendrogram(col_clust),
  
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  labCol = clean_labels,
  ColSideColors = col_side_colors,
  
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = -0.9,
  cexCol = 0.7,
  
  density.info = "none",
  key = TRUE,
  key.title = "",
  key.xlab = ""
)


#---------------------------------------------------------------------------------------------------------------------------
                       # Correlation between gene expression & Sleep Rebound
#--------------------------------------------------------------------------------------------------------------------------

result_sleep_rebound<-read.csv("Average_Cumulative_Rebound.csv", row.names = 1)


# step 4- prepare sleep rebound data and gene expression data for correlation analysis

sleep_data_reb <- result_sleep_rebound[!grepl("^(ZT)", rownames(result_sleep_rebound)), ]
gene_data_reb <- Mean_mat[, !grepl("^(ZT)", colnames(Mean_mat))]
row_names<-rownames(gene_data_reb)
column_names<-colnames(sleep_data_reb)
gene_data_reb<-t(gene_data_reb)


# Check if all rownames are the same
setequal(rownames(gene_data_reb), rownames(sleep_data_reb))
# Row names in df1 but not in df2
missing_in_gene_data_reb <- setdiff(rownames(gene_data_reb), rownames(sleep_data_reb))
missing_in_gene_data_reb

identical(rownames(gene_data_reb), rownames(sleep_data_reb)) # Returns TRUE → if both have exactly the same row names in the same order

# Make sure both have exactly the same row names in the same order
sleep_data_reb <- sleep_data_reb[rownames(gene_data_reb), , drop = FALSE]

identical(rownames(gene_data_reb), rownames(sleep_data_reb)) # Returns TRUE → if both have exactly the same row names in the same order





#- Step 5 pearson correlation analysis -----------------------------------------------------

# Calculate correlations and p-values for each gene and sleep interval
correlation_matrix <- matrix(nrow = ncol(gene_data_reb), ncol = ncol(sleep_data_reb))
p_value_matrix <- matrix(nrow = ncol(gene_data_reb), ncol = ncol(sleep_data_reb))

for (i in 1:ncol(gene_data_reb)) {
  for (j in 1:ncol(sleep_data_reb)) {
    cor_test <- cor.test(gene_data_reb[,i], sleep_data_reb[,j])
    correlation_matrix[i, j] <- cor_test$estimate
    p_value_matrix[i, j] <- cor_test$p.value
  }
}

# Function to calculate adjusted p-values for each row
adjusted_p_values <- apply(p_value_matrix, 1, function(row) {
  p.adjust(row, method = "BH")
})

adjusted_p_values_t<-round(t(adjusted_p_values), digits = 3)



results_df <- data.frame(ext_gene = rep(row_names, each=12),
                         Sleep = rep(column_names, times=ncol(gene_data_reb)),
                         Cor_val =unlist(lapply(1:nrow(correlation_matrix), function(i) correlation_matrix[i,])),
                         t_P = unlist(lapply(1:nrow(p_value_matrix), function(i) p_value_matrix[i,])),
                         BH_Q=unlist(lapply(1:nrow(adjusted_p_values_t), function(i) adjusted_p_values_t[i,]))
)




Correlaton_results<-results_df[which(results_df$BH_Q<0.05),] # select only significant correlations

positive_cor_results<-Correlaton_results[which(Correlaton_results$Cor_val>0),c("ext_gene", "Sleep")]

Negetive_cor_results<-Correlaton_results[which(Correlaton_results$Cor_val<0),c("ext_gene", "Sleep")]


#
sig_cor_genes<-unique(Correlaton_results$ext_gene)
sig_pos_cor<-unique(positive_cor_results$ext_gene)
sig_neg_cor<-unique(Negetive_cor_results$ext_gene)


rownames(correlation_matrix)<-row_names
sig_correlation_matrix<-correlation_matrix[which(rownames(correlation_matrix) %in% sig_cor_genes),]
sig_correlation_matrix_pos<-correlation_matrix[which(rownames(correlation_matrix) %in% sig_pos_cor),]
sig_correlation_matrix_neg<-correlation_matrix[which(rownames(correlation_matrix) %in% sig_neg_cor),]

#-------------------------------------------------------------------------------------------------------

gene_matrix<-t(gene_data_reb)


sig_gene_matrix<-gene_matrix[which(rownames(gene_matrix) %in% sig_cor_genes),]
sig_gene_matrix_pos<-gene_matrix[which(rownames(gene_matrix) %in% sig_pos_cor),]
sig_gene_matrix_neg<-gene_matrix[which(rownames(gene_matrix) %in% sig_neg_cor),]

#--------------------------------------------------------------------------------------------------------------------------------

library('gplots')
library(RColorBrewer)

x_labels <- as.character(1:12)

#Positive correlation plot
sig_correlation_matrix_df<-as.data.frame(sig_correlation_matrix_pos)
sig_correlation_matrix_df$MaxCol <- unname(apply(sig_correlation_matrix_df, 1, which.max))
sig_correlation_matrix_df <- sig_correlation_matrix_df[order(sig_correlation_matrix_df$MaxCol), ]
sig_correlation_matrix_df<-as.matrix(sig_correlation_matrix_df[-13])

heatmap.2(sig_correlation_matrix_df,Rowv=FALSE,Colv=FALSE,col=rev(brewer.pal(9,"RdBu")), trace='none',
          dendrogram = "none",
          labCol = x_labels,
          key = TRUE,
          density.info = "none")

#Negative correlation plot
sig_correlation_matrix_df<-as.data.frame(sig_correlation_matrix_neg)
sig_correlation_matrix_df$MaxCol <- unname(apply(sig_correlation_matrix_df, 1, which.min))
sig_correlation_matrix_df <- sig_correlation_matrix_df[order(sig_correlation_matrix_df$MaxCol), ]
sig_correlation_matrix_df<-as.matrix(sig_correlation_matrix_df[-13])

heatmap.2(sig_correlation_matrix_df,Rowv=FALSE,Colv=FALSE,col=rev(brewer.pal(9,"RdBu")), trace='none',
          dendrogram = "none",
          labCol = x_labels,
          key = TRUE,
          density.info = "none")

#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------

# Sort in decreasing order of Cor_val and pick top 5 unique ext_gene
corlated_values <- results_df %>%
  arrange(desc(abs(Cor_val))) %>%
  distinct(ext_gene, .keep_all = TRUE) 


mt_corlated_values<-corlated_values[which(corlated_values$ext_gene %in% mt_genes),]



plots <- list()

# Loop over each gene in top5_high
for (i in 1:nrow(mt_corlated_values)) {
  gene <- mt_corlated_values$ext_gene[i]
  sleep_interval <- mt_corlated_values$Sleep[i]
  
  r_val <- round(mt_corlated_values$Cor_val[i], 3)
  q_val <- formatC(mt_corlated_values$BH_Q[i],3)  # scientific notation
  
  
  # Extract sleep values for the given interval
  sleep_vec <- sleep_data_reb[,sleep_interval]
  
  # Extract gene expression values
  expr_vec <- as.numeric(gene_data_reb[,gene])
  
  # Combine into a data frame for plotting
  plot_df <- tibble(
    Name = rownames(gene_data_reb),
    Sleep = sleep_vec,
    Expression = expr_vec
  )
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = Sleep, y = Expression, color = Name)) +
    geom_point(size = 3) +
    #geom_text(aes(label = Name), vjust = -1, size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.7) +
    ggtitle(paste0(gene),
            subtitle = paste0("r = ", r_val, ", q = ", q_val)) +
    theme_minimal() +
    theme(
      legend.position = "none",   # keep legend visible
      panel.grid = element_blank()
    ) +
    xlab("sleep history (min)") +
    ylab("Gene Expression (log2)")
  
  plots[[i]] <- p
}


# Combine all plots into one panel with 3 rows
final_plot <- wrap_plots(plots, nrow = 2)
print(final_plot)


#-----------------------------------------------------------------------------------------------------------------
#Cluster Heat Map Sleep rebound correlted genes Genes
#------------------------------------------------------------------------------------------------------------------




gene_matrix<-TPM_log2

sig_gene_matrix<-gene_matrix[which(rownames(gene_matrix) %in% sig_cor_genes),]
sig_gene_matrix_pos<-gene_matrix[which(rownames(gene_matrix) %in% sig_pos_cor),]
sig_gene_matrix_neg<-gene_matrix[which(rownames(gene_matrix) %in% sig_neg_cor),]


# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")





#---------------------------------------------------------------------------------------


all_genes <- rownames(sig_gene_matrix_z)


gene_sign <- ifelse(all_genes %in% sig_pos_cor, "pos",
                    ifelse(all_genes %in% sig_neg_cor, "neg", NA))


gene_clusters <- ifelse(gene_sign == "pos", 1,
                        ifelse(gene_sign == "neg", 2, NA))

row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z))), method = "complete")
row_dend <- as.dendrogram(row_clust)

row_dend_colored <- color_branches(
  row_dend,
  clusters = gene_clusters,
  col = c("red", "blue"),  # 1 = pos, 2 = neg
  height = max(height(row_dend)) * 0   # color only lowest 10%
)


# Remove .number suffix from column names
colnames(sig_gene_matrix_z) <- sub("\\.\\d+$", "", colnames(sig_gene_matrix_z))



# assign colors based on group
col_side_colors <- group_colors[sample_condition$group]

# FIXED (use sample instead of condition)
col_side_colors <- col_side_colors[
  match(colnames(sig_gene_matrix_z), sample_condition$sample)
]



clean_labels <- sub("\\.\\d+$", "", colnames(sig_gene_matrix_z))

P1 <- heatmap.2(
  sig_gene_matrix_z,
  Rowv = row_dend_colored,
  Colv = as.dendrogram(col_clust),
  
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  labCol = clean_labels,
  ColSideColors = col_side_colors,
  
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = -0.9,
  cexCol = 0.7,
  
  density.info = "none",
  key = TRUE,
  key.title = "",
  key.xlab = ""
)
