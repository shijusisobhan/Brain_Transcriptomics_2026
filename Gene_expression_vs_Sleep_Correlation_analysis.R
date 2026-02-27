rm(list=ls())

library(dplyr)
library(readxl)

# Please set the path where all you data are stored
setwd("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/Brain_Transcriptomics_2026/Data")

df<-read.csv('TPM_all_samples.csv')


#Set threshold
threshold <- 1
max_allowed <- floor(0.5 * 72)  # 50% of 72 samples = 36 → floor to 9

#-------------------------------------------------------------------------------
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

# Sort in decreasing order of Cor_val and pick top 5 unique ext_gene
corlated_values <- results_df %>%
  arrange(desc(abs(Cor_val))) %>%
  distinct(ext_gene, .keep_all = TRUE) 

mt_genes<- c('Trh', 'AstA-R2') 

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
sig_gene_matrix_df <- as.matrix(sig_gene_matrix_pos)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")

# Plot heatmap
P1 <- heatmap.2(
  sig_gene_matrix_z,
  
  # --- Clustering ---
  Rowv = as.dendrogram(row_clust),
  Colv = as.dendrogram(col_clust),
  #dendrogram = "none",  # hide dendrogram
  
  # --- Colors and trace ---
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  # --- X-axis labels ---
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = 0.5,
  cexCol = 0.7,
  
  # --- Color key without histogram ---
  key = TRUE,
  density.info = "none",
  key.title = "",
  key.xlab = ""
)

#------------------------------------------------------------------------------

# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix_neg)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")

# Plot heatmap
P1 <- heatmap.2(
  sig_gene_matrix_z,
  
  # --- Clustering ---
  Rowv = as.dendrogram(row_clust),
  Colv = as.dendrogram(col_clust),
  #dendrogram = "none",  # hide dendrogram
  
  # --- Colors and trace ---
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  # --- X-axis labels ---
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = 0.5,
  cexCol = 0.7,
  
  # --- Color key without histogram ---
  key = TRUE,
  density.info = "none",
  key.title = "",
  key.xlab = ""
)


#-----------------------------------------------------------


# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")



library(gplots)
library(dendextend)

# Convert hclust to dendrogram
row_dend <- as.dendrogram(row_clust)

# Color 3 major clusters
row_dend_colored <- color_branches(row_dend, k = 5)

# Optional: choose specific colors
# row_dend_colored <- color_branches(row_dend, k = 3,
#                                    col = c("red", "blue", "darkgreen"))

# Plot heatmap with colored row dendrogram
P1 <- heatmap.2(
  sig_gene_matrix_z,
  
  Rowv = row_dend_colored,
  Colv = as.dendrogram(col_clust),
  
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = 0.5,
  cexCol = 0.7,
  
  key = TRUE,
  density.info = "none",
  key.title = "",
  key.xlab = ""
)


#--------------------------------------------------------


all_genes <- rownames(sig_gene_matrix_z)


gene_sign <- ifelse(all_genes %in% sig_pos_cor, "pos",
                    ifelse(all_genes %in% sig_neg_cor, "neg", NA))


gene_clusters <- ifelse(gene_sign == "pos", 1,
                        ifelse(gene_sign == "neg", 2, NA))



library(dendextend)

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

#----------------------------------------------------------------
#------------------------------------------------------------------




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


# Open a new blank plot for legend in new plot window

# plot.new()
# 
# # Add legend only
# legend(
#   "center",
#   legend = names(group_colors),
#   fill = group_colors,
#   border = NA,
#   bty = "n",
#   cex = 1.2
# )



#---------------------------------------------------------------------------------------------------------------
                       # Correlation between gene expression & Sleep Rebound
#---------------------------------------------------------------------------------------------------------------


# File paths
file_path <- "Sleep_after_sampling.xlsx"
file_path_2 <- "Sleep_baseline.xlsx"

# Get sheet names from both Excel files
sheets_1 <- excel_sheets(file_path)
sheets_2 <- excel_sheets(file_path_2)

# Print the sheet names
cat("Sheets in first file:\n")
print(sheets_1)
cat("\nSheets in second file:\n")
print(sheets_2)

# Check if names and order are identical
if (identical(sheets_1, sheets_2)) {
  cat("\n✅ The two files have identical sheet names in the same order.\n")
} else {
  cat("\n❌ The sheet names or order differ between the two files.\n")
  
  # Show differences if any
  cat("\nSheets only in first file:\n")
  print(setdiff(sheets_1, sheets_2))
  
  cat("\nSheets only in second file:\n")
  print(setdiff(sheets_2, sheets_1))
  
  cat("\nOrder differences:\n")
  print(data.frame(File1 = sheets_1, File2 = sheets_2))
}

#-------------------------------------------------------------------------------------------------------



# To store results for all sheets
all_results <- list()

for (sheet in sheets_1) {
  # Read data (first row is header)
  data_1 <- read_excel(file_path, sheet = sheet)
  data_2 <- read_excel(file_path_2, sheet = sheet)
  
  # Remove first column
  data_trimmed_1 <- data_1[ , -1]
  data_trimmed_2 <- data_2[ , -1]
  
  # Ensure it's numeric (remove factors or characters if any)
  data_trimmed_1 <- as.data.frame(lapply(data_trimmed_1, as.numeric))
  data_trimmed_2 <- as.data.frame(lapply(data_trimmed_2, as.numeric))
  
  # Get number of columns
  n_cols <- ncol(data_trimmed_1)
  
  # Initialize vector to store means
  rebound_means <- numeric()
  
  # Loop for first 2, 4, ..., 24 columns (12 steps)
  for (i in seq(2, 24, by = 2)) {
    # Get the first i columns
    selected_cols_1 <- data_trimmed_1[, 1:i]
    selected_cols_2 <- data_trimmed_2[, 1:i]
    
    # Row-wise sum
    row_sums_1 <- rowSums(selected_cols_1, na.rm = TRUE)
    row_sums_2 <- rowSums(selected_cols_2, na.rm = TRUE)
    rebound<-row_sums_1-row_sums_2 # sleep_after_sampling - baseline sampling
    
    # Average across all rows
    rebound_means <- c(rebound_means, mean(rebound, na.rm = TRUE))
  }
  
  # Store result for this sheet
  all_results[[sheet]] <- rebound_means
}

# Convert results to data frame
result_sleep_rebound <- do.call(rbind, all_results)
colnames(result_sleep_rebound) <- paste0(seq(1,12,by=1), "_hr")
rownames(result_sleep_rebound) <- sheets_1


# step 4- prepare sleep rebound data and gene expression data for correlation analysis

sleep_data <- result_sleep_rebound[!grepl("^(ZT)", rownames(result_sleep_rebound)), ]
gene_data <- Mean_mat[, !grepl("^(ZT)", colnames(Mean_mat))]
row_names<-rownames(gene_data)
column_names<-colnames(sleep_data)
gene_data<-t(gene_data)


# Check if all rownames are the same
setequal(rownames(gene_data), rownames(sleep_data))
# Row names in df1 but not in df2
missing_in_gene_data <- setdiff(rownames(gene_data), rownames(sleep_data))
missing_in_gene_data

identical(rownames(gene_data), rownames(sleep_data)) # Returns TRUE → if both have exactly the same row names in the same order

# Make sure both have exactly the same row names in the same order
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


gene_matrix <- gene_matrix[, !grepl("^(ZT)", colnames(gene_matrix))]


sig_gene_matrix<-gene_matrix[which(rownames(gene_matrix) %in% sig_cor_genes),]
sig_gene_matrix_pos<-gene_matrix[which(rownames(gene_matrix) %in% sig_pos_cor),]
sig_gene_matrix_neg<-gene_matrix[which(rownames(gene_matrix) %in% sig_neg_cor),]


# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix_pos)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")

# Plot heatmap
P1 <- heatmap.2(
  sig_gene_matrix_z,
  
  # --- Clustering ---
  Rowv = as.dendrogram(row_clust),
  Colv = as.dendrogram(col_clust),
  #dendrogram = "none",  # hide dendrogram
  
  # --- Colors and trace ---
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  # --- X-axis labels ---
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = 0.5,
  cexCol = 0.7,
  
  # --- Color key without histogram ---
  key = TRUE,
  density.info = "none",
  key.title = "",
  key.xlab = ""
)

#------------------------------------------------------------------------------

# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix_neg)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")

# Plot heatmap
P1 <- heatmap.2(
  sig_gene_matrix_z,
  
  # --- Clustering ---
  Rowv = as.dendrogram(row_clust),
  Colv = as.dendrogram(col_clust),
  #dendrogram = "none",  # hide dendrogram
  
  # --- Colors and trace ---
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  # --- X-axis labels ---
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = 0.5,
  cexCol = 0.7,
  
  # --- Color key without histogram ---
  key = TRUE,
  density.info = "none",
  key.title = "",
  key.xlab = ""
)


#-----------------------------------------------------------


# Convert to matrix if not already
sig_gene_matrix_df <- as.matrix(sig_gene_matrix)

# Row-wise z-score: (x - mean) / sd
sig_gene_matrix_z <- t(scale(t(sig_gene_matrix_df)))


# Compute row and column clustering (based on z-scores)
row_clust <- hclust(as.dist(1 - cor(t(sig_gene_matrix_z), method = "pearson")),
                    method = "complete")
col_clust <- hclust(as.dist(1 - cor(sig_gene_matrix_z, method = "pearson")),
                    method = "complete")



library(gplots)
library(dendextend)

# Convert hclust to dendrogram
row_dend <- as.dendrogram(row_clust)

# Color 3 major clusters
row_dend_colored <- color_branches(row_dend, k = 5)

# Optional: choose specific colors
# row_dend_colored <- color_branches(row_dend, k = 3,
#                                    col = c("red", "blue", "darkgreen"))

# Plot heatmap with colored row dendrogram
P1 <- heatmap.2(
  sig_gene_matrix_z,
  
  Rowv = row_dend_colored,
  Colv = as.dendrogram(col_clust),
  
  col = rev(brewer.pal(9, "RdBu")),
  trace = "none",
  
  srtCol = 45,
  adjCol = c(1, 1),
  offsetCol = 0.5,
  cexCol = 0.7,
  
  key = TRUE,
  density.info = "none",
  key.title = "",
  key.xlab = ""
)


#--------------------------------------------------------


all_genes <- rownames(sig_gene_matrix_z)


gene_sign <- ifelse(all_genes %in% sig_pos_cor, "pos",
                    ifelse(all_genes %in% sig_neg_cor, "neg", NA))


gene_clusters <- ifelse(gene_sign == "pos", 1,
                        ifelse(gene_sign == "neg", 2, NA))



library(dendextend)

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

#----------------------------------------------------------------
#------------------------------------------------------------------




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


# Open a new blank plot for legend in new plot window

# plot.new()
# 
# # Add legend only
# legend(
#   "center",
#   legend = names(group_colors),
#   fill = group_colors,
#   border = NA,
#   bty = "n",
#   cex = 1.2
# )



