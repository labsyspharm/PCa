# Load required libraries
library(RANN)       # For KDTree
library(dbscan)     # For HDBSCAN
library(dplyr)      # For data wrangling

# 1) Filter T cells
for (i in 1:29){
  LSP<- ldata[[i]]
  index<-which(LSP$tls_id==0)
  LSP1<- subset(LSP, LSP$tls_id == 0)
  LSP<- LSP1
  tcell_adata <- subset(LSP, LSP$coarse_phen_vec == "T cells")
  
  # 2) Extract X/Y coordinates
  xt_t <- tcell_adata$x
  yt_t <- tcell_adata$y
  
  # 3) Create Nx2 matrix of points
  t_points <- cbind(xt_t, yt_t)
  
  # 4) Build a KDTree and compute nearest neighbor distances
  kd_tree <- nn2(t_points, k = 2)  # k = 2 to get the nearest neighbor (excluding itself)
  nn_dists <- kd_tree$nn.dists[, 2]  # Extract second column (first NN excluding self)
  
  # Compute summary statistics
  median_nn <- median(nn_dists)
  mean_nn <- mean(nn_dists)
  min_nn <- min(nn_dists)
  max_nn <- max(nn_dists)
  
  # Print results
  cat("Median NN distance among all T cells =", median_nn, "\n")
  cat("Mean NN distance among all T cells =", mean_nn, "\n")
  cat("Min NN distance =", min_nn, "\n")
  cat("Max NN distance =", max_nn, "\n")
  

  LSP$tcell_cluster_hdbscan <- -999  # Initialize cluster labels
  global_offset <- 0  # Offset to keep cluster IDs unique
  
  
  # If not enough cells, skip
  if (nrow(tcell_adata) < 2) next
  
  # 2) Extract coordinates
  coords_t <- cbind(tcell_adata$x, tcell_adata$y)
  
  # 3) Run HDBSCAN
  hdb_result <- dbscan(coords_t, eps = 50, minPts = 100)  # minPts ≈ min_samples in Python
  
  labels_local <- hdb_result$cluster  # Cluster labels (-1 = noise)
  
  # 4) Offset labels for uniqueness across datasets
  if (any(labels_local != 0)) {
    max_local_label <- max(labels_local)
  } else {
    max_local_label <- -1
  }
  
  offset_labels <- ifelse(labels_local == 0, -1, labels_local + global_offset)  # Convert noise (0) to -1
  
  # 5) Store results
  LSP$tcell_cluster_hdbscan[LSP$coarse_phen_vec == "T cells"] <- offset_labels
  
  # 6) Update offset
  if (max_local_label >= 0) {
    global_offset <- global_offset + max_local_label + 1
  }
  tcell_dbscan <- numeric(dim(ldata[[i]])[1])
  ldata[[i]]$tcell_dbscan <- tcell_dbscan
  ldata[[i]]$tcell_dbscan[index] <- LSP$tcell_cluster_hdbscan 
  # Load required library
  
}


# Required libraries (from the file)
library(ggplot2)
library(vioplot)  # Though not used in the ggplot, included as in file

#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

# Load sheet and define high/low groups (exact from file, with removals)
sheet <- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")
high <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "High")
high <- high[-3]  # Remove as in file
low <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "Low")
low <- low[-c(1,4)]  # Remove as in file

# Calculate TIC counts per specimen for high group (using tdbscan_clusters, excluding -999, -1, 0)
h.ntic <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  valid_clusters <- unique(l$tcell_dbscan[l$tcell_dbscan != -999 & l$tcell_dbscan != -1 & l$tcell_dbscan != 0])
  h.ntic[c] <- length(valid_clusters)
  c <- c + 1
}

# Calculate TIC counts per specimen for low group
l.ntic <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  valid_clusters <- unique(l$tcell_dbscan[l$tcell_dbscan != -999 & l$tcell_dbscan != -1 & l$tcell_dbscan != 0])
  l.ntic[c] <- length(valid_clusters)
  c <- c + 1
}


# Create data frame for counts
gleason <- rep(c("HGG", "LGG"), times = c(length(h.ntic), length(l.ntic)))
value <- c(h.ntic, l.ntic)
data <- data.frame(Gleason = gleason, Value = value)

# Remove outliers using the IQR method (exact from file)
remove_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.1)
  Q3 <- quantile(df[[column]], 0.78)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  df_filtered <- df[df[[column]] >= lower_bound & df[[column]] <= upper_bound, ]
  return(df_filtered)
}
data_cleaned <- remove_outliers(data, "Value")

# Check normality with Kolmogorov-Smirnov test for each group
hgg_values <- data_cleaned$Value[data_cleaned$Gleason == "HGG"]
lgg_values <- data_cleaned$Value[data_cleaned$Gleason == "LGG"]

ks_hgg <- ks.test(hgg_values, "pnorm", mean = mean(hgg_values), sd = sd(hgg_values))
ks_lgg <- ks.test(lgg_values, "pnorm", mean = mean(lgg_values), sd = sd(lgg_values))

# Perform appropriate test: t-test if both groups are normal, else Mann-Whitney U
if (ks_hgg$p.value > 0.05 && ks_lgg$p.value > 0.05) {
  test_result <- t.test(Value ~ Gleason, data = data_cleaned)
} else {
  test_result <- wilcox.test(Value ~ Gleason, data = data_cleaned)
}

# Define y-position for significance bar (near max, adjusted to match image ~60 for max 75)
bar_y <- max(data_cleaned$Value) * 0.8

# Create violin plot without outliers
p <- ggplot(data_cleaned, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2, color = "black") +  # Black jitter points (no segment coloring)
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Suppress outliers in boxplot
  labs(title = "TIC per specimen", y = "Counts per specimen") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +  # Remove the legend
  ylim(0, max(data_cleaned$Value) * 1.1) +  # Slight extension for annotations
  # Add significance bar
  geom_segment(x = 1.1, xend = 1.9, y = bar_y, yend = bar_y, linetype = "dashed", color = "black") +
  # Add p-value text (above bar)
  annotate("text", x = 1.5, y = bar_y, label = paste("P =", formatC(test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1, color = ifelse(test_result$p.value < 0.05, "red", "black")) +
  # Specify colors for groups
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77"))

# Display the plot
print(p)

