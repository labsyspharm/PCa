# Required libraries (from the file)
library(ggplot2)
library(bitops)  # For bitAnd, though not used here


low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
          "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
# Calculate TIC counts per specimen for high group (needed for size vector allocation)
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

# Calculate individual TIC sizes (cells per TIC) for high and low groups
h.size <- numeric(sum(h.ntic))
l.size <- numeric(sum(l.ntic))

cc <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  valid_clusters <- unique(l$tcell_dbscan[l$tcell_dbscan != -999 & l$tcell_dbscan != -1 & l$tcell_dbscan != 0])
  for (j in valid_clusters) {
    l.size[cc] <- sum(l$tcell_dbscan == j)
    cc <- cc + 1
  }
}

cc <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  valid_clusters <- unique(l$tcell_dbscan[l$tcell_dbscan != -999 & l$tcell_dbscan != -1 & l$tcell_dbscan != 0])
  for (j in valid_clusters) {
    h.size[cc] <- sum(l$tcell_dbscan == j)
    cc <- cc + 1
  }
}


# Create data frame for sizes
gleason <- rep(c("HGG", "LGG"), times = c(length(h.size), length(l.size)))
value <- c(h.size, l.size)
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

# Define y-position for significance bar (fixed as in file for this scale)
bar_y <- 800

# Create violin plot without outliers
p <- ggplot(data_cleaned, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2, color = "black") +  # Black jitter points to match image
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Suppress outliers in boxplot
  labs(title = "Cells per TIC", y = "Number of cells per TIC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +  # Remove the legend
  ylim(0, max(data_cleaned$Value) * 1.1) +  # Slight extension for annotations
  # Add significance bar
  geom_segment(x = 1.1, xend = 1.9, y = bar_y, yend = bar_y, linetype = "dashed", color = "black") +
  # Add p-value text (above bar)
  annotate("text", x = 1.5, y = bar_y, label = paste("P =", formatC(test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1, color = "black") +
  # Specify colors for groups
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77"))

# Display the plot
print(p)




# Required libraries (from the file)
library(ggplot2)
library(bitops)  # For bitAnd, though not used here


low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
          "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")



# Calculate average TIC size per specimen for high group
h.mtic <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  valid_clusters <- unique(l$tcell_dbscan[l$tcell_dbscan != -999 & l$tcell_dbscan != -1 & l$tcell_dbscan != 0])
  if (length(valid_clusters) > 0) {
    size.vector <- numeric(length(valid_clusters))
    cc <- 1
    for (j in valid_clusters) {
      size.vector[cc] <- sum(l$tcell_dbscan == j)
      cc <- cc + 1
    }
    h.mtic[c] <- mean(size.vector)
  } else {
    h.mtic[c] <- NA  # No TICs in specimen
  }
  c <- c + 1
}

# Calculate average TIC size per specimen for low group
l.mtic <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  valid_clusters <- unique(l$tcell_dbscan[l$tcell_dbscan != -999 & l$tcell_dbscan != -1 & l$tcell_dbscan != 0])
  if (length(valid_clusters) > 0) {
    size.vector <- numeric(length(valid_clusters))
    cc <- 1
    for (j in valid_clusters) {
      size.vector[cc] <- sum(l$tcell_dbscan == j)
      cc <- cc + 1
    }
    l.mtic[c] <- mean(size.vector)
  } else {
    l.mtic[c] <- NA  # No TICs in specimen
  }
  c <- c + 1
}


# Create data frame for average sizes
gleason <- rep(c("HGG", "LGG"), times = c(length(h.mtic), length(l.mtic)))
value <- c(h.mtic, l.mtic)
data <- data.frame(Gleason = gleason, Value = value)

# Remove NA rows (specimens with no TICs)
data <- na.omit(data)

# Remove outliers using the IQR method (exact from file)
remove_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.1)
  Q3 <- quantile(df[[column]], 0.9)
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

  test_result <- wilcox.test(Value ~ Gleason, data = data_cleaned, alternative = "greater")


# Define y-position for significance bar (fixed as in file for this scale)
bar_y <- 1000
p <- ggplot(data_cleaned, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "Cells per TIC", y = "Average number of cells per TIC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77")) +
  # extend y-axis without clipping
  expand_limits(y = bar_y * 0.6) +
  geom_segment(x = 1.1, xend = 1.9, y = 600, yend = 600,
               linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 600,
           label = paste("P =", formatC(test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1, color = "black")

print(p)




