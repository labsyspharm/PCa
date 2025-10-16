# Required libraries (from the file)
library(ggplot2)
library(vioplot)  # Though not used in the ggplot, included as in file

# Assume ldata is loaded as in the file
#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki.RData")

# Load sheet and define high/low groups (exact from file, with removals)
sheet <- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")
high <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "High")
high <- high[-3]  # Remove as in file
low <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "Low")
low <- low[-c(1,4)]  # Remove as in file

# Calculate TLS counts per specimen for high group (needed for mean.vector allocation)
h.ntls <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  h.ntls[c] <- length(unique(l$tls_id))
  c <- c + 1
}
h.ntls <- h.ntls - 1  # Subtract non-TLS (ID 0)

# Calculate TLS counts per specimen for low group
l.ntls <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  l.ntls[c] <- length(unique(l$tls_id))
  c <- c + 1
}
l.ntls <- l.ntls - 1  # Subtract non-TLS (ID 0)

# Calculate average TLS size per specimen for high group
h.mtls <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  mean.vector <- numeric(h.ntls[c])
  cc <- 1
  for (j in unique(l$tls_id)[-which(unique(l$tls_id) == 0)]) {
    mean.vector[cc] <- sum(l$tls_id == j)
    cc <- cc + 1
  }
  h.mtls[c] <- mean(mean.vector)
  c <- c + 1
}

# Calculate average TLS size per specimen for low group
l.mtls <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  mean.vector <- numeric(l.ntls[c])
  cc <- 1
  for (j in unique(l$tls_id)[-which(unique(l$tls_id) == 0)]) {
    mean.vector[cc] <- sum(l$tls_id == j)
    cc <- cc + 1
  }
  l.mtls[c] <- mean(mean.vector)
  c <- c + 1
}

# Create data frame for average sizes
gleason <- rep(c("HGG", "LGG"), times = c(length(h.mtls), length(l.mtls)))
value <- c(h.mtls, l.mtls)
data <- data.frame(Gleason = gleason, Value = value)

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

# Perform appropriate test: t-test if both groups are normal, else Mann-Whitney U (two-sided)
if (ks_hgg$p.value > 0.05 && ks_lgg$p.value > 0.05) {
  test_result <- t.test(Value ~ Gleason, data = data_cleaned)
} else {
  test_result <- wilcox.test(Value ~ Gleason, data = data_cleaned)
}

# Define y-position for significance bar (near max, adjusted to match image ~7000 for max ~8000)
bar_y <- max(data_cleaned$Value) * 0.9

# Create violin plot without outliers
p <- ggplot(data_cleaned, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2, color = "black") +  # Black jitter points
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +  # Semi-transparent boxes
  labs(title = "BICs size per specimen", y = "Average number of cells in BIC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +  # Remove the legend
  ylim(0, max(data_cleaned$Value) * 1.1) +  # Slight extension for annotations
  # Add significance bar (dashed, black to match style)
  geom_segment(x = 1.1, xend = 1.9, y = bar_y, yend = bar_y, linetype = "dashed", color = "black") +
  # Add p-value text (above bar, black since p > 0.05)
  annotate("text", x = 1.5, y = bar_y, label = paste("P =", formatC(test_result$p.value, digits = 4, format = "f")),
           hjust = 0.5, vjust = -1, color = "black") +
  # Specify colors for groups
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77"))

# Display the plot
print(p)