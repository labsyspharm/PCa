# Required libraries (from the file)
library(ggplot2)
library(bitops)  # For bitAnd, though not directly used here

# Assume ldata is loaded as in the file
#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

# Load sheet and define high/low groups (exact from provided, with removals)
sheet <- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")
high <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "High")
high <- high[-c(6)]
low <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "Low")
low <- low[-c(10,13,14)]

# Calculate TLS counts per specimen for high group (needed for size vector allocation)
h.ntls <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[i]]
  h.ntls[c] <- length(unique(l$tls_id))
  c <- c + 1
}
h.ntls <- h.ntls - 1  # Subtract non-TLS (ID 0)

# Calculate TLS counts per specimen for low group
l.ntls <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[i]]
  l.ntls[c] <- length(unique(l$tls_id))
  c <- c + 1
}
l.ntls <- l.ntls - 1  # Subtract non-TLS (ID 0)

# Calculate individual TLS sizes (cells per BIC) for high and low groups
h.size <- numeric(sum(h.ntls))
l.size <- numeric(sum(l.ntls))

cc <- 1
for (i in low) {
  l <- ldata[[i]]
  for (j in unique(l$tls_id)[-which(unique(l$tls_id) == 0)]) {
    l.size[cc] <- sum(l$tls_id == j)
    cc <- cc + 1
  }
}

cc <- 1
for (i in high) {
  l <- ldata[[i]]
  for (j in unique(l$tls_id)[-which(unique(l$tls_id) == 0)]) {
    h.size[cc] <- sum(l$tls_id == j)
    cc <- cc + 1
  }
}
# Required libraries (from the file)
library(ggplot2)
library(bitops)  # For bitAnd, though not used here

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
  test_result <- wilcox.test(Value ~ Gleason, data = data_cleaned, alternative = "greater")
}

# Compute mean ratio MH/ML
mean_HGG <- median(data_cleaned$Value[data_cleaned$Gleason == "HGG"])
mean_LGG <- median(data_cleaned$Value[data_cleaned$Gleason == "LGG"])
mean_ratio <- mean_HGG / mean_LGG

# Define y-position for significance bar (fixed as in file for this scale)
bar_y <- 3500

# Create violin plot without outliers
p <- ggplot(data_cleaned, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2, color = "black") +  # Black jitter points to match image
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Suppress outliers in boxplot
  labs(title = "Cells per BIC", y = "Number of cells per BIC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +  # Remove the legend
  ylim(0, max(data_cleaned$Value) * 1.1) +  # Slight extension for annotations
  # Add significance bar
  geom_segment(x = 1.1, xend = 1.9, y = bar_y, yend = bar_y, linetype = "dashed", color = "black") +
  # Add p-value text (above bar)
  annotate("text", x = 1.5, y = bar_y, label = paste("P =", formatC(test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1, color = "black") +
  # Add mean ratio text (below bar, with superscript/subscript)
  annotate("text", x = 1.5, y = bar_y, label = bquote(M^H / M[L] == .(formatC(mean_ratio, digits = 2, format = "f"))),
           hjust = 0.5, vjust = 2, color = "black") +
  # Specify colors for groups
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77"))

# Display the plot
print(p)