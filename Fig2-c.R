# Required libraries (from the file)
library(ggplot2)
library(vioplot)  # Though not used in the ggplot, included as in file

# Assume ldata is loaded as in the file
#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki.RData")

# Load sheet and define high/low groups (exact from file, with removals)
sheet <- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")
high <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "High")
high <- high[-c(3)]  # Remove as in file
low <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "Low")
low <- low[-c(1,4)]  # Remove as in file

# Calculate TLS counts per specimen for high group
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

# Create data frame for counts
gleason <- rep(c("HGG", "LGG"), times = c(length(h.ntls), length(l.ntls)))
value <- c(h.ntls, l.ntls)
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

# Perform t-test on cleaned data
# One-sided Wilcoxon test
t_test_result <- wilcox.test(Value ~ Gleason, 
                             data = data_cleaned, 
                             alternative = "greater")
# Compute mean ratio MH/ML
mean_HGG <- median(data_cleaned$Value[data_cleaned$Gleason == "HGG"])
mean_LGG <- median(data_cleaned$Value[data_cleaned$Gleason == "LGG"])
mean_ratio <- mean_HGG / mean_LGG

# Define y-position for significance bar (near max, adjusted to match image ~45 for max 50)
bar_y <- max(data_cleaned$Value) * 0.9

# Create violin plot without outliers
p <- ggplot(data_cleaned, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 2, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +  # make boxes semi-transparent
  labs(title = "BICs per specimen", y = "Counts per specimen") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +
  ylim(0, max(data_cleaned$Value) * 1.1) +
  geom_segment(x = 1.1, xend = 1.9, y = bar_y, yend = bar_y,
               linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = bar_y,
           label = paste("P =", formatC(t_test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1,
           color = ifelse(t_test_result$p.value < 0.05, "red", "black")) +
  annotate("text", x = 1.5, y = bar_y,
           label = paste("M_H / M_L =", formatC(mean_ratio, digits = 2, format = "f")),
           hjust = 0.5, vjust = 2, color = "black") +
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77"))

# Display the plot
print(p)

