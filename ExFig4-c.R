
# Define the dataset
df <- data.frame(
  d = c(26, 0, 1, 3, 0, 13, 8, 1, 7, 9, 1, 0, 20, 1, 33, 2, 9, 6, 3, 4, 8, 2, 7, 2, 22, 15, 6, 11, 37) 
  + c(23, 0, 12, 1, 2, 6, 1, 3, 13, 7, 8, 5, 25, 6, 0, 13, 3, 29, 35, 22, 8, 8, 23, 5, 51, 33, 11, 42, 191),
  Group = c("High", "Other", "Low", "Low", "Other", "Low", "Low", "Low", "Low", "Low",
            "High", "Other", "High", "High", "High", "Low", "Low", "High", "Low", "High",
            "High", "Low", "Low", "Low", "High", "High", "High", "High", "High")
)

# Filter only High and Low groups
df <- subset(df, Group %in% c("High", "Low"))
# Subset to only High and Low groups

# Define custom colors
colors <- c("High" = "#7570b3", "Low" = "#1b9e77")

# Violin + boxplot + jitter, with p-value shown
library(dplyr)
library(ggplot2)
library(ggpubr)

# Compute outlier-free data for plotting points
df_no_outliers <- df %>%
  group_by(Group) %>%
  mutate(
    Q1 = quantile(d, 0.25),
    Q3 = quantile(d, 0.75),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR
  ) %>%
  filter(d >= lower & d <= upper) %>%
  ungroup()

# Plot
ggplot(df_no_outliers, aes(x = Group, y = d, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    data = df_no_outliers,   # only plot non-outlier points
    width = 0.15, size = 1, alpha = 0.6
  ) +
  scale_fill_manual(values = colors) +
  stat_compare_means(
    comparisons = list(c("Low", "High")),
    method = "wilcox.test",
    alternative = "greater",
    label = "p.format"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of TIC and BIC per specimen",
    y = "TIC + BIC", x = "Group"
  )
