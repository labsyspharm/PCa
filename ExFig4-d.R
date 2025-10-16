
# Load required libraries
library(ggplot2)

# Define the dataset
df <- data.frame(
  nB = c(26, 0, 1, 3, 0, 13, 8, 1, 7, 9, 1, 0, 20, 1, 33, 2, 9, 6, 3, 4, 8, 2, 7, 2, 22, 15, 6, 11, 37),
  nT = c(23, 0, 12, 1, 2, 6, 1, 3, 13, 7, 8, 5, 25, 6, 0, 13, 3, 29, 35, 22, 8, 8, 23, 5, 51, 33, 11, 42, 191),
  Group = c("High", "Other", "Low", "Low", "Other", "Low", "Low", "Low", "Low", "Low",
            "High", "Other", "High", "High", "High", "Low", "Low", "High", "Low", "High",
            "High", "Low", "Low", "Low", "High", "High", "High", "High", "High")
)

# Filter only High and Low groups
df <- subset(df, Group %in% c("High", "Low"))

# Separate data into High and Low groups
df_high <- subset(df, Group == "High")
df_low <- subset(df, Group == "Low")

# Perform linear regression for High and Low groups
fit_high <- lm(nT ~ nB, data = df_high)
fit_low <- lm(nT ~ nB, data = df_low)

# Extract model p-values
p_value_high <- summary(fit_high)$fstatistic
p_value_high <- pf(p_value_high[1], p_value_high[2], p_value_high[3], lower.tail = FALSE)

p_value_low <- summary(fit_low)$fstatistic
p_value_low <- pf(p_value_low[1], p_value_low[2], p_value_low[3], lower.tail = FALSE)

# Define colors
colors <- c("High" = "#7570b3", "Low" = "#1b9e77")

# Generate the plot
ggplot(df, aes(x = nB, y = nT, color = Group)) +
  geom_point(size = 3) +  # Scatter plot with color based on High/Low category
  geom_smooth(data = df_high, aes(x = nB, y = nT), method = "lm", color = "#7570b3", se = TRUE) +  # Regression line for High
  geom_smooth(data = df_low, aes(x = nB, y = nT), method = "lm", color = "#1b9e77", se = TRUE) +  # Regression line for Low
  scale_color_manual(values = colors) +  # Assign colors to groups
  theme_minimal() +
  labs(
    x = "nB",
    y = "nT",
    title = "Regression of nT vs nB for High and Low Groups"
  ) +
  # Annotate p-values
  annotate("text", x = min(df$nB) + 2, y = max(df$nT) - 5,   
           label = paste0("p(Model, High) = ", signif(p_value_high, 3)), color = "#7570b3", hjust = 0) +
  annotate("text", x = min(df$nB) + 2, y = max(df$nT) - 10,   
           label = paste0("p(Model, Low) = ", signif(p_value_low, 3)), color = "#1b9e77", hjust = 0) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold title
  )

