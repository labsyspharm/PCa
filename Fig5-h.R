


low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
          "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")


# Function to compute closest_index for one sample
compute_closest_index <- function(data_vector) {
  # Keep only positive values
  data_vector_pos <- data_vector[data_vector > 0]
  
  # Return NA if all values are non-positive
  if (length(data_vector_pos) == 0) {
    return(NA)
  }
  
  # Compute required sums
  total_sum <- sum(data_vector_pos)
  cumulative_sum <- cumsum(data_vector_pos)
  half_sum <- total_sum / 2
  
  # Find the closest index
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  
  return(closest_index)
}

# Initialize vectors to store results
ECRL <- numeric(length(low))  # One value per sample in low
ECRH <- numeric(length(high))  # One value per sample in high

# Compute closest_index for LOW group
for (i in seq_along(low)) {
  sample_name <- low[i]
  sample_data <- ldata[[sample_name]]$BT.CR500[1:50]  # Extract first 50 elements
  ECRL[i] <- compute_closest_index(sample_data)*30*2 +compute_closest_index(sample_data) # Compute closest index
}

# Compute closest_index for HIGH group
for (i in seq_along(high)) {
  sample_name <- high[i]
  sample_data <- ldata[[sample_name]]$BT.CR500[1:50]  # Extract first 50 elements
  ECRH[i] <- compute_closest_index(sample_data)*30*2 +compute_closest_index(sample_data) # Compute closest index
}

# Print results
print(ECRL)  # Vector of closest_index values for LOW group
print(ECRH)  # Vector of closest_index values for HIGH group

# Load required library
library(ggplot2)

# Define colors
colors <- c("High" = "#7570b3", "Low" = "#1b9e77")

# Create a data frame for boxplot
ECR_data <- data.frame(
  Group = rep(c("High", "Low"), c(length(ECRH), length(ECRL))),  # Group labels
  Value = c(ECRH, ECRL)  # Corresponding closest_index values
)

# Perform Wilcoxon test (non-parametric equivalent of t-test)
wilcox_test <- wilcox.test(ECRH, ECRL)

# Generate boxplot
ggplot(ECR_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Jitter points for visibility
  scale_fill_manual(values = colors) +  # Assign colors to groups
  labs(
    title = "Comparison of Closest Index (ECRH vs ECRL)",
    x = "Group",
    y = "Closest Index"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  # Annotate p-value
  annotate("text", x = 1.5, y = max(ECR_data$Value) + 2, 
           label = paste0("p-value = ", signif(wilcox_test$p.value, 3)), 
           color = ifelse(wilcox_test$p.value < 0.05, "red", "black"), hjust = 0.5)

library(car)

# Example with two numeric vectors
df_ltest <- data.frame(
  value  = c(ECRH, ECRL),
  group  = rep(c("ECRH", "ECRL"), c(length(ECRH), length(ECRL)))
)

leveneTest(value ~ group, data = df_ltest)


all <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647", "LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
         "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")


# Compute mean BT.CR500 for LOW group (50 elements)
example_data_all <- rowMeans(sapply(all, function(sample) ldata[[sample]]$BT.CR500[1:50]), na.rm = TRUE)

# Print results
print(example_data_all)  # Mean of low group (50 elements)




example_data<- example_data_all
example_data1<-example_data[example_data > 0]

total_sum <- sum(example_data1)
cumulative_sum <- cumsum(example_data1)
half_sum <- total_sum / 2
closest_index <- which.min(abs(cumulative_sum - half_sum))
colors <- ifelse(example_data < 0, "darkgrey", "lightyellow1")
colors[1:(closest_index*2)] <- "sandybrown"
colors[(closest_index*2)] <- "navy"

# Define custom x-axis labels
x_labels <- seq(30, 1500, by=30)

# Create the barplot
barplot(
  example_data,
  ylim = c(-100, 250),
  col = colors,
  names.arg = x_labels,
  xlab = "Micron",
  ylab = "Clustering strenght",
  main = "B and T cells interaction C&R1mm LGG", las = 2, 
  border = "black" # Add borders to bars
)

# Add horizontal grid lines
#abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
abline(h = mean(example_data[which(example_data > 0)]), col = "sandybrown", lty = 1, lwd = 1.5)
text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="sandybrown", font= 2)# Clustering intensity
text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*31, sep = ": "), cex = 1.2, col="sandybrown", font= 3, srt = 90)# 50% effective clustering radius
text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*30*2 + closest_index, sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius

# Add a legend
legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "sandybrown"), border = "black")

