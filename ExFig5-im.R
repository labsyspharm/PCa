


low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
          "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")


# Compute mean BT.CR500 for LOW group (50 elements)
example_data_low <- rowMeans(sapply(low, function(sample) ldata[[sample]]$BT.CR500[1:50]), na.rm = TRUE)

# Compute mean BT.CR500 for HIGH group (50 elements)
example_data_high <- rowMeans(sapply(high, function(sample) ldata[[sample]]$BT.CR500[1:50]), na.rm = TRUE)

# Print results
print(example_data_low)  # Mean of low group (50 elements)
print(example_data_high)  # Mean of high group (50 elements)



example_data<- example_data_low
example_data1<-example_data[example_data > 0]

total_sum <- sum(example_data1)
cumulative_sum <- cumsum(example_data1)
half_sum <- total_sum / 2
closest_index <- which.min(abs(cumulative_sum - half_sum))
colors <- ifelse(example_data < 0, "darkgrey", "lightgreen")
colors[1:(closest_index*2)] <- "#1b9e77"
#colors[(closest_index*2)] <- "darkgreen"

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
  main = "B and T cells interaction C&R 1.5 mm LGG", las = 2, 
  border = "black" # Add borders to bars
)

# Add horizontal grid lines
#abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
abline(h = mean(example_data[which(example_data > 0)]), col = "#1b9e77", lty = 1, lwd = 1.5)
text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="#1b9e77", font= 2)# Clustering intensity
#text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*20 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="darkred", font= 3, srt = 90)# 50% effective clustering radius
text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*30*2 + closest_index, sep = ": "), cex = 1.2, col="#1b9e77", font= 2, srt = 90)#estimated clustering radius

# Add a legend
legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "#1b9e77"), border = "black")










####m

example_data<- example_data_high
example_data1<-example_data[example_data > 0]

total_sum <- sum(example_data1)
cumulative_sum <- cumsum(example_data1)
half_sum <- total_sum / 2
closest_index <- which.min(abs(cumulative_sum - half_sum))
colors <- ifelse(example_data < 0, "darkgrey", "#f7d7ff")
colors[1:(closest_index*2)] <- "#7570b3"
#colors[(closest_index*2)] <- "darkgreen"

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
  main = "B and T cells interaction C&R1mm HGG" , las = 2, 
  border = "black" # Add borders to bars
)

# Add horizontal grid lines
#abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
abline(h = mean(example_data[which(example_data > 0)]), col = "#7570b3", lty = 1, lwd = 1.5)
text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="#7570b3", font= 2)# Clustering intensity
#text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*20 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="darkred", font= 3, srt = 90)# 50% effective clustering radius
text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*30*2 + closest_index, sep = ": "), cex = 1.2, col="#7570b3", font= 2, srt = 90)#estimated clustering radius

# Add a legend
legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "#7570b3"), border = "black")



