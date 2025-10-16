library(FNN)
library(data.table) 
library(googlesheets4)
library(vioplot)
library(bitops)
library(ggplot2)
library(ggtext)
#
markers<-"CD20"
mpflag<- numeric(length(markers))
mat<- list()
for (i in c(1:15, 17:21, 24:29)){
  element_name <- paste(names(ldata)[i])
  Tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8)
  Strom<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8 & ldata[[i]]$tls_id != 0)
  if (dim(Strom)[1]!=0){
    nStromaTLSs<- length(unique(Strom$tls_id))
    mat[[element_name]] <- as.data.frame(matrix(0, length(markers), nStromaTLSs))
    rownames(mat[[element_name]])<- markers
    colnames(mat[[element_name]])<- as.character(unique(Strom$tls_id))
    for (k in 1:nStromaTLSs){
      for (j in 1:length(markers)){
        strom <- subset(Strom, bitAnd(Strom$pflag, mpflag[j]) == mpflag[j] & Strom$tls_id == unique(Strom$tls_id)[k])
        data_points <- cbind(Tumor$x, Tumor$y)  # Combine x and y into a matrix
        query_point <- matrix(c(strom$x, strom$y), ncol=2)
        knn_result <- get.knnx(data_points, query_point, k=1)
        nearest_neighbor_index = knn_result$nn.index # Extract the index of the nearest neighbor
        nearest_neighbor <- data_points[nearest_neighbor_index, ]  # Get the coordinates of the nearest neighbor
        diffs <- nearest_neighbor - cbind(strom$x, strom$y)  # Calculate the differences in x and y coordinates between corresponding points
        distances <- sqrt(rowSums(diffs^2)) # Calculate the Euclidean distance for each pair
        mat[[element_name]][j, k] <- mean(distances) # Print the pairwise distances
      }
      print(mat)
    }
    #print(mat)
  }
}

for (i in 1:length(names(mat))){
  for (j in 1:length(colnames(mat[[i]]))){
    colnames(mat[[i]])[j]<- paste(names(mat)[i], colnames(mat[[i]])[j], sep="@")
  }
}


matr<- as.data.frame(mat[[1]])
for (i in names(mat)[-1]){
  matr<- cbind(matr, mat[[i]])
}

mean_distances<- colMeans(matr, na.rm = TRUE)*2

ld<- log(mean_distances)
ld[which(ld== min(ld))]<- 0

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

# Extract names from the mean_distances vector
mean_distances_names <- names(mean_distances)

# Filter names that contain any of the low identifiers
filtered_names <- mean_distances_names[sapply(mean_distances_names, function(name) any(sapply(low, function(h) grepl(h, name))))]

# Subset the mean_distances vector to get only the values whose names match the filtered names
filtered_values <- mean_distances[filtered_names]
filtered_values_low<- filtered_values

# Calculate the mean of these filtered values
median_value.low <- median(filtered_values)

# Extract names from the mean_distances vector
mean_distances_names <- names(mean_distances)

# Filter names that contain any of the high identifiers
filtered_names <- mean_distances_names[sapply(mean_distances_names, function(name) any(sapply(high, function(h) grepl(h, name))))]

# Subset the mean_distances vector to get only the values whose names match the filtered names
filtered_values <- mean_distances[filtered_names]
filtered_values_high<-filtered_values

# Calculate the mean of these filtered values
median_value.high <- median(filtered_values)

data_high <- data.frame(values = filtered_values_high, group = "High")
data_low <- data.frame(values = filtered_values_low, group = "Low")

# Combine both datasets
data_combined <- rbind(data_high, data_low)

# Compute summary statistics
summary_high <- summary(filtered_values_high)
summary_low <- summary(filtered_values_low)

# Define colors for the groups
high_color <- "#7570b3"  # Purple
low_color <- "#1b9e77"   # Green

# Create a formatted text annotation with bold and color
summary_text <- paste0(
  "**High:**\n\n",
  "<span style='color:", high_color, "'><b>1st Qu:</b> ", round(summary_high["1st Qu."], 2), "</span>\n",
  "<span style='color:", high_color, "'><b>3rd Qu:</b> ", round(summary_high["3rd Qu."], 2), "</span>\n\n","<span style='color:", high_color, "'><b>Median:</b> ", round(summary_high["Median"], 2), "</span>\n",
  "<span style='color:", high_color, "'><b>Mean:</b> ", round(summary_high["Mean"], 2), "</span>\n\n",
  
  "**Low:**\n\n",
  "<span style='color:", low_color, "'><b>1st Qu:</b> ", round(summary_low["1st Qu."], 2), "</span>\n",
  "<span style='color:", low_color, "'><b>3rd Qu:</b> ", round(summary_low["3rd Qu."], 2), "</span>\n\n","<span style='color:", low_color, "'><b>Median:</b> ", round(summary_low["Median"], 2), "</span>\n",
  "<span style='color:", low_color, "'><b>Mean:</b> ", round(summary_low["Mean"], 2), "</span>\n"
  
)

# Create the plot with histogram and density curves
ggplot(data_combined, aes(x = values, fill = group, color = group)) +
  geom_histogram(aes(y = ..density..), alpha = 0.3, position = "identity", bins = 12) +  # Transparent histograms
  geom_density(size = 1.5, adjust = 1.2, alpha= 0.01) +  # Density curves
  scale_fill_manual(values = c(high_color, low_color)) +  # Fill colors for histogram
  scale_color_manual(values = c(high_color, low_color)) +  # Line colors for density curves
  labs(title = "**Stromal ICs Distances to Tumor**",
       x = "Value",
       y = "Density") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_markdown(size = 14, face = "bold")  # Make the title bold
  ) +
  annotate("richtext", x = max(data_combined$values) * 0.6, y = 0.0003, 
           label = summary_text, hjust = 0, size = 4, fill = NA, label.color = NA)



