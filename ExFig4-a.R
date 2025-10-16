
library(FNN)
library(data.table) 
library(googlesheets4)
library(vioplot)
library(bitops)


markers<- c("CD3d")
mpflag<- c(2^11)
#marker equabalvazne tof the pflat


mat<- as.data.frame(matrix(0, length(markers), length(names(ldata))))
rownames(mat)<- markers
colnames(mat)<- names(ldata)


#
mat<- list()
for (i in c(1:15, 17:21, 24:29)){
  element_name <- paste(names(ldata)[i])
  Tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8)
  Strom<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8 & !ldata[[i]]$tcell_dbscan %in% c(-999, 0, -1))
  if (dim(Strom)[1]!=0){
    nStromaTLSs<- length(unique(Strom$tcell_dbscan))
    mat[[element_name]] <- as.data.frame(matrix(0, length(markers), nStromaTLSs))
    rownames(mat[[element_name]])<- markers
    colnames(mat[[element_name]])<- as.character(unique(Strom$tcell_dbscan))
    for (k in 1:nStromaTLSs){
      for (j in 1:length(markers)){
        strom <- subset(Strom, bitAnd(Strom$pflag, mpflag[j]) == mpflag[j] & Strom$tcell_dbscan == unique(Strom$tcell_dbscan)[k])
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



sD<- list()
for (i in 1:length(names(ldata))){
  if (length(unique(ldata[[i]]$tcell_dbscan)[which(!unique(ldata[[i]]$tcell_dbscan) %in% c(-999, 0, -1))]) != 0){
    sD[[i]] <- unique(ldata[[i]]$tcell_dbscan)[which(!unique(ldata[[i]]$tcell_dbscan) %in% c(-999, 0, -1))]
    c<- numeric(length(sD[[i]]))
    h<- sD[[i]]
    for (j in 1:length(c)){
      c[j] <- sum(ldata[[i]]$tcell_dbscan == h[j])
    }
    sD[[i]] <- c
  }
}

names(sD) <- names(ldata)


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

par(mar = c(2, 2, 2, 2))
plot(0, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", axes = FALSE, asp = 1)
title(main="TLS Proximity Chart", font.main=2, col.main="navy")

# Define the data - 26 equal sections
values <- rep(1, 26) # Equal values for equal sections

# Generate random numbers for each section to display inside the pie chart
random_numbers <- sample(1:100, 26, replace = TRUE)

# Frequencies to display outside the pie chart
frequencies <- sample(100:1000, 26, replace = TRUE)

# Define 26 colors
colors <- brewer.pal(n = min(26, 8), name = "Dark2")
colors <- colorRampPalette(colors)(26)

colors[1:length(high)] <- "#7570b3"  # Color for high
colors[(length(high) + 1):(length(high) + length(low))] <- "#1b9e77"  # Color for low
#colors_transparent <- adjustcolor(colors, alpha.f = 0.5)
#colors[(length(high) + length(low) + 1):26] <- "grey"  # Remaining colors as grey

# Calculate the starting and ending angles for each slice
angles_deg <- cumsum(c(0, values / sum(values) * 360))
angles_rad <- angles_deg * pi / 180

# Draw each pie slice as a polygon
for (i in 1:length(values)) {
  x <- c(0, 1.4*cos(angles_rad[i]), 1.4*cos(angles_rad[i + 1]), 0)
  y <- c(0, 1.4*sin(angles_rad[i]), 1.4*sin(angles_rad[i + 1]), 0)
  polygon(x, y, col ="white", border = TRUE, lty = 3, lwd = 0.6)
}




radii <- seq(1, 1.4, by = 0.1) # Define radii for the circles
for (r in radii) {
  symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1)
  if (r == 1.2){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 0.7)}
  if (r == 1.2){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1.6)}
  if (r == 1.3){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 0.7)}
  if (r == 1.4){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 2.2)}
}
radiiii <- seq(1.39, 1.401, by = 0.001) # Define radii for the circles
for (r in radiiii) {
  symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1.5)
}
#add ther median values for each group
r.high <- 1+ log(median_value.high)/25
r.low  <- 1+ log(median_value.low)/25
symbols(0, 0, circles = r.high, inches = FALSE, add = TRUE, fg = "#7570b3", lwd = 1.8, lty= 5)
symbols(0, 0, circles = r.low, inches = FALSE, add = TRUE, fg = "#1b9e77", lwd = 2, lty= 5)

# Loop through each section to add the random numbers and frequencies


##soreting the data
tumor.portion<- numeric(length(values))
names(tumor.portion)<- c(high,low)
for (i in 1:13){
  tumor.portion[i]<- sum(bitAnd(ldata[[which(names(ldata)==c(high,low)[i])]]$cflag , 8) == 8)/dim(ldata[[which(names(ldata)==c(high,low)[i])]])[1]
}
hi.portion<- names(tumor.portion)[order(tumor.portion[1:13])]
for (i in 14:26){
  tumor.portion[i]<- sum(bitAnd(ldata[[which(names(ldata)==c(high,low)[i])]]$cflag , 8) == 8)/dim(ldata[[which(names(ldata)==c(high,low)[i])]])[1]
}
lo.portion<- names(tumor.portion[14:26])[order(tumor.portion[14:26])]
high<- hi.portion
low<- lo.portion

for (i in 1:length(values)) {
  x <- c(0, cos(angles_rad[i]), cos(angles_rad[i + 1]), 0)
  y <- c(0, sin(angles_rad[i]), sin(angles_rad[i + 1]), 0)
  polygon(x, y, col = adjustcolor(colors[i], alpha.f = 0.6))
  #polygon(x, y, col = adjustcolor(colors[i], alpha.f = 0.05 + sqrt(sum(bitAnd(ldata[[which(names(ldata)==c(high,low)[i])]]$cflag , 8) == 8)/dim(ldata[[which(names(ldata)==c(high,low)[i])]])[1])), border = TRUE)
}
radiii <- seq(0.99, 1.001, by = 0.001) # Define radii for the circles
for (r in radiii) {
  symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1.5)
}
for (i in 1:length(values)) {
  # Convert angles to radians for text positioning
  angle_rad_inner <- (angles_deg[i] + angles_deg[i+1]) / 2 * pi / 180
  
  # Calculate sqrt of the sum of squares of sin and cos values
  
  # Coordinates for the random numbers inside the pie
  x_pos_inner1 <- 0.38 * cos(angle_rad_inner)
  y_pos_inner1 <- 0.38 * sin(angle_rad_inner)
  
  x_pos_inner2 <- 0.65 * cos(angle_rad_inner)
  y_pos_inner2 <- 0.65 * sin(angle_rad_inner)
  
  x_pos_inner3 <- 0.9 * cos(angle_rad_inner)
  y_pos_inner3 <- 0.9 * sin(angle_rad_inner)
  
  # Coordinates for the frequencies just outside the pie
  x_pos_outer <- 1.5 * cos(angle_rad_inner)
  y_pos_outer <- 1.5 * sin(angle_rad_inner)
  
  x_pos_outer1 <- 1.1 * cos(angle_rad_inner)
  y_pos_outer1 <- 1.1 * sin(angle_rad_inner)
  
  x_pos_outer2 <- 1.2 * cos(angle_rad_inner)
  y_pos_outer2 <- 1.2 * sin(angle_rad_inner)
  
  x_pos_outer3 <- 1.3 * cos(angle_rad_inner)
  y_pos_outer3 <- 1.3 * sin(angle_rad_inner)
  
  x_pos_outer4 <- 1.4 * cos(angle_rad_inner)
  y_pos_outer4 <- 1.4 * sin(angle_rad_inner)
  # Add the random numbers inside the pie
  #text(x_pos_outer, y_pos_outer, labels = random_numbers[i], cex = 0.7)
  jj<-i
  i<- c(which(names(ldata)%in%c(high, low)[jj]))
  
  #text(x_pos_inner2, y_pos_inner2, labels = length(unique(ldata[[i]]$tls_id))-1, cex = (0.65 + (length(unique(ldata[[i]]$tls_id))-1)/33), font=3 )
  
  # Add the frequencies outside the pie; adjust the alignment based on the angle
  if (angle_rad_inner > pi/2 && angle_rad_inner < 3*pi/2) {
    text(x_pos_inner1, y_pos_inner1, labels = names(ldata)[i], cex = 1.1, font= 2, srt = 180 * angle_rad_inner / pi - 180, adj = 1)
  } else {
    text(x_pos_inner1, y_pos_inner1, labels = names(ldata)[i], cex = 1.1, font= 2, srt = 180 * angle_rad_inner / pi, adj = 0)
  }
  if (names(ldata)[i] %in% names(mat)) {
    #text(x_pos_outer, y_pos_outer, labels = dim(mat[[which(names(mat)==names(ldata)[i])]])[2], font= 2, cex = (0.65 + (dim(mat[[which(names(mat)==names(ldata)[i])]])[2])/33))
    text(x_pos_inner3, y_pos_inner3, labels = (length(unique(ldata[[i]]$tcell_dbscan))-3) - (dim(mat[[which(names(mat)==names(ldata)[i])]])[2]), cex = 1.3)
    err<- runif(1, -0.1, 0.1)
    for (j in 1:length(as.numeric(ld[which(unlist(strsplit(names(ld), "@"))[c(TRUE, FALSE)] == names(ldata)[i])])/25)){
      err<- runif(1, -0.08, 0.08)
      x_pos_data <- (1+ (as.numeric(ld[which(unlist(strsplit(names(ld), "@"))[c(TRUE, FALSE)] == names(ldata)[i])])/25))[j] * cos(angle_rad_inner+err)
      y_pos_data <- (1+ (as.numeric(ld[which(unlist(strsplit(names(ld), "@"))[c(TRUE, FALSE)] == names(ldata)[i])])/25))[j] * sin(angle_rad_inner+err)
      opacity<- 1 - sqrt((x_pos_data)^2 + (y_pos_data)^2) +1
      points(x_pos_data, y_pos_data, pch=16, col=adjustcolor(colors[jj], alpha.f = 0.75), cex = 1.6)}
    #points(x_pos_data, y_pos_data, pch=20, col=adjustcolor(colors[jj], alpha.f= opacity^2.5))}
  }
  else{
    #text(x_pos_outer, y_pos_outer, 0, cex = 0.65, font= 2)
    text(x_pos_inner3, y_pos_inner3, labels = length(unique(ldata[[i]]$tls_id))-1, cex = 1.3)
  }
  if (i == 23) {
    text(x_pos_outer1, y_pos_outer1, paste(12, "um", " "), cex = 0.7, col="black", font= 2)
    text(x_pos_outer2, y_pos_outer2, paste(150, "um", " "), cex = 0.75, col="black", font= 2)
    text(x_pos_outer3, y_pos_outer3, paste(1.8, "mm", " "), cex = 0.75, col="black", font= 2)
    text(x_pos_outer4, y_pos_outer4, paste(22, "mm", " "), cex = 0.8, col="black", font= 2)
  }
  i<-jj
}


# Normality tests (Shapiro-Wilk)
sh_high <- shapiro.test(filtered_values_high)
sh_low  <- shapiro.test(filtered_values_low)

if (sh_high$p.value > 0.05 & sh_low$p.value > 0.05) {
  # Both are normal → t-test
  test_result <- t.test(filtered_values_high, filtered_values_low)
} else {
  # Non-normal → Mann-Whitney U test
  test_result <- wilcox.test(filtered_values_high, filtered_values_low)
}

print(test_result)



