#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

i<-27  
  #pdf(paste(names(ldata)[i], ".pdf", sep ="CR1mm"), width=10.5, height=10.5)
  par(mfrow = c(1,2), oma = c(0,0,0,0), mar = c(1,1,1,1))
  
  #imm.cruse(500, LSP = names(ldata)[i], index = "K-integral", phen = "Both", fig = TRUE, creep = 1)
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  example_data<- ldata[[i]]$B.CR500[1:50]
  example_data1<-example_data[example_data > 0]
  
  total_sum <- sum(example_data1)
  cumulative_sum <- cumsum(example_data1)
  half_sum <- total_sum / 2
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  colors <- ifelse(example_data < 0, "darkgrey", "red")
  colors[1:closest_index] <- "darkred"
  colors[(closest_index*2)] <- "navy"
  
  # Define custom x-axis labels
  x_labels <- seq(30, 1500, by=30)
  
  # Create the barplot
  barplot(
    example_data,
    ylim = c(-100, 250),
    col = colors,
    names.arg = x_labels,
    xlab = "Radius (in micron)",
    ylab = "Clustering index",
    main = paste("B cells C&R 1.5 mm", names(ldata)[i], sep=" for "), las = 2, 
    border = "black" # Add borders to bars
  )
  
  # Add horizontal grid lines
  abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
  abline(h = mean(example_data[which(example_data > 0)]), col = "sienna", lty = 1, lwd = 1.5)
  text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="sienna3", font= 2)# Clustering intensity
  text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*31, sep = ": "), cex = 1.2, col="darkred", font= 3, srt = 90)# 50% effective clustering radius
  text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*30*2+ closest_index, sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
  
  # Add a legend
  legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "red"), border = "black")
  
  #imm.cruse(500, LSP = names(ldata)[i], index = "K-integral", phen = "T cells", fig = TRUE, creep = 1)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  example_data<- ldata[[i]]$T.CR500[1:50]
  example_data1<-example_data[example_data > 0]
  total_sum <- sum(example_data1)
  cumulative_sum <- cumsum(example_data1)
  half_sum <- total_sum / 2
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  colors <- ifelse(example_data < 0, "darkgrey", "green")
  colors[1:closest_index] <- "darkgreen"
  colors[(closest_index*2)] <- "navy"
  # Define custom x-axis labels
  x_labels <- seq(30, 1500, by=30)
  
  # Create the barplot
  barplot(
    example_data,
    ylim = c(-100, 250),
    col = colors,
    names.arg = x_labels,
    xlab = "Radius (in micron)",
    ylab = "Clustering index",
    main = paste("T cells C&R 1.5 mm", names(ldata)[i], sep=" for "), las = 2, 
    border = "black" # Add borders to bars
  )
  
  # Add horizontal grid lines
  abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
  abline(h = mean(example_data[which(example_data > 0)]), col = "olivedrab", lty = 1, lwd = 1.5)
  text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="olivedrab4", font= 2)
  text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*31, sep = ": "), cex = 1.2, col="darkgreen", font= 3, srt = 90)
  text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*30*2 + closest_index, sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
  
  #text(, y_pos_outer1, paste(7, "um", " "), cex = 0.7, col="darkgrey", font= 2)
  # Add a legend
  legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "green"), border = "black")
  
 