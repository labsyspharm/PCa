rm(list = setdiff(ls(), c("ldata", "df")))

F_tls_w <- function(patientID, tlsID, W, draw = FALSE) {
  d  <- ldata[[patientID]]
  df <- subset(d, tls_id == tlsID)
  
  ## centre of TLS window
  centre <- c(mean(range(df$x)), mean(range(df$y)))
  
  ## window filter
  keep  <- with(d,
                x > centre[1] - W & x < centre[1] + W &
                  y > centre[2] - W & y < centre[2] + W)
  d_ret <- d[keep, ]
  
  if (draw) {
    
    ## ------------- background: all cells -------------
    plot(d_ret$x, d_ret$y,
         xlab = "X-coordinate (µm)",
         ylab = "Y-coordinate (µm)",
         main = paste(tlsID, patientID, sep = "@"),
         pch  = 16,                                 # filled circle
         cex  = 0.8,                                # a bit larger than "."
         col  = adjustcolor("grey60", alpha.f = 0.5))
    
    ## ------------- B- and T-cells with black outline -------------
    idx_B <- bitwAnd(d_ret$pflag, 64) == 64
    idx_T <- d_ret$coarse_phen_vec == "T cells"
    
    points(d_ret$x[idx_B], d_ret$y[idx_B],
           pch = 21, bg = adjustcolor("red",   0.9),
           col = "black", cex = 1.3, lwd = 0.4)
    
    points(d_ret$x[idx_T], d_ret$y[idx_T],
           pch = 21, bg = adjustcolor("green", 0.9),
           col = "black", cex = 1.3, lwd = 0.4)
  }
  
  invisible(d_ret)
}

radius_half_integral <- function(values, x_vals = NULL) {
  # Only use positive values
  values_pos <- values[na.omit(values) > 0]
  
  # Create x-values if not provided (assume unit spacing)
  if (is.null(x_vals)) {
    x_vals <- seq_along(values_pos)
  } else {
    x_vals <- x_vals[values > 0]
  }
  
  # Compute cumulative sum
  cumulative_sum <- cumsum(values_pos)
  total_sum <- cumulative_sum[length(cumulative_sum)]
  half_sum <- total_sum / 2
  
  # Find the interval where the half point lies
  idx <- which(cumulative_sum >= half_sum)[1]
  
  if (idx == 1 || is.na(idx)) {
    return(x_vals[1])
  }
  
  # Interpolate within the interval [idx-1, idx]
  x0 <- x_vals[idx - 1]
  x1 <- x_vals[idx]
  y0 <- cumulative_sum[idx - 1]
  y1 <- cumulative_sum[idx]
  
  x_half <- x0 + (half_sum - y0) * (x1 - x0) / (y1 - y0)
  return(x_half)
}


F_BT<- function(patientID, tlsID, plot = FALSE){
  ret <- numeric(6)
  #first the IC of the B 
  W <- 500 # 500 micron to each side
  d <- ldata[[patientID]]
  d$x <- 2*d$x
  d$y <- 2*d$y
  df <- subset(d, d$tls_id == tlsID)
  c <- c((min(df$x) + max(df$x)) / 2, (min(df$y) + max(df$y)) / 2) #tls center 
  dfff<- subset(d, d$x < c[1] + W & d$x > c[1] - W & d$y < c[2] + W & d$y > c[2] - W)
  
  ##for the B cells 
  dff<- subset(dfff, bitAnd(dfff$pflag, 64) == 64)
  ppp <- as.ppp(dff[, c("x", "y")], W = owin(c(c[1] - W, c[1] + W), c(c[2] - W, c[2] + W))) 
  r_vals <- seq(0, 500, 10)
  L <- Lest(ppp, r = r_vals[-length(r_vals)])
  index.holder <- mean(L$border - L$theo, na.rm = TRUE)
  differences <- L$border - L$theo
  ret[1]<- index.holder ## the clustering Intensity 
  ret[2]<- 2*radius_half_integral(differences, r_vals) ## the Estimated clustering radius 
  
  differencesB <<- differences
  
  ###for the T cells
  dff<- subset(dfff, bitAnd(dfff$pflag, 2048) == 2048)
  ppp <- as.ppp(dff[, c("x", "y")], W = owin(c(c[1] - W, c[1] + W), c(c[2] - W, c[2] + W))) 
  r_vals <- seq(0, 500, 10)
  L <- Lest(ppp, r = r_vals[-length(r_vals)])
  index.holder <- mean(L$border - L$theo, na.rm = TRUE)
  differences <- L$border - L$theo
  ret[3]<- index.holder ## the clustering Intensity 
  ret[4]<- 2*radius_half_integral(differences, r_vals[-length(r_vals)]) ## the Estimated clustering radius 
  
  differencesT <<- differences
  
  ##for the B and T interaction 
  phen.data <- subset(dfff, bitAnd(dfff$pflag, 2048) == 2048 | bitAnd(dfff$pflag, 64) == 64)
  if (sum(phen.data$coarse_phen_vec == "T cells") > 10 & sum(phen.data$coarse_phen_vec == "B cells") > 10) {
    # Assign type based on bitwise flags
    type <- ifelse(bitAnd(phen.data$pflag, 2048) == 2048, "T cell", "B cell")
    # Build the data frame
    cell_data <- data.frame(
      x = phen.data$x,
      y = phen.data$y,
      type = type  # newly constructed character vector
    )
    
    win <- owin(c(c[1] - W, c[1] + W), c(c[2] - W, c[2] + W))  # Adjust to actual tissue dimensions
    ppp_cells <- ppp(
      x = cell_data$x,
      y = cell_data$y,
      window = win,
      marks = factor(cell_data$type)  # Marks define cell types
    )
    L <- Lcross(ppp_cells, i = "B cell", j = "T cell", r = r_vals[-length(r_vals)])
    index.holder<- mean(L$border - L$theo, na.rm = TRUE)
    differences <- L$border - L$theo
    ret[5]<- index.holder ## the clustering Intensity 
    ret[6]<- 2*radius_half_integral(differences, r_vals) ## the Estimated clustering radius 
    
    differencesBT <<- differences 
  }
  if (plot){
    par(mfrow = c(2, 2))
    
    ##1
    F_tls_w(patientID, tlsID, W/2, TRUE)
    
    
    #2
    
    example_data<- differencesB
    example_data1<-example_data[example_data > 0]
    
    if(length(example_data1) > 0){
      total_sum <- sum(example_data1, na.rm = TRUE)
      cumulative_sum <- cumsum(example_data1)
      half_sum <- total_sum / 2
      closest_index <- which.min(abs(cumulative_sum - half_sum))
      colors <- ifelse(example_data < 0, "darkgrey", "red")
      colors[1:(closest_index*2)] <- "darkred"
      colors[(closest_index*2)] <- "navy"
      
      
      barplot(
        example_data,
        ylim = c(-200, 200),
        col = colors,
        names.arg = seq(10,500, by=10),
        xlab = "Micron",
        ylab = "Clustering intensity",
        main = "B clustering", las = 2, 
        border = "black" # Add borders to bars
      )
      # Add horizontal grid lines
      #abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
      abline(h = mean(example_data[which(example_data > 0)]), col = "darkred", lty = 1, lwd = 1.5)
      text(x= 5, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="darkred", font= 2)# Clustering intensity
      text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*10 , sep = ": "), cex = 1.2, col="darkred", font= 3, srt = 90)# 50% effective clustering radius
      text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*10*2, sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
      
      # Add a legend
      legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "darkred"), border = "black")
    }
    
    else {# Create points forming an "X" shape
      n <- 100
      x1 <- seq(0, 1, length.out = n)
      y1 <- x1
      
      x2 <- seq(0, 1, length.out = n)
      y2 <- rev(x2)
      
      # Combine both diagonal lines
      x_vals <- c(x1, x2)
      y_vals <- c(y1, y2)
      
      # Plot as black points shaped like "X"
      plot(x_vals, y_vals,
           pch = 16, col = "black",
           axes = FALSE, xlab = "", ylab = "",
           frame.plot = FALSE)}
    
    
    #3
    
    example_data<- differencesT
    example_data1<-example_data[example_data > 0]
    
    if(length(na.omit(example_data1)) > 0){
      total_sum <- sum(example_data1, na.rm = TRUE)
      cumulative_sum <- cumsum(example_data1)
      half_sum <- total_sum / 2
      closest_index <- which.min(abs(cumulative_sum - half_sum))
      colors <- ifelse(example_data < 0, "darkgrey", "lightgreen")
      colors[1:(closest_index*2)] <- "darkgreen"
      colors[(closest_index*2)] <- "navy"
      
      
      barplot(
        example_data,
        ylim = c(-200, 200),
        col = colors,
        names.arg = seq(10,500, by=10),
        xlab = "Micron",
        ylab = "Clustering intensity",
        main = "T clustering", las = 2, 
        border = "black" # Add borders to bars
      )
      # Add horizontal grid lines
      #abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
      abline(h = mean(example_data[which(example_data > 0)]), col = "darkgreen", lty = 1, lwd = 1.5)
      text(x= 5, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="darkgreen", font= 2)# Clustering intensity
      text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*10, sep = ": "), cex = 1.2, col="darkgreen", font= 3, srt = 90)# 50% effective clustering radius
      text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*10*2, sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
      
      # Add a legend
      legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "darkgreen"), border = "black")
    }
    
    else {# Create points forming an "X" shape
      n <- 100
      x1 <- seq(0, 1, length.out = n)
      y1 <- x1
      
      x2 <- seq(0, 1, length.out = n)
      y2 <- rev(x2)
      
      # Combine both diagonal lines
      x_vals <- c(x1, x2)
      y_vals <- c(y1, y2)
      
      # Plot as black points shaped like "X"
      plot(x_vals, y_vals,
           pch = 16, col = "black",
           axes = FALSE, xlab = "", ylab = "",
           frame.plot = FALSE)}
    
    
    ##4
    
    
    
    example_data<- differencesBT
    example_data1<-example_data[example_data > 0]
    
    if(length(example_data1) > 0){
      total_sum <- sum(example_data1, na.rm = TRUE)
      cumulative_sum <- cumsum(example_data1)
      half_sum <- total_sum / 2
      closest_index <- which.min(abs(cumulative_sum - half_sum))
      colors <- ifelse(example_data < 0, "darkgrey", "lightyellow1")
      colors[1:(closest_index*2)] <- "sandybrown"
      colors[(closest_index*2)] <- "navy"
      
      
      barplot(
        example_data,
        ylim = c(-200, 200),
        col = colors,
        names.arg = seq(10,500, by=10),
        xlab = "Micron",
        ylab = "Clustering intensity",
        main = "B and T interaction", las = 2, 
        border = "black" # Add borders to bars
      )
      # Add horizontal grid lines
      #abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
      abline(h = mean(example_data[which(example_data > 0)]), col = "sandybrown", lty = 1, lwd = 1.5)
      text(x= 5, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="sandybrown", font= 2)# Clustering intensity
      text(x= closest_index + closest_index*0.19, y = -50,las = 2,  labels = paste("R50", closest_index*10, sep = ": "), cex = 1.2, col="sandybrown", font= 3, srt = 90)# 50% effective clustering radius
      text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -50,las = 2,  labels = paste("ECR", closest_index*10*2, sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
      
      # Add a legend
      legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "sandybrown"), border = "black")
    }
    
    else {# Create points forming an "X" shape
      n <- 100
      x1 <- seq(0, 1, length.out = n)
      y1 <- x1
      
      x2 <- seq(0, 1, length.out = n)
      y2 <- rev(x2)
      
      # Combine both diagonal lines
      x_vals <- c(x1, x2)
      y_vals <- c(y1, y2)
      
      # Plot as black points shaped like "X"
      plot(x_vals, y_vals,
           pch = 16, col = "black",
           axes = FALSE, xlab = "", ylab = "",
           frame.plot = FALSE)}
    
    
  }
  return(ret)
} 

# Extended Fig 5e
F_BT("LSP12611",4, plot = TRUE)
