

#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")


library(fastICA)
library(ggplot2)
library(dplyr)
library(bitops)
library(fastICA)
library(readr)


markers<- c("Hoechst", "Ki67", "AMCAR", "HMWCK", "CD19", "SMA", "CD20", "CD11b", "CD68", "CD163", "CD4", "CD3d", "CD8a", "TCF1", "FOXP3", "PD1", "CD57", "CD11c", "GranzymeB", "CD15", "HLADR", "CD103", "CD31", "pTBK1", "HLAA", "CD24", "CD44", "CD206")
markers<- markers[length(markers):1]
mpflag<- numeric(length(markers))
for (i in 0:(length(markers)-1)){
  mpflag[i+1]<- 2^i
}
mpflag<- mpflag[length(mpflag):1]

mpflag<- mpflag[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]
markers<- markers[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]
immune.markers<- c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")




F_ki67<- function(patientID, tlsID){#assumes the latest ldata read into memory
  d<- ldata[[patientID]]
  return(dim(subset(d, d$tls_id==tlsID & d$KiGate == 1))[1]/dim(subset(d, d$tls_id==tlsID))[1])
} 
F_Bki67<- function(patientID, tlsID){#assumes the latest ldata read into memory
  d<- ldata[[patientID]]
  return(dim(subset(d, d$tls_id==tlsID & d$KiGate == 1 & d$coarse_phen_vec == "B cells"))[1]/dim(subset(d, d$tls_id==tlsID & d$coarse_phen_vec == "B cells"))[1])
} 
F_TTCF1<- function(patientID, tlsID){#assumes the latest ldata read into memory
  d<- ldata[[patientID]]
  return(dim(subset(d, d$tls_id==tlsID & d$TCF1gate == 1 & d$coarse_phen_vec == "T cells"))[1]/dim(subset(d, d$tls_id==tlsID & d$coarse_phen_vec == "T cells"))[1])
} 
F_T4TCF1<- function(patientID, tlsID){#assumes the latest ldata read into memory
  d<- ldata[[patientID]]
  return(dim(subset(d, d$tls_id==tlsID & d$TCF1gate == 1 & d$phen_vec == "CD4 T cells"))[1]/dim(subset(d, d$tls_id==tlsID & d$phen_vec == "CD4 T cells"))[1])
} 
F_T8TCF1<- function(patientID, tlsID){#assumes the latest ldata read into memory
  d<- ldata[[patientID]]
  return(dim(subset(d, d$tls_id==tlsID & d$TCF1gate == 1 & d$phen_vec == "CD8 T cells"))[1]/dim(subset(d, d$tls_id==tlsID & d$phen_vec == "CD8 T cells"))[1])
} 
F_CD31<- function(patientID, tlsID){
  d<- ldata[[patientID]]
  return(dim(subset(d, d$tls_id==tlsID & bitAnd(d$pflag, 4194304) == 4194304))[1]/dim(subset(d, d$tls_id==tlsID))[1])
}
F_BT<- function(patientID, tlsID, plot = FALSE){
  ret <- numeric(6)
  #first the IC of the B 
  W <- 1500 # 500 micron to each side
  d <- ldata[[patientID]]
  d$x <- 2*d$x
  d$y <- 2*d$y
  df <- subset(d, d$tls_id == tlsID)
  c <- c((min(df$x) + max(df$x)) / 2, (min(df$y) + max(df$y)) / 2) #tls center 
  dfff<- subset(d, d$x < c[1] + W & d$x > c[1] - W & d$y < c[2] + W & d$y > c[2] - W)
  
  ##for the B cells 
  dff<- subset(dfff, bitAnd(dfff$pflag, 64) == 64)
  ppp <- as.ppp(dff[, c("x", "y")], W = owin(c(c[1] - W, c[1] + W), c(c[2] - W, c[2] + W))) 
  r_vals <- seq(0, 1500, 30)
  L <- Lest(ppp, r = r_vals[-length(r_vals)])
  index.holder <- mean(L$border - L$theo, na.rm = TRUE)
  differences <- L$border - L$theo
  ret[1]<- index.holder ## the clustering Intensity 
  ret[2]<- 5*radius_half_integral(differences, r_vals) ## the Estimated clustering radius 
  
  differencesB <<- differences
  
  ###for the T cells
  dff<- subset(dfff, bitAnd(dfff$pflag, 2048) == 2048)
  ppp <- as.ppp(dff[, c("x", "y")], W = owin(c(c[1] - W, c[1] + W), c(c[2] - W, c[2] + W))) 
  r_vals <- seq(0, 1500, 30)
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
        ylim = c(-500, 1500),
        col = colors,
        names.arg = seq(10,1500, by=30),
        xlab = "Micron",
        ylab = "Clustering intensity",
        main = "B clustering", las = 2, 
        border = "black" # Add borders to bars
      )
      # Add horizontal grid lines
      #abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
      abline(h = mean(example_data[which(example_data > 0)]), col = "darkred", lty = 1, lwd = 1.5)
      text(x= 50, y = mean(example_data[which(example_data > 0)]) + 50, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="darkred", font= 2)# Clustering intensity
      text(x= closest_index + closest_index*0.19, y = -250,las = 2,  labels = paste("R50", closest_index*30 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="darkred", font= 3, srt = 90)# 50% effective clustering radius
      text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -250,las = 2,  labels = paste("ECR", closest_index*30*2 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
      
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
        ylim = c(-500, 1500),
        col = colors,
        names.arg = seq(10,1500, by=30),
        xlab = "Micron",
        ylab = "Clustering intensity",
        main = "T clustering", las = 2, 
        border = "black" # Add borders to bars
      )
      # Add horizontal grid lines
      #abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
      abline(h = mean(example_data[which(example_data > 0)]), col = "darkgreen", lty = 1, lwd = 1.5)
      text(x= 50, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="darkgreen", font= 2)# Clustering intensity
      text(x= closest_index + closest_index*0.19, y = -250,las = 2,  labels = paste("R50", closest_index*30 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="darkgreen", font= 3, srt = 90)# 50% effective clustering radius
      text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -250,las = 2,  labels = paste("ECR", closest_index*30*2 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
      
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
        ylim = c(-500, 1500),
        col = colors,
        names.arg = seq(10,1500, by=30),
        xlab = "Micron",
        ylab = "Clustering intensity",
        main = "B and T interaction", las = 2, 
        border = "black" # Add borders to bars
      )
      # Add horizontal grid lines
      #abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
      abline(h = mean(example_data[which(example_data > 0)]), col = "sandybrown", lty = 1, lwd = 1.5)
      text(x= 50, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("CI", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="sandybrown", font= 2)# Clustering intensity
      text(x= closest_index + closest_index*0.19, y = -250,las = 2,  labels = paste("R50", closest_index*30 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="sandybrown", font= 3, srt = 90)# 50% effective clustering radius
      text(x= ((closest_index *2) + (closest_index *2)*0.19) , y = -250,las = 2,  labels = paste("ECR", closest_index*30*2 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="navy", font= 2, srt = 90)#estimated clustering radius
      
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
F_icat <- function(patientID, tlsID){
  d<- ldata[[patientID]]
  df_xy<- subset(d, d$tls_id==tlsID)
  X  <- as.matrix(df_xy[, c("x","y")]*2)
  mu <- colMeans(X, na.rm = TRUE)         # mean of original coords
  ica <- fastICA(X, n.comp = 2)           # centers internally
  S <- ica$S                              # n x 2
  A <- ica$A                              # 2 x 2 (mixing)
  # Reconstruct centered data, then add mean to go back to original location
  Xhat_centered <- S %*% t(A)             # n x 2
  Xhat <- sweep(Xhat_centered, 2, mu, "+")# add back mean
  #par(mfrow = c(2,1))
  #plot(X[,1],   X[,2],   main = "Original (x,y)",  xlab = "x", ylab = "y")
  #plot(Xhat[,1], Xhat[,2], main = "Reconstructed from ICA", xlab = "x", ylab = "y")
  #if(!sum(diag(sqrt(cov(Xhat)))) == sqrt(diag((cov(Xhat)))[1]+diag((cov(Xhat)))[2] +2*sqrt(diag((cov(Xhat)))[1]*diag((cov(Xhat)))[2]))){
  #  stop("Check your input! ")
  #}
  traceSD<- sqrt(diag((cov(Xhat)))[1]+diag((cov(Xhat)))[2] +2*sqrt(diag((cov(Xhat)))[1]*diag((cov(Xhat)))[2]))
  ICAT <- 100*traceSD/dim(X)[1]
  return(ICAT)
}
F_tls_c <- function(patientID, tlsID){## tls size 
  d<- ldata[[patientID]]
  df<- subset(d, d$tls_id==tlsID)
  c <- dim(df)[1]
  return(c)
}

for (i in 1:length(names(ldata))){
  LSPid<- as.character(rep(names(ldata)[i], dim(ldata[[i]])[1]))
  ldata[[i]] <- cbind(ldata[[i]], LSPid)
}

patients_with_tls<- names(ldata)[-c(2,5,12)] ##patients without BICs

t<- data.frame()
for (i in patients_with_tls){
  tlsp<- subset(ldata[[i]], ldata[[i]]$tls_id != 0)
  tlsp<- tlsp[, c(1:36, 135)]
  #lsptlsid<- as.character(rep(i, dim(tlsp)[1]))
  #tlspp<- cbind(tlsp, lsptlsid)
  t<- rbind(t, tlsp)
  print(length(unique(tlsp$tls_id)))
}

tlsP<- character()
for (i in 1:dim(t)[1]){
  tlsP[i]<- paste(t$tls_id[i], t$LSPid[i], sep= "@")
}
t<- cbind(t, tlsP)

df <- data.frame(
  TLS = unique(t$tlsP),
  ICA = numeric(length(unique(t$tlsP))),
  Ki67 = numeric(length(unique(t$tlsP))),
  CD31 = numeric(length(unique(t$tlsP))),
  BKi67= numeric(length(unique(t$tlsP))),
  TTCF1= numeric(length(unique(t$tlsP))),
  T4TCF1= numeric(length(unique(t$tlsP))),
  T8TCF1= numeric(length(unique(t$tlsP))),
  count = numeric(length(unique(t$tlsP)))
)

for (i in 1:dim(df)[1]){
  a<- strsplit(df[i,1], "@")[[1]][2]
  b<- as.numeric(strsplit(df[i,1], "@")[[1]][1])
  df[i,2] <- F_icat(a,b)
  df[i,3] <- F_ki67(a,b)
  df[i,4] <- F_CD31(a,b)
  df[i,5] <- F_Bki67(a,b)
  df[i,6] <- F_TTCF1(a,b)
  df[i,7] <- F_T4TCF1(a,b)
  df[i,8] <- F_T8TCF1(a,b)
  df[i,9] <- F_tls_c(a,b)
}


## ---- cases: use the first four only ----
cases <- c("49@LSP12625", "21@LSP12657",
           "17@LSP12655", "8@LSP12653")
selected_cases <- cases[1:4]

## ---- set up 2x2 layout & book-keeping ----
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
ica_total <- 0

## ---- helper: safe index for marker -> mpflag ----
flag_of <- function(marker_name) {
  i <- which(markers == marker_name)
  if (length(i) == 0) stop(sprintf("Marker '%s' not found in 'markers'.", marker_name))
  mpflag[i]
}

## ---- prefetch bitmasks for speed/readability ----
flag_CD20  <- flag_of(immune.markers[1])        # as in your code
flag_CD4   <- flag_of("CD4")
flag_CD11c <- flag_of(immune.markers[8])
flag_CD8a  <- flag_of(immune.markers[5])

## ---- plotting loop over 4 cases ----
for (ID in selected_cases) {
  
  ## subset one TLS/case
  subt <- subset(t, t$tlsP == ID)
  if (nrow(subt) == 0) {
    plot.new(); title(main = paste("No data for", ID)); next
  }
  df <- data.frame(subt)
  
  ## categories via bitwise flags
  df$category <- ifelse(
    bitAnd(df$pflag, flag_CD20) == flag_CD20, "CD20",
    ifelse(
      bitAnd(df$pflag, flag_CD4) == flag_CD4, "CD4",
      ifelse(
        bitAnd(df$pflag, flag_CD11c) == flag_CD11c, "CD11c",
        ifelse(
          bitAnd(df$pflag, flag_CD8a) == flag_CD8a, "CD8a",
          "Others"
        )
      )
    )
  )
  
  ## ICA on (x, y)
  mean_x <- mean(df$x)
  mean_y <- mean(df$y)
  
  ica_result <- fastICA(as.matrix(df[, c("x", "y")]), n.comp = 2)
  ica_trace  <- sum(apply(ica_result$S, 2, var))
  ica_total  <- ica_total + round(10000 * (ica_trace - 2), 0)
  a<- strsplit(ID, "@")[[1]][2]
  b<- as.numeric(strsplit(ID, "@")[[1]][1])
  ica_total <- F_icat(a,b)
  
  arrow_data <- data.frame(
    mean_x    = mean_x,
    mean_y    = mean_y,
    ica_end_x1 = mean_x + ica_result$A[1, 1],
    ica_end_y1 = mean_y + ica_result$A[2, 1],
    ica_end_x2 = mean_x + ica_result$A[1, 2],
    ica_end_y2 = mean_y + ica_result$A[2, 2]
  )
  
  ## aspect ratio from data ranges
  x_range <- range(df$x)
  y_range <- range(df$y)
  aspect_ratio <- diff(y_range) / diff(x_range)
  if (!is.finite(aspect_ratio) || aspect_ratio <= 0) aspect_ratio <- 1
  
  ## empty plot with correct aspect
  plot(df$x, df$y, type = "n", asp = aspect_ratio,
       xlab = "x", ylab = "y",
       main = paste0("ICA Trace: ", round(ica_total, 1),
                     " / ", ID))
  
  ## layered points by category (with alpha via rgb)
  with(subset(df, !(category %in% c("CD20","CD4","CD8a","CD11c"))),
       points(x, y, col = rgb(0.5, 0.5, 0.5, alpha = 0.8), pch = 16, cex = 0.65))  # Others (grey)
  
  with(subset(df, category == "CD20"),
       points(x, y, col = rgb(1, 0, 0, alpha = 0.5), pch = 16, cex = 0.45))        # red
  
  with(subset(df, category == "CD4"),
       points(x, y, col = rgb(0, 0, 1, alpha = 0.8), pch = 16, cex = 0.65))        # blue
  
  with(subset(df, category == "CD8a"),
       points(x, y, col = rgb(0, 1, 0, alpha = 0.8), pch = 16, cex = 0.65))        # green
  
  with(subset(df, category == "CD11c"),
       points(x, y, col = rgb(1, 0.5, 0, alpha = 0.8), pch = 16, cex = 0.65))      # orange
  
  ## ICA arrows
  arrows(arrow_data$mean_x, arrow_data$mean_y,
         arrow_data$ica_end_x1, arrow_data$ica_end_y1,
         col = "#00AEEF", length = 0.1, lwd = 6)
  
  arrows(arrow_data$mean_x, arrow_data$mean_y,
         arrow_data$ica_end_x2, arrow_data$ica_end_y2,
         col = "#00AEEF", length = 0.1, lwd = 6)
  
  ## legend
  legend("bottomleft",
         legend = c("CD20", "CD8a", "CD4", "CD11c", "Others"),
         col    = c("red", "green", "blue", "orange", "grey"),
         pt.cex = 1.5, pch = 16, title = "Cell Type", bty = "n")
}

## (plots render in the 2x2 grid)