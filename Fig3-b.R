

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


ID<- "4@LSP12611"

a<- strsplit(ID, "@")[[1]][2]
b<- as.numeric(strsplit(ID, "@")[[1]][1])
F_icat(a,b)

ID<- "5@LSP12639"

a<- strsplit(ID, "@")[[1]][2]
b<- as.numeric(strsplit(ID, "@")[[1]][1])
F_icat(a,b)
