
#load("listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

imm.cruse <- function(ws, LSP, index, phen, fig = TRUE, creep = 1) {
  # Window size (ws) in micron, LSP identifies in the list of ldata object,
  # index to be calculated, phenotype for the indexing e.g B cells (reading off the coarse_phen_vec)
  L.models<- list()
  pws <- ws / 0.325 # pixel window size
  d <- ldata[[LSP]]
  xstart <- 1
  ystart <- 1
  if (phen %in% c("T cells", "B cells")){
    if (fig==TRUE){
      a<- max(d$x); b<- max(d$y)
      par(mar = c((10.66-(b/1500)), (10.66-(a/1500)), (10.66-(b/1500)), (10.66-(a/1500))))
      plot(d$x, d$y, pch=19, cex=0.01, cex.axis=1,  # Axis tick labels larger
           cex.lab=1,   # Axis title labels larger
           cex.main=1.4,    # Main title larger 
           col="lightgrey", main = LSP, xlab = paste("Window Size", ws, sep = "="), ylab = "", ylim = range(d$y), xlim =range(d$x), col.main = "navy"
           #,panel.first = grid()
      )  
      if (sum(bitAnd(d$cflag, 8) == 8) > 0){
        ss<- subset(d, bitAnd(d$cflag, 8) == 8)
        #s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
        #ss<- ss[s1, ]
        points(ss$x, ss$y, col=adjustcolor("grey57", alpha.f = 0.2), pch=19, cex=0.005)
      }
      if (sum(d$coarse_phen_vec == phen) > 0){
        if (phen == "T cells"){co<- "green"}
        if (phen == "B cells"){co<- "red"}
        ss<- subset(d, d$coarse_phen_vec == phen)
        #s1<-sample(dim(ss)[1], round(dim(ss)[1]/4))
        #ss<- ss[s1, ]
        points(ss$x, ss$y, col=adjustcolor(co, alpha.f = 0.4), pch=19, cex=0.005)
      }
      for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
      for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i + pws/2, col = adjustcolor("darkolivegreen", alpha.f = 0.5), lty = 3, lwd = 1)}
      for (j in 0:ceiling(max(d$x) / pws)){abline(v = xstart + pws * j, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
    }
    
    ny <- ceiling(max(d$y) / pws)
    nx <- ceiling(max(d$x) / pws)
    index.holder <- matrix(0, ny, nx)
    xstart <- 1
    ystart <- 1
    
    # Initialize progress bar
    total_steps <- nx * ny * creep^2
    pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    step <- 0
    c<-1
    for(k in 1:creep){
      xstart<- k*round(pws/creep)-round(pws/creep)
      for (l in 1:creep){
        ystart<- l*round(pws/creep)-round(pws/creep)
        for (i in 1:nx) {
          for (j in 1:ny) {
            data <- subset(d, d$x < xstart + pws * i & d$x > xstart + pws * (i - 1) & d$y < ystart + pws * j & d$y > ystart + pws * (j - 1))
            if (nrow(data) > 50) { # if the box has more than 10 cells
              phen.data <- subset(data, data$coarse_phen_vec == phen)
              if (nrow(phen.data) > 5) {
                if (index == "K-integral") {
                  ppp <- as.ppp(phen.data[, c("x", "y")], W = owin(c(min(phen.data$x), max(phen.data$x)), c(min(phen.data$y), max(phen.data$y))))
                  L <- Lest(ppp, rmax = ws)
                  index.holder[j, i] <- mean(L$border - L$theo, na.rm = TRUE)
                  differences <- L$border - L$theo
                  L.models[[c]]<- L
                  if(fig){
                    if(k + l < 3){
                      smoothed <- loess(differences ~ seq_along(differences), span = 0.3)
                      smoothed_values <- predict(smoothed, seq_along(differences))
                      lines(seq(xstart + pws * (i-1), xstart + pws * (i), pws/513)[-1], c(smoothed_values + pws/2 + pws*(j-1)), col='plum1', lwd = 2*ws/500)
                    }
                    D<- L$border - L$theo
                    Dp<- D[D>0]
                    text(x = xstart + pws * (i-0.5), y = ystart + pws/2 + pws*(j-1), labels = round(mean(Dp, na.rm = TRUE)), col = "plum4", cex = ws/700 + ws/700 * (round(mean(Dp, na.rm = TRUE))/700), font = 3)
                    # Add the smoothed line to the plot
                    #lines(smoothed_values, col = "red", lwd = 2)
                  }
                  c<- c+1
                }
              }
            }
            # Update progress bar
            step <- step + 1
            setTxtProgressBar(pb, step)
          }
        }
      }
    }
  }
  if(phen == "Both"){#only plotting option for both
    if (fig==TRUE){
      plot(d$x, d$y, pch = 19, cex = 0.01, 
           cex.axis = 1,   # Axis tick labels larger
           cex.lab  = 1,   # Axis title labels larger
           cex.main = 1.4, # Main title larger
           col      = "lightgrey", 
           main     = LSP, 
           xlab     = paste("Window Size", ws*6, sep = "="), 
           ylab     = "", 
           ylim     = range(d$y), 
           xlim     = range(d$x), 
           col.main = "navy",
           axes     = FALSE)   # suppress default axes
      
      # Redraw axes with ticks multiplied by 6
      #axis(1, at = axTicks(1), labels = axTicks(1) * 2)  # X-axis
      #axis(2, at = axTicks(2), labels = axTicks(2) * 2)  # Y-axis
      box()   # redraw box around plot
      if (sum(bitAnd(d$cflag, 8) == 8) > 0){
        ss<- subset(d, bitAnd(d$cflag, 8) == 8)
        #s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
        #ss<- ss[s1, ]
        points(ss$x, ss$y, col=adjustcolor("grey57", alpha.f = 0.2), pch=19, cex=0.005)
      }
      if (sum(d$coarse_phen_vec %in% c("B cells", "T cells")) > 0){
        coT<- "green"
        coB<- "red"
        ssB<- subset(d, d$coarse_phen_vec == "B cells")
        ssT<- subset(d, d$coarse_phen_vec == "T cells")
        #s1<-sample(dim(ss)[1], round(dim(ss)[1]/4))
        #ss<- ss[s1, ]
        points(ssT$x, ssT$y, col=adjustcolor(coT, alpha.f = 0.4), pch=19, cex=0.005)
        points(ssB$x, ssB$y, col=adjustcolor(coB, alpha.f = 0.4), pch=19, cex=0.005)
      }
      for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
      for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i + pws/2, col = adjustcolor("darkolivegreen", alpha.f = 0.5), lty = 3, lwd = 1)}
      for (j in 0:ceiling(max(d$x) / pws)){abline(v = xstart + pws * j, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
    }
    ny <- ceiling(max(d$y) / pws)
    nx <- ceiling(max(d$x) / pws)
    index.holder <- matrix(0, ny, nx)
    xstart <- 1
    ystart <- 1
    
    total_steps <- nx * ny * creep^2
    pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    step <- 0
    c<-1
    for(k in 1:creep){
      xstart<- k*round(pws/creep)-round(pws/creep)
      for (l in 1:creep){
        ystart<- l*round(pws/creep)-round(pws/creep)
        for (i in 1:nx) {
          for (j in 1:ny) {
            data <- subset(d, d$x < xstart + pws * i & d$x > xstart + pws * (i - 1) & d$y < ystart + pws * j & d$y > ystart + pws * (j - 1))
            if (nrow(data) > 50) { # if the box has more than 10 cells
              phen.data <- subset(data, data$coarse_phen_vec == "T cells")
              if (nrow(phen.data) > 5) {
                if (index == "K-integral") {
                  ppp <- as.ppp(phen.data[, c("x", "y")], W = owin(c(min(phen.data$x), max(phen.data$x)), c(min(phen.data$y), max(phen.data$y))))
                  L <- Lest(ppp, rmax = ws)
                  index.holder[j, i] <- mean(L$border - L$theo, na.rm = TRUE)
                  differences <- L$border - L$theo
                  L.models[[c]]<- L
                  if(fig){
                    if(k + l < 3){
                      smoothed <- loess(differences ~ seq_along(differences), span = 0.3)
                      smoothed_values <- predict(smoothed, seq_along(differences))
                      lines(seq(xstart + pws * (i-1), xstart + pws * (i), pws/513)[-1], c(smoothed_values + pws/2 + pws*(j-1)), col='darkgreen', lwd = 2*ws/500, lty =2)
                    }
                    D<- L$border - L$theo
                    Dp<- D[D>0]
                    text(x = xstart + pws * (i-0.85), y = ystart + pws/10 + pws*(j-1), labels = round(mean(Dp, na.rm = TRUE)), col = "darkgreen", cex = ws/700 + ws/700 * (round(mean(Dp, na.rm = TRUE))/700), font = 3)
                    # Add the smoothed line to the plot
                    #lines(smoothed_values, col = "red", lwd = 2)
                  }
                  c<- c+1
                }
              }
            }
            # Update progress bar
            step <- step + 1
            setTxtProgressBar(pb, step)
          }
        }
      }
    }
    for(k in 1:creep){
      xstart<- k*round(pws/creep)-round(pws/creep)
      for (l in 1:creep){
        ystart<- l*round(pws/creep)-round(pws/creep)
        for (i in 1:nx) {
          for (j in 1:ny) {
            data <- subset(d, d$x < xstart + pws * i & d$x > xstart + pws * (i - 1) & d$y < ystart + pws * j & d$y > ystart + pws * (j - 1))
            if (nrow(data) > 50) { # if the box has more than 10 cells
              phen.data <- subset(data, data$coarse_phen_vec == "B cells")
              if (nrow(phen.data) > 5) {
                if (index == "K-integral") {
                  ppp <- as.ppp(phen.data[, c("x", "y")], W = owin(c(min(phen.data$x), max(phen.data$x)), c(min(phen.data$y), max(phen.data$y))))
                  L <- Lest(ppp, rmax = ws)
                  index.holder[j, i] <- mean(L$border - L$theo, na.rm = TRUE)
                  differences <- L$border - L$theo
                  L.models[[c]]<- L
                  if(fig){
                    if(k + l < 3){
                      smoothed <- loess(differences ~ seq_along(differences), span = 0.3)
                      smoothed_values <- predict(smoothed, seq_along(differences))
                      lines(seq(xstart + pws * (i-1), xstart + pws * (i), pws/513)[-1], c(smoothed_values + pws/2 + pws*(j-1)), col='maroon', lwd = 2*ws/500, lty=2)
                    }
                    D<- L$border - L$theo
                    Dp<- D[D>0]
                    text(x = xstart + pws * (i-0.85), y = ystart + pws/3 + pws*(j-1), labels = round(mean(Dp, na.rm = TRUE)), col = "maroon", cex = ws/700 + ws/700 * (round(mean(Dp, na.rm = TRUE))/700), font = 3)
                    # Add the smoothed line to the plot
                    #lines(smoothed_values, col = "red", lwd = 2)
                  }
                  c<- c+1
                }
              }
            }
          }
        }
      }
    }
  }
  # Close progress bar
  close(pb)
  
  return(L.models)
}

imm.cruse(500, names(ldata)[1], "K-integral", phen="Both", fig = TRUE, creep = 1)
