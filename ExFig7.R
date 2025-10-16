
#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

# Define High and Low groups
# low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
#          "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")
# 
# high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
#           "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
# 


sheet <- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")
high <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "High")

low <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "Low")


pdf(file = "ExtFig7a.pdf", width=10.5, height=13.125)
par(mfrow = c(5, 4), mar = c(1, 1, 1, 1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
#13000 is the max range of the points
for (j in 1:20){
  f<- ldata[[j]]
  f$x <- 2*ldata[[j]]$x
  f$y <- 2*ldata[[j]]$y
  ff<-f
  so<-sample(dim(f)[1], round(dim(f)[1]/30))
  f<- f[so,]
  #if (j == 15){f<- ldata[[j]]; f$x <- 2*ldata[[j]]$x; f$y <- 2*ldata[[j]]$y;ff<-f}
  name<- names(ldata)[j]
  determine_color <- function(nam) {
    if (nam %in% high) {
      return("#7570b3")
    } else if (nam %in% low) {
      return("#1b9e77")
    } else {
      return("grey")
    }
  }
  main_color <- determine_color(name)
  ltls<- length(unique(ff$tls_id))-1
  #ldbs<- length(unique(ff$dbscan_cluster))-3
  #uni1<- unique(f$dbscan_cluster)[-c(1,2,4)]
  uni<- unique(f$tls_id)[-1]
  a<- max(ff$x)/2; b<- max(ff$y)/2
  par(mar = c((10.66-(b/1500)), (10.66-(a/1500)), (10.66-(b/1500)), (10.66-(a/1500))))
  plot(f$x, f$y, pch=19, cex=0.01, cex.axis=1,  # Axis tick labels larger
       cex.lab=1,   # Axis title labels larger
       cex.main=1.4,    # Main title larger 
       col="lightgrey", main = paste(names(ldata)[[j]]), xlab = paste("BIC count", length(uni), sep = "="), ylab = "", ylim = range(f$y), xlim =range(f$x), col.main = main_color 
       #,panel.first = grid()
  )  
  if (sum(bitAnd(ff$cflag, 8) == 8) > 0){
    ss<- subset(ff, bitAnd(ff$cflag, 8) == 8)
    s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
    ss<- ss[s1, ]
    points(ss$x, ss$y, col="lightblue", pch=19, cex=0.01)
  }
  for (i in 1:ltls){
    s<- subset(ff, ff$tls_id==uni[i])
    points(s$x, s$y, col="orange", pch=19, cex=0.01)#rainbow(ltls)[i], pch=19, cex=0.01)
    if (dim(s)[1]>0){text(s$x[1], s$y[1], label = uni[i], cex = 0.7)}
  }
  #for (i in 1:ldbs){
  #  s<- subset(ff, ff$dbscan_cluster==uni1[i])
  #  points(s$x, s$y, col="lightgreen", pch=19, cex=0.01)#rainbow(ltls)[i], pch=19, cex=0.01)
  #  if (dim(s)[1]>0){text(s$x[1], s$y[1], label = uni1[i], cex = 0.7)}
  #}
  
}
dev.off()
dev.off()




##for BIC and TIC together 
for (j in 17:29){
  f<- ldata[[j]]
  f$x <- 2*ldata[[j]]$x
  f$y <- 2*ldata[[j]]$y
  ff<-f
  so<-sample(dim(f)[1], round(dim(f)[1]/30))
  f<- f[so,]
  if (j == 15){f<- ldata[[j]]}
  name<- names(ldata)[j]
  determine_color <- function(nam) {
    if (nam %in% high) {
      return("#7570b3")
    } else if (nam %in% low) {
      return("#1b9e77")
    } else {
      return("grey")
    }
  }
  main_color <- determine_color(name)
  ltls<- length(unique(ff$tls_id))-1
  ldbs<- length(unique(ff$dbscan_cluster))-3
  uni1<- unique(f$dbscan_cluster)[-c(1,2,4)]
  uni<- unique(f$tls_id)[-1]
  a<- max(ff$x)/2; b<- max(ff$y)/2
  par(mar = c((10.66-(b/1500)), (10.66-(a/1500)), (10.66-(b/1500)), (10.66-(a/1500))))
  plot(f$x, f$y, pch=19, cex=0.01, cex.axis=1,  # Axis tick labels larger
       cex.lab=1,   # Axis title labels larger
       cex.main=1.4,    # Main title larger 
       col="lightgrey", main = paste(names(ldata)[[j]]), xlab = paste("BIC count", length(uni), sep = "="), ylab = paste("TIC count", length(uni1), sep = "="), ylim = range(f$y), xlim =range(f$x), col.main = main_color 
       #,panel.first = grid()
  )  
  if (sum(bitAnd(ff$cflag, 8) == 8) > 0){
    ss<- subset(ff, bitAnd(ff$cflag, 8) == 8)
    s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
    ss<- ss[s1, ]
    points(ss$x, ss$y, col="lightblue", pch=19, cex=0.01)
  }
  for (i in 1:ltls){
    s<- subset(ff, ff$tls_id==uni[i])
    points(s$x, s$y, col="orange", pch=19, cex=0.01)#rainbow(ltls)[i], pch=19, cex=0.01)
    if (dim(s)[1]>0){text(s$x[1], s$y[1], label = uni[i], cex = 0.7)}
  }
  for (i in 1:ldbs){
    s<- subset(ff, ff$dbscan_cluster==uni1[i])
    points(s$x, s$y, col="lightgreen", pch=19, cex=0.01)#rainbow(ltls)[i], pch=19, cex=0.01)
    if (dim(s)[1]>0){text(s$x[1]+20, s$y[1]+20, label = uni1[i], cex = 0.7)}
  }
  
}

dev.off()