
###external sample 

csv<- read.csv("/Users/aliamiryousefi/Desktop/LSP30197_P198.csv")

##communities 198, 79, and 132

## Packages
library(dbscan)
library(fastICA)
library(RColorBrewer)

## 0) Threshold for CD20+
hi <- log(csv$CD20) > 6.5

## 1) Get the BIC clusters on CD20+ coordinates
# Tune eps/minPts to your pixel/µm scale
eps_val  <- 100    
min_pts  <- 10    # These were visually optimized 

cd20_xy  <- csv[hi, c("X_centroid", "Y_centroid")]
db       <- dbscan(cd20_xy, eps = eps_val, minPts = min_pts)## using dbscan as we are missing the other immune markers for the KNN

## 2) Attach community labels back to full data
csv$community <- NA_integer_
csv$community[hi] <- db$cluster  # 0 = noise

## 3) ICAT function (generalized to work on any community)
# Your original function assumed BT_community==1; we generalize:
F_icat_extra <- function(d_comm, scale_xy = 0.325) {
  # d_comm: data.frame with columns X_centroid, Y_centroid for ONE community
  if (nrow(d_comm) < 3) return(NA_real_)   # guard for tiny clusters
  
  X  <- as.matrix(d_comm[, c("X_centroid", "Y_centroid")] * scale_xy)
  mu <- colMeans(X, na.rm = TRUE)
  
  ica <- fastICA(X, n.comp = 2)
  S   <- ica$S          # n x 2
  A   <- ica$A          # 2 x 2
  
  Xhat_centered <- S %*% t(A)
  Xhat <- sweep(Xhat_centered, 2, mu, "+")
  
  # trace of SD (your formula)
  vx <- var(Xhat[,1]); vy <- var(Xhat[,2])
  traceSD <- sqrt(vx + vy + 2*sqrt(vx*vy))
  ICAT <- 100 * traceSD / nrow(X)  # per your definition
  return(ICAT)
}

## 4) Compute ICAT per community (>0)
clusters <- sort(unique(csv$community[csv$community > 0]))
icat_by_cluster <- setNames(numeric(length(clusters)), clusters)

for (k in clusters) {
  d_k <- subset(csv, community == k, select = c("X_centroid","Y_centroid"))
  icat_by_cluster[as.character(k)] <- F_icat_extra(d_k)
}

## 5) Add ICAT values back to each CD20+ row (NA for others/noise)
csv$ICAT_community <- NA_real_
idx <- which(!is.na(csv$community) & csv$community > 0)
csv$ICAT_community[idx] <- icat_by_cluster[as.character(csv$community[idx])]

## 6) Plot (fast base R): all cells faint, then communities colored & labeled
plot(csv$X_centroid, csv$Y_centroid,
     pch = ".", col = "lightgrey",
     xlab = "X", ylab = "Y", main = "CD20⁺ Communities (DBSCAN, ICAT labels)")

# palette for clusters (max 12 distinct; recycle if more)
n_clust <- max(1, length(clusters))
cols <- brewer.pal(max(3, min(12, n_clust)), "Set3")
if (n_clust > length(cols)) cols <- rep(cols, length.out = n_clust)
names(cols) <- clusters

# plot noise (optional): CD20+ but unclustered
noise_idx <- which(csv$community == 0)
if (length(noise_idx) > 0) {
  points(csv$X_centroid[noise_idx], csv$Y_centroid[noise_idx],
         pch = ".", col = rgb(0.4, 0.4, 0.4, 0.6))
}

# plot each cluster and label with ICAT at its centroid
for (k in c(198, 79, 195)) {
  j <- which(csv$community == k)
  points(csv$X_centroid[j], csv$Y_centroid[j],
         pch = ".", col = cols[as.character(k)])
  # label position at cluster mean
  cx <- mean(csv$X_centroid[j]); cy <- mean(csv$Y_centroid[j])
  lab <- sprintf("ID = %d ICAT= %.1f",k, icat_by_cluster[as.character(k)])
  text(cx, cy, labels = lab, cex = 0.8, col = "black", pos = 3)
}

