
high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

#patients for umap
up<- c(1,4,6, 7, 8, 9,10,17, 18, 19,20, 21, 24, 25, 26, 27, 28, 29)
#markers for umap
um<- c(9, 12, 13, 14, 15, 16, 17, 18, 19, 21, 24, 29, 3, 4, 129, 131, 135) 

set.seed(12345)

sldata<- list()
for (i in up){
  lsp<- ldata[[i]]
  Tlsp<- subset(lsp, bitAnd(lsp$cflag, 8) == 8) 
  Slsp<- subset(lsp, bitAnd(lsp$cflag, 8) != 8)
  
  tlsp<- Tlsp[sample(1:dim(Tlsp)[1], 10000, replace = TRUE), um]
  slsp<- Slsp[sample(1:dim(Slsp)[1], 10000, replace = TRUE), um]
  
  sldata[[i]] <- rbind(tlsp, slsp)
}

####removing the top and lowever percent 
clean_sldata <- list()

for (i in up) {
  ddf <- sldata[[i]]
  df<-ddf[, -c(13,14, 15, 16,17)]
  # For each column (marker), calculate 1% and 99% quantiles
  lower_bounds <- apply(df, 2, quantile, probs = 0.01, na.rm = TRUE)
  upper_bounds <- apply(df, 2, quantile, probs = 0.99, na.rm = TRUE)
  
  # Keep rows that are within bounds for all markers
  in_range <- apply(df, 1, function(row) {
    all(row >= lower_bounds & row <= upper_bounds)
  })
  
  clean_sldata[[i]] <- sldata[[i]][in_range, ]
}

##make the matrix for the umap
udata<- clean_sldata[[up[1]]]
for (i in up[-1]){
  udata<- rbind(udata, clean_sldata[[i]])
}

udata[which(udata$coarse_phen_vec=="Others" & bitAnd(udata$cflag, 8) == 8), which(colnames(udata)=="coarse_phen_vec")] <- "Tumor"

udata[which(udata$phen_vec=="Others" & bitAnd(udata$cflag, 8) == 8), which(colnames(udata)=="phen_vec")] <- "Tumor"

udata$Gleason <- numeric(dim(udata)[1])

udata$Gleason <- ifelse(udata$LSPid %in% high, "HGG", "LGG")

udata[, -c(13:18)] <- scale(udata[, -c(13:18)])

####normalizing with the gates valurs as 0.5

gscale <- function(udata){# a function to scale the values to 0 to 1 with gate being 0.5
  mg<- colnames(udata)[1:12]
  pcode<- c(4, 32,64, 128, 256, 512, 1024, 2048, 4096, 16384, 131072, 4194304)
  sp<- unique(udata$LSPid)
  sgdata<- udata
  
  rescale_to_half_one <- function(v) {
    v_min <- min(v, na.rm = TRUE)
    v_max <- max(v, na.rm = TRUE)
    0.5 + 0.5 * (v - v_min) / (v_max - v_min)
  }
  
  rescale_to_zero_half <- function(v) {
    v_min <- min(v, na.rm = TRUE)
    v_max <- max(v, na.rm = TRUE)
    0.5 * (v - v_min) / (v_max - v_min)
  }
  
  
  for (j in sp){
    for (i in 1:length(mg)){
      gate<- min(sgdata[which(bitAnd(udata$pflag, pcode[i]) == pcode[i] & udata$LSPid == j) , i])
      v_up<-sgdata[which(bitAnd(udata$pflag, pcode[i]) == pcode[i] & udata$LSPid == j) , i]
      v_lo<-sgdata[which(bitAnd(udata$pflag, pcode[i]) != pcode[i] & udata$LSPid == j) , i]
      sgdata[which(bitAnd(udata$pflag, pcode[i]) == pcode[i] & udata$LSPid == j) , i] <- rescale_to_half_one(v_up)
      sgdata[which(bitAnd(udata$pflag, pcode[i]) != pcode[i] & udata$LSPid == j) , i] <- rescale_to_zero_half(v_lo)
    }
  }
  return(sgdata)
}

gs<- gscale(udata = udata)

library(uwot)
set.seed(42)  # For reproducibility
gs_clean <- gs[complete.cases(gs[, -c(13, 14, 15, 16, 17, 18)]), ]

f3 <- uwot::umap(gs_clean[, -c(13, 14, 15, 16, 17, 18)],
                 n_neighbors = 50, min_dist = 0.05, metric = "euclidean")


umap_result <- f3
# Prepare a new data frame with UMAP coordinates
umap_df <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  phen_vec = gs_clean$phen_vec,
  coarse_phen_vec = gs_clean$coarse_phen_vec,
  LSPid = gs_clean$LSPid,
  Gleason = gs_clean$Gleason
)



library(ggplot2)
library(patchwork)


#pdf("umap.pdf",  height = 10, width = 12)

shared_layers <- list(
  coord_fixed(ratio = 1),     # square aspect
  theme_void(),               # empty background
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5)
  )
)

# -----------------------------
p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = phen_vec)) +
  geom_point(size = 0.005, alpha = 0.07) +
  labs(title = "Phenotype") +
  shared_layers +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = coarse_phen_vec)) +
  geom_point(size = 0.005, alpha = 0.08) +
  labs(title = "Cell type") +
  shared_layers +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

p3 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = LSPid)) +
  geom_point(size = 0.005, alpha = 0.07) +
  labs(title = "Specimen") +
  shared_layers +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

p4 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Gleason)) +
  geom_point(size = 0.005, alpha = 0.05) +
  labs(title = "Gleason") +
  shared_layers +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# -----------------------------
# Combine panels and collect legends
# -----------------------------
(p1 | p2) / (p3 | p4) + plot_layout(guides = "collect")
#dev.off()
