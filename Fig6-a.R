
rm(list = setdiff(ls(), c("ldata", "df")))


library(RANN)             # nn2
library(dbscan)           # frNN
library(spatstat.geom)    # ppp, owin
library(spatstat.explore) # Lcross
library(ggplot2)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# YOUR INPUT VECTORS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
samples             <- names(ldata)
high                <- c("LSP12601","LSP12621","LSP12623", "LSP12625","LSP12627",
                         "LSP12629","LSP12635","LSP12639","LSP12641",
                         "LSP12649","LSP12651","LSP12653","LSP12655","LSP12657")
low                 <- c("LSP12603", "LSP12605","LSP12607","LSP12609","LSP12611","LSP12613",
                         "LSP12615","LSP12617","LSP12619","LSP12631",
                         "LSP12633","LSP12637","LSP12643","LSP12645","LSP12647")
high_granzymeB_gate <- c(7,7,6.9,7,7.04,6.8,6.8,6.84,6.95,6.954,6.87,6.88,6.95,6.94)
low_granzymeB_gate  <- c(7,6.94,7.10,7.3,7.09,7.01,7.21,6.82,7.11,7.26,7.04,6.86,7.51,6.96,7.1)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1. BUILD Gr_Tu_ind
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
group <- ifelse(samples %in% high, "high", "low")
gate  <- numeric(length(samples))

# fill gate by matching sample in high/low
for(i in seq_along(samples)) {
  if(group[i]=="high")
    gate[i] <- high_granzymeB_gate[ match(samples[i], high) ]
  else
    gate[i] <- low_granzymeB_gate [ match(samples[i], low ) ]
}

# … assume ldata, high, low, high_granzymeB_gate, low_granzymeB_gate already defined

# Restrict to exactly those 26 samples
samples <- intersect(names(ldata), union(high, low))  # length 26

# Re-build Gr_Tu_ind from just those
group <- ifelse(samples %in% high, "high", "low")
gate  <- numeric(length(samples))
for (i in seq_along(samples)) {
  if (group[i]=="high")
    gate[i] <- high_granzymeB_gate[ match(samples[i], high) ]
  else
    gate[i] <- low_granzymeB_gate[  match(samples[i], low)  ]
}
Gr_Tu_ind <- data.frame(sample=samples, group=group, gate=gate, 
                        stringsAsFactors=FALSE)

# And initialize results_df likewise:
results_df <- data.frame(
  sample       = samples,
  group        = group,
  gate         = gate,
  # … the same 13 summary columns as before
  tumor_count  = NA_integer_,
  ctl_count    = NA_integer_,
  dist_q1      = NA_real_,
  contact_frac = NA_real_,
  stringsAsFactors = FALSE
)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3. MAIN LOOP OVER SAMPLES
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
radii <- c(10, 20, 50)
contact_cut <- 3
kcross_r    <- seq(0, 100, by=5)
for(i in seq_len(nrow(Gr_Tu_ind))) {
  samp  <- Gr_Tu_ind$sample[i]
  gate  <- Gr_Tu_ind$gate [i]
  lp    <- ldata[[samp]]
  
  # --- build masks & coords as before ---
  tumor_mask <- bitwAnd(lp$cflag,8)==8
  ctl_mask   <- lp$phen_vec=="CD8 T cells" & log(lp$GranzymeB)>gate
  tumor_mask <- tumor_mask & !ctl_mask
  tumor_xy   <- as.matrix(lp[tumor_mask, c("x","y")])
  ctl_xy     <- as.matrix(lp[   ctl_mask, c("x","y")])
  
  # --- EARLY GUARD: if either set is empty, fill NA/0 and skip ---
  if(nrow(tumor_xy)==0 || nrow(ctl_xy)==0) {
    results_df[i, c("tumor_count","ctl_count",
                    "dist_q1","dist_med","dist_q3",
                    "prop_10","prop_20","prop_50",
                    "dens_q1","dens_med","dens_q3",
                    "kcross_avg","contact_frac")] <- NA
    next
  }
  
  # --- otherwise proceed as before ---
  results_df$tumor_count[i] <- nrow(tumor_xy)
  results_df$ctl_count  [i] <- nrow(ctl_xy)
  
  # 1) nearest distances
  nn    <- nn2(tumor_xy, ctl_xy, k=1)
  dists <- nn$nn.dists[,1]
  # etc…
}
for(i in seq_len(nrow(Gr_Tu_ind))[c(-16,-22, -23)]) {#samples without tumor!
  samp  <- Gr_Tu_ind$sample[i]
  gate  <- Gr_Tu_ind$gate [i]
  lp    <- ldata[[samp]]
  
  # a) identify masks
  tumor_mask <- bitwAnd(lp$cflag,8)==8
  ctl_mask   <- lp$phen_vec=="CD8 T cells" & log(lp$GranzymeB) > gate
  tumor_mask <- tumor_mask & !ctl_mask
  
  # b) coords
  tumor_xy <- as.matrix(lp[tumor_mask, c("x","y")])
  ctl_xy   <- as.matrix(lp[   ctl_mask, c("x","y")])
  
  # store counts
  results_df$tumor_count[i] <- nrow(tumor_xy)
  results_df$ctl_count[i]   <- nrow(ctl_xy)
  
  # define window
  win <- owin(range(lp$x), range(lp$y))
  
  #––– 1) DISTANCES
  nn    <- nn2(tumor_xy, ctl_xy, k=1)
  dists <- nn$nn.dists[,1]
  qs    <- quantile(dists, probs=c(0.25,0.5,0.75))
  results_df$dist_q1[i]  <- qs[1]
  results_df$dist_med[i] <- qs[2]
  results_df$dist_q3[i]  <- qs[3]
  
  #––– 2) PROPORTIONS
  pw <- sapply(radii, function(r) mean(dists < r))
  results_df$prop_10[i] <- pw[1]
  results_df$prop_20[i] <- pw[2]
  results_df$prop_50[i] <- pw[3]
  
  #––– 3) LOCAL DENSITY (r = 50)
  R    <- 50
  nctl <- nrow(ctl_xy)
  counts <- integer(nctl)
  chunk  <- 10000
  for(start in seq(1, nctl, by=chunk)) {
    end   <- min(start+chunk-1, nctl)
    sub_q <- ctl_xy[start:end,,drop=FALSE]
    fr    <- frNN(tumor_xy, eps=R, query=sub_q)
    counts[start:end] <- lengths(fr$id)
  }
  dens   <- counts / (pi * R^2)
  qd     <- quantile(dens, probs=c(0.25,0.5,0.75))
  results_df$dens_q1[i]  <- qd[1]
  results_df$dens_med[i] <- qd[2]
  results_df$dens_q3[i]  <- qd[3]
  
  #––– 4) CROSS‐K
  all_xy <- rbind(ctl_xy, tumor_xy)
  marks  <- factor(c(rep("CTL", nrow(ctl_xy)), rep("Tumor", nrow(tumor_xy))))
  pp     <- ppp(all_xy[,1], all_xy[,2], window=win, marks=marks)
  Kcr    <- Lcross(pp, i="CTL", j="Tumor", r=kcross_r)
  results_df$kcross_avg[i] <- mean(Kcr$border - Kcr$theo)
  
  #––– 5) CONTACT FRACTION
  results_df$contact_frac[i] <- mean(dists < contact_cut)
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4. QUICK LOOK
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
print(Gr_Tu_ind)     # 26×3
print(results_df)    # 26×16
GranzymeTumor <-cbind(Gr_Tu_ind, results_df)
# assuming GranzymeTumor, ldata are already in your workspace

GranzymeTumor$ctl_prop <- mapply(function(samp, ctls){
  # total CD8+ cells in that sample
  total_cd8 <- sum(ldata[[samp]]$phen_vec == "CD8 T cells")
  # fraction of those that are GZMB+
  ctls / total_cd8
}, 
GranzymeTumor$sample, GranzymeTumor$ctl_count)

# Quick check
head(GranzymeTumor[, c("sample","ctl_count","ctl_prop")])
GranzymeTumor$dist_med <- log(GranzymeTumor$dist_med)



# Re-filter NAs & top 1%
df0    <- subset(GranzymeTumor, !is.na(ctl_prop))
cutoff <- quantile(df0$ctl_prop, 0.99)
df     <- subset(df0, ctl_prop <= cutoff)

# One-sided Wilcoxon
p_value     <- wilcox.test(ctl_prop ~ group, data = df,
                           alternative = "greater", exact = FALSE)$p.value
overall_frac <- mean(df$ctl_prop)

ggplot(df, aes(x = group, y = ctl_prop, fill = group)) +
  
  # boxplot without outliers
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  
  # larger, semi-transparent jittered points
  geom_jitter(width = 0.15, shape = 19, size = 2.5, alpha = 0.4, color = "black") +
  
  # your purple/green fills
  scale_fill_manual(values = c(high = "#7570b3", low = "#1b9e77")) +
  
  # raw fraction on y-axis, 3 decimal places
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  
  labs(
    title    = "CTL fraction",
    subtitle = sprintf(
      "Avg = %.3f; one-sided Wilcoxon p = %.3f",
      overall_frac, p_value
    ),
    x = NULL,
    y = "CTL Fraction"
  ) +
  
  theme_minimal() +
  theme(
    legend.position  = "none",
    plot.subtitle    = element_text(
      color = if (p_value < 0.05) "red" else "black",
      face  = if (p_value < 0.05) "bold" else "plain"
    ),
    axis.text.x      = element_text(size = 12, face = "bold"),
    axis.text.y      = element_text(size = 10),
    axis.title.y     = element_text(size = 12)
  )
