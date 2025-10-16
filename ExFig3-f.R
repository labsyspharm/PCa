

low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617",
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")
high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635",
          "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")


#high ICAT and high Gleason Q3
tlf_hc_hg<- subset(df, df$count < 20000 & df$ICAT > 32 & df$second_part %in% high) 
i<-1
speci<- strsplit(tlf_hc_hg$TLS[i], "@")[[1]][2]
ID<- strsplit(tlf_hc_hg$TLS[i], "@")[[1]][1]
d_hc_hg<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]

for (i in 2:dim(tlf_hc_hg)[1]){
  speci<- strsplit(tlf_hc_hg$TLS[i], "@")[[1]][2]
  ID<- strsplit(tlf_hc_hg$TLS[i], "@")[[1]][1]
  d<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]
  d_hc_hg<- rbind(d_hc_hg, d)
}
sum(tlf_hc_hg$count) == dim(d_hc_hg)[1]
#high ICAT and low Gleason Q3
tlf_hc_lg<- subset(df, df$count < 20000 & df$ICAT > 32 & df$second_part %in% low) 
i<-1
speci<- strsplit(tlf_hc_lg$TLS[i], "@")[[1]][2]
ID<- strsplit(tlf_hc_lg$TLS[i], "@")[[1]][1]
d_hc_lg<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]

for (i in 2:dim(tlf_hc_lg)[1]){
  speci<- strsplit(tlf_hc_lg$TLS[i], "@")[[1]][2]
  ID<- strsplit(tlf_hc_lg$TLS[i], "@")[[1]][1]
  d<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]
  d_hc_lg<- rbind(d_hc_lg, d)
}
sum(tlf_hc_lg$count) == dim(d_hc_lg)[1]

#low ICAT and high Gleason Q3
tlf_lc_hg<- subset(df, df$count < 20000 & df$ICAT < 13.5 & df$second_part %in% high) 
i<-1
speci<- strsplit(tlf_lc_hg$TLS[i], "@")[[1]][2]
ID<- strsplit(tlf_lc_hg$TLS[i], "@")[[1]][1]
d_lc_hg<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]

for (i in 2:dim(tlf_lc_hg)[1]){
  speci<- strsplit(tlf_lc_hg$TLS[i], "@")[[1]][2]
  ID<- strsplit(tlf_lc_hg$TLS[i], "@")[[1]][1]
  d<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]
  d_lc_hg<- rbind(d_lc_hg, d)
}
sum(tlf_lc_hg$count) == dim(d_lc_hg)[1]

#low ICAT and low Gleason Q3
tlf_lc_lg<- subset(df, df$count < 20000 & df$ICAT < 13.5 & df$second_part %in% low) 
i<-1
speci<- strsplit(tlf_lc_lg$TLS[i], "@")[[1]][2]
ID<- strsplit(tlf_lc_lg$TLS[i], "@")[[1]][1]
d_lc_lg<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]

for (i in 2:dim(tlf_lc_lg)[1]){
  speci<- strsplit(tlf_lc_lg$TLS[i], "@")[[1]][2]
  ID<- strsplit(tlf_lc_lg$TLS[i], "@")[[1]][1]
  d<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]
  d_lc_lg<- rbind(d_lc_lg, d)
}
sum(tlf_lc_lg$count) == dim(d_lc_lg)[1]
# install.packages("ggvenn")  # if you haven’t already
library(ggvenn)
# install.packages("patchwork")
library(ggvenn)
library(ggplot2)
library(patchwork)

library(ggvenn)
library(ggplot2)
library(patchwork)

make_venn <- function(df, title){
  sets <- list(
    `T cells`   = which(df$coarse_phen_vec == "T cells" | bitAnd(df$pflag, 1024) == 1024 | bitAnd(df$pflag, 4096) == 4096 | bitAnd(df$pflag, 2048) == 2048),
    `TCF1gate⁺` = which(df$TCF1gate   == TRUE),
    `KiGate⁺`   = which(df$KiGate     == 1),
    `B cells`   = which(df$coarse_phen_vec == "B cells" & !bitAnd(df$pflag, 1024) == 1024 & !bitAnd(df$pflag, 4096) == 4096 & !bitAnd(df$pflag, 2048) == 2048)
  )
  ggvenn(
    sets,
    fill_color      = c("seagreen","orchid","skyblue","tomato"),
    stroke_size     = 0.5,
    set_name_size   = 5,
    show_percentage = TRUE,
    text_size       = 4
  ) +
    theme_void() +                                  # remove axes & grid
    theme(
      plot.title    = element_text(face="bold"),   # re-add title styling
      plot.subtitle = element_text(size=8),
      legend.position = "right"                    # keep the legend
    ) +
    labs(
      title    = title,
      subtitle = "Numbers are % of total cells"
    )
}

p1 <- make_venn(d_hc_hg, "HcHg")
p2 <- make_venn(d_hc_lg, "HcLg")
p3 <- make_venn(d_lc_hg, "LcHg")
p4 <- make_venn(d_lc_lg, "LcLg")

# 2×2 layout
(p1 | p2) /
  (p3 | p4)


make_venn_dist <- function(df) {
  # four logical vectors
  isT  <- df$coarse_phen_vec == "T cells" | bitAnd(df$pflag, 1024) == 1024 | bitAnd(df$pflag, 4096) == 4096 | bitAnd(df$pflag, 2048) == 2048
  isB  <- df$coarse_phen_vec == "B cells" & !bitAnd(df$pflag, 1024) == 1024 & !bitAnd(df$pflag, 4096) == 4096 & !bitAnd(df$pflag, 2048) == 2048
  isK  <- df$KiGate   == 1
  isT1 <- df$TCF1gate == TRUE
  
  # all 16 patterns of membership
  patterns <- expand.grid(T= c(FALSE,TRUE),
                          B= c(FALSE,TRUE),
                          K= c(FALSE,TRUE),
                          T1=c(FALSE,TRUE))
  
  # count rows matching each
  counts <- sapply(seq_len(nrow(patterns)), function(i){
    pat <- patterns[i,]
    sum(isT  == pat$T &
          isB  == pat$B &
          isK  == pat$K &
          isT1 == pat$T1)
  })
  
  # turn into probability
  probs <- counts / sum(counts)
  
  # attach the pattern labels
  lab <- apply(patterns, 1, function(r) paste0(
    if(r["T"])  "T" else "",
    if(r["B"])  "B" else "",
    if(r["K"])  "K" else "",
    if(r["T1"]) "1" else ""
  ))
  # rename the zero‐set as “none”
  lab[lab==""] <- "none"
  names(probs) <- lab
  probs
}
# --- 1) Build the four full distributions as before ---
venn_dists <- list(
  HcHg = make_venn_dist(d_hc_hg),
  HcLg = make_venn_dist(d_hc_lg),
  LcHg = make_venn_dist(d_lc_hg),
  LcLg = make_venn_dist(d_lc_lg)
)

# --- 2) Prune away any pattern that's zero in *all* four  ---
# bind into a matrix: rows = patterns, cols = your cohorts
mat <- do.call(cbind, venn_dists)
# keep only rows whose sum > 0 (i.e. at least one cohort has that pattern)
keep <- rowSums(mat) > 0

# trim each vector
venn_dists_pruned <- lapply(venn_dists, function(v) v[keep])

# --- 3) Re‐compute all pairwise tests on the pruned distributions ---
pairs <- combn(names(venn_dists_pruned), 2, simplify=FALSE)
results_pruned <- lapply(pairs, function(pp){
  a <- pp[1]; b <- pp[2]
  pA <- venn_dists_pruned[[a]]
  pB <- venn_dists_pruned[[b]]
  
  # build counts for chi2
  N   <- 1e5
  obs <- rbind(round(pA*N), round(pB*N))
  
  # 1) Monte‐Carlo χ²
  chi <- suppressWarnings(
    chisq.test(obs, simulate.p.value=TRUE, B= 10000)
  )
  # 2) JSD
  jsd <- philentropy::JSD(rbind(pA, pB), unit="log2")
  
  data.frame(
    group1 = a, group2 = b,
    chi2_p = chi$p.value,
    JSD    = jsd,
    stringsAsFactors = FALSE
  )
})

res_tbl_pruned <- dplyr::bind_rows(results_pruned) %>%
  dplyr::arrange(chi2_p)

print(res_tbl_pruned)

