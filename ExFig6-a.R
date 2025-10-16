####HLAA analysis 

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

#get the AMACR postive and section them into tumor and none-tumor and check the mean log expression of the HLAA to see if there is a difference in the cohort

markers<- c("Hoechst", "Ki67", "AMCAR", "HMWCK", "CD19", "SMA", "CD20", "CD11b", "CD68", "CD163", "CD4", "CD3d", "CD8a", "TCF1", "FOXP3", "PD1", "CD57", "CD11c", "GranzymeB", "CD15", "HLADR", "CD103", "CD31", "pTBK1", "HLAA", "CD24", "CD44", "CD206")
markers<- markers[length(markers):1]
mpflag<- numeric(length(markers))
for (i in 0:(length(markers)-1)){
  mpflag[i+1]<- 2^i
}
mpflag<- mpflag[length(mpflag):1]

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")


tHLA<- numeric(length(names(ldata)))
sHLA<- numeric(length(names(ldata)))
for(i in 1:length(names(ldata))){
  lsp<- ldata[[i]]
  t<- lsp[which(bitAnd(lsp$pflag, mpflag[26]) == mpflag[26] & bitAnd(lsp$cflag, 8) == 8), "HLAA"]
  s<- lsp[which(bitAnd(lsp$pflag, mpflag[26]) == mpflag[26] & bitAnd(lsp$cflag, 8) != 8), "HLAA"]
  tHLA[i]<- mean(t) / sd(t) 
  sHLA[i]<- mean(s) / sd(s)
}


# 1) assemble & clean
df <- data.frame(
  group = factor(rep(c("tHLA","sHLA"),
                     times = c(length(tHLA), length(sHLA)))),
  value = c(tHLA, sHLA)
)
df <- na.omit(df)

# 2) t‑test
tt  <- t.test(value ~ group, data = df)
p   <- signif(tt$p.value, 3)
# 


paired_df <- data.frame(
  sample = names(ldata),
  tHLA   = tHLA,
  sHLA   = sHLA
) %>%
  drop_na()

df_long <- paired_df %>%
  pivot_longer(cols = c(tHLA, sHLA),
               names_to = "Test",
               values_to = "Value")

## 3) Normality test: check differences (tHLA - sHLA)
diffs <- paired_df$tHLA - paired_df$sHLA
sh <- shapiro.test(diffs)

if (sh$p.value > 0.05) {
  test_res <- t.test(paired_df$tHLA, paired_df$sHLA, paired = TRUE)
  method_used <- "Paired t-test"
} else {
  test_res <- wilcox.test(paired_df$tHLA, paired_df$sHLA, paired = TRUE)
  method_used <- "Wilcoxon signed-rank test"
}

pval <- signif(test_res$p.value, 3)

## 4) Plot


# Define fill and point colours
fill_colors  <- c(tHLA = "#4e79a7", sHLA = "#f28e2b")

# Define the high and low sample vectors

# Add a Group column to df_long based on the sample belonging to high or low
df_long$Group <- ifelse(df_long$sample %in% high, "High", ifelse(df_long$sample %in% low, "Low", NA))

point_colors <- c(High = "#7570b3", Low = "#1b9e77")
ggplot(df_long, aes(x = Test, y = Value, fill = Test)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, size = 2, alpha = 0.6) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    title = paste0("tHLA vs sHLA (", method_used, ", p = ", pval, ")"),
    subtitle = paste("n =", nrow(paired_df), "paired samples"),
    x = NULL,
    y = "Mean log HLAA expression"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
