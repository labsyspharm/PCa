

rm(list = setdiff(ls(), c("ldata")))


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


# Compute fractions for high group
high_fractions <- numeric(length(high))
for (i in seq_along(high)) {
  sample <- high[i]
  gate <- high_granzymeB_gate[i]
  df <- ldata[[sample]]
  numerator <- sum(bitAnd(df$pflag, 32768) == 32768 & log(df$GranzymeB) > gate & df$phen_vec == "CD8 T cells" & bitAnd(df$cflag, 8) == 8)
  denominator <- sum(df$coarse_phen_vec == "T cells" & bitAnd(df$cflag, 8) == 8)
  high_fractions[i] <- numerator / denominator
}

# Compute fractions for low group
low_fractions <- numeric(length(low))
for (i in seq_along(low)) {
  sample <- low[i]
  gate <- low_granzymeB_gate[i]
  df <- ldata[[sample]]
  numerator <- sum(bitAnd(df$pflag, 32768) == 32768 & log(df$GranzymeB) > gate & df$phen_vec == "CD8 T cells" & bitAnd(df$cflag, 8) == 8)
  denominator <- sum(df$coarse_phen_vec == "T cells" & bitAnd(df$cflag, 8) == 8)
  low_fractions[i] <- numerator / denominator
}

# Create data frame for plotting and testing
df <- data.frame(
  Group = factor(c(rep("HGG", length(high)), rep("LGG", length(low))),
                 levels = c("HGG", "LGG")),
  Fraction = c(high_fractions, low_fractions)
)

# Perform Mann-Whitney U test (Wilcoxon rank sum test)
test_result <- wilcox.test(Fraction ~ Group, data = df)
pval <- round(test_result$p.value, 3)  # Adjust rounding as needed

# Define colors
fill_colors <- c(HGG = "#7570b3", LGG = "#1b9e77")
point_colors <- c(HGG = "black", LGG = "black")  # Darker shades for points; adjust if needed

# Plot
ggplot(df, aes(x = Group, y = Fraction, fill = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, size = 2, alpha = 0.6) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    title = paste0("PD1+GZB+CD8+ fraction of T cells (Mann-Whitney, p = ", pval, ")"),
    x = "Gleason status",
    y = "Cells in tumor"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")