# Required libraries (from the file)
library(ggplot2)
library(bitops)  # For bitAnd

# Assume ldata is loaded as in the file

#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki.RData")

# Load sheet and define high/low groups (exact from file, with removals)
sheet <- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")
high <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "High")

high <- high[-c(6)]
low <- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow == "Low")
low <- low[-c(10,13,14)]

# high
# [1] "LSP12601" "LSP12621" "LSP12623" "LSP12625" "LSP12627" "LSP12635" "LSP12639" "LSP12641" "LSP12649" "LSP12651" "LSP12653" "LSP12655" "LSP12657"
# > low
# [1] "LSP12603" "LSP12605" "LSP12607" "LSP12609" "LSP12611" "LSP12613" "LSP12615" "LSP12617" "LSP12619" "LSP12633" "LSP12637" "LSP12647"

# Calculate % B cells and % T cells in tumor for high group
h_percent_B <- numeric(length(high))
h_percent_T <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[i]]  # Use direct name indexing as in file
  tumor <- subset(l, bitAnd(cflag, 8) == 8)
  if (nrow(tumor) > 0) {  # Guard against empty tumor subsets
    h_percent_B[c] <- sum(tumor$coarse_phen_vec == "B cells") / nrow(tumor) * 100
    h_percent_T[c] <- sum(tumor$coarse_phen_vec == "T cells") / nrow(tumor) * 100
  } else {
    h_percent_B[c] <- NA
    h_percent_T[c] <- NA
  }
  c <- c + 1
}

# Calculate % B cells and % T cells in tumor for low group
l_percent_B <- numeric(length(low))
l_percent_T <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[i]]
  tumor <- subset(l, bitAnd(cflag, 8) == 8)
  if (nrow(tumor) > 0) {
    l_percent_B[c] <- sum(tumor$coarse_phen_vec == "B cells") / nrow(tumor) * 100
    l_percent_T[c] <- sum(tumor$coarse_phen_vec == "T cells") / nrow(tumor) * 100
  } else {
    l_percent_B[c] <- NA
    l_percent_T[c] <- NA
  }
  c <- c + 1
}

# Create data frame, removing any NA rows (similar to filtering in dff)
data <- data.frame(
  percent_B = c(h_percent_B, l_percent_B),
  percent_T = c(h_percent_T, l_percent_T),
  Group = rep(c("HGG", "LGG"), c(length(high), length(low)))
)
data <- na.omit(data)  # Remove specimens with no tumor data

# Perform linear regressions for each group (exact pattern from fit_high/fit_low in file)
fit_HGG <- lm(percent_T ~ percent_B, data = data[data$Group == "HGG", ])
fit_LGG <- lm(percent_T ~ percent_B, data = data[data$Group == "LGG", ])

# Extract p-values for slopes (exact pattern from summary(fit)$fstatistic and pf in file)
p_value_HGG <- summary(fit_HGG)$coefficients[2, 4]  # p for slope
p_value_LGG <- summary(fit_LGG)$coefficients[2, 4]

# Plot with matching aesthetics: scatter, regression lines, 95% CI shading, colors, annotations
# (Pattern from violin plot and end-of-file boxplots, with geom_smooth for lm + CI)
ggplot(data, aes(x = percent_B, y = percent_T, color = Group)) +
  geom_point(size = 2, alpha = 0.8) +  # Points for specimens
  geom_smooth(method = "lm", fill = "grey70", alpha = 0.3, level = 0.95) +  # Regression + gray 95% CI
  scale_color_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77"),
                     labels = c("HGG specimen", "LGG specimen")) +
  labs(title = "B and T cells regression",
       x = "% B cells on the tumor",
       y = "% T cells on the tumor") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  annotate("text", x = max(data$percent_B) * 0.7, y = max(data$percent_T) * 0.9,
           label = paste("HGG P =", formatC(p_value_HGG, digits = 4, format = "f")),
           color = "#7570b3") +
  annotate("text", x = max(data$percent_B) * 0.7, y = max(data$percent_T) * 0.8,
           label = paste("LGG P =", formatC(p_value_LGG, digits = 4, format = "f")),
           color = "#1b9e77") +
  annotate("text", x = max(data$percent_B) * 0.9, y = min(data$percent_T) + 1,
           label = "95% CI", color = "grey50", angle = 30, size = 3.5)