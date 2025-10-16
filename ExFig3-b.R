# Required libraries (from the file)
library(ggplot2)
library(bitops)  # For bitAnd on KiGate
rm(list = setdiff(ls(), c("ldata", "df")))

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657") ##also removing the hugely infiltrated LSP12629
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")


# Calculate TLS counts per specimen for high group (needed for per-BIC loop)
h.ntls <- numeric(length(high))
c <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  h.ntls[c] <- length(unique(l$tls_id))
  c <- c + 1
}
h.ntls <- h.ntls - 1  # Subtract non-TLS (ID 0)

# Calculate TLS counts per specimen for low group
l.ntls <- numeric(length(low))
c <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  l.ntls[c] <- length(unique(l$tls_id))
  c <- c + 1
}
l.ntls <- l.ntls - 1  # Subtract non-TLS (ID 0)

# Calculate % Ki67+ per BIC for high and low groups
h.ki <- numeric(sum(h.ntls))
h.spec <- rep(high, h.ntls)  # Track specimen for coloring
cc <- 1
for (i in high) {
  l <- ldata[[which(names(ldata) == i)]]
  for (j in unique(l$tls_id)[-which(unique(l$tls_id) == 0)]) {
    s <- subset(l, l$tls_id == j)
    if (nrow(s) > 0) {
      # Assume Ki67 positivity is bitAnd(KiGate, 1) == 1 (adjust bit based on gating)
      h.ki[cc] <- sum(s$KiGate == 1) / nrow(s) * 100
    } else {
      h.ki[cc] <- NA
    }
    cc <- cc + 1
  }
}

l.ki <- numeric(sum(l.ntls))
l.spec <- rep(low, l.ntls)  # Track specimen for coloring
cc <- 1
for (i in low) {
  l <- ldata[[which(names(ldata) == i)]]
  for (j in unique(l$tls_id)[-which(unique(l$tls_id) == 0)]) {
    s <- subset(l, l$tls_id == j)
    if (nrow(s) > 0) {
      l.ki[cc] <- sum(s$KiGate == 1) / nrow(s) * 100
    } else {
      l.ki[cc] <- NA
    }
    cc <- cc + 1
  }
}

# Create data frame for % Ki67+
gleason <- rep(c("HGG", "LGG"), times = c(length(h.ki), length(l.ki)))
value <- c(h.ki, l.ki)
specimen <- c(h.spec, l.spec)
data <- data.frame(Gleason = gleason, Value = value, Specimen = specimen)
data <- na.omit(data)  # Remove any NA BICs

# Random colors per specimen (rainbow palette)
unique_specs <- unique(data$Specimen)
spec_colors <- rainbow(length(unique_specs))
names(spec_colors) <- unique_specs

# Check normality with Kolmogorov-Smirnov test for each group
hgg_values <- data$Value[data$Gleason == "HGG"]
lgg_values <- data$Value[data$Gleason == "LGG"]

ks_hgg <- ks.test(hgg_values, "pnorm", mean = mean(hgg_values), sd = sd(hgg_values))
ks_lgg <- ks.test(lgg_values, "pnorm", mean = mean(lgg_values), sd = sd(lgg_values))

# Perform appropriate test: t-test if both groups are normal, else Mann-Whitney U
if (ks_hgg$p.value > 0.05 && ks_lgg$p.value > 0.05) {
  test_result <- t.test(Value ~ Gleason, data = data)
} else {
  test_result <- wilcox.test(Value ~ Gleason, data = data)
}

# Define y-position for significance bar (near max, adjusted to match image ~25)
bar_y <- max(data$Value) * 0.9

# Create boxplot with colored jitter
p <- ggplot(data, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_boxplot(width = 0.5, alpha = 0.5) +
  geom_jitter(aes(color = Specimen), width = 0.2, size = 2, alpha = 0.8) +  # Colored by specimen
  scale_color_manual(values = spec_colors) +  # Random/rainbow colors
  scale_fill_manual(values = c("HGG" = "#7570b3", "LGG" = "#1b9e77")) +
  labs(title = "Proliferation in BICs", y = "% Ki67+ cells per BIC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") +  # No legend for specimens (as in image)
  ylim(0, 30) +  # Slight extension
  # Add significance bar
  geom_segment(x = 1.1, xend = 1.9, y = 30, yend = 30, linetype = "dashed", color = "black") +
  # Add p-value text
  annotate("text", x = 1.5, y = 30, label = paste("P =", formatC(test_result$p.value, digits = 4, format = "f")),
           hjust = 0.5, vjust = -1, color = ifelse(test_result$p.value < 0.05, "red", "black"))

# Display the plot
print(p)