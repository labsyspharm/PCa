# packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)


low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
         "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
          "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")

###saving the M+ and R50 for B and T cells for each of the samples then plot them in a bar plot!
B.low_high.M <- numeric(26)
T.low_high.M <- numeric(26)

B.low_high.R <- numeric(26)
T.low_high.R <- numeric(26)
c<-1
for (i in c(low,high)){
  example_data<- ldata[[i]]$B.CR500[1:50]
  example_data1<-example_data[example_data > 0]
  total_sum <- sum(example_data1)
  cumulative_sum <- cumsum(example_data1)
  half_sum <- total_sum / 2
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  B.low_high.M[c] <- mean(example_data[which(example_data > 0)])
  B.low_high.R[c] <- closest_index*10
  
  example_data<- ldata[[i]]$T.CR500[1:50]
  example_data1<-example_data[example_data > 0]
  total_sum <- sum(example_data1)
  cumulative_sum <- cumsum(example_data1)
  half_sum <- total_sum / 2
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  T.low_high.M[c] <- mean(example_data[which(example_data > 0)])
  T.low_high.R[c] <- closest_index*10
  
  c<- c+1
}


# ------------------------------------------------------------------
# Inputs assumed from your code: B.low_high.R, T.low_high.R (ECR)
# and B.low_high.M, T.low_high.M (CI-like metric),
# plus the ordered sample ids you looped over:
sample_ids <- c(low, high)   # length 26
stopifnot(length(sample_ids) == length(B.low_high.R),
          length(sample_ids) == length(T.low_high.R),
          length(sample_ids) == length(B.low_high.M),
          length(sample_ids) == length(T.low_high.M))

# helper for placing labels above the top points
top_n_by_group <- function(df, n = 4) {
  df %>% group_by(Cell) %>% slice_max(Value, n = n)
}

# ---------- Panel A: ECR (B vs T) ----------
df_ecr <- tibble(
  Name = rep(sample_ids, times = 2),
  Cell = rep(c("B", "T"), each = length(sample_ids)),
  Value = c(B.low_high.R, T.low_high.R)
) %>% drop_na()

# normality per group
sw_B_ecr <- shapiro.test(df_ecr$Value[df_ecr$Cell == "B"])
sw_T_ecr <- shapiro.test(df_ecr$Value[df_ecr$Cell == "T"])

if (sw_B_ecr$p.value > 0.05 && sw_T_ecr$p.value > 0.05) {
  test_ecr <- t.test(Value ~ Cell, data = df_ecr)
  method_ecr <- "t-test"
} else {
  test_ecr <- wilcox.test(Value ~ Cell, data = df_ecr)
  method_ecr <- "Wilcoxon"
}
p_ecr <- signif(test_ecr$p.value, 3)

top_ecr <- top_n_by_group(df_ecr, n = 4)

pA <- ggplot(df_ecr, aes(Cell, Value)) +
  geom_boxplot(fill = "grey85", width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Cell), width = 0.18, size = 2.6, alpha = 0.8) +
  geom_text_repel(data = top_ecr, aes(label = Name),
                  size = 3.8, min.segment.length = 0,
                  segment.color = "grey60") +
  scale_color_manual(values = c(B = "#c23b22", T = "#1b8c3a")) +
  labs(title = "ECR of B and T cells",
       subtitle = paste("P =", format(p_ecr, digits = 3), "(", method_ecr, ")"),
       x = "Cell type", y = "ECR (µm)") +
  scale_y_continuous(
    limits = c(min(df_ecr$Value, na.rm = TRUE),
               max(df_ecr$Value, na.rm = TRUE) * 1.1),
    labels = function(y) y * 6
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# ---------- Panel B: CI (B vs T) ----------
df_ci <- tibble(
  Name = rep(sample_ids, times = 2),
  Cell = rep(c("B", "T"), each = length(sample_ids)),
  Value = c(B.low_high.M, T.low_high.M)
) %>% drop_na()

sw_B_ci <- shapiro.test(df_ci$Value[df_ci$Cell == "B"])
sw_T_ci <- shapiro.test(df_ci$Value[df_ci$Cell == "T"])

if (sw_B_ci$p.value > 0.05 && sw_T_ci$p.value > 0.05) {
  test_ci <- t.test(Value ~ Cell, data = df_ci)
  method_ci <- "t-test"
} else {
  test_ci <- wilcox.test(Value ~ Cell, data = df_ci)
  method_ci <- "Wilcoxon"
}
p_ci <- signif(test_ci$p.value, 3)

top_ci <- top_n_by_group(df_ci, n = 4)

pB <- ggplot(df_ci, aes(Cell, Value)) +
  geom_boxplot(fill = "grey85", width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Cell), width = 0.18, size = 2.6, alpha = 0.8) +
  geom_text_repel(data = top_ci, aes(label = Name),
                  size = 3.8, min.segment.length = 0,
                  segment.color = "grey60") +
  scale_color_manual(values = c(B = "#c23b22", T = "#1b8c3a")) +
  labs(title = "CI of B and T cells",
       subtitle = paste("P =", format(p_ci, digits = 3), "(", method_ci, ")"),
       x = "Cell type", y = "CI") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# ---------- Combine 1 × 2 ----------
(pA | pB)