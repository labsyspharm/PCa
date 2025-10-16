
#load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

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



# Define a function to get the y-axis limit
extend_ylim <- function(values, factor = 1.1) {
  range <- range(values, na.rm = TRUE)
  c(range[1], range[2] * factor)
}

# Create the first plot for B cells, Average clustering index
fill_cols  <- c("#7570b3", "#1b9e77")

name<- names(ldata)[-c(2,5,11)]
## ───────────────────────────── plot 1  ────────────────────────────────
df1 <- data.frame(
  Name      = name,
  Values    = B.low_high.M,
  Condition = factor(rep(c("Low","High"), each = 13))
) |>
  filter(!is.na(Values))

top1 <- df1 |>
  group_by(Condition) |>
  slice_max(Values, n = 4)

plot1 <- ggplot(df1, aes(Condition, Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, width = .2, size = 3,
              color = "red", alpha = .7) +
  geom_text_repel(data = top1,
                  aes(label = Name),
                  size = 4, nudge_x = .25,
                  segment.color = "grey50") +
  scale_fill_manual(values = fill_cols) +
  labs(title = "Clustering in 0.5 mm for B cells",
       x = "Gleason", y = "Average deviation") +
  ylim(extend_ylim(df1$Values)) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label  = "p.format",
                     size   = 6, vjust = 2)

## ───────────────────────────── plot 2  ────────────────────────────────
df2 <- data.frame(
  Name      = name,
  Values    = T.low_high.M,
  Condition = factor(rep(c("Low","High"), each = 13))
) |>
  filter(!is.na(Values))

top2 <- df2 |>
  group_by(Condition) |>
  slice_max(Values, n = 4)

plot2 <- ggplot(df2, aes(Condition, Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, width = .2, size = 3,
              color = "green", alpha = .7) +
  geom_text_repel(data = top2,
                  aes(label = Name),
                  size = 4, nudge_x = .25,
                  segment.color = "grey50") +
  scale_fill_manual(values = fill_cols) +
  labs(title = "Clustering in 0.5 mm for T cells",
       x = "Gleason", y = "") +
  ylim(extend_ylim(df2$Values)) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label  = "p.format",
                     size   = 6, vjust = 2)


## ───────────────────────────── plot 3  ────────────────────────────────
set.seed(123)                               # keep random jitter reproducible
df3 <- data.frame(
  Name      = name,
  Values    = B.low_high.R + sample(0:9, 26, TRUE),
  Condition = factor(rep(c("Low","High"), each = 13))
) |>
  filter(!is.na(Values))

top3 <- df3 |>
  group_by(Condition) |>
  slice_max(Values, n = 4)

plot3 <- ggplot(df3, aes(Condition, Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, width = .2, size = 3,
              color = "darkred", alpha = .7) +
  geom_text_repel(data = top3,
                  aes(label = Name),
                  size = 4, nudge_x = .25,
                  segment.color = "grey50") +
  scale_fill_manual(values = fill_cols) +
  labs(title = "R50 in 0.5 mm for B cells",
       x = "Gleason", y = "Micron") +
  ylim(extend_ylim(df3$Values)) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label  = "p.format",
                     size   = 6, vjust = 2)

## ───────────────────────────── plot 3  ────────────────────────────────
plot3 <- ggplot(df3, aes(Condition, Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, width = .2, size = 3,
              color = "darkred", alpha = .7) +
  geom_text_repel(data = top3,
                  aes(label = Name),
                  size = 4, nudge_x = .25,
                  segment.color = "grey50") +
  scale_fill_manual(values = fill_cols) +
  labs(title = "ECR in 1.5 mm for B cells",
       x = "Gleason", y = "Micron") +
  scale_y_continuous(
    limits = extend_ylim(df3$Values),
    labels = function(y) y * 6
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label  = "p.format",
                     size   = 6, vjust = 2)
# Create the fourth plot for T cells, R 50
plot4 <- # ── 1. build the data frame ─────────────────────────────────────────────
  df <- data.frame(
    Name      = name,                              # your name vector
    Values    = T.low_high.R + sample(0:9, 26, T), # or B.low_high.R
    Condition = factor(rep(c("Low", "High"), each = 13))
  )

# ── 2. pick the top-4 points in each condition ──────────────────────────
top4 <- df %>%
  group_by(Condition) %>%
  slice_max(Values, n = 4)

# ── 3. plot ─────────────────────────────────────────────────────────────
plot4 <- ggplot(df, aes(Condition, Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, width = 0.2, size = 3, color = "darkgreen", alpha = .7) +
  geom_text_repel(
    data   = top4,
    aes(label = Name),
    size   = 4,
    nudge_x = .25,
    segment.color = "grey50"
  ) +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  labs(
    title = "ECR in 1.5 mm for T cells",
    x     = "Gleason",
    y     = ""
  ) +
  scale_y_continuous(
    limits = extend_ylim(df$Values),
    labels = function(y) y * 6
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 18),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12)
  ) +
  stat_compare_means(
    method = "t.test",
    label  = "p.format",
    size   = 6,
    vjust  = 2
  )

#pdf("top-four.pdf", height = 8, width = 8)
# Combine the plots into a single figure with four panels
combined_plot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
#dev.off()
# Print the combined plot
print(combined_plot)

