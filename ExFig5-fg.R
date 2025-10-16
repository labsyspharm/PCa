

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



name<- c(low,high)

###f
df_scatter <- tibble(
  Name = name,
  B    = B.low_high.M + rnorm(26, 0,1),
  T    = T.low_high.M + rnorm(26, 0,1),
  Group = if_else(name %in% high, "High", "Low")
) %>%
  filter(!is.na(B) & !is.na(T))          # drop incomplete pairs

# 2) colour palette for High / Low
cols <- c(High = "#7570b3", Low = "#1b9e77")

# 3) scatter-plot with labels
ggplot(df_scatter, aes(x = B, y = T, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = Name),
                  size = 4, max.overlaps = 30,
                  segment.color = "grey50") +
  scale_color_manual(values = cols) +
  labs(
    title  = "Clustering intensity over 1.5 mm",
    x      = "B cells",
    y      = "T cells",
    color  = "Condition"
  ) +
  theme_minimal(base_size = 14)

#############

df_scatter <- tibble(
  Name = name,
  B    = 4*B.low_high.R + rnorm(26, 0,5),
  T    = 4*T.low_high.R + rnorm(26, 0,5),
  Group = if_else(name %in% high, "High", "Low")
) %>%
  filter(!is.na(B) & !is.na(T))          # drop incomplete pairs

# 2) colour palette for High / Low
cols <- c(High = "#7570b3", Low = "#1b9e77")

# 3) scatter-plot with labels
library(scales)  # for pretty_breaks if needed

ggplot(df_scatter, aes(x = B, y = T, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = Name),
                  size = 4, max.overlaps = 30,
                  segment.color = "grey50") +
  scale_color_manual(values = cols) +
  labs(
    title  = "Clustering radius over 1.5 mm",
    x      = "B cells",
    y      = "T cells",
    color  = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    breaks = pretty_breaks(),
    labels = function(x) x * 1.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(),
    labels = function(y) y * 1.5
  )

