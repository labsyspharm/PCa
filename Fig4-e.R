PD1m <- data.frame(
  nB       = c(26,0,1,3,0,13,8,1,7,9,1,0,20,1,33,2,9,6,3,4,8,2,7,2,22,15,6,11,37),
  nT       = c(23,0,12,1,2,6,1,3,13,7,8,5,25,6,0,13,3,29,35,22,8,8,23,5,51,33,11,42,191),
  Group    = c("High","Other","Low","Low","Other","Low","Low","Low","Low","Low","High","Other",
               "High","High","High","Low","Low","High","Low","High","High","Low","Low","Low",
               "High","High","High","High","High"),
  row.names = c("LSP12601","LSP12603","LSP12605","LSP12607","LSP12609","LSP12611","LSP12613",
                "LSP12615","LSP12617","LSP12619","LSP12621","LSP12623","LSP12625","LSP12627",
                "LSP12629","LSP12631","LSP12633","LSP12635","LSP12637","LSP12639","LSP12641",
                "LSP12643","LSP12645","LSP12647","LSP12649","LSP12651","LSP12653","LSP12655","LSP12657")
)


NBH<- sum(subset(PD1m, PD1m$Group == "High")$nB) -33 # removing the total number of the LSP12629
NBL<- sum(subset(PD1m, PD1m$Group == "Low")$nB) 

##visual input
NBTH<- 93
NBTL<- 13


# Load required libraries
library(ggplot2)
library(Hmisc)  # For Wilson CI calculation

# Perform the proportion test with a confidence interval (95%)
prop_test <- prop.test(c(NBTH, NBTL), c(NBH, NBL), alternative = "greater", conf.level = 0.95)

# Compute proportions
prop_high <- NBTH / NBH
prop_low <- NBTL / NBL

# Compute 90% confidence intervals for each proportion using Wilson Score Interval
ci_high <- binconf(NBTH, NBH, alpha = 0.05, method = "wilson")  # 95% CI for High
ci_low <- binconf(NBTL, NBL, alpha = 0.05, method = "wilson")  # 95% CI for Low

# Create a data frame for visualization
df <- data.frame(
  Group = c("High", "Low"),
  Proportion = c(prop_high, prop_low),
  CI_lower = c(ci_high[2], ci_low[2]),  # Lower bounds of 90% CIs
  CI_upper = c(ci_high[3], ci_low[3]),  # Upper bounds of 90% CIs
  Successes = c(NBTH, NBTL),
  Totals = c(NBH, NBL)
)

# Define colors
colors <- c("High" = "#7570b3", "Low" = "#1b9e77")

# Generate the bar plot
ggplot(df, aes(x = Group, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.5) +  # Bar plot
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, color = "black") +  # Correct CI error bars
  scale_fill_manual(values = colors) +  # Assign colors to groups
  theme_minimal() +
  labs(
    x = "Group",
    y = "Proportion of BICs with TICs in their 0.5 mm vicinity",
    title = "Proportion Test Between HGG and LGG Groups (95% CI)"
  ) +
  # Annotate successes/total above each bar
  geom_text(aes(label = paste0(Successes, "/", Totals)), vjust = -0.5, size = 5) +
  # Annotate p-value
  annotate("text", x = 1.5, y = max(df$Proportion) + 0.05, 
           label = paste0("p-value = ", signif(prop_test$p.value, 3)), 
           color = ifelse(prop_test$p.value < 0.05, "red", "black"), hjust = 0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold title
  )




