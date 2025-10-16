


rm(list = setdiff(ls(), c("ldata", "df")))

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")



# Create data frame
dff <- data.frame(
  Var0<- df[, 10],
  Var1 = df[,2],  # Independent variable ICAT
  Var2 = df[,3], #ki-67
  Var3 = df[,4],# Dependent variable (before -log transformation) CD31
  Var4 = df[,5], #bKi67
  Var5 = df[,6], #T TCF1
  Var6 = df[,7], #T4 TCF1
  Var7 = df[,8],#T8 TCF1 
  Var8 = as.numeric(df[,9])# size of tls
)

# Assign group based on high/low
dff$Group <- ifelse(dff$Var0 %in% high, "High", 
                    ifelse(dff$Var0 %in% low, "Low", NA)) # NA if not in either

dff_res<- dff ##reserving the dff

dff<- subset(dff_res, dff_res$Var8 < 20000) ##removing the giant 2 BICs with more than 20000 cells 

# Perform -log transformation of Var2
dff$log_Var2 <- -log(dff$Var2)
dff$log_Var3 <- -log(dff$Var3)
dff$log_Var4 <- -log(dff$Var4*dff$Var8)
dff$log_Var4<- (-dff$log_Var4)
dff$log_Var5 <- log(dff$Var5*dff$Var8)
dff$log_Var6 <- -log(dff$Var6*dff$Var8)
dff$log_Var7 <- -log(dff$Var7*dff$Var8)


# Subset data for regression
dff_high <- subset(dff, Group == "High")
dff_low <- subset(dff, Group == "Low")

# Parse sample id after "@"
df$sample_id <- sub(".*@", "", df$TLS)

# Keep only High/Low samples and finite ICAT
df2 <- subset(df, sample_id %in% c(high, low) & is.finite(ICAT) & ICAT < 80)

# OPTIONAL: drop zeros if they are known errors
# df2 <- subset(df2, ICAT > 0)

# Group factor
df2$Group <- factor(ifelse(df2$sample_id %in% high, "High", "Low"),
                    levels = c("Low", "High"))

# Welch t-test (unequal variances)
# Split the data by group
g1 <- df2$ICAT[df2$Group == unique(df2$Group)[1]]
g2 <- df2$ICAT[df2$Group == unique(df2$Group)[2]]

# Normality tests
sh1 <- shapiro.test(g1)
sh2 <- shapiro.test(g2)

# Choose test depending on normality
if (sh1$p.value > 0.05 & sh2$p.value > 0.05) {
  test_result <- t.test(ICAT ~ Group, data = df2)
  cat("\nBoth groups normal so using t-test\n")
} else {
  test_result <- wilcox.test(ICAT ~ Group, data = df2)
  cat("\nAt least one group not normal so using Mann-Whitney U test\n")
}

print(test_result)
tt<- test_result


# Boxplot with p-value label
ymax <- max(df2$ICAT, na.rm = TRUE)
p_lab <- paste0("Welch t-test p = ", format.pval(tt$p.value, digits = 3))

ggplot(df2, aes(x = Group, y = ICAT, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, color = "black") +
  scale_fill_manual(values = c(Low = "#1b9e77", High = "#7570b3")) +
  labs(title = "ICAT by High/Low (parsed after '@')", x = NULL, y = "ICAT") +
  annotate("text", x = 1.5, y = ymax * 1.05, label = p_lab, size = 4) +
  expand_limits(y = ymax * 1.1) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
