


rm(list = setdiff(ls(), c("ldata", "df")))

load("df.RData")
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


# 0) Prep (same as before)
df$flag <- droplevels(factor(df$flag))
df2 <- subset(df, !is.na(flag) & is.finite(ICAT)) 
# If zeros are errors, uncomment:
# df2 <- subset(df2, ICAT > 0)

# Ensure the needed levels exist
need <- c("T","P","N")
if (!all(need %in% levels(df2$flag))) {
  stop("Missing one of required flag levels: T, P, N")
}

# 1) Pairwise Welch t-tests
tp <- wilcox.test(ICAT ~ flag, data = subset(df2, flag %in% c("T","P")))
pn <- wilcox.test(ICAT ~ flag, data = subset(df2, flag %in% c("P","N")))

cat(sprintf("T vs P: MW p = %.4g\n", tp$p.value))
cat(sprintf("P vs N: MW p = %.4g\n", pn$p.value))

# (Optional) adjust the two p-values (BH == same as none here, but shown for practice)


# 2) Three-group boxplot (same as before)
ymax <- max(df2$ICAT, na.rm = TRUE)
p_overall <- kruskal.test(ICAT ~ flag, data = df2)$p.value
p_lab <- paste0("Kruskal test p = ", format.pval(p_overall, digits = 5))

p <- ggplot(df2, aes(x = flag, y = ICAT, fill = flag)) +
  geom_boxplot(width = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "ICAT by flag", x = NULL, y = "ICAT") +
  annotate("text", x = 2, y = ymax * 1.05, label = p_lab, size = 4) +
  expand_limits(y = ymax * 1.1) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
print(p)