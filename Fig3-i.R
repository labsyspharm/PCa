




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

## Clean data (finite Var1/log_Var5 and non-missing Group)
dff_clean <- dff %>%
  filter(is.finite(Var1), is.finite(log_Var5), !is.na(Group),  Var1 < 90)

## Split groups (for modeling)
dff_high <- filter(dff_clean, Group == "High")
dff_low  <- filter(dff_clean,  Group == "Low")

## Fit simple linear models
fit5_high <- lm(log_Var5 ~ Var1, data = dff_high)
fit5_low  <- lm(log_Var5 ~ Var1, data = dff_low)

## Extract p-values from F-statistics
f_pval <- function(fit) {
  fs <- summary(fit)$fstatistic
  pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
}
p_value_high <- f_pval(fit5_high)
p_value_low  <- f_pval(fit5_low)

## R-squared
r2_high <- summary(fit5_high)$r.squared
r2_low  <- summary(fit5_low)$r.squared

## Annotation coordinates
x_annot  <- min(dff_clean$Var1, na.rm = TRUE) + 2
y_max    <- max(dff_clean$log_Var5, na.rm = TRUE)
y_annot1 <- y_max - 0.5
y_annot2 <- y_max - 1.0

## Plotting data (remove y = 0 only for display)
df_all_nz  <- filter(dff_clean, log_Var5 != 0)
df_high_nz <- filter(dff_high,  log_Var5 != 0)
df_low_nz  <- filter(dff_low,   log_Var5 != 0)

## Plot
ggplot(df_all_nz, aes(x = Var1, y = log_Var5)) +
  geom_point(aes(color = Group), size = 3, na.rm = TRUE) +
  geom_smooth(data = df_high_nz, aes(x = Var1, y = log_Var5),
              method = "lm", color = "#7570b3", se = TRUE, na.rm = TRUE) +
  geom_smooth(data = df_low_nz,  aes(x = Var1, y = log_Var5),
              method = "lm", color = "#1b9e77", se = TRUE, na.rm = TRUE) +
  scale_color_manual(values = c(High = "#7570b3", Low = "#1b9e77")) +
  labs(
    x = "ICAT",
    y = "log(CD3d+TCF1+ prop)",
    title = "CD3d+TCF1+ and ICAT Regressions"
  ) +
  annotate("text", x = x_annot, y = y_annot1,
           label = paste0("p(High) = ", signif(p_value_high, 3),
                          ", R² = ", signif(r2_high, 3)),
           color = "#7570b3", hjust = 0) +
  annotate("text", x = x_annot, y = y_annot2,
           label = paste0("p(Low) = ", signif(p_value_low, 3),
                          ", R² = ", signif(r2_low, 3)),
           color = "#1b9e77", hjust = 0) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")
##stratify with Gleason! and colour for patients!

