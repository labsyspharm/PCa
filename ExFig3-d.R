
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

## Start from dff_res; drop rows where Var4==0 or Var5==0 (as intended)

dff <- subset(dff, dff$Var1 < 90) %>% filter(!(Var4 == 0 | Var5 == 0))
q1 <- quantile(dff$Var1, 0.25, na.rm = TRUE)
dff_clean <- dff %>%
  filter(is.finite(Var1), is.finite(Var5), is.finite(Var4), !is.na(Group),  Var1 < 90)
stars_df <- dff_clean %>% filter(Var1 <= q1)

## Clean data: finite logs and non-missing Group; build x,y once
dff_clean <- dff %>%
  filter(is.finite(Var4), is.finite(Var5), is.finite(Var8), !is.na(Group)) %>%
  mutate(
    x = log(Var4 * Var8),
    y = log(Var5 * Var8)
  ) %>%
  filter(is.finite(x), is.finite(y))

## Split groups for modeling and lines
dff_high <- filter(dff_clean, Group == "High")
dff_low  <- filter(dff_clean,  Group == "Low")

## Fits: y ~ x within each group
fitBkvsTTCF1_high <- lm(y ~ x, data = dff_high)
fitBkvsTTCF1_low  <- lm(y ~ x, data = dff_low)

## p-values from F-statistics
f_pval <- function(fit) { fs <- summary(fit)$fstatistic; pf(fs[1], fs[2], fs[3], lower.tail = FALSE) }
p_high <- f_pval(fitBkvsTTCF1_high)
p_low  <- f_pval(fitBkvsTTCF1_low)

## R-squared
r2_high <- summary(fitBkvsTTCF1_high)$r.squared
r2_low  <- summary(fitBkvsTTCF1_low)$r.squared

## Annotation coordinates
x_annot  <- min(dff_clean$x, na.rm = TRUE) + 2
y_max    <- max(dff_clean$y, na.rm = TRUE)
y_annot1 <- y_max - 0.5
y_annot2 <- y_max - 1.0

## Plot
ggplot(dff_clean, aes(x = x, y = y, color = Group)) +
  # base points
  geom_point(size = 3) +
  
  # asterisk markers (expects stars_df with Var4, Var5, Var8, Group)
  geom_point(
    data  = stars_df,
    aes(x = log(Var4 * Var8), y = log(Var5 * Var8), color = Group),
    shape = 8, size = 4, inherit.aes = FALSE
  ) +
  
  # regression lines by group
  geom_smooth(data = dff_high, method = "lm", se = TRUE, color = "#7570b3") +
  geom_smooth(data = dff_low,  method = "lm", se = TRUE, color = "#1b9e77") +
  
  scale_color_manual(values = c(High = "#7570b3", Low = "#1b9e77")) +
  labs(
    x = "log(CD20+ Ki67+ cells)",      # x = log(Var4*Var8)
    y = "log(CD3d+ TCF1+ cells)",      # y = log(Var5*Var8)
    title = "Progenitor T cells vs Proliferating B cells with identified Q1 ICAT"
  ) +
  annotate("text", x = x_annot, y = y_annot1,
           label = sprintf("p(High)=%.3g, R²=%.3g", p_high, r2_high),
           color = "#7570b3", hjust = 0) +
  annotate("text", x = x_annot, y = y_annot2,
           label = sprintf("p(Low)=%.3g, R²=%.3g", p_low, r2_low),
           color = "#1b9e77", hjust = 0) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

