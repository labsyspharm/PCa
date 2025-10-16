
library(dplyr)

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

###Fig 3D

# Parse specimen ID after "@"
df$specimen <- sub(".*@", "", df$TLS)

# Keep only High/Low specimens and finite ICAT
df2 <- df %>%
  filter(specimen %in% c(high, low),
         is.finite(ICAT), ICAT < 80)

# EXCLUDE specimens ending in 05, 15, 27, or 21
#df2 <- df2 %>% filter(!grepl("(05|15|27|21)$", specimen))

# Group label and colors
df2 <- df2 %>%
  mutate(Group = if_else(specimen %in% high, "High", "Low"),
         Group = factor(Group, levels = c("Low", "High")))

# Compute 10th percentile per specimen and order
p10_tbl <- df2 %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarize(p10 = quantile(ICAT, 0.10, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(p10)

df2$specimen <- factor(df2$specimen, levels = p10_tbl$specimen)

# Plot
ggplot(df2, aes(x = specimen, y = ICAT, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1.8, color = "black") +
  scale_fill_manual(values = c(Low = "#1b9e77", High = "#7570b3")) +
  labs(title = "ICAT by specimen (ordered by 10th percentile; filtered suffixes)",
       x = "Specimen (after '@')", y = "ICAT") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 60, hjust = 1))
# If labels are crowded, add: + coord_flip()
