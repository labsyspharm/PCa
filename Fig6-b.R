
#######granzymbeB gate


rm(list = setdiff(ls(), c("ldata", "df")))

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

high_granzymeB_gate<- c(7, 7, 7, 7.04, 6.8, 6.8, 6.84, 6.95, 6.954, 6.87, 6.88, 6.95, 6.94)
low_granzymeB_gate<- c(6.98, 7.14, 7.13, 7.05, 6.98, 6.86, 6.93, 7.3, 7.08, 6.9, 7.5, 7, 6.88)

#####run it for one sample!
na<-"LSP12601"
lp<-ldata[[na]]
####interactions

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 0. PACKAGES
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library(RANN)             # nn2
library(dbscan)           # frNN
library(spatstat.geom)    # ppp
library(spatstat.explore) # Kcross
library(ggplot2)
library(dplyr)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1. SETUP & FILTERING
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# assume lp is already in memory
N <- nrow(lp)
message("Data has ", N, " cells × ", ncol(lp), " features")

## identify tumor & CTLs
tumor_mask <- bitwAnd(lp$cflag, 8) == 8
ctl_mask   <- lp$phen_vec == "CD8 T cells" & log(lp$GranzymeB) > 6.95

## remove CTLs from tumor set
tumor_mask <- tumor_mask & !ctl_mask

tumor_xy <- as.matrix(lp[tumor_mask, c("x","y")])
ctl_xy   <- as.matrix(lp[   ctl_mask, c("x","y")])

message("Tumor: ", nrow(tumor_xy), " cells; CTLs: ", nrow(ctl_xy))

# define analysis window
xrange <- range(lp$x); yrange <- range(lp$y)
win <- owin(xrange, yrange)

#[CTL per tumors; CTL per CD8+ cells]

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 2. DISTANCE TO NEAREST TUMOR CELL
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
nn <- nn2(tumor_xy, ctl_xy, k=1)
dist_nearest <- nn$nn.dists[,1]


#[first quantile, median, mean, and third quantile of the dist_nearest]
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3. PROPORTION WITHIN RADIUS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
radii <- c(10,20,50)
prop_within <- sapply(radii, function(r) mean(dist_nearest < r))
df_prop <- data.frame(radius=factor(radii), proportion=prop_within)

#[fraction of the CTLs in 10, 20 and 50 microns] 

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4. LOCAL TUMOR DENSITY (cells per µm²)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
radius <- 50
n_ctl  <- nrow(ctl_xy)
counts <- integer(n_ctl)

# process in chunks to avoid single huge query
chunk <- 10000
for (start in seq(1, n_ctl, by=chunk)) {
  end     <- min(start+chunk-1, n_ctl)
  sub_q   <- ctl_xy[start:end, , drop=FALSE]
  fr      <- frNN(tumor_xy, eps=radius, query=sub_q) 
  counts[start:end] <- lengths(fr$id)
}

local_density <- counts / (pi * radius^2)

# [1st quantile, median, mean, and 3rd quantile of the local_density]


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 5. RIPLEY’S CROSS‐K FUNCTION
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# build a single multitype pattern
all_xy <- rbind(ctl_xy, tumor_xy)
marks  <- factor(c(rep("CTL", nrow(ctl_xy)), rep("Tumor", nrow(tumor_xy))))
pp_multi <- ppp(all_xy[,1], all_xy[,2], window=win, marks=marks)

# cross‐K from CTL → Tumor
Kcross_res <- Lcross(pp_multi, i="CTL", j="Tumor",
                     r = seq(0, 100, by=5))


#[average Kcross_res$border-Kcross_res$theo over radii]
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 6. DIRECT CONTACT FRACTION
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
contact_cut   <- 2  # µm
frac_contact  <- mean(dist_nearest < contact_cut)
ctl_df <- data.frame(
  x       = ctl_xy[,1],
  y       = ctl_xy[,2],
  contact = dist_nearest < contact_cut
)

# 1) Prepare each subset with a class label
bg_df    <- data.frame(x=lp$x, y=lp$y) %>%    mutate(class="All cells")
tumor_df <- subset(lp, bitwAnd(cflag,8)==8)[,c("x","y")] %>% mutate(class="Tumor")
cd8_df   <- subset(lp, phen_vec=="CD8 T cells")[,c("x","y")] %>% mutate(class="CD8+ T cells")
ctl_not  <- ctl_df %>% filter(!contact) %>%       mutate(class="CTL (no contact)")
ctl_yes  <- ctl_df %>% filter(contact)  %>%       mutate(class="CTL (in contact)")

# 2) Combine
all_df <- bind_rows(bg_df, tumor_df, cd8_df, ctl_not, ctl_yes)

# 4) Plot with single legend (color) and override key size

# … assume bg_df, tumor_df, cd8_df, ctl_df, contact_cut, frac_contact, colors, shapes, sizes, alphas, na are all defined …

# 1) Prepare each subset with an on-the-fly “class” aesthetic
bg_df    <- data.frame(x=lp$x, y=lp$y)
tumor_df <- subset(lp, bitwAnd(cflag,8)==8)[,c("x","y")]
cd8_df   <- subset(lp, phen_vec=="CD8 T cells")[,c("x","y")]
ctl_not  <- filter(ctl_df, !contact)
ctl_yes  <- filter(ctl_df,  contact)

ggplot() +
  geom_point(data=bg_df,    aes(x,y, color="All cells"),
             shape=46, size=0.5, alpha=0.4) +
  geom_point(data=tumor_df, aes(x,y, color="Tumor"),
             shape=46, size=0.5, alpha=0.4) +
  geom_point(data=cd8_df,   aes(x,y, color="CD8+ T cells"),
             shape=46, size=0.5, alpha=0.9) +
  stat_density_2d(data=ctl_df, aes(x=x,y=y),
                  colour="black", size=0.8, alpha=0.5,
                  bins=10, geom="density_2d") +
  geom_point(data = filter(ctl_df, !contact),
             aes(x,y, color="CTL (no contact)"),
             shape=19, size=0.7, alpha=0.8) +
  # <— this is the updated layer:
  geom_point(data = filter(ctl_df, contact),
             aes(x, y, color="CTL (in contact)"),
             shape=21, fill="yellow", color="black",
             size=5.0, alpha=0.9) +
  
  scale_color_manual(
    name   = "Class",
    values = c(
      "All cells"         = "grey80",
      "Tumor"             = "lightblue",
      "CD8+ T cells"      = "darkgreen",
      "CTL (no contact)"  = "tomato",
      "CTL (in contact)"  = "yellow"
    ),
    guide  = guide_legend(override.aes = list(size=3, shape=19, alpha=1))
  ) +
  labs(
    title = sprintf(
      "CTLs in Direct Contact (<%g µm): %.1f%% of CTLs in %s",
      contact_cut, 100 * frac_contact, na
    ),
    x = "X position (µm)",
    y = "Y position (µm)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title      = element_text(size=14, face="bold")
  )

