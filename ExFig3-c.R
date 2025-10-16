
# Extended Fig 3C

#high ICAT and high Gleason Q3
tlf<- subset(df, df$count < 20000) 
i<-1
speci<- strsplit(tlf$TLS[i], "@")[[1]][2]
ID<- strsplit(tlf$TLS[i], "@")[[1]][1]
da<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]

for (i in 2:dim(tlf)[1]){
  speci<- strsplit(tlf$TLS[i], "@")[[1]][2]
  ID<- strsplit(tlf$TLS[i], "@")[[1]][1]
  d<- subset(ldata[[speci]], ldata[[speci]]$tls_id == ID)[,c(3,4,5,6,36,129:135)]
  da<- rbind(da, d)
}
# install.packages("ggvenn")  # if you haven’t already
library(ggvenn)

# build your sets as named lists of row-indices
sets <- list(
  `T cells`   = which(da$coarse_phen_vec == "T cells" | bitAnd(da$pflag, 1024) == 1024 | bitAnd(da$pflag, 4096) == 4096 | bitAnd(da$pflag, 2048) == 2048),
  `TCF1gate` = which(da$TCF1gate   == TRUE),
  `KiGate`    = which(da$KiGate      == 1),
  `B cells`   = which(da$coarse_phen_vec == "B cells" & !bitAnd(da$pflag, 1024) == 1024 & !bitAnd(da$pflag, 4096) == 4096 & !bitAnd(da$pflag, 2048) == 2048)
)

# now plot with percentages

ggvenn(
  sets,
  fill_color      = c("seagreen","orchid","skyblue","tomato"),
  stroke_size     = 0.5,
  set_name_size   = 5,
  show_percentage = TRUE,
  text_size       = 4
) +
  theme_void() +
  theme(
    plot.title     = element_text(face="bold", size=14),
    plot.subtitle  = element_text(size=10),
    legend.position = "right"
  ) +
  labs(
    title    = "Overlap of T vs B cells and KiGate vs TCF1gate",
    subtitle = "Numbers in each region are % of total cells"
  )



# build your sets as named lists of row-indices
sets <- list(
  `T cells`   = which(da$coarse_phen_vec == "T cells" | bitAnd(da$pflag, 1024) == 1024 | bitAnd(da$pflag, 4096) == 4096 | bitAnd(da$pflag, 2048) == 2048),
  `TCF1gate` = which(da$TCF1gate   == TRUE)
)

# now plot with percentages

ggvenn(
  sets,
  fill_color      = c("seagreen","orchid","skyblue","tomato"),
  stroke_size     = 0.5,
  set_name_size   = 5,
  show_percentage = TRUE,
  text_size       = 4
) +
  theme_void() +
  theme(
    plot.title     = element_text(face="bold", size=14),
    plot.subtitle  = element_text(size=10),
    legend.position = "right"
  ) +
  labs(
    title    = "Overlap of T vs B cells and KiGate vs TCF1gate",
    subtitle = "Numbers in each region are % of total cells"
  )



# build your sets as named lists of row-indices
sets <- list(
  `KiGate`    = which(da$KiGate      == 1),
  `B cells`   = which(da$coarse_phen_vec == "B cells" & !bitAnd(da$pflag, 1024) == 1024 & !bitAnd(da$pflag, 4096) == 4096 & !bitAnd(da$pflag, 2048) == 2048)
)

# now plot with percentages

ggvenn(
  sets,
  fill_color      = c("skyblue","tomato"),
  stroke_size     = 0.5,
  set_name_size   = 5,
  show_percentage = TRUE,
  text_size       = 4
) +
  theme_void() +
  theme(
    plot.title     = element_text(face="bold", size=14),
    plot.subtitle  = element_text(size=10),
    legend.position = "right"
  ) +
  labs(
    title    = "Overlap of T vs B cells and KiGate vs TCF1gate",
    subtitle = "Numbers in each region are % of total cells"
  )


