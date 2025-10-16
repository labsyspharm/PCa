SpaCont <- function(pID, res=0.9){
  d <- ldata[[pID]]
  cont <- d[, c(3, 4, 5, 6, 13, 18, 129)]
  
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(bitops)
  
  # Define High and Low groups
  low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", 
           "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")
  
  high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", 
            "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
  
  # Set title color based on pID category
  title_color <- ifelse(pID %in% low, "#1b9e77", ifelse(pID %in% high, "#7570b3", "black"))
  
  # Sample the data (without replacement)
  set.seed(123)
  cont_sample <- cont %>%
    sample_n(min(nrow(cont), 200000 * res), replace = FALSE)
  
  # Create subsets based on conditions
  cont_background <- subset(cont_sample, bitAnd(cflag, 8) != 8)  # Background (Stromal)
  cont_highlight <- subset(cont_sample, bitAnd(cflag, 8) == 8)  # Tumor overlay
  cont_t8cells <- subset(cont_sample, phen_vec == "CD8 T cells")  # CD8 T cells
  cont_bcells <- subset(cont_sample, phen_vec == "B cells")  # B cells
  cont_t4cells <- subset(cont_sample, phen_vec == "CD4 T cells")  # CD4 T cells
  cont_CD31 <- subset(cont_sample, bitAnd(pflag, 4194304) == 4194304)  # CD31+ cells
  
  # Add a label column to assign legend categories
  cont_background$Type <- "Stromal"
  cont_highlight$Type <- "Tumor"
  cont_t8cells$Type <- "CD8 T cells"
  cont_t4cells$Type <- "CD4 T cells"
  cont_bcells$Type <- "B cells"
  cont_CD31$Type <- "CD31+"
  
  # Flip the y-coordinates to adjust to image space
  y_min <- min(cont_sample$y)
  y_max <- max(cont_sample$y)
  
  flip_y <- function(df) {
    df$y <- y_max - df$y
    return(df)
  }
  
  cont_background <- flip_y(cont_background)
  cont_highlight <- flip_y(cont_highlight)
  cont_t8cells <- flip_y(cont_t8cells)
  cont_t4cells <- flip_y(cont_t4cells)
  cont_bcells <- flip_y(cont_bcells)
  cont_CD31 <- flip_y(cont_CD31)
  
  # Combine data
  cont_all <- bind_rows(cont_background, cont_highlight)
  
  # Create the plot with legend
  contour_plot <- ggplot() +
    # Plot background points (Stromal and Tumor) with transparency
    geom_point(data = cont_all, aes(x = x, y = y, color = Type), alpha = 0.3, size = 0.3) +
    #geom_point(data = cont_CD31, aes(x = x, y = y, color = Type), alpha = 0.5, size = 0.5) +
    
    # Overlay finer contour lines for cell types
    geom_density_2d(data = cont_bcells, aes(x = x, y = y, linetype = "B cells"), color = "red", bins = 20, h = c(1000, 1000), alpha = 0.65) +
    geom_density_2d(data = cont_t8cells, aes(x = x, y = y, linetype = "CD8 T cells"), color = "green", bins = 20, h = c(1000, 1000), alpha = 0.9) +
    geom_density_2d(data = cont_t4cells, aes(x = x, y = y, linetype = "CD4 T cells"), color = "blue", bins = 20, h = c(1000, 1000), alpha = 0.6) +
    
    # Manually set colors for categories
    scale_color_manual(
      values = c("Stromal" = "lightgrey", "Tumor" = "grey57", "CD31+" = "maroon4"),
      name = "Cell Type"
    ) +
    
    # Adjust legend for line types
    scale_linetype_manual(
      values = c("B cells" = "solid", "CD4 T cells" = "solid", "CD8 T cells" = "solid"),
      name = "Contour Type"
    ) +
    
    # Set y-axis limits
    scale_y_continuous(limits = c(0, y_max)) +
    
    # Theme and labels
    theme_minimal() +
    labs(
      title = paste(pID),
      x = "x",
      y = "y"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = title_color),  # Set dynamic color
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "right"  # Place the legend to the right
    )
  
  # Display the plot
  print(contour_plot)
}
SpaCont("LSP12651")
SpaCont("LSP12619")
