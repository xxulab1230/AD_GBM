#This R script includes a custom function for creating a modified violin plot, referred to as geom_flat_violin, 
#and uses it in conjunction with other ggplot2 elements to create a comparative plot of gene expression data. 

# Load necessary libraries for data manipulation and visualization
library(ggplot2)
library(dplyr)
library(ggpubr)

# Define a custom operator to provide a fallback value if the first is NULL
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

# Define a new geom for creating flat violin plots
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

# Define the GeomFlatViolin as an extension of the Geom class
GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
                          # Define the setup for the data used in the geom
                          setup_data = function(data, params) {
                            # Calculate the width of the violin and define the bounding box for each group
                            # ...
                          },
                          
                          # Define how to draw a group of data in the plot
                          draw_group = function(data, panel_scales, coord) {
                            # Calculate points for the outline of the violin and sort them
                            # Close the polygon to ensure it's a complete shape
                            # Use the draw_panel method from GeomPolygon to draw the shape
                            # ...
                          },
                          
                          # Define how to draw the legend key for this geom
                          draw_key = draw_key_polygon,
                          
                          # Define default aesthetics for the geom
                          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                                            alpha = NA, linetype = "solid"),
                          
                          # Define which aesthetics are required by the geom
                          required_aes = c("x", "y")
)

# Load or define the gene expression data (expr_ad for ADNI, datExpr0 for ROSMAP)
data = expr_ad[,genelist]
datExpr0 = data.frame(datExpr0)
expr_rush <- datExpr0[,colnames(datExpr0) %in% genelist]

# Ensure the expression data is numeric
expr_rush[] <- lapply(expr_rush, as.numeric)
#expr_rush$mean <- rowMeans(expr_rush)

# Assign a group variable based on diagnosis data
expr_rush$Group <- clin_mr$DIAGNOSIS 

# Convert the data to a tidy format for plotting
expr_gather <- tidyr::gather(data, key = Genenames, value = Geneexpr, -Group)

# Create a PDF file to save the plot
pdf("./ADNI-4gene-compare.pdf", height = 4, width = 10)

# Define the comparisons to be shown in the plot
my_comparisons <- list(c("Dementia", "CN"), c("MCI", "CN"), c("Dementia", "MCI")) # Example for ADNI
my_comparisons <- list(c("AD", "CONTROL"), c("MCI", "CONTROL"), c("AD", "MCI")) # Example for ROSMAP

# Create the plot using ggplot2 with the custom flat violin geom and other elements
ggplot(expr_gather, aes(x=Group, y=Geneexpr)) +
  # Define color and fill scales
  scale_color_manual(values = c("#DED6F2", "#F5AEBF", "#C0D9A3", "#FFE9B3")) +
  scale_fill_manual(values = c("#DED6F2", "#F5AEBF", "#C0D9A3", "#FFE9B3")) +
  # Add the custom flat violin plot, jitter, and boxplot geoms
  geom_flat_violin(aes(fill=Group, color=Group), position=position_nudge(x=.2), width=0.6) +
  geom_jitter(aes(color=Group), width=0.1, alpha=0.5) +
  geom_boxplot(width=.07, position=position_nudge(x=0.2), fill="white", size=0.2) +
  # Apply a theme and remove gridlines
  theme_bw() + theme(panel.grid = element_blank()) +
  # Customize axis text and title
  theme(axis.text = element_text(size=12, family="serif", colour = "black"),
        axis.title = element_text(size=14, family="serif", colour = "black"),
        legend.position="none") +
  # Label the axes and add comparisons using stat_compare_means
  labs(x = NULL, y = 'ROSMAP') +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.format",
                     method = "wilcox.test") +
  # Wrap the plot by Genenames with a single row
  facet_wrap(~Genenames, nrow = 1)

# Close the PDF device
dev.off()