# Loading necessary R packages for data manipulation and visualization
# ---------------------------
library(ggrepel)  # For adding text annotations to the plot
library(ggfun)    # Utility functions for ggplot2
library(grid)     # Grid graphics

# Plotting the volcano plot to visualize gene expression changes
# ---------------------------
# Specify genes of interest to be displayed on the plot

# The following script selects the top 10 up-regulated and top 10 down-regulated genes based on their significance
# The ggplot function is used to create a volcano plot with points colored and sized according to their log-fold change and adjusted p-values

# Define the ggplot object with layers for points, text annotations, and additional plot elements
{p <- ggplot(data = df) + 
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = logFC,
                 size = -log10(adj.P.Val))) + 
  # Add text annotations for the top 10 up-regulated genes using ggrepel
  geom_text_repel(data = df %>%
                    dplyr::filter(significant != "None") %>%
                    dplyr::arrange(desc(-log10(adj.P.Val))) %>%
                    dplyr::slice(1:10) %>%
                    dplyr::filter(significant == "Up"),
                  aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                  nudge_x = 0.5, nudge_y = 0.2, segment.curvature = -0.1, direction = "y", hjust = "left", max.overlaps = 200) +
  # Add text annotations for the top 10 down-regulated genes using ggrepel
  geom_text_repel(data = df %>%
                    dplyr::filter(significant != "None") %>%
                    dplyr::filter(significant != "Up") %>%
                    dplyr::arrange(desc(-log10(adj.P.Val))) %>%
                    dplyr::slice(1:10) %>%
                    dplyr::filter(significant == "Down"),
                  aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
                  nudge_x = -0.2, nudge_y = 0.2, segment.curvature = -0.1, segment.angle = 20, direction = "y", hjust = "left", max.overlaps = 200) + 
  # Add color and fill gradients based on log-fold change values
  scale_color_gradientn(colours = c("#DCADDF", "#FFB1CC", "#FCD0BE", "#E4F9BE", "#A4D9F9"), values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#DCADDF", "#FFB1CC", "#FCD0BE", "#E4F9BE", "#A4D9F9"), values = seq(0, 1, 0.2)) +
  # Add vertical and horizontal lines for visual guides
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  # Set size scaling for points
  scale_size(range = c(1,7)) + 
  # Add plot title and adjust plot limits
  ggtitle(label = "Volcano Plot") + 
  xlim(c(-3, 3)) + ylim(c(-0.5, 30)) +
  # Apply a theme to the plot
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)
  ) 
# Additional annotation and custom graphical elements can be added here if needed

# Save the plot to a PDF file with specified dimensions
pdf("Volcano Plot.pdf", height = 8, width = 10)
print(p)
dev.off()
}
# Export the processed data to a CSV file for further analysis or sharing
write.csv(df, "./deg.csv")