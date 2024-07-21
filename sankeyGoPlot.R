# Define global variables to avoid duplication in the function
globalVariables(c("gene_ratio", 'next_node', "next_x", "node", "term_ratio"))

# Load the ggsankey library for creating Sankey diagrams
library(ggsankey)

# Define the sankeyGoPlot function for generating a Sankey diagram and a GO plot
sankeyGoPlot <- function(goData = NULL,
                         topGenes = 5,
                         sankeyExpand = c(0.5,1),
                         nodeSize = 2.5,
                         nodeColor = NULL,
                         goCol = NULL,
                         xShift = 0.05,
                         downShift = 4.5,
                         upShift = 0.25){
  
  # =========== Data Preparation ============
  # Group the GO data by Description and calculate the gene ratio and term ratio
  ego_df <- goData %>%
    dplyr::group_by(Description) %>%
    dplyr::mutate(gene_ratio = eval(parse(text = GeneRatio))) %>%
    dplyr::arrange(pvalue)
  
  # Set the order of terms based on the Description factor
  ego_df$Description <- factor(ego_df$Description, levels = rev(ego_df$Description))
  
  # Calculate the term_ratio as the difference in genes / total genes in the pathway
  ego_df$term_gene <- sapply(strsplit(ego_df$BgRatio, split = "\\/"), "[", 1) %>% as.numeric()
  ego_df$term_ratio <- ego_df$Count / ego_df$term_gene
  
  # Select top genes for each term and prepare data for the Sankey diagram
  sankey_df <- ego_df %>%
    dplyr::select(Description, geneID) %>%
    tidyr::separate_longer_delim(geneID, delim = "/") %>%
    dplyr::group_by(Description) %>%
    dplyr::slice_head(n = topGenes) %>%
    as.data.frame() %>%
    dplyr::mutate(Description = as.character(Description))
  
  # Convert the data to long format for the Sankey diagram
  sankey_df_long <- sankey_df %>% ggsankey::make_long(geneID, Description)
  
  # Set the order of nodes in the Sankey diagram
  sankey_df_long$node <- factor(sankey_df_long$node,
                                levels = c(unique(sankey_df$Description), unique(sankey_df$geneID)))
  
  # =========== Plotting ============
  # Define color palette for nodes
  if(is.null(nodeColor)){
    mycol <- cols4all::c4a('rainbow_wh_rd', length(unique(sankey_df_long$node)))
  }else{
    mycol <- grDevices::colorRampPalette(colors = nodeColor)(length(unique(sankey_df_long$node)))
  }
  
  # Create the Sankey plot with specified aesthetics
  ps <-
    ggplot(data = sankey_df_long,
           mapping = aes(x = x,
                         next_x = next_x,
                         node = node,
                         next_node = next_node,
                         fill = factor(node),
                         label = node)) +
    geom_sankey(flow.alpha = 0.5,
                flow.fill = 'grey',
                flow.color = 'grey',
                width = 0.1,
                node.fill = mycol) +
    theme_void() +
    scale_x_discrete(expand = expansion(mult = sankeyExpand)) +
    theme(legend.position = 'none')
  
  # Extract data from the Sankey plot for further manipulation
  ps_data <- ggplot_build(ps)
  ps_data_info <- ps_data$data[[2]]
  
  # Add text labels to the Sankey plot
  ps2 <-
    ps +
    geom_text(data = ps_data_info,
              mapping = aes(next_x = 0, next_node = 0,
                            x = xmin, y = (ymin + ymax) / 2,
                            label = label, hjust = 1),
              fontface = "bold", size = nodeSize)
  
  # =========== GO Plot ============
  # Define color scale for the GO plot
  if(is.null(goCol)){
    pcol <- scale_fill_viridis_c(option = "plasma", direction = -1)
  }else{
    pcol <- scale_fill_gradient(low = goCol[1], high = goCol[2])
  }
  
  # Create the GO plot with points representing GO terms and their significance
  pp <-
    ggplot(ego_df) +
    geom_point(aes(x = -log10(pvalue), y = Description,
                   size = term_ratio, fill = gene_ratio),
               color = "black", shape = 21) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(colour = "black"),
          plot.background = element_blank()) +
    ylab("") + pcol
  
  # =========== Position Calculation ============
  # Calculate the relative positions for the GO plot within the Sankey plot
  xmin <- max(ps_data_info$xmax)
  y_range <- subset(ps_data_info, x == "2")
  ymin <- min(c(y_range$ymin, y_range$ymax))
  ymax <- max(c(y_range$ymin, y_range$ymax))
  
  # =========== Final Plot Assembly ============
  # Combine the Sankey plot with the GO plot using the calculated positions
  cb <-
    ps2 +
    annotation_custom(grob = ggplotGrob(pp),
                      xmin = xmin - xShift, xmax = 2 + sankeyExpand[2],
                      ymin = ymin - downShift, ymax = ymax + upShift)
  
  # Return the combined plot as a ggplot object
  return(cb)
}