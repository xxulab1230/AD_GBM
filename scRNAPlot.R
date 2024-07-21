# Load the Single-Cell Pipeline (SCP) library for single-cell data analysis
library(SCP)

# Check if the 'renv' package is installed and install it if necessary for managing R environments
if (!require("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Create a directory for the SCP environment and initialize it with renv
dir.create("~/SCP_env", recursive = TRUE) # The directory cannot be the home directory "~"
renv::init(project = "~/SCP_env", bare = TRUE, restart = TRUE)

# Subset the Seurat object 'sce.big' to include only cells with non-NA cell types
sce.sub <- subset(sce.big, cells = rownames(sce.big@meta.data[!is.na(sce.big@meta.data$cell_type),]))
# Set the default assay for the Seurat object to "RNA"
DefaultAssay(sce.sub) <- "RNA"

# Perform integration using the CCAIntegration method and create a new reduction for the integrated data
sce.sub <- IntegrateLayers(
  object = sce.sub, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# Perform integration using the HarmonyIntegration method and create a new reduction for the integrated data
sce.sub <- IntegrateLayers(
  object = sce.sub, 
  method = HarmonyIntegration,
  orig.reduction = "tsne", 
  new.reduction = "harmony",
  verbose = FALSE
)

# Plot the dimensionality reduction (harmony) with cell types labeled
DimPlot(sce.sub , label = TRUE, reduction = "harmony", group.by = "orig.ident")

# Join the integrated layers into the Seurat object
sce.sub <- JoinLayers(sce.sub)

# Calculate the percentage of features matching a pattern (e.g., mitochondrial genes) and create a new assay
sce.sub <- PercentageFeatureSet(sce.sub, pattern = "^MT-", col.name = "percent.mt")

# Run SCTransform on the Seurat object, regressing out the percentage of mitochondrial genes
sce.sub <- SCTransform(sce.sub, vars.to.regress = "percent.mt", verbose = FALSE)

# Run PCA, UMAP, and t-SNE on the Seurat object for dimensionality reduction
sce.sub <- RunPCA(sce.sub, verbose = FALSE)
sce.sub <- RunUMAP(sce.sub , dims = 1:30)
sce.sub <- RunTSNE(sce.sub, check_duplicates = FALSE)

# Run the standard SCP pipeline
sce.sub <- Standard_SCP(srt = sce.sub)

# Plot the dimensionality reduction (t-SNE) with cell types grouped and using a blank theme
CellDimPlot(
  srt = sce.sub, group.by = c("cell_type"),
  reduction = "tsne", theme_use = "theme_blank"
)

# Create a 3D plot of the cell types using the Seurat object's dimensionality reduction
CellDimPlot3D(srt = sce.sub, group.by = "cell_type")

# Convert the RNA and SCT assays in the Seurat object to the Assay class
sce.sub[["RNA"]] <- as(sce.sub[["RNA"]], "Assay")
sce.sub[["SCT"]] <- as(sce.sub[["SCT"]], "Assay")

# Run PAGA analysis on the Seurat object, grouping by cell type
sce.sub <- RunPAGA(
  srt = sce.sub, group_by = "cell_type",
  linear_reduction = "pca", nonlinear_reduction = "umap"
)

# Plot the PAGA analysis results with various labeling options
PAGAPlot(srt = pancreas_sub, reduction = "UMAP", label = TRUE, label_insitu = TRUE, label_repel = TRUE)

# Run RNA velocity analysis on the Seurat object, grouping by cell type
sce.sub <- RunSCVELO(
  srt = sce.sub, group_by = "cell_type",
  linear_reduction = "pca", nonlinear_reduction = "umap"
)

# Run differential expression testing, filtering for significant markers
sce.sub <- RunDEtest(srt = sce.sub, group_by = "cell_type", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sce.sub, group_by = "cell_type")

# Run enrichment analysis on the Seurat object, filtering for significant GO biological processes
wb.normalized <- RunEnrichment(
  srt = wb.normalized, group_by = "celltype.1", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)

# Set the working directory for outputting plots and results
setwd("~/AD/GBM/scRNA/SCP/intergrated")

# Plot enrichment analysis results as bar charts for specified cell types
EnrichmentPlot(
  srt = wb.normalized, group_by = "celltype.1", group_use = c("Astrocyte", "Microglia", "Oligodendrocyte", "Neoplastic Cells", "Macrophage", "Myeloid"),
  plot_type = "bar"
)

# Generate and save various plots as PDF files, including word clouds and network diagrams for enrichment analysis

pdf("EnrichmentPlot-wordcloud.pdf",height = 15,width = 20)
EnrichmentPlot(
  srt = wb.normalized, group_by = "celltype", group_use = c("Astrocyte","Microglia","Oligodendrocyte","Neoplastic Cells","Macrophage","Myeloid"),
  plot_type = "wordcloud"
)
dev.off()

pdf("EnrichmentPlot-wordcloud-feature.pdf",height = 15,width = 20)
EnrichmentPlot(
  srt = wb.normalized, group_by = "celltype", group_use = c("Astrocyte","Microglia","Oligodendrocyte","Neoplastic Cells","Macrophage","Myeloid"),
  plot_type = "wordcloud", word_type = "feature"
)
dev.off()


pdf("EnrichmentPlot-net.pdf",height = 10,width = 12)
EnrichmentPlot(
  srt = sce.sub, group_by = "celltype", group_use = "Myeloid",
  plot_type = "network"
)
dev.off()

pdf("EnrichmentPlot-netmap.pdf",height = 10,width = 12)
EnrichmentPlot(
  srt = sce.sub, group_by = "celltype", group_use = "Myeloid",
  plot_type = "enrichmap"
)
dev.off()

# "Endothelial_cell","Neoplastic_cell","Myeloid"
# "Oligodendrocyte"，"Astrocyte"

pdf("EnrichmentPlot-compare.pdf",height = 18,width = 15)
EnrichmentPlot(srt = sce.sub, group_by = "celltype", plot_type = "comparison")
dev.off()

##Enrichment analysis(GSEA)
sce.sub <- RunGSEA(
  srt = sce.sub, group_by = "celltype", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)

pdf("GSEAPlot-Neoplastic_cell.pdf",height = 10,width = 12)
GSEAPlot(
  srt = sce.sub, group_by = "cell_type", group_use = "Neoplastic_cell", plot_type = "bar",
  direction = "both", topTerm = 20
)
dev.off()

GSEAPlot(srt = sce.sub, group_by = "celltype", group_use = "Neoplastic_cell", id_use = "GO:0007186")

pdf("GSEAPlot-comparison.pdf",height = 10,width = 12)
GSEAPlot(srt = sce.sub, group_by = "celltype", plot_type = "comparison")
dev.off()

# Run the Slingshot trajectory analysis, using the Seurat object 'sce.sub' and grouping by 'cell_type'
sce.sub <- RunSlingshot(srt = sce.sub, group.by = "cell_type", reduction = "tsne")

# Create a PDF file to store the plot with specified height and width
pdf("TJ-FeatureDimPlot-tsne.pdf", height = 8, width = 15)

# Generate a feature dimensionality plot for lineage-specific genes using the t-SNE reduction and a blank theme
FeatureDimPlot(sce.sub, features = paste0("Lineage", 1:3), reduction = "tsne", theme_use = "theme_blank")

# Close the PDF file device to finish writing the plot to the PDF
dev.off()

# Generate a cell dimensionality plot for the Seurat object 'sce.sub', using the t-SNE reduction
# and displaying lineages with a specified span
CellDimPlot(sce.sub, group.by = "cell_type", reduction = "tsne", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)

# Run the DynamicFeatures function to identify dynamic features across specified lineages in the Seurat object 'sce.sub'
sce.sub <- RunDynamicFeatures(srt = sce.sub, lineages = c("Lineage1", "Lineage2", "Lineage3"), n_candidates = 200)

# Generate a dynamic heatmap using the results of the dynamic features analysis
# with various annotation and display options
ht <- DynamicHeatmap(
  srt = sce.sub, lineages = c("Lineage1", "Lineage2", "Lineage3"),
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
  species = "Homo_sapiens", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis", cell_annotation = "cell_type",
  separate_annotation = list("cell_type", c("NOTCH1", "TNF")), separate_annotation_palette = c("Paired", "Set1"),
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)

# Print the heatmap plot
print(ht$plot)

# Generate a dynamic plot for the Seurat object 'sce.sub' with specified lineages and features
# comparing lineages and not comparing features
DynamicPlot(
  srt = sce.sub, lineages = c("Lineage1", "Lineage2", "Lineage3"), group.by = "celltype",
  features = c("NOTCH1", "TNF", "AKT1", "ABCA2", "NRP1"),
  compare_lineages = TRUE, compare_features = FALSE
)

# Generate a feature statistics plot for the Seurat object 'sce.sub' with specified grouping and backgrounding
# adding boxplots and making specific comparisons between cell types
FeatureStatPlot(
  srt = sce.sub, group.by = "celltype", bg.by = "celltype",
  stat.by = c("NOTCH1", "TNF", "AKT1", "ABCA2", "NRP1"), add_box = TRUE,
  comparisons = list(
    c("Neoplastic_cell", "Myeloid"),
    c("Oligodendrocyte", "Astrocyte")
  )
)
