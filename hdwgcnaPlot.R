# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 5)

library(future)
plan(multisession, workers = 8)

library(scCustomize)

pbmc_sub <- pbmc_sub %>% 
  RunHarmony("projid", plot_convergence = T)

Idents(pbmc_sub) <- "broad.cell.type"

seurat_obj <- SetupForWGCNA(
  pbmc_sub,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("broad.cell.type"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'broad.cell.type' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
name = "AD-Ast-"
groupname = "Astrocyte"

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = groupname, # the name of the group of interest in the group.by column
  group.by='broad.cell.type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
pdf(paste0(name,"PlotSoftPowers.pdf"),width = 8,height = 6)
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,overwrite_tom = TRUE,
  tom_name = 'Excitatory Neurons' # name of the topoligical overlap matrix written to disk
)

pdf(paste0(name,"PlotDendrogram.pdf"),width = 8,height = 6)
PlotDendrogram(seurat_obj, main=paste0(groupname,'hdWGCNA Dendrogram'))
dev.off()

TOM <- GetTOM(seurat_obj)

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="projid"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'broad.cell.type', group_name = groupname
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = groupname
)

pdf(paste0(name,"PlotKMEs.pdf"),width = 10,height = 5)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=4)

p
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

write.csv(modules,file = paste0(name,"modules.csv"))


# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)
write.csv(hub_df,file = paste0(name,"hub_df.csv"))

saveRDS(seurat_obj, file=paste0(name,'hdWGCNA.rds'))

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
pdf(paste0(name,"ModuleFeaturePlot.pdf"),width = 10,height = 5)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
dev.off()

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
pdf(paste0(name,"ModuleFeaturePlot-scores.pdf"),width = 10,height = 5)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
dev.off()

# plot module correlagram
pdf(paste0(name,"ModuleCorrelogram.pdf"),width = 10,height = 10)
ModuleCorrelogram(seurat_obj)
dev.off()

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'broad.cell.type')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='#EC83C5', mid='#FFF8CC', low='#79BBEA')

# plot output
pdf(paste0(name,"DotPlot.pdf"),width = 10,height = 5)
p
dev.off()

####Network Visualization#####
ModuleNetworkPlot(
  seurat_obj, 
  outdir= paste0(name,"ModuleNetworks"), # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)


# hubgene network
pdf(paste0(name,"HubGeneNetworkPlot.pdf"))
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 5, n_other=10,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

#Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=5 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)

##### get modules and TOM from the seurat obj####
modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey') %>% 
  mutate(module = droplevels(module))
mods <- levels(modules$module)
TOM <- GetTOM(seurat_obj)

# get module colors for plotting 
mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct()
mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)

# load the GO Biological Pathways file (donwloaded from EnrichR website)
pathways <- fgsea::gmtPathways("./GO_Biological_Process_2021.txt")

# remove GO Term ID for simplicity:
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")

# selected pathway 
cur_pathway <- 'mitochondrial ATP synthesis coupled electron transport'

# get genes in this pathway 
cur_genes <- pathways[[cur_pathway]]
cur_genes <- cur_genes[cur_genes %in% modules$gene_name]

# subset the TOM 
cur_TOM <- TOM[cur_genes,cur_genes] 

# set up the graph object with igraph & tidygraph
graph <- cur_TOM %>% 
  igraph::graph_from_adjacency_matrix(mode='undirected', weighted=TRUE) %>% 
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes) 

library(ggraph)
# make the plot with ggraph
p <- ggraph(graph) + 
  geom_edge_link(color='grey', alpha=0.2) + 
  geom_node_point(color='black') +
  geom_node_label(aes(label=name), repel=TRUE, max.overlaps=Inf, fontface='italic') 

pdf("./ggraph.pdf")
p
dev.off()

# set up the graph object with igraph & tidygraph
graph <- cur_TOM %>% 
  igraph::graph_from_adjacency_matrix(mode='undirected', weighted=TRUE) %>% 
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes) 

# add the module name to the graph:
V(graph)$module <- modules[V(graph)$name,'module']

# make the plot with gggraph
p <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight), color='grey') + 
  geom_node_point(aes(color=module)) +
  geom_node_label(aes(label=name), repel=TRUE, max.overlaps=Inf, fontface='italic') +
  scale_colour_manual(values=mod_cp)  
pdf("./ggraph2.pdf")
p
dev.off()


# only keep the upper triangular part of the TOM:
cur_TOM[upper.tri(cur_TOM)] <- NA

# cast the network from wide to long format
cur_network <- cur_TOM %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene1 = Var1, gene2 = Var2, weight=value) %>%
  subset(!is.na(weight))

# get the module & color info for gene1
temp1 <- dplyr::inner_join(
  cur_network,
  modules %>% 
    dplyr::select(c(gene_name, module, color)) %>% 
    dplyr::rename(gene1 = gene_name, module1=module, color1=color),
  by = 'gene1'
) %>% dplyr::select(c(module1, color1))

# get the module & color info for gene2
temp2 <- dplyr::inner_join(
  cur_network,
  modules %>% 
    dplyr::select(c(gene_name, module, color)) %>% 
    dplyr::rename(gene2 = gene_name, module2=module, color2=color),
  by = 'gene2'
) %>% dplyr::select(c(module2, color2))

# add the module & color info 
cur_network <- cbind(cur_network, temp1, temp2)

# set the edge color to the module's color if they are the two genes are in the same module 
cur_network$edge_color <- ifelse(
  cur_network$module1 == cur_network$module2, 
  as.character(cur_network$module1),
  'grey'
)

# keep this network before subsetting
cur_network_full <- cur_network 

# keep the top 10% of edges 
edge_percent <- 0.1
cur_network <- cur_network_full %>% 
  dplyr::slice_max(
    order_by = weight, 
    n = round(nrow(cur_network)*edge_percent)
  )

# make the graph object with tidygraph
graph <- cur_network %>% 
  igraph::graph_from_data_frame() %>%
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes)

# add the module name to the graph:
V(graph)$module <- modules[V(graph)$name,'module']

# get the top 25 hub genes for each module
hub_genes <- GetHubGenes(seurat_obj, n_hubs=50) %>% .$gene_name
V(graph)$hub <- ifelse(V(graph)$name %in% hub_genes, V(graph)$name, "")

# make the plot with gggraph
p <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  geom_node_label(aes(label=hub), repel=TRUE, max.overlaps=Inf, fontface='italic') +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) 

pdf("./AD-Ast-ggraph3.pdf")
p
dev.off()


# subset to only keep edges between genes in the same module
cur_network <- cur_network_full %>% 
  subset(module1 == module2)

# make the graph object with tidygraph
graph <- cur_network %>% 
  igraph::graph_from_data_frame() %>%
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes)

# add the module name to the graph:
V(graph)$module <- modules[V(graph)$name,'module']

# get the top 25 hub genes for each module
hub_genes <- GetHubGenes(seurat_obj, n_hubs=25) %>% .$gene_name
V(graph)$hub <- ifelse(V(graph)$name %in% hub_genes, V(graph)$name, "")

# make the plot with gggraph
p <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  geom_node_label(aes(label=hub), repel=TRUE, max.overlaps=Inf, fontface='italic') +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  NoLegend()
pdf("./ggraph4.pdf")
p
dev.off()


# randomly sample 50% of the edges within the same module
cur_network1 <- cur_network_full %>% 
  subset(module1 == module2) %>%
  group_by(module1) %>%
  sample_frac(0.5) %>% 
  ungroup()

# keep the top 10% of other edges 
edge_percent <- 0.50
cur_network2 <- cur_network_full %>% 
  subset(module1 != module2) %>%
  dplyr::slice_max(
    order_by = weight, 
    n = round(nrow(cur_network)*edge_percent)
  )

cur_network <- rbind(cur_network1, cur_network2)

# set factor levels for edges:
cur_network$edge_color <- factor(
  as.character(cur_network$edge_color),
  levels = c(mods, 'grey')
)

# rearrange so grey edges are on the bottom:
cur_network %<>% arrange(rev(edge_color))

# make the graph object with tidygraph
graph <- cur_network %>% 
  igraph::graph_from_data_frame() %>%
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes)

# add the module name to the graph:
V(graph)$module <- modules[V(graph)$name,'module']

# get the top 25 hub genes for each module
hub_genes <- GetHubGenes(seurat_obj, n_hubs=25) %>% .$gene_name
V(graph)$hub <- ifelse(V(graph)$name %in% hub_genes, V(graph)$name, "")

# 1. default layout
p1 <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  ggtitle("layout = 'stress' (auto)") +
  NoLegend()

# 2. Kamada Kawai (kk) layout
graph2 <- graph; E(graph)$weight <- E(graph)$weight + 0.0001
p2 <- ggraph(graph, layout='kk', maxiter=100) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  ggtitle("layout = 'kk'") +
  NoLegend()

# 3. igraph layout_with_fr
p3 <- ggraph(graph, layout=layout_with_fr(graph)) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  ggtitle("layout_with_fr()") +
  NoLegend()

# 4. igraph layout_as_tree
p4 <- ggraph(graph, layout=layout_as_tree(graph)) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  ggtitle("layout_as_tree()") +
  NoLegend()

# 5. igraph layout_nicely
p5 <- ggraph(graph, layout=layout_nicely(graph)) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  ggtitle("layout_nicely()") +
  NoLegend()

# 6. igraph layout_in_circle
p6 <- ggraph(graph, layout=layout_in_circle(graph)) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) + 
  geom_node_point(aes(color=module)) +
  scale_colour_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  ggtitle("layout_in_circle()") +
  NoLegend()


# make a combined plot
pdf("./ggraph5.pdf",width = 8,height = 5)
(p1 | p2 | p3) / (p4 | p5 | p6) 
dev.off()


# get the UMAP df and subset by genes that are in our graph
umap_df <- GetModuleUMAP(seurat_obj)
umap_layout <- umap_df[names(V(graph)),] %>% dplyr::rename(c(x=UMAP1, y = UMAP2, name=gene))
rownames(umap_layout) <- 1:nrow(umap_layout)

# create the layout
lay <- ggraph::create_layout(graph, umap_layout)
lay$hub <- V(graph)$hub

p <- ggraph(lay) + 
  geom_edge_link(aes(alpha=weight, color=edge_color)) +
  geom_node_point(data=subset(lay, hub == ''), aes(color=module, size=kME)) + 
  geom_node_point(data=subset(lay, hub != ''), aes(fill=module, size=kME), color='black', shape=21) +
  scale_colour_manual(values=mod_cp) +
  scale_fill_manual(values=mod_cp) +
  scale_edge_colour_manual(values=mod_cp) +
  geom_node_label(aes(label=hub), repel=TRUE, max.overlaps=Inf, fontface='italic') +
  NoLegend()

#####Enrichment analysis######
# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
write.csv(enrich_df,paste0(name,"enrich_df.csv"))

#
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 5, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
pdf(paste0(name,"EnrichrDotPlot.pdf"),width = 10,height = 15)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=5 # number of terms for each module
)
dev.off()


###Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$broad.cell.type
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf(paste0(name,"OverlapModulesDEGs.pdf"),width = 12,height = 10)
wrap_plots(plot_list, ncol=4)
dev.off()

# plot odds ratio of the overlap as a dot plot
pdf(paste0(name,"OverlapDotPlot.pdf"),width = 10,height = 10)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()