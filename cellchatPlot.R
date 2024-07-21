library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

data.input <- sce.sub[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
meta <- sce.sub@meta.data# create a dataframe of the cell labels
meta$samples <- meta$orig.ident

data.input.T <- GetAssayData(subset(x = wb.normalized, group=="T"), assay = "RNA", slot = "data") # normalized data matrix

# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
meta.T <- subset(x = wb.normalized, group=="T")@meta.data# create a dataframe of the cell labels
meta.T$samples <- meta.T$orig.ident

cellchat <- createCellChat(object = data.input.T, meta = meta.T, group.by = "celltype.1")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communicatiocn analysis
CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#计算通信概率并推断细胞通信网络
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 5)

cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()

groupSize <- as.numeric(table(cellchat@idents))
pdf("./netVisual_circle.pdf",height = 5,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("./netVisual_circle-single.pdf",height =10,width = 13)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i])
}
dev.off()

cellchat@netP$pathways
# [1] "ApoE"       "APP"        "SPP1"       "PTN"        "MIF"        "CypA"       "PSAP"       "COMPLEMENT" "GALECTIN"   "CD99"      
# [11] "MK"         "MHC-II"     "IL1"        "NCAM"       "GRN"        "ANNEXIN"    "EGF"        "COLLAGEN"   "MAG"        "GAP"       
# [21] "CLDN"       "NRXN"       "CADM"       "FN1"        "TNF"        "CNTN"       "LAIR1"      "OSM"        "ADGRE"      "PLAU"      
# [31] "ICAM"       "GAS"        "CDH"        "JAM"        "MHC-I"      "THY1"       "LAMININ"    "NOTCH"      "PTPR"       "CD39"      
# [41] "FGF"        "KLK"        "MPZ"        "TENASCIN"   "ESAM"       "PDGF"       "PECAM1"     "SEMA4"      "PCDH"
pathways.show <- c("ApoE") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,3:5), remove.isolate = FALSE)
#> Comparing communications on a single object
#> 
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 3, targets.use = c(1:2,3:5), lab.cex = 0.5,legend.pos.y = 40)

plotGeneExpression(cellchat, signaling = "ApoE", enriched.only = TRUE, type = "violin")

ptm = Sys.time()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("ApoE", "APP"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

pdf("./netAnalysis_signalingRole_heatmap.pdf",width =10,height = 8)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 12,color.use = brewer.pal(n = 11, name = "Set3"),color.heatmap = "PuBu")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 12,color.use = brewer.pal(n = 11, name = "Set3"),color.heatmap = "PuBu")
ht1 + ht2
dev.off()

# Load the NMF library for Non-negative Matrix Factorization and ggalluvial for visualization
library(NMF)
library(ggalluvial)

# Create a PDF file to store the plot for selecting the optimal number of communication patterns (outgoing)
pdf("./outgoing-selectK.pdf", width = 8, height = 4)
# Use the selectK function from the cellchat package to identify and visualize outgoing communication patterns
selectK(cellchat, pattern = "outgoing")
# Close the PDF file device to finish writing the plot to the PDF
dev.off()

# Set the number of communication patterns to identify
nPatterns = 4

# Use the identifyCommunicationPatterns function to recognize outgoing communication patterns
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, height = 12)

# Create a PDF file for the river plot of outgoing communication patterns
pdf("./outgoing-netAnalysis_river.pdf", width = 10, height = 6)
# Perform network analysis using the netAnalysis_river function and specify colors
netAnalysis_river(cellchat, pattern = "outgoing", color.use = brewer.pal(n = 11, name = "Set3"),
                  color.use.pattern = c("#e3bec6","#efdad7","#9ad0ec","#1572a1","#AEDFFE"),
                  color.use.signaling = colorRampPalette(brewer.pal(11, "Set3"))(53))
# Close the PDF file device
dev.off()

# Create a PDF file for the dot plot of outgoing communication patterns
pdf("./outgoing-netAnalysis_dot.pdf", width = 10, height = 6)
# Perform network analysis using the netAnalysis_dot function and specify colors
netAnalysis_dot(cellchat, pattern = "outgoing", color.use = brewer.pal(n = 11, name = "Set3"))
# Close the PDF file device
dev.off()

# Repeat the process for incoming communication patterns
# ...

# Compute the similarity of signaling networks based on functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
# Perform embedding of signaling networks in a low-dimensional space based on the functional similarity
cellchat <- netEmbedding(cellchat, type = "functional")

# Cluster signaling networks based on functional similarity
cellchat <- netClustering(cellchat, type = "functional")

# Create a PDF file for visualizing the embedded signaling networks in 2D space (functional)
pdf("./functional-netVisual_embedding.pdf", width = 8, height = 6)
# Visualize the embedded signaling networks
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# Close the PDF file device
dev.off()

# Compute the similarity of signaling networks based on structural similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
# Perform embedding of signaling networks in a low-dimensional space based on the structural similarity
cellchat <- netEmbedding(cellchat, type = "structural")

# Cluster signaling networks based on structural similarity
cellchat <- netClustering(cellchat, type = "structural")

# Create a PDF file for visualizing the embedded signaling networks in 2D space (structural)
pdf("./structural-netVisual_embedding.pdf", width = 8, height = 6)
# Visualize the embedded signaling networks
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
# Close the PDF file device
dev.off()

# Create a PDF file for a signaling role network plot focusing on the 'CypA' signaling molecule
pdf("cellchat-CypA.pdf")
# Perform network analysis to visualize the signaling role of 'CypA' in the cellchat object
netAnalysis_signalingRole_network(cellchat, signaling = "CypA", 
                                  width = 10, height = 4, font.size = 10, 
                                  color.use = brewer.pal(n = 11, name = "Set3"),
                                  color.heatmap = "Blues")
# Close the PDF file device to finish writing the plot to the PDF
dev.off()

# Create a PDF file for a scatter plot visualization of the signaling role
pdf("cellchat-signalingRole_scatter.pdf")
# Perform a scatter plot analysis to visualize the signaling roles in the cellchat object
netAnalysis_signalingRole_scatter(cellchat)
# Close the PDF file device to finish writing the scatter plot to the PDF
dev.off()

# Define a vector of signaling molecules or pathways to visualize, in this case, only 'APP' is included
pathways.show <- c("APP") # Optionally include "PSAP", "CypA"

# Aggregate and visualize the signaling networks for the specified pathways in a circular layout
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Save the current state of the cellchat object to an RDS file for future use or sharing
saveRDS(cellchat, file = "cellchat.rds")
