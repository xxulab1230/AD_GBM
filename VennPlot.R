##### Model Important Variable Venn Diagram
# Read gene lists from the clipboard
ADNI_MT <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
RUSH_MT <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
ADNI_INTER <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
RUSH_INTER <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
TCGA_MT <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
CGGA_MT <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
TCGA_INTER <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)
CGGA_INTER <- read.table(file = "clipboard", header = FALSE, stringsAsFactors = FALSE)

# Create a list of gene lists
gene_list <- list("ADNI_MT" = ADNI_MT$V1,
                  "RUSH_MT" = RUSH_MT$V1,
                  "ADNI_INTER" = ADNI_INTER$V1,
                  "RUSH_INTER" = RUSH_INTER$V1,
                  "TCGA_MT" = TCGA_MT$V1,
                  "CGGA_MT" = CGGA_MT$V1,
                  "TCGA_INTER" = TCGA_INTER$V1,
                  "CGGA_INTER" = CGGA_INTER$V1)

# Load the UpSetR library for visualization of set intersections
library(UpSetR)
library(RColorBrewer)

# Create an UpSet plot with specified parameters
upset(data.frame(gene_list), nsets = 8, 
      nintersects = 30, 
      keep.order = T, 
      mb.ratio = c(0.6, 0.4), 
      point.size = 3, 
      line.size = 1, 
      matrix.color = "#AFAFAF", 
      shade.color = "#64FFFF", 
      matrix.dot.alpha = 0.5, 
      shade.alpha = 0.25, 
      main.bar.color = "#808080", 
      mainbar.y.label = "Number of Sample Overlaps", 
      mainbar.y.max = 50, 
      show.numbers = "yes", 
      sets.bar.color = brewer.pal(5, "Pastel1"), 
      sets.x.label = "Number of Samples", 
      set_size.scale_max = 800, 
      text.scale = c(1.5,1,1.5,1,1.5,1),  ### Adjustments for sorting and text size in the plot
      group.by = "degree", 
      order.by = "freq")

# Initialize total sum to 0
total_sum <- 0

# Iterate over each sublist and calculate the sum of elements
for (sub_list in gene_list) {
  sub_list_sum <- length(sub_list)
  total_sum <- total_sum + sub_list_sum
}

# Output the total sum
print(total_sum)

# Load the SuperExactTest library for testing and visualization of set intersections
library(SuperExactTest)
# Perform a super exact test on the gene lists
Result <- supertest(gene_list, n = total_sum)

# Plot the results with a landscape layout, sorted by size
plot(Result, Layout = "landscape", sort.by = "size", keep = FALSE,
     bar.split = c(70, 180), show.elements = TRUE, elements.cex = 0.6,
     elements.list = subset(summary(Result)$Table, Observed.Overlap <= 20),
     show.expected.overlap = TRUE, expected.overlap.style = "hatchedBox",
     color.expected.overlap = 'red')

# Save the plot to a PDF file
pdf('examples/ex5.png', width = 4000, height = 2000, res = 300)
grid.newpage()
# Create viewports for a grid layout
vp0 <- viewport(layout = grid.layout(1, 2))
vp1 <- viewport(layout.pos.col = 1, layout.pos.row = 1, name = "plot_left")
vp2 <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "plot_right")
vps <- vpTree(vp0, vpList(vp1, vp2))

pushViewport(vps)

# Plot the results in the left viewport with a landscape layout
seekViewport("plot_left")
plot(Result, Layout = "landscape", sort.by = 'size', keep = FALSE,
     bar.split = c(70, 180), show.elements = TRUE, elements.cex = 0.5,
     elements.list = subset(summary(Result)$Table, Observed.Overlap <= 5),
     show.expected.overlap = TRUE, expected.overlap.style = "hatchedBox",
     color.expected.overlap = '#F35185', color.on = c("#6BBE01", "#F94F35", "#A2E9EC", "#266DD7", "#A68BEC", "#F6D839", "#EC69E8", "#F47328"),
     title = 'Figure A. Landscape layout', new.gridPage = FALSE)

# Plot the results in the right viewport with a circular layout
seekViewport("plot_right")
plot(Result, Layout = "circular", sort.by = 'size', keep = FALSE,
     show.expected.overlap = TRUE, expected.overlap.style = "hatchedBox",
     color.expected.overlap = '#F35185', color.on = c("#6BBE01", "#F94F35", "#A2E9EC", "#266DD7", "#A68BEC", "#F6D839", "#EC69E8", "#F47328"),
     title = 'Figure B. Circular layout', new.gridPage = FALSE)