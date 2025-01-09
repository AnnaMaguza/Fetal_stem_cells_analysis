#                                          PART 1 - Data loading
### Load required modules

library(NMF)
library(dplyr)
library(igraph)
library(Matrix)
library(ggplot2)
library(CellChat) 
library(patchwork)
library(ggalluvial)
library(reticulate)
options(stringsAsFactors = FALSE)
library(png)
library(Cairo)
library(ggalluvial)

use_python("/Users/annamaguza/miniconda3/envs/scanpy_env/bin", required = TRUE)

### Read in data

ad <- import("anndata", convert = FALSE)
pd <- import("pandas", convert = FALSE)
scipy <- import("scipy", convert = FALSE)
ad_object <- ad$read_h5ad("/Users/annamaguza/Desktop/CellChat/FetalSC_prepared_for_CellChat.h5ad")

### Access expression matrix

data.input <- t(py_to_r(ad_object$X))

### Add metadata

rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))

meta.data <- py_to_r(ad_object$obs)
meta <- meta.data

#                                          PART 2 - CellChat Object Creation and Preprocessing 

### Create `cellchat` object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cluster")

rm(ad_object)
rm(data.input)
rm(meta)
rm(meta.data)

### Set up ligand-receptor interaction database for `cellchat`

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

### Process expression data

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

rm(CellChatDB.use)

cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)

### Infer cell-cell communication
cellchat <- computeCommunProbPathway(cellchat)

### Calculate aggregated cell-cell communication

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))


#                                          PART 3 - Visualise number of interactions and interaction strength  


# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Number_of_interactions.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
# Close the device (this will save the file)
dev.off()


# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Interaction_weights_strength_epi.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
# Close the device (this will save the file)
dev.off()


#                                          PART 4 - Visualise top pathways
### Export results as dataframe

df.net <- subsetCommunication(cellchat)
head(df.net)

# Filter the results to include only Club cell interactions
cells <- c("ASS1+_SLC40A1+_SC", 'RPS10+_RPS17+_SC', 'FXYD3+_CKB+_SC')

### Print unique pathway names found in the df.net dataframe.
unique(df.net$pathway_name)

# Define the signaling pathways
signaling_pathways <- c("EGF", "VISFATIN", "MK", "PARs", "LAMININ", "COLLAGEN", "THBS", "APP", "CD99", "CDH", "DESMOSOME", "EPHA", "JAM", "NECTIN", "OCLN")

# Define the directory where the plots will be saved
output_dir <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Pathways"

# Generate and save the plots
for (pathway in signaling_pathways) {
  # Open the PNG device
  png(filename = file.path(output_dir, paste0(pathway, ".png")), width = 3000, height = 3500, res = 300)
  
  # Generate the plot
  netVisual_aggregate(cellchat, signaling = pathway, layout = "chord", sources.use = cells, targets.use = NULL, pt.title = 20, title.space = 12, vertex.label.cex = 1.5, show.legend = FALSE)
  
  # Close the PNG device
  dev.off()
}

#                                          PART 5 - Sankey Plot - outgoing patterns

selectK(cellchat, pattern = "outgoing")

options(repr.plot.width = 15, repr.plot.height = 15)
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,  width = 10, height = 10)

# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Sankey_plot_outgoing.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
netAnalysis_river(cellchat, pattern = "outgoing", font.size = 2.5, font.size.title = 20)
# Close the device (this will save the file)
dev.off()

# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Dot_plot_outgoing.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
netAnalysis_dot(cellchat, pattern = "outgoing")
# Close the device (this will save the file)
dev.off()

# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Enriched_outgoing_pathways.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
pairLR.pathway <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.pathway[1,] # show one ligand-receptor pair
# Close the device (this will save the file)
dev.off()

#                                          PART 5 - Sankey Plot - incoming patterns

### Identify global communication patterns

selectK(cellchat, pattern = "incoming")

options(repr.plot.width = 15, repr.plot.height = 15)
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,  width = 10, height = 10)

# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Sankey_plot_incoming.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
netAnalysis_river(cellchat, pattern = "incoming", font.size = 2.5, font.size.title = 20)
# Close the device (this will save the file)
dev.off()

# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Dot_plot_incoming.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
netAnalysis_dot(cellchat, pattern = "incoming")
# Close the device (this will save the file)
dev.off()


#                                          PART 7 - Hierarchy plot

# File path where the image will be saved
file_path <- "/Users/annamaguza/Desktop/Stem_cells_trajectory_analysis/3_fetalSC_characteristics/figures/CellChat/Hierarchy_plot.png"
# Start the png device
png(filename = file_path, width = 3000, height = 3000, res = 300)
# Generate your plot
vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# Close the device (this will save the file)
dev.off()
