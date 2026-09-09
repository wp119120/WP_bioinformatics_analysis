# ==============================================================================
# Script: 03_Target_Venn_Intersection.R
# Purpose: Intersect PFOA-predicted targets with lung fibrosis DEGs to identify
#          common candidate genes.
# Relevant Figures: Figure 2E, Figure 2F
# ==============================================================================

suppressPackageStartupMessages({
  library(ggvenn)
})

compoundName <- "PFOA"
diseaseName <- "Lung_Fibrosis"

# 1. Read Target Lists ---------------------------------------------------------
# Compound.txt: 363 targets predicted by SwissTargetPrediction/PharmMapper/ChEMBL
pfoa_targets <- as.vector(read.table("Compound.txt", header = FALSE, sep = "\t")[, 1])
# Disease.txt: DEGs from lung fibrosis differential expression analysis
fibrosis_degs <- as.vector(read.table("Disease.txt", header = FALSE, sep = "\t")[, 1])

geneList <- list(
  PFOA = pfoa_targets,
  Pulmonary_Fibrosis = fibrosis_degs
)

# 2. Venn Diagram (Figure 2E) --------------------------------------------------
pdf(file = "Figure2E_venn.pdf", width = 6, height = 6)
ggvenn(geneList,
       show_percentage = TRUE,
       stroke_color = "white",
       stroke_size = 0.5,
       fill_color = c("#E41A1C", "#1E90FF"),
       set_name_color = c("#E41A1C", "#1E90FF"),
       set_name_size = 6,
       text_size = 4.5)
dev.off()

# 3. Intersect Targets and Network Attributes (Figure 2F) ----------------------
interGenes <- Reduce(intersect, geneList)
write.table(interGenes, file = "interGenes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Generate interaction and node tables for Cytoscape network
networkTab <- rbind(
  cbind(compoundName, interGenes, "Compound"),
  cbind(diseaseName, interGenes, "Disease")
)
colnames(networkTab) <- c("Node1", "Node2", "Type")
write.table(networkTab, file = "net.network.txt", sep = "\t", quote = FALSE, row.names = FALSE)

nodeTab <- rbind(
  cbind(compoundName, "Compound"),
  cbind(diseaseName, "Disease"),
  cbind(interGenes, "Gene")
)
colnames(nodeTab) <- c("Node", "Type")
write.table(nodeTab, file = "net.node.txt", sep = "\t", quote = FALSE, row.names = FALSE)