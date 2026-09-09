# ==============================================================================
# Script: 04_GO_Enrichment_Analysis.R
# Purpose: Perform GO pathway enrichment analysis on candidate targets.
# Relevant Figures: Figure 2G
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(ggpubr)
  library(ComplexHeatmap)
})

pvalueFilter <- 0.05
p.adjustFilter <- 0.05
colorSel <- ifelse(p.adjustFilter > 0.05, "pvalue", "p.adjust")

# 1. Gene ID Conversion --------------------------------------------------------
rt <- read.table("interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- na.omit(as.character(mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)))

# 2. GO Enrichment -------------------------------------------------------------
kk <- enrichGO(gene = entrezIDs, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
GO <- as.data.frame(kk)
GO <- GO[(GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter), ]
write.table(GO, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 3. Barplot & Dotplot ---------------------------------------------------------
pdf(file = "Figure2G_GO_barplot.pdf", width = 11, height = 7)
bar <- barplot(kk, drop = TRUE, showCategory = 10, label_format = 100, split = "ONTOLOGY", color = colorSel) +
       facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

pdf(file = "Figure2G_GO_bubble.pdf", width = 11, height = 7)
bub <- dotplot(kk, showCategory = 10, orderBy = "GeneRatio", label_format = 100, split = "ONTOLOGY", color = colorSel) +
       facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

# 4. GO Circlize Plot ----------------------------------------------------------
ontology.col <- c("#00CC33FF", "#FFC20AFF", "#CC33FFFF")
data <- GO[order(GO$p.adjust), ]
datasig <- data[data$pvalue < 0.05, , drop = FALSE]

BP <- head(datasig[datasig$ONTOLOGY == "BP", , drop = FALSE], 6)
CC <- head(datasig[datasig$ONTOLOGY == "CC", , drop = FALSE], 6)
MF <- head(datasig[datasig$ONTOLOGY == "MF", , drop = FALSE], 6)
data <- rbind(BP, CC, MF)
main.col <- ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio <- Gene / BgGene
logpvalue <- -log10(data$pvalue)
logpvalue.col <- brewer.pal(n = 6, name = "Reds")
f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10), colors = logpvalue.col)
BgGene.col <- f(logpvalue)

df <- data.frame(GO = data$ID, start = 1, end = max(BgGene))
rownames(df) <- df$GO
bed2 <- data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 <- data.frame(GO = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 <- data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5

pdf(file = "Figure2G_GO_circlize.pdf", width = 10, height = 10)
par(omi = c(0.1, 0.1, 0.1, 1.5))
circos.par(track.margin = c(0.01, 0.01))
circos.genomicInitialize(df, plotType = "none")

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
}

circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[, 2], border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], adj = 0, cex = 0.8, ...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], cex = 0.9, adj = 0, ...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white", bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim <- get.cell.meta.data("cell.xlim")
                      cell.ylim <- get.cell.meta.data("cell.ylim")
                      for (j in 1:9) {
                        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[, 1], col = value[, 2], border = NA, ...)
                    })
circos.clear()

# Legends
middle.legend <- Legend(
  labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
  type = "points", pch = c(15, 15, 17), legend_gp = gpar(col = c('pink', '#BA55D3', ontology.col[1])),
  title = "", nrow = 3, size = unit(3, "mm")
)
circle_size <- unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend <- Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),
  type = "points", pch = 15, legend_gp = gpar(col = ontology.col),
  title_position = "topcenter", title = "ONTOLOGY", nrow = 3,
  size = unit(3, "mm"), grid_height = unit(5, "mm"), grid_width = unit(5, "mm")
)

logp.legend <- Legend(
  labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '>=10'),
  type = "points", pch = 16, legend_gp = gpar(col = logpvalue.col),
  title = "-log10(Pvalue)", title_position = "topcenter",
  grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), size = unit(3, "mm")
)

lgd <- packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
dev.off()