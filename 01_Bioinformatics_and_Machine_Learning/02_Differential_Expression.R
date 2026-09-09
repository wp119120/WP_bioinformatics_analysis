# ==============================================================================
# Script: 02_Differential_Expression.R
# Purpose: Identify DEGs in lung fibrosis using limma; generate volcano and heatmaps.
# Relevant Figures: Figure 2C, Figure 2D
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
})

# Thresholds set according to the manuscript (|log2FC| > 0.14, FDR < 0.05)
logFCfilter <- 0.14
adj.P.Val.Filter <- 0.05
inputFile <- "merge.normalize.txt"

# 1. Data Preparation ----------------------------------------------------------
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

# Extract phenotypic groups and sort
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]
Project <- gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

# 2. Linear Modeling (limma) ----------------------------------------------------
design <- model.matrix(~0 + factor(Type))
colnames(design) <- c("Control", "Treat")
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(Treat - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Export all genes and statistically significant DEGs
allDiff <- topTable(fit2, adjust = 'fdr', number = Inf)
allDiffOut <- rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = "all.txt", sep = "\t", quote = FALSE, col.names = FALSE)

diffSig <- allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut <- rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = "diff.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Export DEG expression values
diffGeneExp <- data[rownames(diffSig), ]
diffGeneExpOut <- rbind(id = paste0(colnames(diffGeneExp), "_", Type), diffGeneExp)
write.table(diffGeneExpOut, file = "diffGeneExp.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# 3. Heatmap of Top DEGs (Figure 2D) -------------------------------------------
geneNum <- 50
diffUp <- diffSig[diffSig$logFC > 0, ]
diffDown <- diffSig[diffSig$logFC < 0, ]
geneUp <- rownames(diffUp)[1:min(nrow(diffUp), geneNum)]
geneDown <- rownames(diffDown)[1:min(nrow(diffDown), geneNum)]

hmExp <- data[c(geneUp, geneDown), ]
annotation_df <- data.frame(Project = Project, Type = Type, row.names = colnames(data))

pdf(file = "Figure2D_heatmap.pdf", width = 10, height = 7)
pheatmap(hmExp,
         annotation_col = annotation_df,
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         fontsize = 8,
         fontsize_row = 5.5,
         fontsize_col = 8)
dev.off()

# 4. Volcano Plot (Figure 2C) --------------------------------------------------
rt_vol <- read.table("all.txt", header = TRUE, sep = "\t", check.names = FALSE)
rt_vol$Sig <- ifelse((rt_vol$adj.P.Val < adj.P.Val.Filter) & (abs(rt_vol$logFC) > logFCfilter),
                     ifelse(rt_vol$logFC > logFCfilter, "Up", "Down"), "Not")

p_vol <- ggplot(rt_vol, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(col = Sig), size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("Up" = "red2", "Down" = "green2", "Not" = "grey")) +
  theme_bw() +
  labs(title = "Differentially Expressed Genes in Lung Fibrosis", x = "log2(Fold Change)", y = "-log10(FDR)") +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

pdf(file = "Figure2C_volcano.pdf", width = 5.5, height = 4.5)
print(p_vol)
dev.off()