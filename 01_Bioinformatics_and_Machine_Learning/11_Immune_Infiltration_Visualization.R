# ==============================================================================
# Script: 11_Immune_Infiltration_Visualization.R
# Purpose: Visualize immune landscape (stacked barplot, boxplot, and correlation).
# Relevant Figures: Figure 5A, Figure 5B, Figure 5C, Figure 5D, Figure 5E
# ==============================================================================

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)
  library(corrplot)
})

ciber_res <- read.csv("CIBERSORT_Results.csv", row.names = 1, check.names = FALSE)
ciber_res$Group <- ifelse(grepl("Control", rownames(ciber_res), ignore.case = TRUE), "Control", "Fibrosis")
ciber_res$Sample <- rownames(ciber_res)

# 1. Immune Cell Composition Stacked Barplot (Figure 5A) -----------------------
data_long <- melt(ciber_res, id.vars = c("Sample", "Group"), variable.name = "CellType", value.name = "Fraction")
palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(data_long$CellType)))

p_bar <- ggplot(data_long, aes(x = Sample, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Infiltration Proportions of 22 Immune Cell Types", y = "Relative Proportion", fill = "Cell Type")

pdf("Figure5A_Immune_StackedBarplot.pdf", width = 12, height = 6)
print(p_bar)
dev.off()

# 2. Immune Infiltration Boxplot Comparison (Figure 5B) ------------------------
pdf("Figure5B_Immune_Difference_Boxplot.pdf", width = 12, height = 6)
ggplot(data_long, aes(x = CellType, y = Fraction, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "#377EB8", "Fibrosis" = "#E41A1C")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test")
dev.off()

# 3. Correlation between Core Genes and Immune Fractions (Figure 5D-E) ---------
rt_exp <- read.table("merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
gene_tbl <- read.table("important.genes.txt", header = TRUE, sep = "\t")
core_genes <- intersect(gene_tbl$Gene[1:6], rownames(rt_exp))

ciber_cells <- ciber_res[, !(colnames(ciber_res) %in% c("Group", "Sample"))]
common_samples <- intersect(rownames(ciber_cells), colnames(rt_exp))

cor_mat <- cor(t(rt_exp[core_genes, common_samples]), ciber_cells[common_samples, ], method = "spearman")

pdf("Figure5D_Gene_Immune_Correlation.pdf", width = 10, height = 5)
corrplot(cor_mat, method = "circle", tl.col = "black", tl.srt = 45, col = colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(100))
dev.off()