# ==============================================================================
# Script: 12_Immune_Correlation_and_Network.R
# Purpose: Inter-immune cell correlation matrices and correlation profiling between
#          core target genes (LCN2, TYMS, etc.) and immune infiltrates.
# Relevant Figures: Figure 5C, Figure 5D, Figure 5E
# ==============================================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(corrplot)
  library(dplyr)
  library(linkET)
  library(ggplot2)
})

# 1. Load Datasets -------------------------------------------------------------
immune_df <- read.csv("CIBERSORT_Results.csv", header = TRUE, row.names = 1, check.names = FALSE)
immune_mat <- immune_df[, 1:22]
# Remove cell types with zero variation
immune_mat <- immune_mat[, apply(immune_mat, 2, sd) > 0, drop = FALSE]

expr_mat <- read.table("merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
core_genes <- c("LCN2", "TYMS", "BCAT2", "NQO1", "PRDX6", "CLPP")
core_genes <- intersect(core_genes, rownames(expr_mat))

common_samples <- intersect(rownames(immune_mat), colnames(expr_mat))
immune_mat <- immune_mat[common_samples, ]
expr_sub <- t(expr_mat[core_genes, common_samples])

# 2. Inter-Cellular Correlation Heatmap with Significance (Figure 5C) ----------
cor_cell <- cor(immune_mat, method = "spearman")
p_cell <- matrix(NA, nrow = ncol(immune_mat), ncol = ncol(immune_mat))

for (i in 1:ncol(immune_mat)) {
  for (j in 1:ncol(immune_mat)) {
    p_cell[i, j] <- cor.test(immune_mat[, i], immune_mat[, j], method = "spearman")$p.value
  }
}

stars_cell <- matrix("", nrow = nrow(cor_cell), ncol = ncol(cor_cell))
stars_cell[p_cell < 0.05]  <- "*"
stars_cell[p_cell < 0.01]  <- "**"
stars_cell[p_cell < 0.001] <- "***"

pdf("Figure5C_Immune_Interactions_Heatmap.pdf", width = 10, height = 10)
pheatmap(
  cor_cell,
  color = colorRampPalette(c("#3498DB", "#FFFFFF", "#E74C3C"))(100),
  display_numbers = stars_cell,
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Correlation Matrix of Infiltrating Immune Cells"
)
dev.off()

# 3. Core Genes vs. Immune Cells Correlation Heatmap (Figure 5D) ---------------
cor_gene_cell <- matrix(NA, nrow = ncol(expr_sub), ncol = ncol(immune_mat),
                        dimnames = list(colnames(expr_sub), colnames(immune_mat)))
p_gene_cell <- cor_gene_cell

for (g in colnames(expr_sub)) {
  for (c in colnames(immune_mat)) {
    test <- cor.test(expr_sub[, g], immune_mat[, c], method = "spearman")
    cor_gene_cell[g, c] <- test$estimate
    p_gene_cell[g, c] <- test$p.value
  }
}

label_gene_cell <- matrix("", nrow = nrow(cor_gene_cell), ncol = ncol(cor_gene_cell))
for (i in 1:nrow(cor_gene_cell)) {
  for (j in 1:ncol(cor_gene_cell)) {
    val <- sprintf("%.2f", cor_gene_cell[i, j])
    star <- if (p_gene_cell[i, j] < 0.001) "***" else if (p_gene_cell[i, j] < 0.01) "**" else if (p_gene_cell[i, j] < 0.05) "*" else ""
    label_gene_cell[i, j] <- paste0(val, "\n", star)
  }
}

pdf("Figure5D_Gene_Immune_Correlation_Heatmap.pdf", width = 12, height = 4.5)
pheatmap(
  cor_gene_cell,
  color = colorRampPalette(c("#27AE60", "#FFFFFF", "#C0392B"))(100),
  display_numbers = label_gene_cell,
  fontsize_number = 7.5,
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Correlation between Core Biomarkers and Immune Infiltration"
)
dev.off()

# 4. LinkET Correlation Network (Figure 5E) -------------------------------------
corr_list <- data.frame()
for (g in colnames(expr_sub)) {
  for (c in colnames(immune_mat)) {
    test <- cor.test(immune_mat[, c], expr_sub[, g], method = "spearman")
    corr_list <- rbind(corr_list, data.frame(
      spec = g, env = c, r = abs(as.numeric(test$estimate)),
      p = as.numeric(test$p.value),
      sign = ifelse(test$estimate > 0, "Positive", "Negative")
    ))
  }
}

corr_list$pd <- ifelse(corr_list$p < 0.05, corr_list$sign, "Not significant")
corr_list$rd <- cut(corr_list$r, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6"))

write.csv(corr_list, file = "gene_immune_correlation_details.csv", row.names = FALSE)

pdf("Figure5E_Gene_Immune_LinkET_Network.pdf", width = 11, height = 8)
p_link <- qcorrplot(correlate(immune_mat, method = "spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = corr_list, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "Pastel2")), name = "Cell Correlation") +
  scale_size_manual(values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3), name = "|Cor|") +
  scale_colour_manual(values = c("Positive" = "#E74C3C", "Negative" = "#3498DB", "Not significant" = "#BDC3C7"), name = "Association") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_link)
dev.off()