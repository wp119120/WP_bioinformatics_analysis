# ==============================================================================
# Script: 08_Core_Biomarkers_Evaluation.R
# Purpose: Validate 6 core genes via multi-gene ROC curves and differential boxplots.
# Relevant Figures: Figure 3F, Figure 3G
# ==============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})

expFile  <- "merge.normalize.txt"
geneFile <- "important.genes.txt"

# 1. Load Data -----------------------------------------------------------------
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
y <- ifelse(grepl("Control", colnames(rt), ignore.case = TRUE), 0, 1)

gene_tbl <- read.table(geneFile, header = TRUE, sep = "\t", check.names = FALSE)
core_genes <- intersect(as.vector(gene_tbl[, 1])[1:6], rownames(rt))

# 2. Multi-gene ROC Curves (Figure 3F) -----------------------------------------
cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
pdf("Figure3F_CoreGenes_ROC.pdf", width = 6, height = 6)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity", main = "Diagnostic Efficacy of Core Biomarkers")
abline(0, 1, lty = 2, col = "gray70")

legend_text <- c()
for (i in seq_along(core_genes)) {
  g <- core_genes[i]
  roc_res <- roc(y, as.numeric(rt[g, ]), direction = "auto", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_res))
  if (auc_val < 0.5) auc_val <- 1 - auc_val
  
  lines(1 - roc_res$specificities, roc_res$sensitivities, col = cols[i], lwd = 2.5)
  legend_text <- c(legend_text, sprintf("%s (AUC = %.3f)", g, auc_val))
}
legend("bottomright", legend = legend_text, col = cols, lwd = 2.5, bty = "n", cex = 0.85)
dev.off()

# 3. Expression Boxplots (Figure 3G) -------------------------------------------
box_data <- data.frame()
for (g in core_genes) {
  tmp <- data.frame(Gene = g, Expression = as.numeric(rt[g, ]), Type = ifelse(y == 0, "Control", "Fibrosis"))
  box_data <- rbind(box_data, tmp)
}

p_box <- ggplot(box_data, aes(x = Gene, y = Expression, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 0.6, alpha = 0.5) +
  scale_fill_manual(values = c("Control" = "#4682B4", "Fibrosis" = "#FF6347")) +
  theme_minimal(base_size = 13) +
  labs(title = "Core Biomarker Expression", y = "Log2 Normalized Expression") +
  stat_compare_means(aes(group = Type), label = "p.signif", method = "t.test")

pdf("Figure3G_CoreGenes_Boxplot.pdf", width = 8, height = 5)
print(p_box)
dev.off()