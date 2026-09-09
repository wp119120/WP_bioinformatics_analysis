# ==============================================================================
# Script: 09_Nomogram_and_DCA_Analysis.R
# Purpose: Build a multivariable diagnostic nomogram with calibration and DCA validation.
# Relevant Figures: Figure 4A, Figure 4B, Figure 4C, Figure 4D
# ==============================================================================

suppressPackageStartupMessages({
  library(rms)
  library(rmda)
  library(pROC)
  library(ggplot2)
})

# 1. Prepare Model Environment -------------------------------------------------
rt <- read.table("merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
gene_tbl <- read.table("important.genes.txt", header = TRUE, sep = "\t")
core_genes <- intersect(gene_tbl$Gene[1:6], rownames(rt))

df_expr <- as.data.frame(t(rt[core_genes, ]))
df_expr$GroupType <- ifelse(grepl("Control", rownames(df_expr), ignore.case = TRUE), 0, 1)

dd <- datadist(df_expr)
options(datadist = "dd")

# 2. Nomogram Construction (Figure 4A) -----------------------------------------
reg_formula <- as.formula(paste("GroupType ~", paste(core_genes, collapse = " + ")))
lrm_fit <- lrm(reg_formula, data = df_expr, x = TRUE, y = TRUE)

nomo <- nomogram(lrm_fit, fun = plogis, fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9), lp = FALSE, funlabel = "Fibrosis Risk")
pdf("Figure4A_Nomogram.pdf", width = 10, height = 6)
plot(nomo)
dev.off()

# 3. Calibration Curve (Figure 4B) ---------------------------------------------
cal <- calibrate(lrm_fit, method = "boot", B = 1000)
pdf("Figure4B_Calibration.pdf", width = 5.5, height = 5.5)
plot(cal, xlab = "Nomogram Predicted Probability", ylab = "Observed Rate", main = "Calibration Curve")
dev.off()

# 4. Decision Curve Analysis (Figure 4C) ---------------------------------------
dca_fit <- decision_curve(formula = reg_formula, data = df_expr, bootstraps = 100)
pdf("Figure4C_DCA.pdf", width = 6, height = 6)
plot_decision_curve(dca_fit, curve.names = "6-Gene Nomogram", col = "#DE3E54", confidence.intervals = FALSE)
dev.off()

# 5. Diagnostic Model ROC (Figure 4D) ------------------------------------------
pred_probs <- predict(lrm_fit, type = "fitted")
roc_obj <- roc(df_expr$GroupType, pred_probs, quiet = TRUE)

pdf("Figure4D_Nomogram_ROC.pdf", width = 5.5, height = 5.5)
plot(roc_obj, print.auc = TRUE, col = "#DE3E54", lwd = 3, legacy.axes = TRUE,
     main = sprintf("Nomogram ROC (AUC = %.3f)", auc(roc_obj)))
dev.off()