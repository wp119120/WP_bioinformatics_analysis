# ==============================================================================
# Script: 07_SHAP_Interpretability.R
# Purpose: Compute SHAP values for the top ML ensemble to quantify feature importance.
# Relevant Figures: Figure 3B, Figure 3C, Figure 3D, Figure 3E
# ==============================================================================

suppressPackageStartupMessages({
  library(caret)
  library(DALEX)
  library(ggplot2)
  library(randomForest)
  library(kernelshap)
  library(shapviz)
})

inputFile <- "merge.normalize.txt"
geneFile  <- "model.genes.txt"
bestMethod <- "Stepglm[backward]+RF"

# 1. Data Preparation ----------------------------------------------------------
data <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
geneRT <- read.table(geneFile, header = TRUE, sep = "\t", check.names = FALSE)
geneRT <- geneRT[geneRT$algorithm == bestMethod, ]
common_genes <- intersect(as.vector(geneRT[, 1]), rownames(data))

data <- t(data[common_genes, , drop = FALSE])
rownames(data) <- gsub("-", "_", rownames(data))
data <- as.data.frame(data)
data$Type <- ifelse(grepl("Control", rownames(data), ignore.case = TRUE), 0, 1)

# 2. Fit Model & Calculate SHAP ------------------------------------------------
set.seed(12345)
control <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE)
model <- train(Type ~ ., data = data, method = "rf", trControl = control)

X_data <- data[, -ncol(data)]
fit_shap <- permshap(model, X_data)
shp <- shapviz(fit_shap, X_pred = X_data, X = X_data, interactions = TRUE)

# Export Top Ranked Genes
important <- sort(colMeans(abs(shp$S)), decreasing = TRUE)
write.table(data.frame(Gene = names(important), Importance = important), 
            file = "important.genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 3. Visualizations ------------------------------------------------------------
pdf("Figure3B_SHAP_Barplot.pdf", width = 6, height = 5)
sv_importance(shp, kind = "bar", show_numbers = TRUE) + theme_bw()
dev.off()

pdf("Figure3C_SHAP_Beeswarm.pdf", width = 7, height = 6)
sv_importance(shp, kind = "bee", show_numbers = TRUE) + theme_bw()
dev.off()

pdf("Figure3D_SHAP_Dependence.pdf", width = 10, height = 6)
sv_dependence(shp, v = names(important)[1:min(6, length(important))]) + theme_bw()
dev.off()

pdf("Figure3E_SHAP_Force.pdf", width = 9, height = 5)
sv_force(shp, row_id = 1)
dev.off()