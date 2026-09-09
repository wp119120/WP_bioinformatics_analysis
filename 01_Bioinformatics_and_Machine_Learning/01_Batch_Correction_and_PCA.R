# ==============================================================================
# Script: 01_Batch_Correction_and_PCA.R
# Purpose: Merge multi-center lung fibrosis transcriptomic datasets, perform 
#          ComBat batch correction (sva), and assess clustering via PCA and Boxplots.
# Relevant Figures: Figure 2A, Figure 2B
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
})

# 1. Read and Intersect Genes across Datasets -----------------------------------
files <- grep("normalize.txt$", dir(), value = TRUE)
files <- setdiff(files, c("merge.preNorm.txt", "merge.normalize.txt"))

geneList <- list()
for (file in files) {
  rt <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  geneNames <- trimws(rt[, 1])
  header <- unlist(strsplit(file, "\\.|\\-"))[1]
  geneList[[header]] <- unique(geneNames)
}
interGenes <- Reduce(intersect, geneList)

# 2. Harmonize Datasets and Merge -----------------------------------------------
allTab <- data.frame()
batchType <- c()

for (i in seq_along(files)) {
  inputFile <- files[i]
  header <- unlist(strsplit(inputFile, "\\.|\\-"))[1]
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  geneNames <- trimws(rt[, 1])
  expr <- rt[, -1, drop = FALSE]
  
  data <- apply(expr, 2, as.numeric)
  rownames(data) <- geneNames
  
  rt_avg <- avereps(data)
  colnames(rt_avg) <- paste0(header, "_", colnames(rt_avg))
  rt_subset <- rt_avg[interGenes, , drop = FALSE]
  
  if (i == 1) {
    allTab <- rt_subset
  } else {
    allTab <- cbind(allTab, rt_subset)
  }
  batchType <- c(batchType, rep(i, ncol(rt_subset)))
}

# Export pre-normalized expression matrix
outTab_pre <- rbind(geneNames = colnames(allTab), allTab)
write.table(outTab_pre, file = "merge.preNorm.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# 3. ComBat Empirical Bayes Batch Correction -----------------------------------
combat_data <- ComBat(allTab, batchType, par.prior = TRUE)
outTab_combat <- rbind(geneNames = colnames(combat_data), combat_data)
write.table(outTab_combat, file = "merge.normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# 4. Boxplot Quality Assessment ------------------------------------------------
plot_batch_boxplot <- function(inputFile, outFile, titleName) {
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  data <- t(rt)
  Project <- gsub("(.*?)\\_.*", "\\1", rownames(data))
  Sample <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", rownames(data))
  data <- cbind(as.data.frame(data), Sample, Project)
  
  rt1 <- melt(data, id.vars = c("Project", "Sample"))
  colnames(rt1) <- c("Project", "Sample", "Gene", "Expression")

  pdf(file = outFile, width = 10, height = 5)
  p <- ggplot(rt1, aes(x = Sample, y = Expression)) +
    geom_boxplot(aes(fill = Project), notch = TRUE, outlier.shape = NA) +
    ggtitle(titleName) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 2),
          plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

plot_batch_boxplot("merge.preNorm.txt", "boxplot.preNorm.pdf", "Before batch correction")
plot_batch_boxplot("merge.normalize.txt", "boxplot.normalize.pdf", "After batch correction")

# 5. Principal Component Analysis (Figure 2A-B) --------------------------------
plot_pca <- function(inputFile, outFile, titleName) {
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  data <- t(rt)
  Project <- gsub("(.*?)\\_.*", "\\1", rownames(data))
  
  data.pca <- prcomp(data)
  pcaPredict <- predict(data.pca)
  PCA <- data.frame(PC1 = pcaPredict[, 1], PC2 = pcaPredict[, 2], Type = Project)

  pdf(file = outFile, width = 5.5, height = 4.25)
  p <- ggscatter(data = PCA, x = "PC1", y = "PC2", color = "Type", shape = "Type",
                 ellipse = TRUE, ellipse.type = "norm", ellipse.border.remove = FALSE, ellipse.alpha = 0.1,
                 size = 2, main = titleName, legend = "right") +
       theme(plot.margin = unit(rep(1.5, 4), 'lines'), plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

plot_pca("merge.preNorm.txt", "Figure2A_PCA.preNorm.pdf", "Before batch correction")
plot_pca("merge.normalize.txt", "Figure2B_PCA.normalize.pdf", "After batch correction")