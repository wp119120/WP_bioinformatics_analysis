# ==============================================================================
# Script: 05_ML_Data_Preparation.R
# Purpose: Align and batch-correct multi-cohort transcriptomic matrices for ML input.
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(sva)
})

geneFile <- "interGenes.txt"

# 1. Read files and extract intersection genes ---------------------------------
files <- grep("normalize.txt$", dir(), value = TRUE)
files <- setdiff(files, c("merge.preNorm.txt", "merge.normalize.txt"))

geneList <- list()
for (file in files) {
  rt <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  geneNames <- trimws(as.vector(rt[, 1]))
  header <- unlist(strsplit(file, "\\.|\\-"))[1]
  geneList[[header]] <- unique(geneNames)
}
interGenes <- Reduce(intersect, geneList)

# 2. Merge expression matrices -------------------------------------------------
allTab <- NULL
batchType <- c()
geneOrder <- c()

for (i in seq_along(files)) {
  inputFile <- files[i]
  header <- unlist(strsplit(inputFile, "\\.|\\-"))[1]
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  geneNames <- trimws(rt[, 1])
  expr <- as.matrix(rt[, -1])
  expr <- apply(expr, 2, as.numeric)
  rownames(expr) <- geneNames
  
  expr_avg <- avereps(expr)
  commonGenes <- intersect(rownames(expr_avg), interGenes)
  expr_sub <- expr_avg[commonGenes, , drop = FALSE]
  
  if (i == 1) {
    geneOrder <- commonGenes
  } else {
    expr_sub <- expr_sub[geneOrder, , drop = FALSE]
  }
  
  colnames(expr_sub) <- paste0(header, "_", colnames(expr_sub))
  
  if (is.null(allTab)) {
    allTab <- expr_sub
  } else {
    allTab <- cbind(allTab, expr_sub)
  }
  batchType <- c(batchType, rep(i, ncol(expr_sub)))
}

# 3. ComBat batch adjustment ---------------------------------------------------
if (length(unique(batchType)) > 1) {
  svaTab <- ComBat(allTab, batchType, par.prior = TRUE)
} else {
  svaTab <- allTab
}

# 4. Subset candidate targets and split cohorts --------------------------------
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
targetGenes <- trimws(as.vector(geneRT[, 1]))
commonTargets <- intersect(rownames(svaTab), targetGenes)
geneTab <- t(svaTab[commonTargets, ])

train_idx <- grepl("^merge", rownames(geneTab), ignore.case = TRUE)
trainExp <- geneTab[train_idx, , drop = FALSE]
testExp <- geneTab[!train_idx, , drop = FALSE]
rownames(trainExp) <- gsub("merge_", "Train.", rownames(trainExp))

extractType <- function(names) {
  ifelse(grepl("Control", names, ignore.case = TRUE), 0, 1)
}

trainExp <- cbind(as.data.frame(trainExp), Type = extractType(rownames(trainExp)))
testExp  <- cbind(as.data.frame(testExp), Type = extractType(rownames(testExp)))

write.table(cbind(ID = rownames(trainExp), trainExp), file = "data.train.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(ID = rownames(testExp), testExp), file = "data.test.txt", sep = "\t", quote = FALSE, row.names = FALSE)