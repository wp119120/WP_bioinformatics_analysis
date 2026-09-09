# ==============================================================================
# Script: 10_CIBERSORT_Immune_Infiltration.R
# Purpose: Estimate 22 immune cell infiltration fractions using CIBERSORT (SVR).
# ==============================================================================

suppressPackageStartupMessages({
  library(e1071)
  library(parallel)
  library(preprocessCore)
})

coreAlgorithm <- function(refer_matrix, target_vector) {
  steps <- c(0.25, 0.5, 0.75)
  models <- lapply(steps, function(nu) svm(refer_matrix, target_vector, type = "nu-regression", kernel = "linear", nu = nu, scale = FALSE))
  
  rmse_vec <- sapply(models, function(m) {
    w <- t(m$coefs) %*% m$SV
    w[w < 0] <- 0
    w <- w / sum(w)
    sqrt(mean(((refer_matrix %*% t(w)) - target_vector)^2))
  })
  
  best_m <- models[[which.min(rmse_vec)]]
  weights <- t(best_m$coefs) %*% best_m$SV
  weights[weights < 0] <- 0
  list(final_weights = weights / sum(weights), best_rmse = min(rmse_vec))
}

runCIBERSORT <- function(sig_file, mix_file, perm = 1000) {
  sig_mat <- data.matrix(read.table(sig_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  mix_mat <- data.matrix(read.table(mix_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  
  common_genes <- intersect(rownames(sig_mat), rownames(mix_mat))
  sig_mat <- sig_mat[common_genes, ]
  mix_mat <- mix_mat[common_genes, ]
  sig_mat <- scale(sig_mat)
  
  res <- apply(mix_mat, 2, function(y) coreAlgorithm(sig_mat, scale(y))$final_weights)
  res_df <- as.data.frame(t(res))
  colnames(res_df) <- colnames(sig_mat)
  write.csv(res_df, "CIBERSORT_Results.csv")
  return(res_df)
}

# Execute
cibersort_res <- runCIBERSORT("LM22.txt", "merge.normalize.txt", perm = 1000)