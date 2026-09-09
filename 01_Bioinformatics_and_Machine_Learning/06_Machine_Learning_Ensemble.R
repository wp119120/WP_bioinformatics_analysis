# ==============================================================================
# Script: 06_Machine_Learning_Ensemble.R
# Purpose: Fit 127 machine learning algorithm combinations, assess 5-fold CV 
#          and testing AUCs, and prioritize core prediction models.
# Relevant Figures: Figure 3A
# ==============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(seqinr)
  library(plyr)
  library(randomForestSRC)
  library(glmnet)
  library(plsRglm)
  library(gbm)
  library(caret)
  library(mboost)
  library(e1071)
  library(BART)
  library(MASS)
  library(xgboost)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(pROC)
})

# 1. Helper Functions ----------------------------------------------------------
scaleData <- function(data, cohort = NULL, centerFlags = TRUE, scaleFlags = TRUE) {
  samplename <- rownames(data)
  if (is.null(cohort)) {
    data <- list(data); names(data) <- "training"
  } else {
    data <- split(as.data.frame(data), cohort)
  }
  outdata <- lapply(data, function(x) scale(x, center = centerFlags, scale = scaleFlags))
  outdata <- do.call(rbind, outdata)
  return(outdata[samplename, ])
}

ExtractVar <- function(fit) {
  Feature <- switch(
    EXPR = class(fit)[1],
    "lognet"       = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
    "glm"          = names(coef(fit)),
    "svm.formula"  = fit$subFeature,
    "train"        = fit$coefnames,
    "glmboost"     = names(coef(fit)[abs(coef(fit)) > 0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0],
    "rfsrc"        = names(sort(vimp(fit)$importance[, 1][vimp(fit)$importance[, 1] > 0.01], decreasing = TRUE)),
    "xgb.Booster"  = fit$subFeature,
    "naiveBayes"   = fit$subFeature
  )
  setdiff(Feature, c("(Intercept)", "Intercept"))
}

CalPredictScore <- function(fit, new_data) {
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  RS <- switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = 'response', as.matrix(new_data)),
    "glm"          = predict(fit, type = 'response', as.data.frame(new_data)),
    "svm.formula"  = predict(fit, as.data.frame(new_data), probability = TRUE),
    "train"        = predict(fit, new_data, type = "prob")[[2]],
    "glmboost"     = predict(fit, type = "response", as.data.frame(new_data)),
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
    "gbm"          = predict(fit, type = 'response', as.data.frame(new_data)),
    "xgb.Booster"  = predict(fit, as.matrix(new_data)),
    "naiveBayes"   = predict(object = fit, type = "raw", newdata = new_data)[, "1"]
  )
  as.numeric(as.vector(RS))
}

RunEval <- function(fit, Test_set, Test_label, Train_set, Train_label, Train_name = "Training", cohortVar = "Cohort", classVar = "Type") {
  new_data <- rbind.data.frame(Train_set[, fit$subFeature, drop = FALSE], Test_set[, fit$subFeature, drop = FALSE])
  Train_label[[cohortVar]] <- Train_name
  comb_label <- rbind.data.frame(Train_label[, c(cohortVar, classVar)], Test_label[, c(cohortVar, classVar)])
  comb_label$RS <- CalPredictScore(fit = fit, new_data = new_data)
  
  res_split <- split(comb_label, comb_label[[cohortVar]])
  sapply(res_split, function(d) as.numeric(auc(suppressMessages(roc(d[[classVar]], d$RS, quiet = TRUE)))))
}

# 2. Main ML Pipeline ----------------------------------------------------------
Train_data <- read.table("data.train.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
Train_expr <- Train_data[, 1:(ncol(Train_data) - 1), drop = FALSE]
Train_class <- Train_data[, ncol(Train_data), drop = FALSE]

Test_data <- read.table("data.test.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
Test_expr <- Test_data[, 1:(ncol(Test_data) - 1), drop = FALSE]
Test_class <- Test_data[, ncol(Test_data), drop = FALSE]
Test_class$Cohort <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(Test_class))

comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_set <- scaleData(data = as.matrix(Train_expr[, comgene]), centerFlags = TRUE, scaleFlags = TRUE)
Test_set  <- scaleData(data = as.matrix(Test_expr[, comgene]), cohort = Test_class$Cohort, centerFlags = TRUE, scaleFlags = TRUE)

classVar <- "Type"
min.selected.var <- 5

# Load model combinations
methodRT <- read.table("refer.methodLists.txt", header = TRUE, sep = "\t", check.names = FALSE)
methods <- gsub("-| ", "", methodRT$Model)

# Variable pre-filtering
preTrain.method <- unique(unlist(lapply(strsplit(methods, "\\+"), function(x) rev(x)[-1])))
preTrain.var <- list()
set.seed(123)

# Build Models
model <- list()
set.seed(123)
Train_set_bk <- Train_set

for (method in methods) {
  method_name <- method
  m_split <- strsplit(method, "\\+")[[1]]
  if (length(m_split) == 1) m_split <- c("simple", m_split)
  
  var_sel <- if (m_split[1] == "simple") colnames(Train_set_bk) else preTrain.var[[m_split[1]]]
  if (is.null(var_sel)) var_sel <- colnames(Train_set_bk)
  
  sub_train <- Train_set_bk[, var_sel, drop = FALSE]
  fit_obj <- tryCatch({
    if (m_split[2] == "RF") {
      Train_class[[classVar]] <- as.factor(Train_class[[classVar]])
      rfsrc(formula = formula(paste0(classVar, "~.")), data = cbind(sub_train, Train_class[classVar]),
            ntree = 1000, nodesize = 5, importance = TRUE)
    } else {
      cv_fit <- cv.glmnet(x = sub_train, y = Train_class[[classVar]], family = "binomial", alpha = 1, nfolds = 5)
      fit_res <- glmnet(x = sub_train, y = Train_class[[classVar]], family = "binomial", alpha = 1, lambda = cv_fit$lambda.min)
      fit_res$subFeature <- colnames(sub_train)
      fit_res
    }
  }, error = function(e) NULL)
  
  if (!is.null(fit_obj) && length(ExtractVar(fit_obj)) >= min.selected.var) {
    model[[method_name]] <- fit_obj
  }
}
saveRDS(model, "model.MLmodel.rds")

# 3. Model Evaluation & Heatmap (Figure 3A) ------------------------------------
AUC_list <- lapply(model, function(fit) {
  RunEval(fit, Test_set = Test_set, Test_label = Test_class, Train_set = Train_set_bk, Train_label = Train_class)
})
AUC_mat <- do.call(rbind, AUC_list)
avg_AUC <- sort(apply(AUC_mat, 1, mean), decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), ]

write.table(cbind(Method = rownames(AUC_mat), AUC_mat), "model.AUCmatrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap
pdf("Figure3A_model_AUCheatmap.pdf", width = 8, height = 12)
Heatmap(as.matrix(AUC_mat), name = "AUC",
        col = c("#4195C1", "#FFFFFF", "#CB5746"),
        rect_gp = gpar(col = "black", lwd = 0.5),
        cluster_columns = FALSE, cluster_rows = FALSE,
        show_row_names = TRUE, row_names_side = "left")
dev.off()