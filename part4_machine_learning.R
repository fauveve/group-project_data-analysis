# ============================================================
# LSABD group project - Part 4: Machine learning
# Task: classify Schizophrenia vs Control from AnCg gene expression
# Models: elastic net (glmnet) and random forest (ranger),
#         compared with a covariate-only baseline and a permutation baseline
# NOTE: starter script, written but not run. Check every step before you trust it.
# Needs: results/preprocessed_GSE80655.rds from the Part 1 script.
# ============================================================

# install.packages(c("glmnet", "ranger", "pROC", "caret", "tidyverse"))
library(tidyverse)
library(glmnet)
library(ranger)
library(pROC)
library(caret)

set.seed(42)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Data ---------------------------------------------------------------------
obj    <- readRDS("results/preprocessed_GSE80655.rds")
logcpm <- obj$logcpm
meta   <- obj$meta

# One region only: each donor then appears once, so a plain stratified split has no
# donor leakage. If you pool regions, split BY DONOR (caret::groupKFold or similar).
sel <- meta$brain_region == "AnCg" &
       meta$diagnosis %in% c("Control", "Schizophrenia") &
       !is.na(meta$brain_pH)
X <- t(logcpm[, sel])                                     # samples x genes
y <- factor(as.character(meta$diagnosis[sel]), levels = c("Control", "Schizophrenia"))
covs <- meta[sel, c("age", "PMI", "brain_pH", "sex")]
table(y)                                                  # ~48 samples: very small

# ---- 2. Pipeline function -----------------------------------------------------------
# Everything that learns from labels or data is done INSIDE the training set only:
# variable-gene selection, scaling, hyperparameter tuning. The test set is only transformed.
run_enet <- function(Xtr, ytr, Xte, n_var = 2000, alphas = c(0.2, 0.5, 1), seed = 1) {
  top <- names(sort(apply(Xtr, 2, var), decreasing = TRUE))[1:n_var]
  mu  <- colMeans(Xtr[, top]); s <- apply(Xtr[, top], 2, sd)
  Ztr <- scale(Xtr[, top], mu, s); Zte <- scale(Xte[, top], mu, s)
  set.seed(seed)
  foldid <- createFolds(ytr, k = 5, list = FALSE)         # stratified folds
  fits <- lapply(alphas, function(a)
    cv.glmnet(Ztr, ytr, family = "binomial", alpha = a,
              type.measure = "auc", foldid = foldid))
  best <- which.max(sapply(fits, function(f) max(f$cvm)))
  fit  <- fits[[best]]
  list(prob = as.numeric(predict(fit, Zte, s = "lambda.min", type = "response")),
       fit = fit, alpha = alphas[best], genes = top, Ztr = Ztr, Zte = Zte)
}

# ---- 3. One stratified train/test split (70/30) ------------------------------------------
tr <- createDataPartition(y, p = 0.7, list = FALSE)[, 1]
Xtr <- X[tr, ]; Xte <- X[-tr, ]; ytr <- y[tr]; yte <- y[-tr]
cat("Train:", length(ytr), " Test:", length(yte), "\n")

enet <- run_enet(Xtr, ytr, Xte)
cat("Chosen alpha:", enet$alpha, " lambda.min:", enet$fit$lambda.min, "\n")

roc_enet <- roc(yte, enet$prob, levels = c("Control", "Schizophrenia"), direction = "<")
cat("Elastic net test AUC:", round(auc(roc_enet), 3), "\n")
print(ci.auc(roc_enet))                                   # wide CI expected with ~14 test samples

pred_class <- factor(ifelse(enet$prob > 0.5, "Schizophrenia", "Control"),
                     levels = levels(y))
print(confusionMatrix(pred_class, yte, positive = "Schizophrenia"))

# ---- 4. Inspect the fitted model ------------------------------------------------------------
co <- coef(enet$fit, s = "lambda.min")
co <- data.frame(ensembl_id = rownames(co), coef = as.numeric(co)) %>%
  filter(coef != 0, ensembl_id != "(Intercept)") %>% arrange(desc(abs(coef)))
cat("Non-zero genes:", nrow(co), "\n")
head(co, 20)
write.csv(co, "results/enet_coefficients.csv", row.names = FALSE)

# Cross-validation curve (tuning of lambda)
png("figures/enet_cv_curve.png", width = 600, height = 450); plot(enet$fit); dev.off()

# ---- 5. Comparison models ---------------------------------------------------------------------
# (a) Random forest on the same training-selected genes
rf <- ranger(x = enet$Ztr, y = ytr, probability = TRUE, num.trees = 1000,
             importance = "permutation", seed = 42)
p_rf <- predict(rf, enet$Zte)$predictions[, "Schizophrenia"]
roc_rf <- roc(yte, p_rf, levels = c("Control", "Schizophrenia"), direction = "<")
cat("Random forest test AUC:", round(auc(roc_rf), 3), "\n")
head(sort(rf$variable.importance, decreasing = TRUE), 20)

# (b) Covariate-only baseline: do the genes add anything beyond age/PMI/pH/sex?
cov_tr <- covs[tr, ]; cov_te <- covs[-tr, ]
glm_cov <- glm(ytr ~ age + PMI + brain_pH + sex, data = cov_tr, family = binomial)
p_cov <- predict(glm_cov, cov_te, type = "response")
roc_cov <- roc(yte, p_cov, levels = c("Control", "Schizophrenia"), direction = "<")
cat("Covariate-only test AUC:", round(auc(roc_cov), 3), "\n")

# ROC curves
png("figures/roc_curves.png", width = 600, height = 500)
plot(roc_enet, col = "firebrick", legacy.axes = TRUE)
plot(roc_rf, col = "steelblue", add = TRUE)
plot(roc_cov, col = "grey40", add = TRUE)
legend("bottomright", c("Elastic net", "Random forest", "Covariates only"),
       col = c("firebrick", "steelblue", "grey40"), lwd = 2)
dev.off()

# ---- 6. How stable is the result? 50 repeated random splits ----------------------------------------
# A single split of ~14 test samples is very noisy; report the distribution instead.
aucs <- replicate(50, {
  i <- createDataPartition(y, p = 0.7, list = FALSE)[, 1]
  p <- run_enet(X[i, ], y[i], X[-i, ], seed = sample.int(1e6, 1))$prob
  as.numeric(auc(roc(y[-i], p, levels = c("Control", "Schizophrenia"), direction = "<", quiet = TRUE)))
})
cat("Repeated-split AUC: mean", round(mean(aucs), 3), " sd", round(sd(aucs), 3), "\n")

# ---- 7. Permutation baseline: AUC when labels carry no information ----------------------------------
# Shuffle training labels, repeat the WHOLE pipeline (incl. gene selection and tuning).
perm_auc <- replicate(100, {
  ytr_perm <- sample(ytr)
  p <- run_enet(Xtr, ytr_perm, Xte, seed = sample.int(1e6, 1))$prob
  as.numeric(auc(roc(yte, p, levels = c("Control", "Schizophrenia"), direction = "<", quiet = TRUE)))
})
cat("Permutation AUC: mean", round(mean(perm_auc), 3), "\n")
cat("Permutation p-value:", (sum(perm_auc >= auc(roc_enet)) + 1) / (length(perm_auc) + 1), "\n")

ggplot(data.frame(auc = perm_auc), aes(auc)) +
  geom_histogram(bins = 25, colour = "white") +
  geom_vline(xintercept = as.numeric(auc(roc_enet)), colour = "firebrick") +
  labs(x = "Test AUC with permuted labels", y = "Count")
ggsave("figures/permutation_auc.png", width = 5, height = 4)

# ---- 8. Interpretation reminders ------------------------------------------------------------------------
# - Report test-set (held-out) performance only, never training performance.
# - Compare with the covariate-only model and the permutation baseline.
# - With ~48 samples, results are unstable: report the spread and the limitations.
# - Class balance is ~50/50 so accuracy is fair, but still report AUC and the confusion matrix.
# - Small caveat: logCPM/TMM from Part 1 were computed on all samples (unsupervised, label-free),
#   which is a minor, usually accepted, form of leakage. Mention it in the report.
sessionInfo()
