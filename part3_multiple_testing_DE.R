# ============================================================
# LSABD group project - Part 3: Multiple testing study
# Contrast of interest: Schizophrenia vs Control in the anterior cingulate cortex (AnCg)
# Model per gene: limma-voom, ~ diagnosis + age + PMI + brain_pH + sex
# NOTE: starter script, written but not run. Check every step before you trust it.
# Needs: results/preprocessed_GSE80655.rds from the Part 1 script.
# ============================================================

# BiocManager::install(c("edgeR", "limma", "org.Hs.eg.db", "AnnotationDbi"))
library(tidyverse)
library(edgeR)
library(limma)
library(ggrepel)

set.seed(42)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Load and select samples -----------------------------------------------
obj  <- readRDS("results/preprocessed_GSE80655.rds")
meta <- obj$meta
raw  <- obj$dge$counts      # raw counts, genes already pre-filtered in Part 1

# ONE contrast, defined in advance: SZ vs CTL in AnCg.
# Within one region every donor occurs once, so samples are independent
# (this avoids the repeated-donor problem of the full dataset).
# The one donor with missing brain_pH (4463, SZ) drops out because pH is in the model.
sel <- meta$brain_region == "AnCg" &
       meta$diagnosis %in% c("Control", "Schizophrenia") &
       !is.na(meta$brain_pH)
m <- meta[sel, ]
m$diagnosis <- droplevels(factor(m$diagnosis, levels = c("Control", "Schizophrenia")))
m$sex <- droplevels(m$sex)
cat("Samples used:\n"); print(table(m$diagnosis, m$sex))   # check the female counts!

# ---- 2. Filter uninformative genes (WITHOUT using diagnosis) ----------------------
dge <- DGEList(counts = raw[, rownames(m)])
keep <- filterByExpr(dge, design = model.matrix(~1, data = m),
                     min.count = 10, min.total.count = 15)
dge <- dge[keep, , keep.lib.sizes = FALSE]
n_tested <- nrow(dge)
cat("Number of genes tested:", n_tested, "\n")            # REPORT THIS NUMBER

dge <- calcNormFactors(dge, method = "TMM")

# ---- 3. Model: voom + linear model with covariates ----------------------------------
design <- model.matrix(~ diagnosis + age + PMI + brain_pH + sex, data = m)
colnames(design)

# voomWithQualityWeights down-weights low-quality samples (useful for post-mortem tissue)
v   <- voomWithQualityWeights(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# ---- 4. Results table + multiple testing correction ------------------------------------
res <- topTable(fit, coef = "diagnosisSchizophrenia", number = Inf,
                sort.by = "P", adjust.method = "BH") %>%
  rownames_to_column("ensembl_id")

# Optional: gene symbols (the count file only has Ensembl IDs)
res$symbol <- tryCatch(
  AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = res$ensembl_id,
                        keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"),
  error = function(e) NA_character_)

cat("Genes with FDR < 0.05:", sum(res$adj.P.Val < 0.05), "\n")
cat("Genes with raw p < 0.05:", sum(res$P.Value < 0.05), "\n")
cat("Expected by chance at raw p < 0.05:", round(0.05 * n_tested), "\n")
write.csv(res, "results/DE_SZ_vs_CTL_AnCg.csv", row.names = FALSE)

# ---- 5. Visualisations ------------------------------------------------------------------------
# (a) P-value histogram: flat = no signal, spike near 0 = real signal
ggplot(res, aes(P.Value)) + geom_histogram(bins = 40, boundary = 0, colour = "white") +
  labs(x = "Raw p-value", y = "Number of genes")
ggsave("figures/pvalue_histogram.png", width = 5, height = 4)

# (b) Volcano plot
res <- res %>% mutate(sig = adj.P.Val < 0.05,
                      label = ifelse(rank(P.Value) <= 10, symbol, NA))
ggplot(res, aes(logFC, -log10(P.Value), colour = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_text_repel(aes(label = label), size = 3, na.rm = TRUE, colour = "black") +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"),
                      name = "FDR < 0.05") +
  labs(x = "log2 fold change (SZ vs CTL)", y = "-log10(p-value)")
ggsave("figures/volcano_SZ_vs_CTL_AnCg.png", width = 6, height = 5)

# ---- 6. Sanity checks ---------------------------------------------------------------------------
# The paper (Ramaker 2017) reports ~87 genes at FDR < 0.05 (DESeq2) and EGR1 down in SZ.
# Check the Ensembl IDs yourself before relying on them:
# EGR1 ENSG00000120738, GAD1 ENSG00000128683, GAD2 ENSG00000136750
res %>% filter(ensembl_id %in% c("ENSG00000120738", "ENSG00000128683", "ENSG00000136750"))

# Sensitivity analysis: does the result change without covariates?
design0 <- model.matrix(~ diagnosis, data = m)
fit0 <- eBayes(lmFit(voom(dge, design0), design0), robust = TRUE)
cat("FDR < 0.05 without covariates:",
    sum(topTable(fit0, coef = 2, number = Inf)$adj.P.Val < 0.05), "\n")

# ---- 7. Interpretation reminders ------------------------------------------------------------------
# - FDR < 0.05 means: among the genes called significant, we expect about 5% to be false
#   positives. It does NOT mean each gene has a 95% chance of being real.
# - Non-significant is not "no difference": the sample is small (~24 per group).
# - A different pipeline (e.g. DESeq2) will give a similar but not identical gene list.
# - Only association: post-mortem data cannot separate disease from medication or cause of death.
sessionInfo()
