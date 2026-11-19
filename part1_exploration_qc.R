# ============================================================
# LSABD group project - Part 1: Explore, preprocess, QC
# Dataset: GSE80655 (post-mortem brain RNA-seq, SZ / BPD / MDD / CTL)
# NOTE: starter script, written but not run on your machine. Check every step.
# ============================================================

# ---- 0. Packages -------------------------------------------------------
# install.packages(c("tidyverse", "pheatmap", "ggrepel"))
# BiocManager::install(c("edgeR", "limma"))
library(tidyverse)
library(edgeR)
library(limma)
library(pheatmap)

set.seed(42)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Load data ------------------------------------------------------
counts <- read.csv("GSE80655_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("GSE80655_metadata.csv", stringsAsFactors = FALSE)

cat("Counts:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata:", nrow(meta), "samples x", ncol(meta), "variables\n")

# The count columns (SL...) match meta$description, NOT meta$GSM or meta$title
stopifnot(setequal(colnames(counts), meta$description))
meta   <- meta[match(colnames(counts), meta$description), ]
rownames(meta) <- meta$description
stopifnot(identical(colnames(counts), rownames(meta)))

# ---- 2. Variable types and factors --------------------------------------
# nominal: diagnosis, brain_region, sex, ethnicity, donor_id (an identifier!)
# continuous: age, PMI, brain_pH
meta <- meta %>%
  mutate(
    diagnosis    = factor(diagnosis,
                          levels = c("Control", "Schizophrenia",
                                     "Bipolar Disorder", "Major Depression")),
    brain_region = factor(brain_region),
    sex          = factor(sex),
    ethnicity    = factor(ethnicity),
    donor_id     = factor(donor_id)
  )
rownames(meta) <- meta$description
str(meta)

# ---- 3. Design overview (tables, no figures for the report) -------------
table(meta$diagnosis, meta$brain_region)
table(meta$diagnosis, meta$sex)             # sex is strongly unbalanced
table(table(meta$donor_id))                 # samples per donor (3 or 2)
meta %>% group_by(diagnosis) %>%
  summarise(across(c(age, PMI, brain_pH), list(mean = ~mean(.x, na.rm = TRUE),
                                               sd = ~sd(.x, na.rm = TRUE))))

# ---- 4. Missing values ---------------------------------------------------
colSums(is.na(meta))
meta[is.na(meta$brain_pH), c("description", "brain_region", "diagnosis", "donor_id")]
# brain_pH is missing for ONE donor (4463, SZ) in all 3 regions.
# Options to document and justify: (a) drop the donor when pH is in the model,
# (b) impute with the median/regression, (c) check the paper (pH >= 6.5 was required).
sum(is.na(counts)); sum(counts < 0)          # counts: no NA / negatives expected

# ---- 5. Library size and detected genes -----------------------------------
meta$lib_size <- colSums(counts)
meta$n_detected <- colSums(counts > 0)

ggplot(meta, aes(reorder(description, lib_size), lib_size / 1e6, fill = brain_region)) +
  geom_col() + coord_flip() + theme(axis.text.y = element_blank()) +
  labs(x = "Sample", y = "Library size (million reads)")
ggsave("figures/library_sizes.png", width = 6, height = 5)

ggplot(meta, aes(brain_region, lib_size / 1e6, fill = diagnosis)) +
  geom_boxplot() + labs(y = "Library size (million reads)")

# ---- 6. Gene filtering (independent of diagnosis!) -------------------------
# Use design = ~1 / no group argument so the criterion does not use group labels.
dge <- DGEList(counts = counts)
cat("Genes before filtering:", nrow(dge), "\n")
cat("All-zero genes:", sum(rowSums(counts) == 0), "\n")

keep <- filterByExpr(dge, design = model.matrix(~1, data = meta),
                     min.count = 10, min.total.count = 15)
# Alternative, easy to explain: CPM > 1 in at least 20 samples
# keep <- rowSums(cpm(dge) > 1) >= 20
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes after filtering:", nrow(dge), "\n")     # REPORT THIS NUMBER

# ---- 7. Normalisation + log transformation ---------------------------------
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

# Before/after normalisation check
par(mfrow = c(1, 2))
boxplot(log2(counts[rownames(dge), ] + 1), las = 2, outline = FALSE,
        main = "Raw (log2 + 1)", names = NA)
boxplot(logcpm, las = 2, outline = FALSE, main = "TMM + logCPM", names = NA)
par(mfrow = c(1, 1))

# ---- 8. Sample outliers and clustering -----------------------------------------
pca <- prcomp(t(logcpm), scale. = TRUE)
var_expl <- summary(pca)$importance[2, ] * 100
pcs <- cbind(meta, pca$x[, 1:10])

plot_pca <- function(colour, shape = NULL) {
  ggplot(pcs, aes(PC1, PC2, colour = .data[[colour]],
                  shape = if (is.null(shape)) NULL else .data[[shape]])) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = sprintf("PC1 (%.1f%%)", var_expl[1]),
         y = sprintf("PC2 (%.1f%%)", var_expl[2]), colour = colour)
}
plot_pca("brain_region")   # expected: nAcc separates from AnCg/DLPFC
plot_pca("diagnosis")
plot_pca("sex")
plot_pca("brain_pH")
plot_pca("PMI")
plot_pca("lib_size")

# Which technical/biological variable drives which PC?
pc_assoc <- sapply(c("age", "PMI", "brain_pH", "lib_size"), function(v)
  sapply(paste0("PC", 1:10), function(pc) cor(pcs[[pc]], pcs[[v]], use = "complete.obs")))
round(pc_assoc, 2)
pheatmap(pc_assoc, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Correlation of PCs with continuous covariates")

# Since region dominates, repeat PCA WITHIN one region (e.g. AnCg)
anc <- meta$brain_region == "AnCg"
pca_anc <- prcomp(t(logcpm[, anc]), scale. = TRUE)
pcs_anc <- cbind(meta[anc, ], pca_anc$x[, 1:5])
ggplot(pcs_anc, aes(PC1, PC2, colour = diagnosis)) + geom_point(size = 2)

# Sample-to-sample distance heatmap (all samples)
d <- as.matrix(dist(t(logcpm)))
pheatmap(d, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = meta[, c("brain_region", "diagnosis")])

# ---- 9. Sex sanity check (possible sample swaps) ---------------------------------
# XIST = ENSG00000229807 (female-high); RPS4Y1 = ENSG00000129824 (male-high)
sex_genes <- c(XIST = "ENSG00000229807", RPS4Y1 = "ENSG00000129824")
sex_genes <- sex_genes[sex_genes %in% rownames(logcpm)]
sex_df <- data.frame(t(logcpm[sex_genes, , drop = FALSE]), sex = meta$sex,
                     description = meta$description)
if (all(c("XIST", "RPS4Y1") %in% colnames(sex_df))) {
  ggplot(sex_df, aes(XIST, RPS4Y1, colour = sex)) + geom_point()
}
# Look for samples where recorded sex disagrees with expression.

# ---- 10. Donor consistency ---------------------------------------------------------
# Each donor has up to 3 regions => samples are NOT independent.
# This matters for Parts 2-4: model per region, or account for donor_id,
# and split train/test by donor in Part 4.

# ---- 11. Save the preprocessed objects for the whole group ---------------------------
saveRDS(list(logcpm = logcpm, dge = dge, meta = meta[colnames(logcpm), ]),
        "results/preprocessed_GSE80655.rds")
write.csv(logcpm, "results/logCPM_filtered_TMM.csv")
write.csv(meta, "results/metadata_clean.csv", row.names = FALSE)

# Keep a record of the package versions you used (for the Methods section)
sessionInfo()
