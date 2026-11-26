setwd("C:/Users/karen/Documents/Ugent/Academiejaar 26-27/Semester 1/Large scale analysis/project/DATA")
library(tidyverse)
library(edgeR)
library(limma)
library(pheatmap)

counts <- read.csv("GSE80655_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("GSE80655_metadata.csv", stringsAsFactors = FALSE)

# Subset metadata to DLPFC (which() avoids NA rows)
meta_dlpfc <- meta[which(meta$brain_region == "DLPFC"), ]
nrow(meta_dlpfc)
table(meta_dlpfc$diagnosis)
n_distinct(meta_dlpfc$donor_id)

# Check that all sample IDs exist in the counts matrix
missing <- setdiff(meta_dlpfc$description, colnames(counts))
if (length(missing) > 0) stop("Not in counts: ", paste(head(missing), collapse = ", "))

# Subset counts (as a numeric matrix) in metadata order
counts_dlpfc <- as.matrix(counts[, meta_dlpfc$description])
stopifnot(identical(colnames(counts_dlpfc), meta_dlpfc$description))

# Normalising on all genes(library sieze and TMM) afterwards logCPM
dge <- DGEList(counts = counts_dlpfc)
dge <- calcNormFactors(dge)
logcpm_dlpfc <- cpm(dge, log = TRUE, prior.count = 1)

# 3. select 2 genes
genes <- c("ENSG00000128683", "ENSG00000157005")
ontbrekend <- setdiff(genes, rownames(logcpm_dlpfc))
if (length(ontbrekend) > 0) stop("Gene not found: ", paste(ontbrekend, collapse = ", "))

logcpm_4 <- logcpm_dlpfc[genes, , drop = FALSE]

# 4. have a control of the sexes
table(meta_dlpfc$sex, useNA = "ifany")

female_ids <- meta_dlpfc$description[which(meta_dlpfc$sex %in% c("F", "female", "Female"))]
male_ids   <- meta_dlpfc$description[which(meta_dlpfc$sex %in% c("M", "male", "Male"))]

# 5. Exporting the data
export <- function(ids, naam) {
  write.csv(logcpm_4[, ids, drop = FALSE], paste0("logcpm_4genen_", naam, ".csv"))
  write.csv(meta_dlpfc[meta_dlpfc$description %in% ids, ],
            paste0("meta_4genen_", naam, ".csv"), row.names = FALSE)
}

export(female_ids,                 "female")
export(male_ids,                   "male")
export(colnames(logcpm_4),        "alle")   # male + female samen

# Controle
c(female = length(female_ids), male = length(male_ids), alle = ncol(logcpm_4))

