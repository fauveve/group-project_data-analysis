###############################################################
# GAD1 - SST co-expression by diagnosis in the DLPFC
# Dataset: Females only
#
# Hypothesis: the association between GAD1 and SST expression
# differs according to psychiatric diagnosis.
#
# Input files (from the export script):
#   logcpm_4genen_vrouwen.csv
#   meta_4genen_vrouwen.csv
###############################################################

library(tidyverse)

setwd("C:/Users/karen/Documents/Ugent/Academiejaar 26-27/Semester 1/Large scale analysis/project/DATA")

# ---------------- SETTINGS ----------------
gene_x <- "GAD1"      # must match the rownames in the exported file
gene_y <- "SST"       # (use "gen1", "gen2", ... if the rownames were not renamed)
ref_level  <- NULL    # control label, e.g. "control"; NULL = first level alphabetically
covariates <- c()     # optional metadata columns, e.g. c("age", "RIN", "PMI")
label <- "Females"
tag   <- "vrouwen"
# ------------------------------------------

# ---- 1. Load data ----
expr <- read.csv(paste0("logcpm_4genen_female", tag, ".csv"), row.names = 1, check.names = FALSE)
meta <- read.csv(paste0("meta_4genen_female", tag, ".csv"), stringsAsFactors = FALSE)
stopifnot(all(c(gene_x, gene_y) %in% rownames(expr)))

df <- data.frame(description = colnames(expr),
                 x = as.numeric(expr[gene_x, ]),
                 y = as.numeric(expr[gene_y, ])) %>%
  left_join(meta, by = "description") %>%
  mutate(diagnosis = factor(diagnosis))
if (!is.null(ref_level)) df$diagnosis <- relevel(df$diagnosis, ref = ref_level)

cat("Samples per diagnosis (", label, "):\n")
print(table(df$diagnosis))

# ---- 2. Plots ----
p_box <- df %>%
  pivot_longer(c(x, y), names_to = "gene", values_to = "logCPM") %>%
  mutate(gene = recode(gene, x = gene_x, y = gene_y)) %>%
  ggplot(aes(diagnosis, logCPM, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.2) +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() + theme(legend.position = "none") +
  labs(title = paste("DLPFC expression -", label))

p_scatter <- ggplot(df, aes(x, y, colour = diagnosis)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  theme_bw() +
  labs(title = paste(gene_x, "vs", gene_y, "-", label),
       x = paste(gene_x, "(logCPM)"), y = paste(gene_y, "(logCPM)"))

p_scatter_facet <- p_scatter + facet_wrap(~diagnosis) + theme(legend.position = "none")

print(p_box); print(p_scatter); print(p_scatter_facet)
ggsave(paste0("boxplot_", label, ".png"), p_box, width = 8, height = 5, dpi = 300)
ggsave(paste0("scatter_", label, ".png"), p_scatter, width = 7, height = 5, dpi = 300)
ggsave(paste0("scatter_by_diagnosis_", label, ".png"), p_scatter_facet, width = 9, height = 6, dpi = 300)

# ---- 3. Does the GAD1-SST association differ by diagnosis? ----
cov_term <- if (length(covariates) > 0) paste("+", paste(covariates, collapse = "+")) else ""
m0 <- lm(as.formula(paste("y ~ x + diagnosis", cov_term)), data = df)   # common slope
m1 <- lm(as.formula(paste("y ~ x * diagnosis", cov_term)), data = df)   # slope per diagnosis
interaction_test <- tibble(group = label,
                           test = paste(gene_x, "x diagnosis interaction on", gene_y),
                           p_value = anova(m0, m1)$`Pr(>F)`[2])
cat("\nInteraction test:\n"); print(interaction_test)
cat("\nSlope per diagnosis model:\n"); print(summary(m1)$coefficients)

# ---- 4. Correlation per diagnosis ----
cors <- map_dfr(levels(df$diagnosis), function(dg) {
  d  <- filter(df, diagnosis == dg)
  ct <- if (nrow(d) > 3) cor.test(d$x, d$y) else NULL
  tibble(group = label, diagnosis = dg, n = nrow(d),
         pearson_r    = cor(d$x, d$y),
         spearman_rho = cor(d$x, d$y, method = "spearman"),
         p_value = if (!is.null(ct)) ct$p.value else NA,
         ci_low  = if (!is.null(ct)) ct$conf.int[1] else NA,
         ci_high = if (!is.null(ct)) ct$conf.int[2] else NA)
})
cat("\nCorrelations per diagnosis:\n"); print(cors)

p_cor <- ggplot(cors, aes(diagnosis, pearson_r, fill = diagnosis)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, na.rm = TRUE) +
  geom_hline(yintercept = 0) +
  theme_bw() + theme(legend.position = "none") +
  labs(title = paste("Correlation", gene_x, "-", gene_y, "-", label),
       y = "Pearson r (95% CI)")
print(p_cor)
ggsave(paste0("correlation_", label, ".png"), p_cor, width = 6, height = 5, dpi = 300)

# ---- 5. Each diagnosis vs reference: Fisher z-test on correlations ----
fisher_p <- function(r1, n1, r2, n2) {
  z <- (atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  2 * pnorm(-abs(z))
}
ref_dg  <- levels(df$diagnosis)[1]
ref_row <- filter(cors, diagnosis == ref_dg)
fisher <- cors %>%
  filter(diagnosis != ref_dg, n > 3, ref_row$n > 3) %>%
  mutate(vs_reference = ref_dg,
         r_reference  = ref_row$pearson_r,
         p_fisher = fisher_p(ref_row$pearson_r, ref_row$n, pearson_r, n),
         fdr = p.adjust(p_fisher, "BH"))
cat("\nFisher z-test vs reference (", ref_dg, "):\n"); print(fisher)

# ---- 6. Export results ----
write.csv(interaction_test, paste0("result_interaction_", label, ".csv"), row.names = FALSE)
write.csv(cors,   paste0("result_correlations_", label, ".csv"), row.names = FALSE)
write.csv(fisher, paste0("result_fisher_", label, ".csv"), row.names = FALSE)