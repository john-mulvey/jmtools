# Script to generate sample metadata based on real proteomics data
# Derives clinical variables correlated with proteomics patterns
#
# DATA SOURCE:
# Mulvey, J.F. et al. (2026). An unbiased molecular characterisation of
# peripartum cardiomyopathy hearts identifies mast cell chymase as a new
# diagnostic candidate. Molecular & Cellular Proteomics, Volume 0, Issue 0, 101510.
#
# The proteomics data (example_proteome.csv) is from the published study.
# This script generates synthetic sample metadata correlated with the proteomics
# patterns for demonstration purposes.

# ============================================================================
# CONFIGURATION
# ============================================================================

seed <- 42
set.seed(seed)

# ============================================================================
# Load proteomics data
# ============================================================================

cat("Reading proteomics data...\n")

proteomics <- read.csv("example_proteome.csv", check.names = FALSE)

# Extract sample columns (all columns except Gene symbol and Uniprot Identifier)
annotation_cols <- c("Gene symbol", "Uniprot Identifier")
sample_cols <- setdiff(colnames(proteomics), annotation_cols)

cat("Found", length(sample_cols), "samples:", paste(sample_cols, collapse = ", "), "\n")
cat("Found", nrow(proteomics), "proteins\n")

# Extract expression matrix (samples as columns, proteins as rows)
expr_matrix <- as.matrix(proteomics[, sample_cols])
rownames(expr_matrix) <- proteomics$`Gene symbol`

# Define groups based on sample names
groups <- ifelse(grepl("^Ctrl", sample_cols), "control", "heart_failure")
n_samples <- length(sample_cols)

cat("\nGroup distribution:\n")
print(table(groups))

# ============================================================================
# Run PCA to derive correlated metadata
# ============================================================================

cat("\nRunning PCA on proteomics data...\n")

# Transpose so samples are rows for PCA
pca_result <- prcomp(t(expr_matrix), scale. = TRUE, center = TRUE)
pc1_scores <- pca_result$x[, 1]
pc2_scores <- pca_result$x[, 2]

var_explained <- summary(pca_result)$importance["Proportion of Variance", 1:2] * 100
cat("PC1 explains", round(var_explained[1], 1), "% variance\n")
cat("PC2 explains", round(var_explained[2], 1), "% variance\n")

# ============================================================================
# Generate sample metadata correlated with proteomics
# ============================================================================

cat("\nGenerating sample metadata...\n")

# Sex: randomly assigned, balanced overall
sex <- sample(rep(c("Male", "Female"), length.out = n_samples))

# Ejection fraction: primarily determined by group, with within-group correlation to PC1
# Controls: EF ~55-60%, HF: EF ~18-28%
# PC1 correlates negatively with EF (lower PC1 = higher EF within each group)

ejection_fraction <- rep(NA, n_samples)

# Standardise PC1 within each group for consistent correlation
pc1_ctrl <- scale(pc1_scores[groups == "control"])[,1]
pc1_hf <- scale(pc1_scores[groups == "heart_failure"])[,1]

ctrl_idx <- 1
hf_idx <- 1
for (i in 1:n_samples) {
  if (groups[i] == "control") {
    # Controls have high EF (55-60 range)
    base_ef <- 57
    pc_effect <- -pc1_ctrl[ctrl_idx] * 1.5  # Negative correlation
    ejection_fraction[i] <- round(base_ef + pc_effect + rnorm(1, sd = 1))
    ctrl_idx <- ctrl_idx + 1
  } else {
    # HF patients have low EF (18-28 range)
    base_ef <- 23
    pc_effect <- -pc1_hf[hf_idx] * 2  # Negative correlation
    ejection_fraction[i] <- round(base_ef + pc_effect + rnorm(1, sd = 1.5))
    hf_idx <- hf_idx + 1
  }
  # Clamp to realistic range
  ejection_fraction[i] <- max(15, min(62, ejection_fraction[i]))
}

# Age: slight correlation with PC2 + noise (for demonstration of PC-metadata associations)
age <- round(48 + pc2_scores * 2 + rnorm(n_samples, sd = 6))
age <- pmax(35, pmin(65, age))

# BMI: random, no correlation needed
bmi <- round(rnorm(n_samples, mean = 27, sd = 3), 1)
bmi <- pmax(20, pmin(35, bmi))

# ============================================================================
# Create metadata data frame
# ============================================================================

sample_meta <- data.frame(
  sample_id = sample_cols,
  group = groups,
  Sex = sex,
  Age = age,
  BMI = bmi,
  `Ejection fraction (%)` = ejection_fraction,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Save metadata
write.csv(sample_meta, "sample_meta.csv", row.names = FALSE)

cat("\nSaved sample metadata with", nrow(sample_meta), "samples\n")

cat("\nMetadata summary:\n")
print(summary(sample_meta))

cat("\nEjection fraction by group:\n")
print(tapply(sample_meta$`Ejection fraction (%)`, sample_meta$group, summary))

cat("\n=== Complete ===\n")
cat("Created: sample_meta.csv\n")
cat("Use with: example_proteome.csv\n")
