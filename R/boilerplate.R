# Boilerplate Code Generators
# jmtools package
#
# These functions generate template code that can be copied into analysis
# scripts and adapted as needed.

#' Generate PCA boilerplate code
#'
#' Outputs boilerplate R code for performing PCA on a matrix with features
#' (e.g., proteins, genes) as rows and samples as columns. The generated code
#' includes PCA computation, scores plot, screeplot with cumulative variance,
#' and pairwise component scatterplots.
#'
#' @param data_name Character string specifying the name of the data matrix
#'   variable in the generated code.
#' @param metadata_name Character string specifying the name of the metadata
#'   data frame variable.
#' @param colour_var Character string specifying the metadata column to use
#'   for colouring points (default: "group").
#' @param id_col Character string specifying the column in metadata containing
#'   sample identifiers that match the column names of the data matrix
#'   (default: "sample_id").
#' @param n_pcs_scree Number of principal components to show in the screeplot
#'   (default: 10).
#' @param n_pcs_pairs Number of principal components to include in the ggpairs
#'   plot (default: 5).
#' @param scale Logical indicating whether to scale variables to unit variance
#'   in the generated code (default: TRUE).
#' @param center Logical indicating whether to centre variables in the
#'   generated code (default: TRUE).
#'
#' @return Returns `invisible(NULL)`. The function prints code to the console
#'   which can be copied into an analysis script.
#'
#' @details
#' The generated code uses `prcomp()` from base R for PCA computation, avoiding
#' additional package dependencies. The code includes:
#'
#' \enumerate{
#'   \item Data preparation (transposing so samples are rows)
#'   \item PCA computation with `prcomp()`
#'   \item Scores plot (PC1 vs PC2) coloured by a metadata variable
#'   \item Screeplot showing per-component variance (bars) and cumulative
#'     variance (points/line on secondary y-axis)
#'   \item Pairwise scatterplot matrix of the first N components using
#'     `GGally::ggpairs()`
#' }
#'
#' @export
#'
#' @examples
#' # Generate boilerplate for proteomics analysis
#' generate_pca_boilerplate(
#'   data_name = "abundance",
#'   metadata_name = "sample_meta"
#' )
#'
#' # Customise variable names
#' generate_pca_boilerplate(
#'   data_name = "protein_matrix",
#'   metadata_name = "clinical_data",
#'   colour_var = "treatment"
#' )
generate_pca_boilerplate <- function(data_name,
                                     metadata_name,
                                     colour_var = "group",
                                     id_col = "sample_id",
                                     n_pcs_scree = 10,
                                     n_pcs_pairs = 5,
                                     scale = TRUE,
                                     center = TRUE) {

  # Generate column selection string for ggpairs
  pairs_cols <- paste0("PC", seq_len(n_pcs_pairs), collapse = '", "')

  code <- glue::glue('
library(ggplot2)
library(GGally)

# Transpose so samples are rows, run PCA
pca_input <- t({data_name})
pca_result <- prcomp(pca_input, scale. = {toupper(scale)}, center = {toupper(center)})

# Extract scores and combine with metadata
pca_scores <- as.data.frame(pca_result$x)
pca_scores${id_col} <- rownames(pca_scores)
pca_data <- merge(pca_scores, {metadata_name}, by = "{id_col}")

# Calculate variance explained
var_explained <- summary(pca_result)$importance
pct_var <- var_explained["Proportion of Variance", ] * 100
cumulative_var <- var_explained["Cumulative Proportion", ] * 100

# Scores plot (PC1 vs PC2)
ggplot(pca_data, aes(x = PC1, y = PC2, colour = {colour_var})) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(pct_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pct_var[2], 1), "%)"),
    title = "PCA Scores Plot"
  ) +
  theme_jm()

# Screeplot with cumulative variance
n_pcs_plot <- min({n_pcs_scree}, length(pct_var))
scree_data <- data.frame(
  PC = factor(paste0("PC", seq_len(n_pcs_plot)), levels = paste0("PC", seq_len(n_pcs_plot))),
  variance = pct_var[seq_len(n_pcs_plot)],
  cumulative = cumulative_var[seq_len(n_pcs_plot)]
)
scale_factor <- 100 / max(scree_data$variance)

ggplot(scree_data, aes(x = PC)) +
  geom_col(aes(y = variance), fill = "steelblue", alpha = 0.8) +
  geom_point(aes(y = cumulative / scale_factor), colour = "darkred", size = 3) +
  geom_line(aes(y = cumulative / scale_factor, group = 1), colour = "darkred", linewidth = 1) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~ . * scale_factor, name = "Cumulative Variance (%)")
  ) +
  labs(x = "Principal Component", title = "PCA Screeplot") +
  theme_jm() +
  theme(
    axis.title.y.right = element_text(colour = "darkred"),
    axis.text.y.right = element_text(colour = "darkred")
  )

# Pairwise component scatterplots
ggpairs(
  pca_data,
  columns = c("{pairs_cols}"),
  mapping = aes(colour = {colour_var}),
  lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
  upper = list(continuous = wrap("cor", size = 3)),
  progress = FALSE
) +
  theme_jm()
')

  cat(code)
  return(invisible(NULL))
}


#' Generate limma boilerplate code
#'
#' Outputs boilerplate R code for performing differential expression/abundance
#' analysis using limma. The generated code includes design matrix creation,
#' model fitting, contrast specification, and results extraction with
#' snake_case column names.
#'
#' @param data_name Character string specifying the name of the expression/abundance
#'   matrix variable (features as rows, samples as columns).
#' @param metadata_name Character string specifying the name of the metadata
#'   data frame variable.
#' @param id_col Character string specifying the column in metadata containing
#'   sample identifiers that match the column names of the data matrix.
#'   Default: "sample_id".
#' @param formula Character string specifying the model formula for the design
#'   matrix. Use `"~ 0 + group"` for a no-intercept model where coefficients
#'   represent group means (recommended for contrasts). Use `"~ group"` for an
#'   intercept model where coefficients represent differences from the reference
#'   level. Default: `"~ 0 + group"`.
#' @param contrast Character string specifying the contrast to test. Should be
#'   in the form `"groupA - groupB"` matching the column names in the design
#'   matrix. For `~ 0 + group` formulas, use the factor levels directly
#'   (e.g., `"treatment - control"`). For `~ group` formulas with intercept,
#'   column names will be like `grouptreatment`. Default: `"treatment - control"`.
#'   Set to `NULL` to skip pairwise contrast testing (only run ANOVA if
#'   `include_anova = TRUE`).
#' @param output_name Character string specifying the name for the output
#'   results object. If `NULL` (default), auto-generates from the contrast
#'   (e.g., `"heart_failure - control"` becomes
#'   `"limma_results_heart_failure_vs_control"`).
#' @param feature_id_col Character string specifying the name to use for the
#'   feature identifier column in results. Default: `"feature_id"`.
#' @param fc_threshold Default fold-change threshold for the generated code.
#'   Default: 1 (i.e., 2-fold change on linear scale).
#' @param adj_p_threshold Default adjusted p-value threshold for the generated
#'   code. Default: 0.05.
#' @param include_anova Logical. If `TRUE`, includes code for an F-test
#'   (ANOVA) testing for any differences between groups. Default: `FALSE`.
#'
#' @return Returns `invisible(NULL)`. The function prints code to the console
#'   which can be copied into an analysis script.
#'
#' @details
#' The generated code follows a standard limma workflow:
#'
#' \enumerate{
#'   \item Define significance thresholds (`fc_threshold`, `adj_p_threshold`)
#'   \item Create design matrix from the specified formula
#'   \item Fit linear model with `limma::lmFit()`
#'   \item Define and fit contrasts with `limma::makeContrasts()` and
#'     `limma::contrasts.fit()`
#'   \item Compute moderated statistics with `limma::eBayes()`
#'   \item Extract results with `limma::topTable()`
#'   \item Clean column names to snake_case using `janitor::clean_names()`
#'   \item Add `diff_abundant` boolean column based on thresholds
#' }
#'
#' If `include_anova = TRUE`, additional code is generated to perform an F-test
#' comparing all groups against a reference level (the first factor level).
#' This tests the null hypothesis that all group means are equal.
#'
#' Column names in the output are converted to snake_case, with `adj.P.Val`
#' renamed to `adj_p_value` for consistency.
#'
#' @export
#'
#' @examples
#' # Generate boilerplate for basic analysis
#' generate_limma_boilerplate(
#'   data_name = "abundance",
#'   metadata_name = "sample_meta"
#' )
#'
#' # Heart failure vs control example
#' generate_limma_boilerplate(
#'   data_name = "abundance",
#'   metadata_name = "sample_meta",
#'   formula = "~ 0 + group",
#'   contrast = "heart_failure - control"
#' )
#'
#' # Include ANOVA F-test for overall group differences
#' generate_limma_boilerplate(
#'   data_name = "abundance",
#'   metadata_name = "sample_meta",
#'   formula = "~ 0 + treatment",
#'   contrast = "drug - placebo",
#'   include_anova = TRUE
#' )
generate_limma_boilerplate <- function(data_name,
                                       metadata_name,
                                       id_col = "sample_id",
                                       formula = "~ 0 + group",
                                       contrast = "treatment - control",
                                       output_name = NULL,
                                       feature_id_col = "feature_id",
                                       fc_threshold = 1,
                                       adj_p_threshold = 0.05,
                                       include_anova = FALSE) {

  # Validate: need at least one of contrast or anova
  if (is.null(contrast) && !include_anova) {
    stop("At least one of 'contrast' or 'include_anova = TRUE' must be specified")
  }

  # Generate output name from contrast if not provided
  if (is.null(output_name) && !is.null(contrast)) {
    # Convert "heart_failure - control" to "limma_results_heart_failure_vs_control"
    output_name <- contrast |>
      gsub(" ", "", x = _) |>
      gsub("-", "_vs_", x = _) |>
      paste0("limma_results_", x = _)
  } else if (is.null(output_name)) {
    output_name <- "limma_results"
  }

  # Parse formula to extract variable name for factor conversion
  formula_vars <- all.vars(stats::as.formula(formula))
  group_var <- formula_vars[1]  # Primary grouping variable

  # Determine if formula has intercept
  has_intercept <- !grepl("~ *0 *\\+|~ *-1", formula)

  # Generate colnames cleanup code for no-intercept models
  colnames_code <- if (has_intercept) {
    ""
  } else {
    glue::glue('colnames(design) <- gsub("^{group_var}", "", colnames(design))')
  }

  # Build code sections
  header_code <- glue::glue('
library(limma)
library(janitor)
library(dplyr)
library(ggplot2)

# Match samples between data and metadata
stopifnot(
  "Metadata ID column not found" = "{id_col}" %in% colnames({metadata_name}),
  "Sample IDs must match between data and metadata" =
    all(colnames({data_name}) %in% {metadata_name}${id_col})
)
{metadata_name} <- {metadata_name}[match(colnames({data_name}), {metadata_name}${id_col}), ]
stopifnot(all(colnames({data_name}) == {metadata_name}${id_col}))

# Significance thresholds
fc_threshold <- {fc_threshold}
adj_p_threshold <- {adj_p_threshold}

# Create design matrix
{metadata_name}${group_var} <- as.factor({metadata_name}${group_var})
design <- model.matrix({formula}, data = {metadata_name})
{colnames_code}
head(design)

# Fit linear model
fit <- lmFit({data_name}, design)

# SA plot: check if trend = TRUE is needed in eBayes()
fit_eb <- eBayes(fit)
sa_data <- data.frame(average = fit_eb$Amean, sigma = sqrt(fit_eb$s2.post))

ggplot(sa_data, aes(x = average, y = sigma)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", colour = "darkred", se = FALSE) +
  geom_hline(yintercept = sqrt(fit_eb$s2.prior), colour = "steelblue", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Average log-expression",
    y = "Sqrt(posterior variance)",
    title = "SA Plot: Check for Mean-Variance Trend",
    subtitle = "Blue dashed = prior; Red = loess trend. If trend is strong, set use_trend <- TRUE"
  ) +
  theme_jm()

# Set to TRUE if SA plot shows a clear mean-variance trend
use_trend <- FALSE
')

  # Contrast section (optional)
  if (!is.null(contrast)) {
    contrast_code <- glue::glue('

# Pairwise contrast: {contrast}
contrast_matrix <- makeContrasts(contrast_of_interest = {contrast}, levels = design)
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts, trend = use_trend)

# Extract and clean results
{output_name} <- topTable(fit_contrasts, coef = "contrast_of_interest", number = Inf, sort.by = "P") |>
  tibble::rownames_to_column("{feature_id_col}") |>
  janitor::clean_names() |>
  dplyr::rename(adj_p_value = adj_p_val) |>
  dplyr::mutate(diff_abundant = abs(log_fc) >= fc_threshold & adj_p_value < adj_p_threshold)

# Summary
message("Pairwise contrast results ({contrast}):")
message("- Total features tested: ", nrow({output_name}))
message("- Significant (|log2FC| >= ", fc_threshold, " & adj_p < ", adj_p_threshold, "): ", sum({output_name}$diff_abundant))
message("- Up-regulated: ", sum({output_name}$diff_abundant & {output_name}$log_fc > 0))
message("- Down-regulated: ", sum({output_name}$diff_abundant & {output_name}$log_fc < 0))

head({output_name})
')
  } else {
    contrast_code <- ""
  }

  # ANOVA section (optional)
  if (include_anova) {
    anova_code <- glue::glue('

# ANOVA F-test (overall group differences)
group_levels <- levels({metadata_name}${group_var})
reference_level <- group_levels[1]
other_levels <- group_levels[-1]

anova_contrast_exprs <- paste0(other_levels, " - ", reference_level)
names(anova_contrast_exprs) <- paste0(other_levels, "_vs_", reference_level)

contrast_matrix_anova <- makeContrasts(contrasts = anova_contrast_exprs, levels = design)
fit_anova <- contrasts.fit(fit, contrast_matrix_anova)
fit_anova <- eBayes(fit_anova, trend = use_trend)

# F-test results
results_anova <- topTable(fit_anova, number = Inf, sort.by = "F") |>
  tibble::rownames_to_column("{feature_id_col}") |>
  janitor::clean_names() |>
  dplyr::rename(adj_p_value = adj_p_val) |>
  dplyr::mutate(significant = adj_p_value < adj_p_threshold)

message("ANOVA F-test results (any group difference):")
message("- Total features tested: ", nrow(results_anova))
message("- Significant (adj_p < ", adj_p_threshold, "): ", sum(results_anova$significant))

head(results_anova)
')
  } else {
    anova_code <- ""
  }

  # Combine and output
  full_code <- paste0(header_code, contrast_code, anova_code)
  cat(full_code)
  return(invisible(NULL))
}


#' Generate limma results plotting boilerplate code
#'
#' Outputs boilerplate R code for visualising limma differential expression/
#' abundance results. Includes volcano plot, MA plot, and p-value histogram.
#'
#' @param results_name Character string specifying the name of the results
#'   data frame variable (output from limma analysis).
#' @param feature_col Character string specifying the column containing feature
#'   identifiers (e.g., gene symbols, protein IDs). Default: "feature_id".
#' @param logfc_col Character string specifying the column containing log
#'   fold-change values. Default: "log_fc".
#' @param pval_col Character string specifying the column containing raw
#'   p-values. Default: "p_value".
#' @param avg_expr_col Character string specifying the column containing
#'   average expression values. Default: "ave_expr".
#' @param signif_col Character string specifying the column containing the

#'   logical significance flag. Default: "diff_abundant".
#' @param title Character string for plot titles. Default: "Differential Abundance".
#'
#' @return Returns `invisible(NULL)`. The function prints code to the console
#'   which can be copied into an analysis script.
#'
#' @details
#' The generated code creates three diagnostic plots:
#'
#' \enumerate{
#'   \item **Volcano plot**: -log10(p-value) vs log fold-change, with
#'     significant features highlighted and optional gene labels
#'   \item **MA plot**: Log fold-change vs average expression, useful for
#'     checking intensity-dependent bias
#'   \item **P-value histogram**: Distribution of raw p-values, which should
#'     be uniform under the null with enrichment near zero
#' }
#'
#' The code includes a `genes_to_label` vector that can be customised to
#' annotate specific features of interest on the volcano and MA plots.
#'
#' @export
#'
#' @examples
#' # Generate plotting boilerplate
#' generate_limma_plots_boilerplate(
#'   results_name = "limma_results_treatment_vs_control"
#' )
#'
#' # Customise for specific column names
#' generate_limma_plots_boilerplate(
#'   results_name = "de_results",
#'   feature_col = "protein_id",
#'   title = "PPCM vs Control"
#' )
generate_limma_plots_boilerplate <- function(results_name,
                                              feature_col = "feature_id",
                                              logfc_col = "log_fc",
                                              pval_col = "p_value",
                                              avg_expr_col = "ave_expr",
                                              signif_col = "diff_abundant",
                                              title = "Differential Abundance") {

  code <- glue::glue('
library(ggplot2)
library(ggrepel)
library(ggpointdensity)

# Features to label (customise as needed)
genes_to_label <- c(
  # "GENE1",
  # "GENE2"
)

# Volcano plot
sig_pvals <- {results_name} |>
  dplyr::filter({signif_col} == TRUE) |>
  dplyr::pull({pval_col})
highest_sig_pval <- if (length(sig_pvals) > 0) max(sig_pvals, na.rm = TRUE) else 1

ggplot({results_name}, aes(x = {logfc_col}, y = -log10({pval_col}))) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = -log10(highest_sig_pval), linetype = "dashed", colour = "grey60") +
  geom_pointdensity(size = 1, alpha = 0.7, show.legend = FALSE) +
  scale_colour_gradient(low = "#6baed6", high = "#08306b") +
  geom_point(
    data = . %>% dplyr::filter({signif_col} == TRUE),
    colour = "red", size = 1.5, alpha = 0.8
  ) +
  geom_text_repel(
    data = . %>% dplyr::filter({signif_col} == TRUE | {feature_col} %in% genes_to_label),
    aes(label = {feature_col}),
    size = 3, max.overlaps = 20, min.segment.length = 0, show.legend = FALSE
  ) +
  labs(
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~p-value),
    title = "{title}: Volcano Plot"
  ) +
  theme_jm()

# MA plot
ggplot({results_name}, aes(x = {avg_expr_col}, y = {logfc_col})) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", colour = "grey60") +
  geom_pointdensity(size = 1, alpha = 0.7, show.legend = FALSE) +
  scale_colour_gradient(low = "#6baed6", high = "#08306b") +
  geom_point(
    data = . %>% dplyr::filter({signif_col} == TRUE),
    colour = "red", size = 1.5, alpha = 0.8
  ) +
  geom_text_repel(
    data = . %>% dplyr::filter({signif_col} == TRUE | {feature_col} %in% genes_to_label),
    aes(label = {feature_col}),
    size = 3, max.overlaps = 20, min.segment.length = 0, show.legend = FALSE
  ) +
  labs(
    x = expression(log[2]~mean~expression),
    y = expression(log[2]~fold~change),
    title = "{title}: MA Plot"
  ) +
  theme_jm()

# P-value histogram
ggplot({results_name}, aes(x = {pval_col})) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "steelblue", colour = "white") +
  labs(
    x = "P-value",
    y = "Count",
    title = "{title}: P-value Distribution",
    subtitle = "Should be uniform with enrichment near 0; U-shape indicates problems"
  ) +
  theme_jm()
')

  cat(code)
  return(invisible(NULL))
}


#' Generate GSEA boilerplate code
#'
#' Outputs boilerplate R code for performing Gene Set Enrichment Analysis (GSEA)
#' using fgsea with gene sets from msigdbr. Uses limma t-statistics as the
#' ranking metric for preranked GSEA.
#'
#' @param results_name Character string specifying the name of the limma results
#'   data frame variable.
#' @param output_name Character string specifying the name for the output
#'   results data frame.
#' @param feature_col Character string specifying the column containing gene
#'   symbols for matching to gene sets. Default: "gene_symbol".
#' @param stat_col Character string specifying the column containing the
#'   t-statistic (or equivalent) for ranking. Default: "t".
#' @param species Character string specifying the species for msigdbr. Common
#'   values: "Homo sapiens", "Mus musculus". Default: "Homo sapiens".
#' @param category Character string specifying the msigdbr category. Common
#'   values: "H" (hallmark), "C2" (curated), "C5" (ontology), "C7" (immunologic).
#'   Default: "H".
#' @param subcategory Character string specifying the msigdbr subcategory, or
#'   `NULL` for all subcategories. Examples: "CP:KEGG", "CP:REACTOME",
#'   "GO:BP", "GO:MF", "GO:CC". Default: `NULL`.
#' @param min_size Minimum gene set size for fgsea filtering. Default: 15.
#' @param max_size Maximum gene set size for fgsea filtering. Default: 500.
#'
#' @return Returns `invisible(NULL)`. The function prints code to the console
#'   which can be copied into an analysis script.
#'
#' @details
#' The generated code follows a standard fgsea workflow:
#'
#' \enumerate{
#'   \item Retrieve gene sets from msigdbr for the specified species and category
#'   \item Convert to fgsea-compatible list format
#'   \item Create a named vector of ranking statistics (t-values) from limma results
#'   \item Run fgsea with the `fgsea()` function
#'   \item Clean column names to snake_case with consistent naming (p_value,
#'     adj_p_value, etc.)
#'   \item Generate summary statistics
#' }
#'
#' The t-statistic from limma is recommended for ranking as it incorporates
#' both the magnitude of change (log fold-change) and the precision of the
#' estimate (standard error), providing a more robust ranking than log
#' fold-change alone.
#'
#' @section msigdbr Categories:
#' \describe{
#'   \item{H}{Hallmark gene sets (50 sets, broadly applicable)}
#'   \item{C2}{Curated gene sets (includes KEGG, Reactome, BioCarta, etc.)}
#'   \item{C5}{Ontology gene sets (GO Biological Process, Molecular Function,
#'     Cellular Component)}
#'   \item{C7}{Immunologic signature gene sets}
#'   \item{C8}{Cell type signature gene sets}
#' }
#'
#' @export
#'
#' @examples
#' # Generate boilerplate for hallmark gene sets (human)
#' generate_gsea_boilerplate(
#'   results_name = "limma_results_treatment_vs_control",
#'   output_name = "gsea_hallmark"
#' )
#'
#' # Use Reactome pathways
#' generate_gsea_boilerplate(
#'   results_name = "de_results",
#'   output_name = "gsea_reactome",
#'   category = "C2",
#'   subcategory = "CP:REACTOME"
#' )
generate_gsea_boilerplate <- function(results_name,
                                       output_name,
                                       feature_col = "gene_symbol",
                                       stat_col = "t",
                                       species = "Homo sapiens",
                                       category = "H",
                                       subcategory = NULL,
                                       min_size = 15,
                                       max_size = 500) {

  # Build subcategory filter code
  if (is.null(subcategory)) {
    subcategory_code <- ""
    subcategory_desc <- paste0("all subcategories in ", category)
  } else {
    subcategory_code <- glue::glue(',
  subcategory = "{subcategory}"')
    subcategory_desc <- subcategory
  }

  code <- glue::glue('
library(msigdbr)
library(fgsea)
library(dplyr)

# Get gene sets from msigdbr
gene_sets_df <- msigdbr(species = "{species}", category = "{category}"{subcategory_code})
gene_sets <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)

message("Loaded ", length(gene_sets), " gene sets")
message("Gene set sizes: ", min(lengths(gene_sets)), "-", max(lengths(gene_sets)))

# Prepare ranking vector (remove NA and duplicates, keep highest |stat|)
ranks_df <- {results_name} |>
  dplyr::filter(!is.na({stat_col}), !is.na({feature_col})) |>
  dplyr::group_by({feature_col}) |>
  dplyr::slice_max(abs({stat_col}), n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

ranks <- setNames(ranks_df${stat_col}, ranks_df${feature_col})
ranks <- sort(ranks, decreasing = TRUE)

message("Ranking vector: ", length(ranks), " genes")
message("Range: ", round(min(ranks), 2), " to ", round(max(ranks), 2))

# Run fgsea
set.seed(42)
{output_name} <- fgsea(pathways = gene_sets, stats = ranks, minSize = {min_size}, maxSize = {max_size}) |>
  tibble::as_tibble() |>
  dplyr::rename(
    gene_set = pathway, p_value = pval, adj_p_value = padj, log2_err = log2err,
    enrichment_score = ES, nes = NES, n_genes = size, leading_edge = leadingEdge
  ) |>
  dplyr::arrange(p_value)

# Summary
message("GSEA Results Summary:")
message("- Gene sets tested: ", nrow({output_name}))
message("- Significant (adj_p < 0.05): ", sum({output_name}$adj_p_value < 0.05, na.rm = TRUE))
message("- Enriched (NES > 0, adj_p < 0.05): ", sum({output_name}$nes > 0 & {output_name}$adj_p_value < 0.05, na.rm = TRUE))
message("- Depleted (NES < 0, adj_p < 0.05): ", sum({output_name}$nes < 0 & {output_name}$adj_p_value < 0.05, na.rm = TRUE))

head({output_name})
')

  cat(code)
  return(invisible(NULL))
}


#' Generate ORA boilerplate code
#'
#' Outputs boilerplate R code for performing Over-Representation Analysis (ORA)
#' using fgsea::fora() with gene sets from msigdbr. Uses the hypergeometric test
#' to identify enriched gene sets.
#'
#' @param gene_list_name Character string specifying the name of the vector
#'   containing genes of interest (e.g., significant genes).
#' @param background_name Character string specifying the name of the vector
#'   containing background genes (e.g., all measured genes).
#' @param output_name Character string specifying the name for the output
#'   results data frame.
#' @param species Character string specifying the species for msigdbr. Common
#'   values: "Homo sapiens", "Mus musculus". Default: "Homo sapiens".
#' @param category Character string specifying the msigdbr category. Common
#'   values: "H" (hallmark), "C2" (curated), "C5" (ontology), "C7" (immunologic).
#'   Default: "H".
#' @param subcategory Character string specifying the msigdbr subcategory, or
#'   `NULL` for all subcategories. Examples: "CP:KEGG", "CP:REACTOME",
#'   "GO:BP", "GO:MF", "GO:CC". Default: `NULL`.
#' @param min_size Minimum gene set size for filtering. Default: 15.
#' @param max_size Maximum gene set size for filtering. Default: 500.
#'
#' @return Returns `invisible(NULL)`. The function prints code to the console
#'   which can be copied into an analysis script.
#'
#' @details
#' The generated code follows a standard fgsea ORA workflow:
#'
#' \enumerate{
#'   \item Retrieve gene sets from msigdbr for the specified species and category
#'   \item Convert to fgsea-compatible list format
#'   \item Run `fgsea::fora()` with hypergeometric test
#'   \item Clean column names to snake_case with consistent naming (p_value,
#'     adj_p_value, etc.)
#'   \item Add an `enrichment_score` column (log2 odds ratio) for compatibility
#'     with plotting functions
#'   \item Generate summary statistics
#' }
#'
#' The output includes an `enrichment_score` column (log2 of the odds ratio)
#' to allow use with the same clustering and plotting functions as GSEA results
#' ([cluster_enrichment_results()], [plot_enrichment_bars()]).
#'
#' @section msigdbr Categories:
#' \describe{
#'   \item{H}{Hallmark gene sets (50 sets, broadly applicable)}
#'   \item{C2}{Curated gene sets (includes KEGG, Reactome, BioCarta, etc.)}
#'   \item{C5}{Ontology gene sets (GO Biological Process, Molecular Function,
#'     Cellular Component)}
#'   \item{C7}{Immunologic signature gene sets}
#'   \item{C8}{Cell type signature gene sets}
#' }
#'
#' @export
#'
#' @examples
#' # Generate boilerplate for hallmark gene sets (human)
#' generate_ora_boilerplate(
#'   gene_list_name = "upregulated_genes",
#'   background_name = "all_detected_genes",
#'   output_name = "ora_hallmark"
#' )
#'
#' # Use Reactome pathways
#' generate_ora_boilerplate(
#'   gene_list_name = "sig_genes",
#'   background_name = "background_genes",
#'   output_name = "ora_reactome",
#'   category = "C2",
#'   subcategory = "CP:REACTOME"
#' )
generate_ora_boilerplate <- function(gene_list_name,
                                      background_name,
                                      output_name,
                                      species = "Homo sapiens",
                                      category = "H",
                                      subcategory = NULL,
                                      min_size = 15,
                                      max_size = 500) {

  # Build subcategory filter code
  if (is.null(subcategory)) {
    subcategory_code <- ""
    subcategory_desc <- paste0("all subcategories in ", category)
  } else {
    subcategory_code <- glue::glue(',
  subcategory = "{subcategory}"')
    subcategory_desc <- subcategory
  }

  code <- glue::glue('
library(msigdbr)
library(fgsea)
library(dplyr)

# Get gene sets from msigdbr
gene_sets_df <- msigdbr(species = "{species}", category = "{category}"{subcategory_code})
gene_sets <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)

message("Loaded ", length(gene_sets), " gene sets")
message("Gene set sizes: ", min(lengths(gene_sets)), "-", max(lengths(gene_sets)))

# Ensure genes of interest are in the background
{gene_list_name} <- intersect({gene_list_name}, {background_name})

message("Genes of interest: ", length({gene_list_name}))
message("Background genes: ", length({background_name}))

# Run ORA
{output_name} <- fora(
  pathways = gene_sets, genes = {gene_list_name}, universe = {background_name},
  minSize = {min_size}, maxSize = {max_size}
) |>
  tibble::as_tibble() |>
  dplyr::rename(
    gene_set = pathway, p_value = pval, adj_p_value = padj,
    n_overlap = overlap, n_gene_set = size, overlap_genes = overlapGenes
  ) |>
  dplyr::mutate(
    fold_enrichment = (n_overlap / length({gene_list_name})) / (n_gene_set / length({background_name})),
    enrichment_score = log2(fold_enrichment)
  ) |>
  dplyr::arrange(p_value)

# Summary
message("ORA Results Summary:")
message("- Gene sets tested: ", nrow({output_name}))
message("- Significant (adj_p < 0.05): ", sum({output_name}$adj_p_value < 0.05, na.rm = TRUE))

head({output_name})
')

  cat(code)
  return(invisible(NULL))
}
