# Exploratory Data Analysis Functions
# jmtools package

#' Test associations between principal components and metadata
#'
#' Tests associations between PCA coordinates and metadata variables, returning
#' a non-negative effect size in \[0, 1\] for every (variable, PC) pair. Two
#' methods are supported: a univariate per-variable linear model with an F-test
#' p-value, and a multivariate Shapley-value (LMG) decomposition of the joint
#' R-squared with no per-cell p-value.
#'
#' @param pca_result A `prcomp` object from [stats::prcomp()]. The function
#'   extracts scores from `pca_result$x` and variance explained from
#'   `pca_result$sdev`.
#' @param metadata Data frame containing metadata variables to test.
#' @param vars_to_test Character vector of variable names from `metadata` to
#'   test for association with PCs.
#' @param id_col Name of the column containing sample identifiers. Must be
#'   present in `metadata` and match the rownames of `pca_result$x`.
#' @param n_pcs Maximum number of principal components to test (default: 5).
#' @param min_var_explained Minimum percentage of variance a PC must explain
#'   to be included (default: 0). Set to a positive value (e.g., 5) to only
#'   test PCs explaining meaningful variance.
#' @param min_complete_pct Minimum proportion of non-missing values required
#'   for a variable to be tested (default: 0.7, i.e., 70%).
#' @param p_adjust_method Method for p-value adjustment. See [stats::p.adjust()]
#'   for options (default: "BH" for Benjamini-Hochberg). Ignored when
#'   `method = "multivariate_shapley"`.
#' @param method Either `"univariate"` (default) or `"multivariate_shapley"`.
#'   See Details.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{variable}{Name of the metadata variable}
#'   \item{principal_component}{Name of the principal component (e.g., "PC1")}
#'   \item{test_type}{Method used: `"Univariate"` or `"Multivariate Shapley"`}
#'   \item{statistic}{Effect size in \[0, 1\] (see Details)}
#'   \item{p_value}{Raw p-value (univariate only; `NA` for shapley)}
#'   \item{n}{Number of complete observations used}
#'   \item{p_adj}{Adjusted p-value (univariate only; `NA` for shapley)}
#' }
#'
#' @details
#' \strong{Univariate.} For each (variable, PC) pair, fits a linear model
#' `lm(PC ~ variable)` and returns the model R-squared as the effect size
#' and the F-test p-value (from `anova(fit)$"Pr(>F)"[1]`) as the significance.
#' Continuous variables yield squared Pearson correlation; categorical
#' variables (factors or character vectors) yield eta-squared (one-way ANOVA),
#' both on the same \[0, 1\] scale. Pairwise complete cases are used per cell.
#' The F-test p-value assumes approximately normal residuals; for categorical
#' variables with strongly unbalanced groups this can be sensitive at small n.
#'
#' \strong{Multivariate Shapley.} For each PC, fits a single joint model
#' `lm(PC ~ var1 + var2 + ...)` and decomposes the joint R-squared into
#' per-variable contributions via the Shapley value (LMG metric: the average
#' marginal R-squared contribution over all variable orderings), implemented
#' by [relaimpo::calc.relimp()] with `type = "lmg"`. The contributions are
#' non-negative and sum exactly to the joint R-squared. This addresses the
#' confounding problem of the univariate approach (correlated variables both
#' appearing associated with the same PC) by partitioning shared variance
#' fairly between predictors. No per-variable p-value is returned: the LMG
#' statistic is a descriptive decomposition and lacks a natural null
#' distribution, so `p_value` and `p_adj` are `NA`. Listwise-complete cases
#' are used per PC (samples with any NA across `vars_to_test` are dropped),
#' and `vars_to_test` must contain at least two variables. Requires the
#' `relaimpo` package.
#'
#' Variables with more than `1 - min_complete_pct` proportion of missing
#' values are automatically excluded from testing in either mode.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pca_result <- prcomp(t(expr_matrix), scale. = TRUE, center = TRUE)
#'
#' # Default: univariate per-variable tests with p-values
#' associations <- test_pc_metadata_associations(
#'   pca_result = pca_result,
#'   metadata = sample_metadata,
#'   vars_to_test = c("age", "sex", "treatment", "bmi"),
#'   id_col = "sample_id",
#'   min_var_explained = 5
#' )
#' associations[associations$p_adj < 0.05, ]
#'
#' # Multivariate decomposition: handles confounding, no p-values
#' shapley <- test_pc_metadata_associations(
#'   pca_result = pca_result,
#'   metadata = sample_metadata,
#'   vars_to_test = c("age", "sex", "treatment", "bmi"),
#'   id_col = "sample_id",
#'   method = "multivariate_shapley"
#' )
#' }
#'
#' @seealso [plot_pc_metadata_associations()] for visualising results
#'
#' @importFrom stats lm anova p.adjust complete.cases as.formula var
test_pc_metadata_associations <- function(pca_result,
                                          metadata,
                                          vars_to_test,
                                          id_col = "sample_id",
                                          n_pcs = 5,
                                          min_var_explained = 0,
                                          min_complete_pct = 0.7,
                                          p_adjust_method = "BH",
                                          method = c("univariate", "multivariate_shapley")) {

  method <- match.arg(method)

  # Validate inputs
  if (!inherits(pca_result, "prcomp")) {
    stop("pca_result must be a prcomp object from stats::prcomp()")
  }

  if (missing(vars_to_test) || is.null(vars_to_test) || length(vars_to_test) == 0) {
    stop("vars_to_test must be specified and cannot be empty")
  }

  if (method == "multivariate_shapley") {
    if (!requireNamespace("relaimpo", quietly = TRUE)) {
      stop("Package 'relaimpo' is required for method = 'multivariate_shapley'. ",
           "Install it with install.packages('relaimpo')")
    }
  }

  # Extract PC scores
  pca_coords <- pca_result$x

  # Calculate variance explained from sdev
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  # Determine which PCs to test
  if (min_var_explained > 0) {
    pcs_to_use <- which(var_explained >= min_var_explained)

    if (length(pcs_to_use) == 0) {
      warning("No PCs explain >= ", min_var_explained,
              "% variance. Using first ", n_pcs, " PCs.")
      pcs_to_use <- seq_len(min(n_pcs, ncol(pca_coords)))
    } else {
      pcs_to_use <- pcs_to_use[pcs_to_use <= n_pcs]
      message("Testing ", length(pcs_to_use), " PCs explaining >= ",
              min_var_explained, "% variance: ",
              paste0("PC", pcs_to_use, " (", round(var_explained[pcs_to_use], 1), "%)",
                     collapse = ", "))
    }
  } else {
    pcs_to_use <- seq_len(min(n_pcs, ncol(pca_coords)))
  }

  # Prepare PC data with sample IDs
  pc_data <- as.data.frame(pca_coords[, seq_len(max(pcs_to_use)), drop = FALSE])
  pc_data[[id_col]] <- rownames(pca_coords)

  # Merge with metadata
  combined_data <- merge(pc_data, metadata, by = id_col)

  # Report sample drops from the merge so a mismatch between PCA scores and
  # metadata sample sets is visible rather than silent.
  n_pca_samples <- nrow(pc_data)
  n_meta_samples <- nrow(metadata)
  n_kept <- nrow(combined_data)
  if (n_kept < n_pca_samples || n_kept < n_meta_samples) {
    message("Merging PCA scores and metadata: ", n_kept,
            " samples retained (",
            n_pca_samples - n_kept, " PCA samples and ",
            n_meta_samples - n_kept, " metadata samples ",
            "without a match were dropped).")
  }

  # Filter variables by completeness and existence
  vars_to_test <- vars_to_test[vars_to_test %in% names(combined_data)]
  if (length(vars_to_test) == 0) {
    warning("None of the specified variables found in metadata")
    return(data.frame())
  }

  valid_pct <- vapply(vars_to_test, function(v) {
    mean(!is.na(combined_data[[v]]))
  }, numeric(1))

  excluded_vars <- vars_to_test[valid_pct < min_complete_pct]
  if (length(excluded_vars) > 0) {
    message("Excluding variables with <", min_complete_pct * 100,
            "% complete data: ", paste(excluded_vars, collapse = ", "))
  }
  vars_to_test <- vars_to_test[valid_pct >= min_complete_pct]

  if (length(vars_to_test) == 0) {
    warning("No variables meet the minimum completeness threshold")
    return(data.frame())
  }

  if (method == "multivariate_shapley" && length(vars_to_test) < 2) {
    stop("method = 'multivariate_shapley' requires at least 2 variables ",
         "after completeness filtering")
  }

  # Initialise results list
  results_list <- list()

  if (method == "univariate") {
    for (var in vars_to_test) {
      var_data <- combined_data[[var]]

      if (all(is.na(var_data))) next

      is_continuous <- is.numeric(var_data) && !is.factor(var_data)

      for (pc_idx in pcs_to_use) {
        pc_name <- colnames(pca_coords)[pc_idx]
        pc_values <- combined_data[[pc_name]]

        complete_idx <- !is.na(var_data) & !is.na(pc_values)
        n_complete <- sum(complete_idx)

        if (n_complete < 3) next

        if (is_continuous) {
          x <- var_data[complete_idx]
          if (stats::var(x) == 0) next
          fit_data <- data.frame(pc = pc_values[complete_idx], var = x)
        } else {
          group <- droplevels(as.factor(var_data[complete_idx]))
          if (length(levels(group)) < 2) next
          fit_data <- data.frame(pc = pc_values[complete_idx], var = group)
        }

        fit <- tryCatch(
          stats::lm(pc ~ var, data = fit_data),
          error = function(e) NULL
        )
        if (is.null(fit)) next

        r_squared <- summary(fit)$r.squared
        p_value <- stats::anova(fit)$"Pr(>F)"[1]

        if (is.na(r_squared) || is.na(p_value)) next

        results_list[[length(results_list) + 1]] <- data.frame(
          variable = var,
          principal_component = pc_name,
          test_type = "Univariate",
          statistic = r_squared,
          p_value = p_value,
          n = n_complete,
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    # Multivariate Shapley: one joint lm per PC, decomposed via relaimpo LMG.
    # Coerce character variables to factors so lm/calc.relimp handle them
    # consistently across PC iterations.
    for (v in vars_to_test) {
      if (is.character(combined_data[[v]])) {
        combined_data[[v]] <- as.factor(combined_data[[v]])
      }
    }

    for (pc_idx in pcs_to_use) {
      pc_name <- colnames(pca_coords)[pc_idx]

      shapley_data <- combined_data[, c(pc_name, vars_to_test), drop = FALSE]
      complete_idx <- stats::complete.cases(shapley_data)
      n_complete <- sum(complete_idx)

      # Need n > p + 1 at minimum for the joint model; require a margin.
      if (n_complete < length(vars_to_test) + 2) {
        warning("Skipping ", pc_name, ": only ", n_complete,
                " complete cases for ", length(vars_to_test),
                "-variable joint model")
        next
      }

      fit_data <- shapley_data[complete_idx, , drop = FALSE]
      fit_formula <- stats::as.formula(
        paste(pc_name, "~", paste(vars_to_test, collapse = " + "))
      )

      fit <- tryCatch(
        stats::lm(fit_formula, data = fit_data),
        error = function(e) {
          warning("Skipping ", pc_name, ": lm failed (", e$message, ")")
          NULL
        }
      )
      if (is.null(fit)) next

      ri <- tryCatch(
        relaimpo::calc.relimp(fit, type = "lmg"),
        error = function(e) {
          warning("Skipping ", pc_name, ": relaimpo::calc.relimp failed (",
                  e$message, ")")
          NULL
        }
      )
      if (is.null(ri)) next

      lmg_values <- ri@lmg

      for (var in vars_to_test) {
        results_list[[length(results_list) + 1]] <- data.frame(
          variable = var,
          principal_component = pc_name,
          test_type = "Multivariate Shapley",
          statistic = as.numeric(lmg_values[var]),
          p_value = NA_real_,
          n = n_complete,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Combine results
  if (length(results_list) == 0) {
    warning("No associations could be tested")
    return(data.frame())
  }

  results_df <- do.call(rbind, results_list)

  # Adjust p-values (NAs propagate cleanly for shapley rows)
  results_df$p_adj <- stats::p.adjust(results_df$p_value, method = p_adjust_method)

  return(results_df)
}


#' Plot PC-metadata associations as a heatmap
#'
#' Creates a heatmap visualisation of associations between principal components
#' and metadata variables, with significant associations highlighted when
#' p-values are available.
#'
#' @param association_results Data frame of results from
#'   [test_pc_metadata_associations()]. The `statistic` column is expected to
#'   be R-squared in \[0, 1\].
#' @param p_threshold Significance threshold for highlighting associations
#'   (default: 0.05). Ignored when `association_results` lacks p-values.
#' @param use_adjusted_p Logical. If `TRUE` (default), uses adjusted p-values
#'   for determining significance. If `FALSE`, uses raw p-values.
#' @param show_all Logical. If `TRUE` (default), shows all tested metadata
#'   variables. If `FALSE`, restricts the plot to variables that have at
#'   least one significant association (across any PC); cells of those
#'   variables that are themselves non-significant are still shown so the
#'   reader can see which PCs are and are not associated. Silently
#'   overridden to `TRUE` when no p-values are available (e.g.,
#'   `method = "multivariate_shapley"` results), since there is no
#'   significance criterion to filter on.
#' @param order_by_pc Integer specifying which PC to use for ordering variables
#'   by effect size (default: 1 for PC1).
#' @param fill_limits Numeric vector of length 2 specifying the limits for the
#'   colour scale (default: `c(0, 1)`, the natural range of R-squared). Values
#'   outside this range are squished to the nearest limit. Narrow this range
#'   (e.g., `c(0, 0.5)`) to enhance contrast when most observed values are
#'   small.
#' @param title Optional plot title. If `NULL`, no title is added.
#'
#' @return A ggplot2 object. Returns `NULL` with a warning if no associations
#'   meet the criteria for plotting.
#'
#' @details
#' The heatmap uses a sequential viridis colour scale, since the underlying
#' statistic (R-squared) is non-negative. Tile text colour switches between
#' white (on darker, low-R-squared tiles) and black (on brighter,
#' high-R-squared tiles) for legibility across the viridis range.
#'
#' When p-values are present (univariate results), significant associations
#' are highlighted with bold black borders and bold text. When p-values are
#' absent or all `NA` (multivariate Shapley results), every tile is plotted
#' without significance highlighting -- the absence of bold borders is
#' itself a visual cue that the values are a variance decomposition rather
#' than per-cell tests.
#'
#' Effect sizes outside `fill_limits` are squished to the nearest limit.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # First run PCA and association tests
#' pca_result <- prcomp(t(expr_matrix), scale. = TRUE, center = TRUE)
#' associations <- test_pc_metadata_associations(
#'   pca_result = pca_result,
#'   metadata = sample_metadata,
#'   vars_to_test = c("age", "sex", "treatment"),
#'   id_col = "sample_id"
#' )
#'
#' # Plot all associations
#' plot_pc_metadata_associations(associations, show_all = TRUE)
#'
#' # Plot only significant associations
#' plot_pc_metadata_associations(associations, show_all = FALSE)
#'
#' # Narrow the colour scale to enhance contrast for small effect sizes
#' plot_pc_metadata_associations(associations, fill_limits = c(0, 0.5))
#'
#' # Multivariate Shapley results: no p-values, all tiles shown unfiltered
#' shapley <- test_pc_metadata_associations(
#'   pca_result = pca_result,
#'   metadata = sample_metadata,
#'   vars_to_test = c("age", "sex", "treatment"),
#'   id_col = "sample_id",
#'   method = "multivariate_shapley"
#' )
#' plot_pc_metadata_associations(shapley)
#' }
#'
#' @seealso [test_pc_metadata_associations()] for generating the input data
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_viridis_c coord_flip theme_minimal labs theme element_text
#' @importFrom rlang .data
plot_pc_metadata_associations <- function(association_results,
                                          p_threshold = 0.05,
                                          use_adjusted_p = TRUE,
                                          show_all = TRUE,
                                          order_by_pc = 1,
                                          fill_limits = c(0, 1),
                                          title = NULL) {

  if (nrow(association_results) == 0) {
    warning("No association results to plot")
    return(NULL)
  }

  p_col <- if (use_adjusted_p) "p_adj" else "p_value"

  # Detect whether usable p-values are present. If not (e.g. shapley results),
  # fall through to plotting all tiles without significance highlighting.
  has_significance <- p_col %in% names(association_results) &&
                      !all(is.na(association_results[[p_col]]))

  if (!has_significance) {
    show_all <- TRUE
  }

  # Restrict to variables that have at least one significant association if
  # requested. We filter at the variable level, not the cell level, so the
  # null cells of retained variables remain visible (informative for
  # confounding analysis).
  plot_data <- association_results
  if (!show_all) {
    sig_vars <- unique(plot_data$variable[!is.na(plot_data[[p_col]]) &
                                          plot_data[[p_col]] < p_threshold])
    plot_data <- plot_data[plot_data$variable %in% sig_vars, ]
  }

  if (nrow(plot_data) == 0) {
    warning("No associations meet the plotting criteria")
    return(NULL)
  }

  # Order variables by association with specified PC
  pc_name <- paste0("PC", order_by_pc)
  pc_subset <- plot_data[plot_data$principal_component == pc_name, ]

  if (nrow(pc_subset) > 0) {
    var_order <- pc_subset$variable[order(pc_subset$statistic, decreasing = TRUE)]
    all_vars <- unique(c(var_order, plot_data$variable))
  } else {
    all_vars <- unique(plot_data$variable)
  }

  plot_data$variable <- factor(plot_data$variable, levels = all_vars)

  # Mark significant associations (always FALSE when p-values absent)
  if (has_significance) {
    plot_data$significant <- !is.na(plot_data[[p_col]]) &
      plot_data[[p_col]] < p_threshold
  } else {
    plot_data$significant <- FALSE
  }

  # Choose text colour for legibility against viridis: white on darker
  # (low-R-squared) tiles, black on brighter (high-R-squared) tiles.
  fill_midpoint <- mean(fill_limits)
  text_colour <- ifelse(plot_data$statistic >= fill_midpoint, "black", "white")
  plot_data$text_colour <- text_colour

  # Create base plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$variable,
                                    y = .data$principal_component,
                                    fill = .data$statistic)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(
      option = "viridis",
      limits = fill_limits,
      name = "Variance\nexplained\n(R\u00b2)",
      oob = scales::squish
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Metadata variable",
      y = "Principal component",
      title = title
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
    )

  # Add borders for significant associations
  sig_data <- plot_data[plot_data$significant, ]
  if (nrow(sig_data) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = sig_data,
        colour = "black", linewidth = 1.5, fill = NA
      )
  }

  # Add text labels (non-significant: regular weight; significant: bold).
  # Use identity scale to honour per-row text colour without a legend entry.
  nonsig_data <- plot_data[!plot_data$significant, ]
  if (nrow(nonsig_data) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = nonsig_data,
        ggplot2::aes(label = sprintf("%.2f", .data$statistic),
                     colour = .data$text_colour),
        size = 3
      )
  }

  if (nrow(sig_data) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = sig_data,
        ggplot2::aes(label = sprintf("%.2f", .data$statistic),
                     colour = .data$text_colour),
        size = 3, fontface = "bold"
      )
  }

  p <- p + ggplot2::scale_colour_identity()

  return(p)
}


#' Standardised mean differences between a binary group and metadata variables
#'
#' Computes Hedges' g (a small-sample-bias-corrected standardised mean
#' difference) between a binary group variable and each metadata variable.
#' Continuous predictors use the standard d = (mean_case - mean_ctrl) / pooled
#' SD, then Hedges' J correction. Two-level factors and character variables
#' are auto-encoded as 0/1 and use the proportion-difference SMD
#' (p_case - p_ctrl) / sqrt(p_pool * (1 - p_pool)) with the same J correction.
#'
#' @param metadata Data frame containing metadata variables to test.
#' @param vars_to_test Character vector of variable names from `metadata` to
#'   test for association with the group variable.
#' @param group_var Name of the binary group column in `metadata`.
#'   Default: `"group"`.
#' @param positive_level Character string specifying which level of `group_var`
#'   is treated as the "case" class. Effect-size signs are reported as
#'   case minus control.
#' @param positive_levels Optional named list mapping binary predictor variable
#'   names to the level that should be encoded as 1 (e.g.
#'   `list(treatment = "treated", smoker = "yes")`). Variables not listed fall
#'   back to the alphabetical default (the second level of
#'   `sort(unique(x))`). Use this to control the sign of `hedges_g` for binary
#'   character predictors where the alphabetical default produces a
#'   counter-intuitive encoding. Default: `list()`.
#' @param min_complete Minimum number of complete (non-NA) observations
#'   required to test a variable. Default: 5.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{variable}{Name of the metadata variable}
#'   \item{type}{`"continuous"` or `"binary"`}
#'   \item{hedges_g}{Hedges' g (signed; positive = higher in case)}
#'   \item{g_lo}{Lower bound of the 95% confidence interval for g}
#'   \item{g_hi}{Upper bound of the 95% confidence interval for g}
#'   \item{n_complete}{Number of complete observations used}
#' }
#' Rows are ordered by descending `|hedges_g|`.
#'
#' @details
#' For continuous predictors:
#' \deqn{d = (\bar{x}_{case} - \bar{x}_{ctrl}) /
#'           \sqrt{((n_1 - 1) s_1^2 + (n_0 - 1) s_0^2) / (n_1 + n_0 - 2)}}
#' For binary predictors (after 0/1 encoding via `positive_levels`):
#' \deqn{d = (p_{case} - p_{ctrl}) / \sqrt{p_{pool} (1 - p_{pool})}}
#' Hedges' J correction is then applied:
#' \deqn{J = 1 - 3 / (4 \cdot df - 1), \qquad df = n_1 + n_0 - 2}
#' \deqn{g = J \cdot d}
#' The 95% CI uses Hedges' approximate SE
#' \eqn{\sqrt{(n_1+n_0)/(n_1 n_0) + g^2 / (2 (n_1+n_0))}} times 1.96. The
#' J correction is exact for continuous d under normal residuals; applying
#' it to the binary SMD is a pragmatic approximation that gives a unified
#' column at the cost of slight imprecision at very small n.
#'
#' Variables with more than two unique levels (categorical with > 2 levels)
#' are skipped with a message. Variables with zero variance or fewer than
#' `min_complete` complete observations are also skipped.
#'
#' No p-values are returned. The propensity-score literature (Austin 2009
#' and others) argues that p-values are not informative for baseline imbalance
#' assessment because their sensitivity to sample size obscures the question
#' of imbalance magnitude. Use Hedges' g together with Cohen's rules of thumb
#' (small \eqn{|g| \approx 0.2}, medium \eqn{0.5}, large \eqn{0.8}).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- test_variable_group_associations(
#'   metadata       = sample_meta,
#'   vars_to_test   = c("age", "sex", "bmi", "ejection_fraction"),
#'   group_var      = "group",
#'   positive_level = "heart_failure"
#' )
#' }
#'
#' @seealso [plot_variable_group_associations()] for visualising results
#'
#' @importFrom stats var sd
test_variable_group_associations <- function(metadata,
                                              vars_to_test,
                                              group_var       = "group",
                                              positive_level,
                                              positive_levels = list(),
                                              min_complete    = 5) {

  if (!group_var %in% names(metadata)) {
    stop("'", group_var, "' not found in metadata")
  }
  if (!positive_level %in% as.character(metadata[[group_var]])) {
    stop("'", positive_level, "' not found in column '", group_var, "'")
  }

  group_binary <- as.integer(metadata[[group_var]] == positive_level)

  vars_to_test <- vars_to_test[vars_to_test %in% names(metadata)]
  if (length(vars_to_test) == 0) {
    stop("None of the specified variables found in metadata")
  }

  results_list <- list()

  for (var in vars_to_test) {
    x <- metadata[[var]]
    type <- "continuous"

    # Auto-encode two-level factors/characters as 0/1
    if (is.factor(x) || is.character(x)) {
      lvls <- sort(unique(x[!is.na(x)]))
      if (length(lvls) == 2) {
        positive <- positive_levels[[var]]
        if (is.null(positive)) {
          positive <- as.character(lvls[2])
          message("Auto-encoding binary variable '", var, "': '", positive,
                  "' = 1 (alphabetical default; pass via positive_levels to override)")
        } else {
          if (!positive %in% as.character(lvls)) {
            stop("positive_levels[['", var, "']] = '", positive,
                 "' is not a level of '", var, "' (levels: ",
                 paste(lvls, collapse = ", "), ")")
          }
          message("Encoding binary variable '", var, "': '", positive, "' = 1")
        }
        x <- as.integer(x == positive)
        type <- "binary"
      } else {
        message("Skipping '", var, "': non-binary categorical (", length(lvls), " levels)")
        next
      }
    }

    complete   <- !is.na(x) & !is.na(group_binary)
    n_complete <- sum(complete)

    if (n_complete < min_complete) {
      message("Skipping '", var, "': fewer than ", min_complete, " complete observations")
      next
    }
    if (stats::var(x[complete]) == 0) {
      message("Skipping '", var, "': zero variance")
      next
    }

    xc <- x[complete]
    gc <- group_binary[complete]
    n1 <- sum(gc == 1)
    n0 <- sum(gc == 0)
    if (n1 < 2 || n0 < 2) {
      message("Skipping '", var, "': fewer than 2 observations in one group")
      next
    }

    if (type == "continuous") {
      m1 <- mean(xc[gc == 1]); m0 <- mean(xc[gc == 0])
      s1 <- stats::sd(xc[gc == 1]); s0 <- stats::sd(xc[gc == 0])
      pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
      d <- (m1 - m0) / pooled_sd
    } else {
      p1 <- mean(xc[gc == 1]); p0 <- mean(xc[gc == 0]); p_pool <- mean(xc)
      d <- (p1 - p0) / sqrt(p_pool * (1 - p_pool))
    }

    df <- n1 + n0 - 2
    J  <- 1 - 3 / (4 * df - 1)
    g  <- J * d
    se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0)))

    results_list[[length(results_list) + 1]] <- data.frame(
      variable   = var,
      type       = type,
      hedges_g   = g,
      g_lo       = g - 1.96 * se,
      g_hi       = g + 1.96 * se,
      n_complete = n_complete,
      stringsAsFactors = FALSE
    )
  }

  if (length(results_list) == 0) {
    warning("No associations could be computed")
    return(data.frame())
  }

  results_df <- do.call(rbind, results_list)
  results_df <- results_df[order(-abs(results_df$hedges_g)), ]
  rownames(results_df) <- NULL

  return(results_df)
}


#' Plot variable-group associations as a forest plot
#'
#' Creates a forest plot of Hedges' g effect sizes between metadata variables
#' and a binary group, with 95% confidence intervals. Vertical reference lines
#' are drawn at Cohen's small (\eqn{\pm 0.2}) and medium (\eqn{\pm 0.5})
#' thresholds.
#'
#' @param association_results Data frame of results from
#'   [test_variable_group_associations()].
#' @param show_thresholds Logical. If `TRUE` (default), draws Cohen's small
#'   (\eqn{\pm 0.2}) and medium (\eqn{\pm 0.5}) reference lines.
#'
#' @return A ggplot2 object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- test_variable_group_associations(
#'   metadata       = sample_meta,
#'   vars_to_test   = c("age", "bmi", "ejection_fraction"),
#'   group_var      = "group",
#'   positive_level = "heart_failure"
#' )
#'
#' plot_variable_group_associations(results)
#' }
#'
#' @seealso [test_variable_group_associations()] for generating the input data
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_vline labs
#' @importFrom rlang .data
#' @importFrom stats reorder
plot_variable_group_associations <- function(association_results,
                                              show_thresholds = TRUE) {

  if (nrow(association_results) == 0) {
    warning("No association results to plot")
    return(NULL)
  }

  required_cols <- c("variable", "hedges_g", "g_lo", "g_hi")
  missing_cols <- setdiff(required_cols, names(association_results))
  if (length(missing_cols) > 0) {
    stop("association_results is missing columns: ",
         paste(missing_cols, collapse = ", "),
         ". Did you generate it with test_variable_group_associations()?")
  }

  plot_data <- association_results
  plot_data$variable <- reorder(plot_data$variable, plot_data$hedges_g)

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$hedges_g, y = .data$variable)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")

  if (show_thresholds) {
    p <- p +
      ggplot2::geom_vline(xintercept = c(-0.5, 0.5),
                          linetype = "dotted", colour = "grey60") +
      ggplot2::geom_vline(xintercept = c(-0.2, 0.2),
                          linetype = "dotted", colour = "grey75")
  }

  p +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$g_lo, xmax = .data$g_hi),
      width = 0.2, orientation = "y", colour = "grey20"
    ) +
    ggplot2::geom_point(size = 2.5, colour = "grey20") +
    ggplot2::labs(
      x = "Hedges' g",
      y = NULL,
      caption = if (show_thresholds)
        "Reference lines: Cohen's small (±0.2) and medium (±0.5)"
        else NULL
    ) +
    theme_jm()
}


#' Clustered heatmap with NA-tolerant distance computation
#'
#' A thin wrapper around [pheatmap::pheatmap()] for the specific case where the
#' abundance matrix contains NA values and you want to cluster rows and/or
#' columns. In all other cases, use [pheatmap::pheatmap()] directly.
#'
#' Standard `pheatmap` fails when clustering rows or columns that contain NAs
#' because `hclust` cannot handle NA distances. This function computes pairwise
#' Euclidean distances using [cluster::daisy()], which uses only shared non-NA
#' values for each pair, then passes pre-computed `hclust` objects to
#' `pheatmap`. NA cells are rendered in a distinct colour (default white).
#'
#' @param mat Numeric matrix with features as rows and samples as columns.
#'   May contain NA values.
#' @param scale_rows Logical. If `TRUE` (default), rows are z-score scaled
#'   before computing distances and plotting. Scaling uses [base::scale()] which
#'   ignores NAs. The heatmap is then drawn with `scale = "none"` so that
#'   clustering and display are consistent.
#' @param clustering_method Agglomeration method passed to [stats::hclust()].
#'   Default: `"ward.D2"`.
#' @param cluster_rows Logical. Whether to cluster rows. Default: `TRUE`.
#' @param cluster_cols Logical. Whether to cluster columns. Default: `TRUE`.
#' @param na_col Colour used for NA cells. Default: `"white"`.
#' @param ... Additional arguments passed to [pheatmap::pheatmap()], e.g.
#'   `annotation_col`, `show_rownames`, `main`, `fontsize_col`.
#'
#' @return Invisibly returns the `pheatmap` object.
#'
#' @details
#' When two rows (or columns) share no non-NA values, `daisy` returns NA for
#' that pair. These NA distances are replaced with the maximum observed distance
#' so that `hclust` can proceed; such features are placed at the periphery of
#' the dendrogram.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pheatmap_clustered_with_na(
#'   abundance_filt[sig_proteins, ],
#'   annotation_col = sample_annotation,
#'   show_rownames = FALSE,
#'   main = "Differentially abundant proteins"
#' )
#' }
#'
#' @importFrom stats hclust
pheatmap_clustered_with_na <- function(mat,
                                       scale_rows = TRUE,
                                       clustering_method = "ward.D2",
                                       cluster_rows = TRUE,
                                       cluster_cols = TRUE,
                                       na_col = "white",
                                       ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install it with install.packages('pheatmap')")
  }
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required. Install it with install.packages('cluster')")
  }

  if (scale_rows) {
    mat <- t(scale(t(mat)))
    n_constant <- sum(rowSums(!is.nan(mat) & !is.na(mat)) == 0)
    if (n_constant > 0) {
      warning(n_constant, " row(s) had zero variance and became all-NaN ",
              "after scaling; they will attach via max-distance replacement ",
              "in the dendrogram. Consider removing them before plotting if ",
              "unexpected.")
    }
  }

  na_tolerant_hclust <- function(m, method) {
    d <- cluster::daisy(m, metric = "euclidean")
    if (any(is.na(d))) {
      d[is.na(d)] <- max(d, na.rm = TRUE)
    }
    stats::hclust(d, method = method)
  }

  if (cluster_rows) {
    cluster_rows <- na_tolerant_hclust(mat, clustering_method)
  }
  if (cluster_cols) {
    cluster_cols <- na_tolerant_hclust(t(mat), clustering_method)
  }

  pheatmap::pheatmap(
    mat,
    scale = "none",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    na_col = na_col,
    ...
  )
}
