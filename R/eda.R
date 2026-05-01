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


#' Test associations between variables and a binary group
#'
#' Tests associations between metadata variables and a binary group variable
#' using point-biserial correlation (equivalent to Pearson r with a binary
#' outcome). Two-level factors and character variables are automatically
#' encoded as 0/1 integers.
#'
#' @param metadata Data frame containing metadata variables to test.
#' @param vars_to_test Character vector of variable names from `metadata` to
#'   test for association with the group variable.
#' @param group_var Name of the binary group column in `metadata`.
#'   Default: `"group"`.
#' @param positive_level Character string specifying which level of `group_var`
#'   is encoded as 1 (the "positive" or "case" class).
#' @param min_complete Minimum number of complete (non-NA) observations
#'   required to test a variable. Default: 5.
#' @param p_adjust_method Method for p-value adjustment passed to
#'   [stats::p.adjust()]. Default: `"BH"`.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{variable}{Name of the metadata variable}
#'   \item{r}{Pearson correlation coefficient}
#'   \item{r_lo}{Lower bound of the 95% confidence interval for r}
#'   \item{r_hi}{Upper bound of the 95% confidence interval for r}
#'   \item{p_value}{Raw p-value from [stats::cor.test()]}
#'   \item{n_complete}{Number of complete observations used}
#'   \item{adj_p_value}{Adjusted p-value using the specified method}
#' }
#' Rows are ordered by ascending p-value.
#'
#' @details
#' Variables that are factors or character vectors with exactly two unique
#' levels are automatically encoded as 0/1 integers (the second level
#' alphabetically becomes 1). Variables with more than two levels are skipped
#' with a message. Variables with zero variance or fewer than `min_complete`
#' observations are also skipped.
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
#' @importFrom stats cor.test p.adjust var
test_variable_group_associations <- function(metadata,
                                              vars_to_test,
                                              group_var      = "group",
                                              positive_level,
                                              min_complete   = 5,
                                              p_adjust_method = "BH") {

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

    # Auto-encode two-level factors/characters as 0/1
    if (is.factor(x) || is.character(x)) {
      lvls <- sort(unique(x[!is.na(x)]))
      if (length(lvls) == 2) {
        message("Auto-encoding binary variable '", var, "': '", lvls[2], "' = 1")
        x <- as.integer(x == lvls[2])
      } else {
        message("Skipping '", var, "': non-binary categorical (", length(lvls), " levels)")
        next
      }
    }

    complete  <- !is.na(x) & !is.na(group_binary)
    n_complete <- sum(complete)

    if (n_complete < min_complete) {
      message("Skipping '", var, "': fewer than ", min_complete, " complete observations")
      next
    }
    if (stats::var(x[complete]) == 0) {
      message("Skipping '", var, "': zero variance")
      next
    }

    result <- tryCatch({
      ct <- stats::cor.test(x[complete], group_binary[complete], method = "pearson")
      data.frame(
        variable   = var,
        r          = as.numeric(ct$estimate),
        r_lo       = ct$conf.int[1],
        r_hi       = ct$conf.int[2],
        p_value    = ct$p.value,
        n_complete = n_complete,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("Skipping '", var, "': ", e$message)
      NULL
    })

    if (!is.null(result)) {
      results_list[[length(results_list) + 1]] <- result
    }
  }

  if (length(results_list) == 0) {
    warning("No associations could be tested")
    return(data.frame())
  }

  results_df <- do.call(rbind, results_list)
  results_df$adj_p_value <- stats::p.adjust(results_df$p_value, method = p_adjust_method)
  results_df <- results_df[order(results_df$p_value), ]
  rownames(results_df) <- NULL

  return(results_df)
}


#' Plot variable-group associations as a forest plot
#'
#' Creates a forest plot of point-biserial correlations between metadata
#' variables and a binary group, with 95% confidence intervals and optional
#' significance colouring.
#'
#' @param association_results Data frame of results from
#'   [test_variable_group_associations()].
#' @param p_threshold Significance threshold for colouring associations
#'   (default: 0.05).
#' @param use_adjusted_p Logical. If `TRUE` (default), uses adjusted p-values
#'   for determining significance. If `FALSE`, uses raw p-values.
#' @param group_var Character string used in the x-axis label. Should match
#'   the `group_var` argument passed to [test_variable_group_associations()].
#'   Default: `"group"`.
#' @param positive_level Character string for the x-axis label indicating
#'   which group level was encoded as 1. If `NULL`, the label omits this
#'   detail. Default: `NULL`.
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
#' plot_variable_group_associations(
#'   results,
#'   group_var      = "group",
#'   positive_level = "heart_failure"
#' )
#' }
#'
#' @seealso [test_variable_group_associations()] for generating the input data
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_vline scale_colour_manual labs
#' @importFrom rlang .data
#' @importFrom stats reorder
plot_variable_group_associations <- function(association_results,
                                              p_threshold    = 0.05,
                                              use_adjusted_p = TRUE,
                                              group_var      = "group",
                                              positive_level = NULL) {

  if (nrow(association_results) == 0) {
    warning("No association results to plot")
    return(NULL)
  }

  p_col <- if (use_adjusted_p) "adj_p_value" else "p_value"

  plot_data <- association_results
  plot_data$significant <- plot_data[[p_col]] < p_threshold
  plot_data$variable    <- reorder(plot_data$variable, plot_data$r)

  x_label <- if (!is.null(positive_level)) {
    paste0("Pearson r with ", group_var, " (", positive_level, " = 1)")
  } else {
    paste0("Pearson r with ", group_var)
  }

  legend_label <- paste0(if (use_adjusted_p) "adj. " else "", "p < ", p_threshold)

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = .data$r, y = .data$variable,
                               colour = .data$significant)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$r_lo, xmax = .data$r_hi),
      width = 0.2, orientation = "y"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey50", "TRUE" = "red"),
      name   = legend_label
    ) +
    ggplot2::labs(x = x_label, y = NULL) +
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
