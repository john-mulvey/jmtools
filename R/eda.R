# Exploratory Data Analysis Functions
# jmtools package

#' Test associations between principal components and metadata
#'
#' Tests associations between PCA coordinates and metadata variables using
#' appropriate statistical tests: Spearman correlation for continuous variables
#' and Kruskal-Wallis test for categorical variables.
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
#'   for options (default: "BH" for Benjamini-Hochberg).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{variable}{Name of the metadata variable}
#'   \item{principal_component}{Name of the principal component (e.g., "PC1")}
#'   \item{test_type}{Statistical test used ("Spearman" or "Kruskal-Wallis")}
#'   \item{statistic}{Effect size: Spearman's rho for continuous variables,
#'     eta-squared for categorical variables}
#'   \item{p_value}{Raw p-value from the statistical test}
#'   \item{n}{Number of complete observations used}
#'   \item{p_adj}{Adjusted p-value using the specified method}
#' }
#'
#' @details
#' For continuous variables, Spearman's rank correlation is used, which is
#' robust to non-normality and outliers. For categorical variables (factors
#' or character vectors), the Kruskal-Wallis test is used with eta-squared
#' as the effect size measure.
#'
#' Variables with more than `1 - min_complete_pct` proportion of missing
#' values are automatically excluded from testing.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run PCA with prcomp
#' pca_result <- prcomp(t(expr_matrix), scale. = TRUE, center = TRUE)
#'
#' # Test associations
#' associations <- test_pc_metadata_associations(
#'   pca_result = pca_result,
#'   metadata = sample_metadata,
#'   vars_to_test = c("age", "sex", "treatment", "bmi"),
#'   id_col = "sample_id",
#'   min_var_explained = 5
#' )
#'
#' # View significant associations
#' associations[associations$p_adj < 0.05, ]
#' }
#'
#' @seealso [plot_pc_metadata_associations()] for visualising results
#'
#' @importFrom stats cor.test kruskal.test p.adjust
test_pc_metadata_associations <- function(pca_result,
                                          metadata,
                                          vars_to_test,
                                          id_col = "sample_id",
                                          n_pcs = 5,
                                          min_var_explained = 0,
                                          min_complete_pct = 0.7,
                                          p_adjust_method = "BH") {

  # Validate inputs
  if (!inherits(pca_result, "prcomp")) {
    stop("pca_result must be a prcomp object from stats::prcomp()")
  }

  if (missing(vars_to_test) || is.null(vars_to_test) || length(vars_to_test) == 0) {
    stop("vars_to_test must be specified and cannot be empty")
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

  # Initialise results list
  results_list <- list()

  # Test each variable against each PC
  for (var in vars_to_test) {
    var_data <- combined_data[[var]]

    # Skip if all NA
    if (all(is.na(var_data))) next

    # Determine variable type
    is_continuous <- is.numeric(var_data) && !is.factor(var_data)

    for (pc_idx in pcs_to_use) {
      pc_name <- colnames(pca_coords)[pc_idx]
      pc_values <- combined_data[[pc_name]]

      # Get complete cases
      complete_idx <- !is.na(var_data) & !is.na(pc_values)
      n_complete <- sum(complete_idx)

      if (n_complete < 3) next

      if (is_continuous) {
        # Spearman correlation for continuous variables
        test_result <- stats::cor.test(
          var_data[complete_idx],
          pc_values[complete_idx],
          method = "spearman",
          exact = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          variable = var,
          principal_component = pc_name,
          test_type = "Spearman",
          statistic = as.numeric(test_result$estimate),
          p_value = test_result$p.value,
          n = n_complete,
          stringsAsFactors = FALSE
        )

      } else {
        # Kruskal-Wallis for categorical variables
        test_data <- data.frame(
          pc = pc_values[complete_idx],
          group = droplevels(as.factor(var_data[complete_idx]))
        )

        # Need at least 2 groups
        if (length(unique(test_data$group)) < 2) next

        kw_test <- stats::kruskal.test(pc ~ group, data = test_data)

        # Calculate eta-squared as effect size
        # Formula: (H - k + 1) / (n - k), where H is the test statistic and k is number of groups
        k <- length(unique(test_data$group))
        eta_sq <- (as.numeric(kw_test$statistic) - k + 1) / (n_complete - k)

        results_list[[length(results_list) + 1]] <- data.frame(
          variable = var,
          principal_component = pc_name,
          test_type = "Kruskal-Wallis",
          statistic = eta_sq,
          p_value = kw_test$p.value,
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

  # Adjust p-values
  results_df$p_adj <- stats::p.adjust(results_df$p_value, method = p_adjust_method)

  return(results_df)
}


#' Plot PC-metadata associations as a heatmap
#'
#' Creates a heatmap visualisation of associations between principal components
#' and metadata variables, with significant associations highlighted.
#'
#' @param association_results Data frame of results from
#'   [test_pc_metadata_associations()].
#' @param p_threshold Significance threshold for highlighting associations
#'   (default: 0.05).
#' @param use_adjusted_p Logical. If `TRUE` (default), uses adjusted p-values
#'   for determining significance. If `FALSE`, uses raw p-values.
#' @param show_all Logical. If `TRUE`, shows all tested associations. If
#'   `FALSE` (default), shows only significant associations.
#' @param order_by_pc Integer specifying which PC to use for ordering variables
#'   by effect size (default: 1 for PC1).
#' @param fill_limits Numeric vector of length 2 specifying the limits for the
#'   colour scale (default: `c(-1, 1)`). Values outside this range are squished
#'   to the nearest limit. For Spearman correlations, \[-1, 1\] is appropriate.
#'   For eta-squared values (Kruskal-Wallis), consider `c(0, 1)` or a narrower
#'   range based on observed values.
#' @param title Optional plot title. If `NULL`, no title is added.
#'
#' @return A ggplot2 object. Returns `NULL` with a warning if no associations
#'   meet the criteria for plotting.
#'
#' @details
#' The heatmap uses a diverging colour scale from blue (negative association)
#' through white (no association) to red (positive association). Significant
#' associations are highlighted with bold black borders and bold text.
#'
#' Effect sizes outside `fill_limits` are squished to the nearest limit.
#' The default \[-1, 1\] is appropriate for Spearman correlations. For
#' Kruskal-Wallis eta-squared values (which are always positive and typically
#' small), consider adjusting the limits or using `c(0, 1)`.
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
#' # Use raw p-values instead of adjusted
#' plot_pc_metadata_associations(associations, use_adjusted_p = FALSE)
#' }
#'
#' @seealso [test_pc_metadata_associations()] for generating the input data
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 coord_flip theme_minimal labs theme element_text
#' @importFrom rlang .data
plot_pc_metadata_associations <- function(association_results,
                                          p_threshold = 0.05,
                                          use_adjusted_p = TRUE,
                                          show_all = FALSE,
                                          order_by_pc = 1,
                                          fill_limits = c(-1, 1),
                                          title = NULL) {

  if (nrow(association_results) == 0) {
    warning("No association results to plot")
    return(NULL)
  }

  p_col <- if (use_adjusted_p) "p_adj" else "p_value"

  # Filter to significant associations if requested
  plot_data <- association_results
  if (!show_all) {
    plot_data <- plot_data[plot_data[[p_col]] < p_threshold, ]
  }

  if (nrow(plot_data) == 0) {
    warning("No associations meet the plotting criteria")
    return(NULL)
  }

  # Order variables by association with specified PC
  pc_name <- paste0("PC", order_by_pc)
  pc_subset <- plot_data[plot_data$principal_component == pc_name, ]

  if (nrow(pc_subset) > 0) {
    var_order <- pc_subset$variable[order(abs(pc_subset$statistic), decreasing = TRUE)]
    all_vars <- unique(c(var_order, plot_data$variable))
  } else {
    all_vars <- unique(plot_data$variable)
  }

  plot_data$variable <- factor(plot_data$variable, levels = all_vars)

  # Mark significant associations
  plot_data$significant <- plot_data[[p_col]] < p_threshold

  # Create base plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$variable,
                                    y = .data$principal_component,
                                    fill = .data$statistic)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limits = fill_limits, space = "Lab",
      name = "Effect size\n(\u03c1 or \u03b7\u00b2)",
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

  # Add text labels
  nonsig_data <- plot_data[!plot_data$significant, ]
  if (nrow(nonsig_data) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = nonsig_data,
        ggplot2::aes(label = sprintf("%.2f", .data$statistic)),
        colour = "grey40", size = 3
      )
  }

  if (nrow(sig_data) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = sig_data,
        ggplot2::aes(label = sprintf("%.2f", .data$statistic)),
        colour = "black", size = 3, fontface = "bold"
      )
  }

  return(p)
}
