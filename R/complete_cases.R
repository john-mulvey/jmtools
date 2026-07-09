# Joint complete-cases subsetting via greedy sample drop with 2-opt
# refinement at each operating point. jmtools package.


# Internal helper: 2-opt local search starting from a kept set K. At
# each iteration, find the (drop sample, keep sample) swap that
# maximally increases the count of jointly-complete proteins; iterate
# to a local optimum. Vectorised over candidate swap-in samples for
# each candidate swap-out sample, using:
#   new_complete(j, j') = T1(j) + T2(j, j')
# T1(j)     = currently-complete proteins where j is not missing
# T2(j, j') = proteins with exactly one miss in K, that miss being
#             at j', and where j is also not missing.
# Returns list(K, complete).
two_opt_refine_internal <- function(detected, K_initial, max_iter = 50) {
  miss <- !detected
  K <- K_initial
  D <- setdiff(colnames(detected), K)
  current_n_miss   <- rowSums(miss[, K, drop = FALSE])
  current_complete <- sum(current_n_miss == 0L)
  for (iter in seq_len(max_iter)) {
    best_delta <- 0L; best_j <- NULL; best_jp <- NULL
    for (j in D) {
      not_missing_j   <- !miss[, j]
      one_miss_subset <- (current_n_miss == 1L) & not_missing_j
      complete_subset <- (current_n_miss == 0L) & not_missing_j
      T1 <- sum(complete_subset)
      T2_vec <- if (any(one_miss_subset)) {
        colSums(miss[one_miss_subset, K, drop = FALSE])
      } else rep(0L, length(K))
      new_complete_per_jp <- T1 + T2_vec
      best_jp_idx <- which.max(new_complete_per_jp)
      delta <- new_complete_per_jp[best_jp_idx] - current_complete
      if (delta > best_delta) {
        best_delta <- delta; best_j <- j; best_jp <- K[best_jp_idx]
      }
    }
    if (best_delta <= 0L) break
    K <- c(setdiff(K, best_jp), best_j)
    D <- c(setdiff(D, best_j), best_jp)
    current_n_miss   <- current_n_miss - miss[, best_jp] + miss[, best_j]
    current_complete <- current_complete + best_delta
  }
  list(K = K, complete = current_complete)
}


# Internal: greedy sample drop with 2-opt refinement at each step.
# Returns a tibble with one row per operating point, columns:
#   n_samples_kept       : samples remaining
#   n_complete_proteins  : proteins jointly observed in the kept set
#   kept_samples         : list column - character vector of kept
#                          sample names at this step
#
# The greedy itself is unaffected by 2-opt. The next greedy drop is
# taken from the greedy K, not from the 2-opt K. 2-opt is applied at
# each operating point as a refinement decorating the trace; it does
# not feed back into the greedy trajectory. The kept set must be
# stored explicitly because 2-opt may swap back a sample that greedy
# had previously dropped, so the kept set is not recoverable from a
# sequential drop list.
greedy_complete_cases_trace <- function(abundance) {
  if (is.null(colnames(abundance))) {
    stop("abundance must have column names (used to record kept samples)")
  }
  detected <- !is.na(abundance)
  greedy_K <- colnames(detected)
  N <- length(greedy_K)

  trace_rows <- vector("list", N)

  trace_rows[[1]] <- list(
    n_samples_kept      = N,
    n_complete_proteins = sum(rowSums(!detected[, greedy_K, drop = FALSE]) == 0),
    kept_samples        = greedy_K
  )

  step <- 1L
  while (length(greedy_K) > 1L) {
    current  <- detected[, greedy_K, drop = FALSE]
    one_miss <- rowSums(!current) == 1L
    gains <- if (any(one_miss)) {
      colSums(!current[one_miss, , drop = FALSE])
    } else {
      colSums(!current)
    }
    drop     <- names(which.max(gains))
    greedy_K <- setdiff(greedy_K, drop)

    K_step <- if (length(greedy_K) >= 1L && length(greedy_K) < N) {
      two_opt_refine_internal(detected, greedy_K)$K
    } else {
      greedy_K
    }

    step <- step + 1L
    trace_rows[[step]] <- list(
      n_samples_kept      = length(K_step),
      n_complete_proteins = sum(rowSums(!detected[, K_step, drop = FALSE]) == 0),
      kept_samples        = K_step
    )
  }

  tibble::tibble(
    n_samples_kept      = vapply(trace_rows, `[[`, integer(1), "n_samples_kept"),
    n_complete_proteins = vapply(trace_rows, `[[`, integer(1), "n_complete_proteins"),
    kept_samples        = lapply(trace_rows, `[[`, "kept_samples")
  )
}


#' Plot the joint complete-cases distribution from greedy + 2-opt
#'
#' Plots the (samples kept, complete proteins) curve produced by
#' greedy sample drop refined with 2-opt local search at each
#' operating point, optionally with the per-protein detection rate
#' curve overlaid as the optimistic upper bound on protein retention.
#'
#' The two curves answer different questions and are easy to confuse:
#'
#' \itemize{
#'   \item \emph{Per-protein detection rate}: proteins where each one
#'     individually clears an X\% detection threshold. The X\% of samples
#'     in which protein A is detected need not overlap with protein B's,
#'     so proteins are treated as fungible across samples. This is the
#'     optimistic upper bound on protein retention - \strong{not} a
#'     complete-cases curve.
#'   \item \emph{Joint complete cases (greedy + 2-opt)}: proteins fully
#'     observed across the *same* subset of samples. This is the
#'     achievable no-imputation sub-matrix.
#' }
#'
#' Finding the best *k* samples to drop for each *k* is NP-hard
#' (equivalent to maximum biclique on the bipartite detection graph),
#' so this function is a heuristic and provides no optimality
#' guarantee. It uses a greedy approach to seed each operating point -
#' at each step drop the sample whose removal unlocks the most
#' newly-complete proteins (or, when no protein is one-drop-away, the
#' sample with the most missing entries as a fallback) - and then
#' refines that kept set with 2-opt local search, accepting any (drop,
#' keep) swap that increases the complete-protein count until no swap
#' improves. On typical proteomics matrices this agrees with a HiGHS
#' MIP solve at most operating points and is within a few percent in
#' the steepest region of the curve.
#'
#' @param abundance Numeric matrix or data frame with proteins (or
#'   other features) as rows and samples as columns. Column names are
#'   required.
#' @param show_per_protein_curve Logical. If `TRUE` (default), overlays
#'   the per-protein detection rate curve as a comparison upper bound.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @seealso [get_greedy_complete_cases()] to extract the corresponding
#'   sub-matrix at a chosen operating point.
#'
#' @examples
#' \dontrun{
#' plot_greedy_complete_cases_distribution(abundance)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 expansion scale_colour_manual labs
plot_greedy_complete_cases_distribution <- function(abundance,
                                                    show_per_protein_curve = TRUE) {
  n_samples <- ncol(abundance)
  trace     <- greedy_complete_cases_trace(abundance)

  greedy_curve <- data.frame(
    curve            = "Joint complete cases (greedy + 2-opt)",
    pct_samples      = trace$n_samples_kept / n_samples * 100,
    n_proteins       = trace$n_complete_proteins,
    stringsAsFactors = FALSE
  )

  if (show_per_protein_curve) {
    per_protein_curve <- data.frame(
      curve            = "Per-protein detection rate",
      pct_samples      = (0:n_samples) / n_samples * 100,
      n_proteins       = sapply(0:n_samples, function(m) {
        sum(rowSums(!is.na(abundance)) >= m)
      }),
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(per_protein_curve, greedy_curve)
    colour_values <- c(
      "Per-protein detection rate"            = "#2166AC",
      "Joint complete cases (greedy + 2-opt)" = "#B2182B"
    )
  } else {
    plot_data <- greedy_curve
    colour_values <- c("Joint complete cases (greedy + 2-opt)" = "#B2182B")
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$pct_samples,
                 y = .data$n_proteins,
                 colour = .data$curve)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 25)) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::labs(
      x      = "% of samples",
      y      = "Number of proteins",
      colour = NULL
    ) +
    theme_jm()
}


#' Extract a joint complete-cases sub-matrix at a chosen operating point
#'
#' Runs the greedy + 2-opt joint complete-cases trace (see
#' [plot_greedy_complete_cases_distribution()]) and returns the
#' sub-matrix at the operating point matching the requested target.
#' Specify exactly one of `n_proteins`, `n_samples`, or `pct_samples`.
#'
#' Finding the best *k* samples to drop is NP-hard so this function
#' is a heuristic and provides no optimality guarantee. The greedy +
#' 2-opt approach is described in
#' [plot_greedy_complete_cases_distribution()].
#'
#' @param abundance Numeric matrix or data frame with proteins as rows
#'   and samples as columns. Column names are required.
#' @param n_proteins Target minimum number of jointly-complete
#'   proteins. The function returns the operating point with the
#'   *largest* sample retention that achieves at least this protein
#'   count.
#' @param n_samples Target number of samples to retain. The trace
#'   evaluates one operating point per `n_samples_kept` value, so this
#'   maps directly to a single row of the trace.
#' @param pct_samples Target percentage of samples to retain.
#'   Converted to `round(ncol(abundance) * pct_samples / 100)` and
#'   treated as `n_samples`.
#'
#' @return A numeric matrix or data frame with proteins as rows and
#'   samples as columns. Rows are restricted to proteins with no
#'   missing values across the kept samples; columns are the kept
#'   samples chosen by greedy + 2-opt.
#'
#' @export
#'
#' @seealso [plot_greedy_complete_cases_distribution()].
#'
#' @examples
#' \dontrun{
#' # Largest sample subset achieving >= 1000 jointly-complete proteins
#' sub <- get_greedy_complete_cases(abundance, n_proteins = 1000)
#'
#' # Drop down to a fixed number of samples
#' sub <- get_greedy_complete_cases(abundance, n_samples = 200)
#'
#' # ... or to a percentage
#' sub <- get_greedy_complete_cases(abundance, pct_samples = 80)
#' }
get_greedy_complete_cases <- function(abundance,
                                      n_proteins  = NULL,
                                      n_samples   = NULL,
                                      pct_samples = NULL) {
  specified <- c(!is.null(n_proteins), !is.null(n_samples), !is.null(pct_samples))
  if (sum(specified) != 1L) {
    stop("Specify exactly one of n_proteins, n_samples, or pct_samples")
  }

  trace <- greedy_complete_cases_trace(abundance)

  if (!is.null(n_proteins)) {
    candidates <- trace[trace$n_complete_proteins >= n_proteins, ]
    if (nrow(candidates) == 0L) {
      stop("Target n_proteins (", n_proteins, ") is unreachable; ",
           "maximum jointly-complete proteins in trace is ",
           max(trace$n_complete_proteins))
    }
    op_row <- candidates[which.max(candidates$n_samples_kept), ]
  } else {
    target <- if (!is.null(n_samples)) {
      n_samples
    } else {
      round(ncol(abundance) * pct_samples / 100)
    }
    if (target < 1L || target > ncol(abundance)) {
      stop("Target sample count (", target, ") must be between 1 and ",
           ncol(abundance))
    }
    op_row <- trace[trace$n_samples_kept == target, ][1L, ]
  }

  samples_kept <- op_row$kept_samples[[1]]

  sub <- abundance[, samples_kept, drop = FALSE]
  sub[rowSums(is.na(sub)) == 0, , drop = FALSE]
}
