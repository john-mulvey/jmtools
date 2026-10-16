# Sample normalisation by pairwise median log-ratios. jmtools package.


#' Normalise samples by pairwise median log-ratios
#'
#' Estimates one log2 scale factor per sample from the median log-ratio
#' of the features each pair of samples shares, reconciled across all
#' pairs by least squares, and subtracts it from that sample. This is
#' the pairwise-ratio step of MaxLFQ (Cox et al. 2014,
#' \doi{10.1074/mcp.M113.031591}) applied to all features at once
#' rather than to the peptides of one protein. It generalises the
#' median-of-ratios scale factor to data with missing values, without
#' restricting to complete cases.
#'
#' For each pair of samples \eqn{(j, k)} sharing at least
#' `min_shared_features` features, \eqn{r_{jk}} is the median of
#' \eqn{x_{ik} - x_{ij}} over the shared features \eqn{i}. The factors
#' \eqn{f} minimise \eqn{\sum_{j<k} (f_k - f_j - r_{jk})^2} subject to
#' \eqn{\sum f = 0}, so normalisation leaves the overall log2 level
#' unchanged.
#'
#' The pairwise ratio is the plain median, as in the iq package. The
#' fallback in Cox et al. 2014 to the ratio of summed intensities when
#' two samples share few features is not applied: between samples
#' measured to different depths that ratio tracks depth rather than
#' loading. The factors were matched to those of `iq::fast_MaxLFQ()`
#' (iq version 2.0.1) run on all features as a single protein, centred
#' on zero, at the time of writing.
#'
#' Factors are only on a common scale when every sample is linked to
#' every other through a chain of pairs sharing at least
#' `min_shared_features` features. If not, the function warns, naming
#' the samples outside the largest linked set, and returns `NULL`.
#'
#' Run time grows with the square of the number of samples times the
#' number of features: about 2 minutes for 800 samples x 40,000
#' features and about 9 minutes for 1,600 x 40,000 on one Apple silicon
#' core. Peak memory is a few times the size of `log2_abundance` (about
#' 2.5 GB at 1,600 x 40,000).
#'
#' @param log2_abundance Numeric matrix of log2 intensities with
#'   features (e.g. precursors) as rows and samples as columns, `NA`
#'   for missing values. Column names are required.
#' @param min_shared_features Minimum number of shared features for a
#'   pair of samples to contribute a ratio. Defaults to 1, as in iq.
#' @param output `"matrix"` (default) to return the normalised matrix,
#'   or `"scale_factors"` to return the scale factors.
#'
#' @return For `output = "matrix"`, `log2_abundance` with each sample's
#'   scale factor subtracted. For `output = "scale_factors"`, a named
#'   numeric vector of log2 scale factors, one per sample, summing to
#'   zero. `NULL` if the samples are not all linked.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' precursor_normalised <- normalise_pairwise_median_ratio(precursor_log2)
#' scale_factors <- normalise_pairwise_median_ratio(precursor_log2,
#'                                                  output = "scale_factors")
#' }
normalise_pairwise_median_ratio <- function(log2_abundance,
                                            min_shared_features = 1,
                                            output = c("matrix", "scale_factors")) {
  output <- match.arg(output)
  if (!is.matrix(log2_abundance) || !is.numeric(log2_abundance)) {
    stop("log2_abundance must be a numeric matrix")
  }
  if (is.null(colnames(log2_abundance))) {
    stop("log2_abundance must have column names")
  }
  if (ncol(log2_abundance) < 2L) {
    stop("log2_abundance must have at least two samples")
  }
  if (any(is.infinite(log2_abundance))) {
    stop("log2_abundance contains infinite values; ",
         "remove zero intensities before the log2 transform")
  }

  n_samples <- ncol(log2_abundance)

  # One row per pair of samples: indices, median log-ratio of the
  # later sample over the earlier one, and the number of shared features
  pair_ratios <- do.call(rbind, lapply(seq_len(n_samples - 1L), function(j) {
    later <- (j + 1L):n_samples
    log_ratios <- log2_abundance[, later, drop = FALSE] - log2_abundance[, j]
    data.frame(sample_j = j,
               sample_k = later,
               log_ratio = matrixStats::colMedians(log_ratios, na.rm = TRUE),
               n_shared = colSums(!is.na(log_ratios)))
  }))
  pair_ratios <- pair_ratios[pair_ratios$n_shared >= min_shared_features, ]

  sample_components <- igraph::components(igraph::make_graph(
    rbind(pair_ratios$sample_j, pair_ratios$sample_k),
    n = n_samples, directed = FALSE))
  if (sample_components$no > 1L) {
    largest <- which.max(sample_components$csize)
    warning("Samples are not all linked by pairs sharing at least ",
            min_shared_features, " features, so their scale factors are ",
            "not on a common scale. Not linked to the largest set: ",
            paste(colnames(log2_abundance)[sample_components$membership != largest],
                  collapse = ", "),
            ". Returning NULL.")
    return(NULL)
  }

  # Normal equations of the least-squares fit (graph Laplacian), bordered
  # by a Lagrange row and column for the sum-to-zero constraint
  normal_matrix <- matrix(0, n_samples + 1L, n_samples + 1L)
  normal_matrix[cbind(pair_ratios$sample_j, pair_ratios$sample_k)] <- -1
  normal_matrix[cbind(pair_ratios$sample_k, pair_ratios$sample_j)] <- -1
  diag(normal_matrix)[seq_len(n_samples)] <-
    tabulate(c(pair_ratios$sample_j, pair_ratios$sample_k), n_samples)
  normal_matrix[n_samples + 1L, seq_len(n_samples)] <- 1
  normal_matrix[seq_len(n_samples), n_samples + 1L] <- 1

  sum_by_sample <- function(sample_index) {
    as.vector(tapply(pair_ratios$log_ratio,
                     factor(sample_index, levels = seq_len(n_samples)),
                     sum, default = 0))
  }
  normal_rhs <- c(sum_by_sample(pair_ratios$sample_k) -
                    sum_by_sample(pair_ratios$sample_j), 0)

  scale_factors <- stats::setNames(
    solve(normal_matrix, normal_rhs)[seq_len(n_samples)],
    colnames(log2_abundance))

  switch(output,
         matrix = sweep(log2_abundance, 2, scale_factors),
         scale_factors = scale_factors)
}
