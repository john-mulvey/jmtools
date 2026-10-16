test_that("noise-free per-sample offsets are recovered exactly", {
  # Without noise every shared feature gives the same log-ratio for a pair,
  # so each pairwise median is the true offset difference and the fit
  # returns the true offsets, centred
  withr::local_seed(1)
  sample_offsets <- c(-1.5, -0.2, 0, 0.4, 0.9, 2.1, -0.7, 1.3)
  log2_abundance <- outer(rnorm(200, 10, 2), sample_offsets, `+`)
  log2_abundance[runif(length(log2_abundance)) < 0.5] <- NA
  colnames(log2_abundance) <- paste0("sample_", seq_along(sample_offsets))

  scale_factors <- normalise_pairwise_median_ratio(log2_abundance,
                                                   output = "scale_factors")

  expect_equal(scale_factors,
               stats::setNames(sample_offsets - mean(sample_offsets),
                               colnames(log2_abundance)))
})

test_that("samples not all linked give a warning and NULL", {
  withr::local_seed(1)
  log2_abundance <- matrix(rnorm(100 * 6, 10, 2), 100, 6,
                           dimnames = list(NULL, paste0("sample_", 1:6)))

  # Two blocks of samples sharing no features
  blocks <- log2_abundance
  blocks[1:50, 4:6] <- NA
  blocks[51:100, 1:3] <- NA
  expect_warning(blocks_normalised <- normalise_pairwise_median_ratio(blocks),
                 "not all linked")
  expect_null(blocks_normalised)

  # Every pair shares all 100 features, below the required 101
  expect_warning(
    threshold_normalised <- normalise_pairwise_median_ratio(
      log2_abundance, min_shared_features = 101),
    "not all linked")
  expect_null(threshold_normalised)
})
