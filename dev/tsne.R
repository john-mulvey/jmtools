# t-SNE Functions
# jmtools package
#
# STATUS: UNDER TESTING - DO NOT USE IN PRODUCTION
# Testing location: /Users/John/Documents/Projects/staged/20260130_plot_for_egfr/analysis/
# Date: 2026-02-02

#' Run t-SNE using Kobak & Berens (2019) recommendations
#'
#' Computes t-SNE embeddings following the best practices described in
#' Kobak & Berens (2019) "The art of using t-SNE for single-cell transcriptomics"
#' (Nature Communications). Uses FIt-SNE with:
#' \itemize{
#'   \item Multi-scale perplexity (combining 30 and n/100)
#'   \item PCA initialisation (scaled to SD = 0.0001)
#'   \item High learning rate (n/12)
#' }
#'
#' @section Dataset size (n > 100,000):
#' For datasets larger than 100,000 cells, set `large_dataset = TRUE`. This
#' enables late exaggeration (factor of 4) as recommended by Kobak & Berens
#' to prevent cluster expansion. Note that for very large datasets (n > 1M),
#' downsampling-based initialisation may also be beneficial but is not
#' implemented here.
#'
#' @param object A Seurat object. Must have PCA already computed.
#' @param dims Dimensions to use from PCA embeddings. Default is `1:50`.
#' @param reduction Name of the dimensional reduction to use. Default is `"pca"`.
#' @param reduction.name Name for the new t-SNE reduction. Default is `"tsne"`.
#' @param reduction.key Key prefix for the new reduction. Default is `"tSNE_"`.
#' @param perplexity_small Small perplexity value for multi-scale. Default is 30.
#' @param perplexity_large Large perplexity value for multi-scale. If `NULL`
#'   (default), uses n/100 where n is the number of cells.
#' @param max_iter Maximum number of iterations. Default is 1000.
#' @param large_dataset Logical. If `TRUE`, enables late exaggeration (factor 4)
#'   for datasets > 100,000 cells. Default is `FALSE`; set to `TRUE` for large
#'   datasets to prevent cluster expansion.
#' @param seed Random seed for reproducibility. Default is 42.
#' @param nthreads Number of threads for parallelisation. Default is 0 (auto).
#' @param fast_tsne_path Path to the FIt-SNE executable. If `NULL` (default),
#'   uses [get_fast_tsne_path()] to find it automatically.
#' @param verbose Logical. Print progress messages? Default is `TRUE`.
#'
#' @return The Seurat object with a new dimensional reduction added (named by
#'   `reduction.name`).
#'
#' @details
#' This function implements the Kobak & Berens (2019) recommendations exactly:
#'
#' \itemize{
#'   \item **Multi-scale perplexity**: Combines perplexity 30 (local structure)
#'     with n/100 (global structure) by averaging their similarity matrices.
#'   \item **PCA initialisation**: The first two PCs are scaled to SD = 0.0001
#'     and used as starting positions, preserving global structure.
#'   \item **High learning rate**: n/12 ensures good convergence and prevents
#'     cluster fragmentation.
#'   \item **Late exaggeration**: For large datasets, a factor of 4 prevents
#'     cluster expansion (enabled via `large_dataset = TRUE`).
#' }
#'
#' @references
#' Kobak D, Berens P (2019). "The art of using t-SNE for single-cell
#' transcriptomics." Nature Communications, 10:5416.
#' \doi{10.1038/s41467-019-13056-x}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard usage
#' seurat_obj <- run_tsne_kobak(seurat_obj)
#' DimPlot(seurat_obj, reduction = "tsne")
#'
#' # For large datasets (> 100,000 cells)
#' seurat_obj <- run_tsne_kobak(seurat_obj, large_dataset = TRUE)
#'
#' # With custom reduction name
#' seurat_obj <- run_tsne_kobak(seurat_obj, reduction.name = "tsne_kobak")
#' DimPlot(seurat_obj, reduction = "tsne_kobak")
#' }
run_tsne_kobak <- function(object,
                         dims = 1:50,
                         reduction = "pca",
                         reduction.name = "tsne",
                         reduction.key = "tSNE_",
                         perplexity_small = 30,
                         perplexity_large = NULL,
                         max_iter = 1000,
                         large_dataset = FALSE,
                         seed = 42,
                         nthreads = 0,
                         fast_tsne_path = NULL,
                         verbose = TRUE) {

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required")
  }

  # Get path to FIt-SNE executable
  if (is.null(fast_tsne_path)) {
    fast_tsne_path <- get_fast_tsne_path()
  }

  # Source the FIt-SNE R wrapper
  fit_sne_wrapper <- file.path(dirname(dirname(fast_tsne_path)), "fast_tsne.R")
  if (!file.exists(fit_sne_wrapper)) {
    stop("Cannot find FIt-SNE R wrapper at: ", fit_sne_wrapper)
  }
  source(fit_sne_wrapper, chdir = TRUE, local = TRUE)

  # Extract PCA embeddings
  input_matrix <- Seurat::Embeddings(object, reduction)[, dims]
  n_cells <- nrow(input_matrix)

  # Calculate large perplexity (n/100) if not specified
  if (is.null(perplexity_large)) {
    perplexity_large <- n_cells / 100
  }

  # Ensure perplexities are valid
  max_perplexity <- (n_cells - 1) / 3
  if (perplexity_large > max_perplexity) {
    if (verbose) {
      message("Capping large perplexity from ", round(perplexity_large, 1),
              " to ", round(max_perplexity, 1))
    }
    perplexity_large <- max_perplexity
  }

  # Multi-scale perplexity list
  perplexity_list <- c(perplexity_small, perplexity_large)

  # Warn if dataset is large but large_dataset not set
  if (n_cells > 100000 && !large_dataset) {
    warning(
      "Dataset has ", n_cells, " cells (> 100,000). ",
      "Consider setting large_dataset = TRUE to enable late exaggeration."
    )
  }

  if (verbose) {
    message("Running FIt-SNE with Kobak & Berens (2019) parameters")
    message("  Number of cells: ", n_cells)
    message("  Multi-scale perplexity: ", perplexity_small, " + ", round(perplexity_large, 1))
    message("  Learning rate: auto (n/12 = ", round(n_cells / 12, 1), ")")
    if (large_dataset) {
      message("  Late exaggeration: 4")
    }
  }

  # Run FIt-SNE
  # perplexity = 0 triggers multi-scale mode using perplexity_list
  # learning_rate = 'auto' gives n/12 (Kobak & Berens recommendation)
  # initialization = 'pca' scales to 0.0001 * SD (Kobak & Berens recommendation)
  tsne_result <- fftRtsne(
    X = input_matrix,
    dims = 2,
    perplexity = 0,  # Use perplexity_list for multi-scale
    perplexity_list = perplexity_list,
    max_iter = max_iter,
    learning_rate = "auto",  # n/12, minimum 200
    initialization = "pca",  # Scaled to 0.0001 * SD
    late_exag_coeff = if (large_dataset) 4 else -1,
    start_late_exag_iter = if (large_dataset) 250 else -1,
    rand_seed = seed,
    nthreads = nthreads,
    fast_tsne_path = fast_tsne_path
  )

  # Create embedding matrix with proper row/column names
  rownames(tsne_result) <- colnames(object)
  colnames(tsne_result) <- paste0(reduction.key, 1:2)

  # Add to Seurat object
  object[[reduction.name]] <- Seurat::CreateDimReducObject(
    embeddings = tsne_result,
    key = reduction.key,
    assay = Seurat::DefaultAssay(object)
  )

  object
}
