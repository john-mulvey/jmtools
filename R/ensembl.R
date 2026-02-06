# Ensembl BioMart Local Database
# jmtools package
#
# Load gene annotations from a locally downloaded Ensembl BioMart TSV file.
# See ~/Documents/databases/ensembl_biomart/download_ensembl_biomart.sh for
# the download script.

#' Load Ensembl BioMart gene annotations from local TSV
#'
#' Reads a locally downloaded Ensembl BioMart TSV file and returns a tibble
#' with standardised column names. The file is expected to live at
#' `~/Documents/databases/ensembl_biomart/` and follow the naming convention
#' `ensembl_<species>_biomart_<date>.tsv`.
#'
#' @param species Character string specifying the species. Must be one of
#'   `"human"` or `"mouse"`. Default: `"human"`.
#' @param file Character string specifying a particular filename (not full
#'   path) to load, e.g. `"ensembl_human_biomart_20260206.tsv"`. When `NULL`
#'   (default), the most recently modified file matching the species is used.
#' @param database_dir Character string specifying the directory containing
#'   the BioMart TSV files. Default: `"~/Documents/databases/ensembl_biomart"`.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene stable ID (e.g., ENSG00000141510)}
#'   \item{gene_name}{HGNC gene symbol (e.g., TP53)}
#'   \item{description}{Gene description from Ensembl}
#'   \item{uniprotswissprot}{UniProtKB/Swiss-Prot accession (e.g., P04637)}
#'   \item{entrezgene_id}{NCBI Entrez gene ID (e.g., 7157)}
#'   \item{gene_biotype}{Gene biotype (e.g., protein_coding, lncRNA)}
#'   \item{chromosome_name}{Chromosome or scaffold name}
#'   \item{uniprot_gn_symbol}{UniProtKB gene name symbol}
#' }
#'
#' @details
#' This function provides a stable alternative to the `biomaRt` R package,
#' which can be unreliable due to dependency conflicts and session management
#' issues. The underlying data are downloaded directly from the Ensembl
#' BioMart web service using curl (see the download script in the database
#' directory).
#'
#' When `file` is not specified, the most recently modified file matching the
#' species is used and its filename is printed so the notebook output records
#' which version was loaded. To pin a specific version for reproducibility,
#' pass the filename via the `file` parameter.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the latest human annotations
#' ensembl <- load_ensembl_mapping()
#'
#' # Pin a specific version
#' ensembl <- load_ensembl_mapping(file = "ensembl_human_biomart_20260206.tsv")
#'
#' # Load mouse annotations
#' ensembl_mouse <- load_ensembl_mapping("mouse")
#'
#' # Filter to UniProt accessions in your dataset
#' ensembl %>%
#'   dplyr::filter(uniprotswissprot %in% my_protein_ids) %>%
#'   dplyr::select(uniprotswissprot, entrezgene_id, ensembl_gene_id) %>%
#'   dplyr::distinct()
#'
#' # Map Ensembl IDs to gene names
#' ensembl %>%
#'   dplyr::filter(ensembl_gene_id %in% my_ensembl_ids) %>%
#'   dplyr::select(ensembl_gene_id, gene_name) %>%
#'   dplyr::distinct()
#' }
load_ensembl_mapping <- function(species = "human",
                                 file = NULL,
                                 database_dir = "~/Documents/databases/ensembl_biomart") {

  species <- match.arg(species, choices = c("human", "mouse"))
  database_dir <- path.expand(database_dir)

  if (!dir.exists(database_dir)) {
    stop("Database directory not found: ", database_dir, "\n",
         "Run the download script at ",
         "~/Documents/databases/ensembl_biomart/download_ensembl_biomart.sh")
  }

  if (!is.null(file)) {
    tsv_path <- file.path(database_dir, file)
    if (!file.exists(tsv_path)) {
      stop("File not found: ", tsv_path)
    }
  } else {
    pattern <- paste0("ensembl_", species, "_biomart_.*\\.tsv$")
    files <- list.files(database_dir, pattern = pattern, full.names = TRUE)

    if (length(files) == 0) {
      stop("No BioMart TSV found for species '", species,
           "' in ", database_dir, "\n",
           "Expected filename pattern: ensembl_", species,
           "_biomart_<date>.tsv")
    }

    tsv_path <- files[order(file.mtime(files), decreasing = TRUE)[1]]
  }

  message("Loading: ", basename(tsv_path))

  col_names <- c("ensembl_gene_id", "gene_name", "description",
                 "uniprotswissprot", "entrezgene_id", "gene_biotype",
                 "chromosome_name", "uniprot_gn_symbol")

  result <- utils::read.delim(tsv_path, header = FALSE, skip = 1,
                              col.names = col_names,
                              stringsAsFactors = FALSE,
                              quote = "")

  tibble::as_tibble(result)
}
