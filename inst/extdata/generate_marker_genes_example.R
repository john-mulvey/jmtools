# Generate example marker genes dataset for jmtools package
#
# Trims the full Heart Cell Atlas v2 Level 1 marker gene CSV to the top 20
# markers per cell type by AUC, retaining all cell type rows for each selected
# gene (so avg_expr and pct_1 are available for all clusters in dot plots).
#
# Run this script manually when the group drive is mounted. The output is
# committed to inst/extdata/ and used by the package vignette.

full_path <- paste0(
  "/Volumes/groupdir/SUN-BMI-CardiacProteomics/scRNAseq_datasets/JM/",
  "atlas_snRNAseq/results/heart_cell_atlas_v2/",
  "heart_cell_atlas_v2_all_cells_annotated_predicted_celltype_l1_marker_genes.csv"
)

marker_genes_raw <- read.csv(full_path, stringsAsFactors = FALSE)

cat("Full dataset:", nrow(marker_genes_raw), "rows,",
    length(unique(marker_genes_raw$cluster)), "cell types\n")

# Select top 20 markers per cell type by AUC (stored as power = AUC - 0.5)
top_genes_per_cluster <- tapply(
  seq_len(nrow(marker_genes_raw)),
  marker_genes_raw$cluster,
  function(idx) {
    cluster_rows <- marker_genes_raw[idx, ]
    cluster_rows <- cluster_rows[order(-cluster_rows$power), ]
    head(cluster_rows$gene, 20)
  }
)
top_genes <- unique(unlist(top_genes_per_cluster))

cat("Top 20 markers per cell type:", length(top_genes), "unique genes\n")

# Retain all rows for those genes so expression is available for all clusters
example_marker_genes <- marker_genes_raw[marker_genes_raw$gene %in% top_genes, ]

cat("Example dataset:", nrow(example_marker_genes), "rows\n")

out_path <- "example_marker_genes.csv"
write.csv(example_marker_genes, out_path, row.names = FALSE)
cat("Written to:", out_path, "\n")
