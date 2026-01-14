# Gene Set Enrichment Analysis Functions
# jmtools package

#' Calculate Jaccard similarity between two sets
#'
#' @param set1 First set (vector)
#' @param set2 Second set (vector)
#' @return Numeric Jaccard similarity coefficient (0-1)
#' @keywords internal
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  if (union_size == 0) return(0)
  intersection / union_size
}

#' Calculate overlap coefficient between two sets
#'
#' Also known as the Szymkiewicz-Simpson coefficient. Less sensitive to
#' differences in set size than Jaccard, making it better suited for
#' hierarchically-related gene sets where one may be a subset of another.
#'
#' @param set1 First set (vector)
#' @param set2 Second set (vector)
#' @return Numeric overlap coefficient (0-1)
#' @keywords internal
overlap_coefficient <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  min_size <- min(length(set1), length(set2))
  if (min_size == 0) return(0)
  intersection / min_size
}

#' Cluster enrichment results by gene overlap similarity
#'
#' Clusters significant gene sets based on similarity of their
#' overlapping genes. Works with both GSEA results (using leading edge genes)
#' and ORA results (using overlap genes). Returns the input data frame with
#' added cluster information.
#'
#' @param enrichment_results Data frame of enrichment results from GSEA or ORA,
#'   typically with cleaned column names from the boilerplate functions.
#' @param gene_set_col Character string specifying the column containing gene
#'   set names. Default: `"gene_set"`.
#' @param score_col Character string specifying the column containing the
#'   enrichment score (NES for GSEA, enrichment_score for ORA). Default: `"nes"`.
#'   For ORA results, use `"enrichment_score"`.
#' @param adj_pval_col Character string specifying the column containing
#'   adjusted p-values for filtering. Default: `"adj_p_value"`.
#' @param genes_col Character string specifying the column containing the
#'   genes driving enrichment (as a list column). For GSEA use `"leading_edge"`,
#'   for ORA use `"overlap_genes"`. Default: `"leading_edge"`.
#' @param n_genes_col Character string specifying the column containing the
#'   number of genes in each gene set. Default: `"n_genes"`. For ORA results,
#'   use `"n_gene_set"`.
#' @param adj_p_threshold Adjusted p-value threshold for filtering gene sets.
#'   Default: 0.05.
#' @param similarity_threshold Minimum similarity to consider two gene
#'   sets as belonging to the same cluster. Default: 0.3.
#' @param similarity_metric Character string specifying which similarity metric
#'   to use. One of `"jaccard"` (Jaccard index) or `"overlap"` (overlap
#'   coefficient / Szymkiewicz-Simpson). The overlap coefficient is less
#'   sensitive to differences in set size, making it better suited for
#'   hierarchically-related gene sets. Default: `"jaccard"`.
#'
#' @return A data frame containing all significant gene sets (one row per gene
#'   set) with added columns:
#' \describe{
#'   \item{cluster}{Integer cluster ID}
#'   \item{cluster_label}{Representative name for the cluster (the gene set
#'     with the largest number of genes in that cluster)}
#' }
#'
#' @details
#' Clustering is performed using the Walktrap algorithm, which detects
#' communities based on random walks through the similarity network. Gene sets
#' that are densely connected (high Jaccard similarity) tend to end up in the
#' same cluster. This method can identify meaningful subcommunities even within
#' connected components.
#'
#' The cluster label is the name of the gene set with the largest number of
#' genes in that cluster.
#'
#' Only gene sets meeting the significance threshold are returned. Non-significant
#' gene sets are excluded from the output.
#'
#' This function works with both GSEA and ORA results:
#' \itemize{
#'   \item **GSEA**: Use `genes_col = "leading_edge"`, `score_col = "nes"`,
#'     `n_genes_col = "n_genes"`
#'   \item **ORA**: Use `genes_col = "overlap_genes"`, `score_col = "enrichment_score"`,
#'     `n_genes_col = "n_gene_set"`
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Cluster GSEA results
#' clustered_gsea <- cluster_enrichment_results(gsea_results)
#'
#' # Cluster ORA results (note different column names)
#' clustered_ora <- cluster_enrichment_results(
#'   ora_results,
#'   score_col = "enrichment_score",
#'   genes_col = "overlap_genes",
#'   n_genes_col = "n_gene_set"
#' )
#'
#' # View cluster assignments
#' head(clustered_gsea[, c("gene_set", "nes", "cluster", "cluster_label")])
#'
#' # Plot as bar chart
#' plot_enrichment_bars(clustered_gsea)
#' }
#'
#' @seealso [plot_enrichment_bars()] for plotting clustered results
cluster_enrichment_results <- function(enrichment_results,
                                  gene_set_col = "gene_set",
                                  score_col = "nes",
                                  adj_pval_col = "adj_p_value",
                                  genes_col = "leading_edge",
                                  n_genes_col = "n_genes",
                                  adj_p_threshold = 0.05,
                                  similarity_threshold = 0.3,
                                  similarity_metric = c("jaccard", "overlap")) {

  similarity_metric <- match.arg(similarity_metric)

  # Select similarity function
 similarity_fn <- switch(similarity_metric,
    jaccard = jaccard_similarity,
    overlap = overlap_coefficient
  )

  # Validate inputs
  required_cols <- c(gene_set_col, score_col, adj_pval_col, genes_col, n_genes_col)
  missing_cols <- setdiff(required_cols, names(enrichment_results))
  if (length(missing_cols) > 0) {
    stop("Missing columns in enrichment_results: ", paste(missing_cols, collapse = ", "))
  }

  # Filter by significance
  sig_data <- enrichment_results[enrichment_results[[adj_pval_col]] < adj_p_threshold, ]

  if (nrow(sig_data) == 0) {
    warning("No gene sets meet the significance threshold (adj_p < ", adj_p_threshold, ")")
    return(NULL)
  }

  if (nrow(sig_data) == 1) {
    # Single gene set - no clustering needed
    sig_data$cluster <- 1L
    sig_data$cluster_label <- sig_data[[gene_set_col]]
    return(sig_data)
  }

  # Extract genes driving enrichment
  gene_sets <- sig_data[[genes_col]]
  names(gene_sets) <- sig_data[[gene_set_col]]

  # Calculate pairwise similarities and build adjacency matrix
  n_sets <- nrow(sig_data)
  adj_matrix <- matrix(0, nrow = n_sets, ncol = n_sets)

  for (i in seq_len(n_sets - 1)) {
    for (j in (i + 1):n_sets) {
      sim <- similarity_fn(gene_sets[[i]], gene_sets[[j]])
      if (sim >= similarity_threshold) {
        adj_matrix[i, j] <- sim
        adj_matrix[j, i] <- sim
      }
    }
  }

  # Clustering via Walktrap algorithm (community detection based on random walks)
  graph <- igraph::graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  walktrap_result <- igraph::cluster_walktrap(graph, weights = igraph::E(graph)$weight)
  cluster_ids <- as.integer(igraph::membership(walktrap_result))

  sig_data$cluster <- cluster_ids

  # Determine cluster labels (gene set with largest n_genes in each cluster)
  cluster_labels <- vapply(sort(unique(cluster_ids)), function(cl) {
    cluster_members <- which(cluster_ids == cl)
    n_genes_values <- sig_data[[n_genes_col]][cluster_members]
    sig_data[[gene_set_col]][cluster_members[which.max(n_genes_values)]]
  }, character(1))

  names(cluster_labels) <- as.character(sort(unique(cluster_ids)))
  sig_data$cluster_label <- cluster_labels[as.character(cluster_ids)]

  # Sort by cluster, then by score within cluster
  sig_data <- sig_data[order(sig_data$cluster, -abs(sig_data[[score_col]])), ]

  return(sig_data)
}


#' Plot clustered enrichment results as a bar plot
#'
#' Creates a bar plot of clustered enrichment results (GSEA or ORA) showing
#' mean enrichment score with individual gene set values overlaid as points.
#'
#' @param clustered_results Data frame from [cluster_enrichment_results()],
#'   containing gene sets with `cluster` and `cluster_label` columns.
#' @param score_col Character string specifying the column containing enrichment
#'   scores. Default: `"nes"` for GSEA; use `"enrichment_score"` for ORA.
#' @param n_top Number of top clusters to show (split between enriched and
#'   depleted). Default: 20.
#' @param title Optional plot title. Default: `"Enrichment Results (Clustered)"`.
#' @param x_label Label for the x-axis. Default: `"Enrichment Score"`.
#' @param colours Named vector of colours for enriched and depleted clusters.
#'   Default: `c(Enriched = "firebrick", Depleted = "steelblue")`.
#' @param wrap_labels Number of characters at which to wrap long cluster labels.
#'   Set to `NULL` to disable wrapping. Default: 30.
#' @param show_n Logical. If `TRUE`, appends the number of gene sets in each
#'   cluster to the label. Default: `TRUE`.
#' @param point_size Size of the individual gene set points. Default: 2.
#' @param point_alpha Alpha transparency of individual points. Default: 0.6.
#'
#' @return A ggplot2 object.
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item **Bars**: Mean enrichment score for each cluster
#'   \item **Points**: Individual gene set score values, jittered vertically
#'   \item **Y-axis**: Cluster labels, ordered by score (most positive at top)
#'   \item **Colour**: Direction of enrichment (Enriched or Depleted)
#' }
#'
#' This function works with both GSEA and ORA results:
#' \itemize{
#'   \item **GSEA**: Use `score_col = "nes"`, `x_label = "Normalised Enrichment Score (NES)"`
#'   \item **ORA**: Use `score_col = "enrichment_score"`, `x_label = "Enrichment Score (log2 fold enrichment)"`
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # GSEA workflow
#' clustered_gsea <- cluster_enrichment_results(gsea_results)
#' plot_enrichment_bars(clustered_gsea)
#'
#' # ORA workflow
#' clustered_ora <- cluster_enrichment_results(
#'   ora_results,
#'   score_col = "enrichment_score",
#'   genes_col = "overlap_genes",
#'   n_genes_col = "n_gene_set"
#' )
#' plot_enrichment_bars(
#'   clustered_ora,
#'   score_col = "enrichment_score",
#'   x_label = "Enrichment Score (log2 fold enrichment)"
#' )
#' }
#'
#' @seealso [cluster_enrichment_results()] for clustering enrichment results,
#'   [plot_gsea_dots()] for plotting unclustered GSEA results
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_vline scale_fill_manual position_jitter facet_grid labs theme element_text
#' @importFrom rlang .data
#' @importFrom stats sd setNames
plot_enrichment_bars <- function(clustered_results,
                            score_col = "nes",
                            n_top = 20,
                            title = "Enrichment Results (Clustered)",
                            x_label = "Enrichment Score",
                            colours = c(Enriched = "firebrick", Depleted = "steelblue"),
                            wrap_labels = 30,
                            show_n = TRUE,
                            point_size = 2,
                            point_alpha = 0.6) {

  if (is.null(clustered_results)) {
    warning("No clustered results to plot")
    return(NULL)
  }

  if (!is.data.frame(clustered_results)) {
    stop("clustered_results must be a data frame from cluster_enrichment_results()")
  }

  if (!all(c("cluster", "cluster_label", score_col) %in% names(clustered_results))) {
    stop("clustered_results must contain 'cluster', 'cluster_label', and '",
         score_col, "' columns")
  }

  if (nrow(clustered_results) == 0) {
    warning("No data to plot")
    return(NULL)
  }

  # Calculate summary statistics per cluster
  clusters <- unique(clustered_results$cluster)

  summary_list <- lapply(clusters, function(cl) {
    cluster_data <- clustered_results[clustered_results$cluster == cl, ]
    score_values <- cluster_data[[score_col]]
    n <- length(score_values)

    data.frame(
      cluster = cl,
      cluster_label = cluster_data$cluster_label[1],
      n_gene_sets = n,
      mean_score = mean(score_values),
      direction = ifelse(mean(score_values) > 0, "Enriched", "Depleted"),
      stringsAsFactors = FALSE
    )
  })

  summary_data <- do.call(rbind, summary_list)

  # Ensure direction is a factor
  summary_data$direction <- factor(summary_data$direction, levels = c("Depleted", "Enriched"))

  # Select top clusters (split between enriched and depleted)
  enriched <- summary_data[summary_data$direction == "Enriched", ]
  depleted <- summary_data[summary_data$direction == "Depleted", ]

  n_each <- ceiling(n_top / 2)

  if (nrow(enriched) > 0) {
    enriched <- enriched[order(-abs(enriched$mean_score)), ]
    enriched <- head(enriched, n_each)
  }

  if (nrow(depleted) > 0) {
    depleted <- depleted[order(-abs(depleted$mean_score)), ]
    depleted <- head(depleted, n_each)
  }

  summary_data <- rbind(enriched, depleted)

  if (nrow(summary_data) == 0) {
    warning("No clusters to plot after filtering")
    return(NULL)
  }

  # Filter individual points to only include selected clusters
  selected_clusters <- summary_data$cluster
  points_data <- clustered_results[clustered_results$cluster %in% selected_clusters, ]

  # Add direction to points data
  points_data$direction <- ifelse(
    points_data$cluster %in% summary_data$cluster[summary_data$direction == "Enriched"],
    "Enriched",
    "Depleted"
  )
  points_data$direction <- factor(points_data$direction, levels = c("Depleted", "Enriched"))

  # Add n to labels if requested
  if (show_n) {
    summary_data$display_label <- paste0(summary_data$cluster_label, " (n=", summary_data$n_gene_sets, ")")
  } else {
    summary_data$display_label <- summary_data$cluster_label
  }

  # Create mapping from cluster to display_label for points
  label_map <- setNames(summary_data$display_label, summary_data$cluster)
  points_data$display_label <- label_map[as.character(points_data$cluster)]

  # Wrap long labels if requested
  if (!is.null(wrap_labels)) {
    summary_data$display_label <- sapply(
      summary_data$display_label,
      function(x) paste(strwrap(x, width = wrap_labels), collapse = "\n")
    )
    points_data$display_label <- sapply(
      points_data$display_label,
      function(x) paste(strwrap(x, width = wrap_labels), collapse = "\n")
    )
  }


  # Order by score (descending so most positive at top of plot)
  summary_data <- summary_data[order(-summary_data$mean_score), ]
  label_levels <- rev(unique(summary_data$display_label))

  summary_data$display_label <- factor(summary_data$display_label, levels = label_levels)
  points_data$display_label <- factor(points_data$display_label, levels = label_levels)

  # Create plot
  p <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(
      x = .data$mean_score,
      y = .data$display_label,
      fill = .data$direction
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_col(width = 0.7, alpha = 0.7) +
    ggplot2::geom_point(
      data = points_data,
      ggplot2::aes(x = .data[[score_col]], y = .data$display_label),
      inherit.aes = FALSE,
      size = point_size,
      alpha = point_alpha,
      colour = "grey20",
      position = ggplot2::position_jitter(width = 0, height = 0.15, seed = 42)
    ) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::labs(
      x = x_label,
      y = NULL,
      title = title
    ) +
    theme_jm()

  return(p)
}


#' Plot GSEA results as a dot plot
#'
#' Creates a dot plot visualisation of Gene Set Enrichment Analysis results,
#' showing top enriched and depleted gene sets with NES on the x-axis and
#' significance indicated by point size.
#'
#' @param gsea_results Data frame of GSEA results, typically from fgsea with
#'   cleaned column names (e.g., output from `generate_gsea_boilerplate()`).
#' @param gene_set_col Character string specifying the column containing gene
#'   set names. Default: "gene_set".
#' @param nes_col Character string specifying the column containing Normalised
#'   Enrichment Scores. Default: "nes".
#' @param pval_col Character string specifying the column containing p-values
#'   for the size aesthetic. Default: "p_value".
#' @param adj_pval_col Character string specifying the column containing
#'   adjusted p-values for filtering. Default: "adj_p_value".
#' @param n_top Number of top gene sets to show (split between enriched and
#'   depleted). Default: 20 (i.e., up to 10 enriched and 10 depleted).
#' @param adj_p_threshold Adjusted p-value threshold for filtering gene sets.
#'   Default: 0.05.
#' @param title Optional plot title. Default: "GSEA Results".
#' @param colours Named vector of colours for enriched and depleted gene sets.
#'   Default: `c(Enriched = "firebrick", Depleted = "steelblue")`.
#' @param wrap_labels Number of characters at which to wrap long gene set names.
#'   Set to `NULL` to disable wrapping. Default: 30.
#'
#' @return A ggplot2 object.
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item **X-axis**: Normalised Enrichment Score (NES)
#'   \item **Y-axis**: Gene set names, ordered by NES (most positive at top)
#'   \item **Point size**: `-log10(p-value)`, so larger points are more significant
#'   \item **Point colour**: Direction of enrichment (Enriched or Depleted)
#' }
#'
#' Gene sets are filtered by adjusted p-value, then the top N are selected
#' based on absolute NES (split evenly between enriched and depleted where
#' possible).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot GSEA results with default settings
#' plot_gsea_dots(gsea_results)
#'
#' # Customise number of gene sets and colours
#' plot_gsea_dots(
#'   gsea_results,
#'   n_top = 30,
#'   colours = c(Enriched = "darkred", Depleted = "darkblue")
#' )
#'
#' # Use different column names
#' plot_gsea_dots(
#'   gsea_results,
#'   gene_set_col = "pathway",
#'   nes_col = "NES",
#'   pval_col = "pval",
#'   adj_pval_col = "padj"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline scale_colour_manual scale_size_continuous facet_grid labs theme element_text
#' @importFrom rlang .data
#' @importFrom utils head
plot_gsea_dots <- function(gsea_results,
                               gene_set_col = "gene_set",
                               nes_col = "nes",
                               pval_col = "p_value",
                               adj_pval_col = "adj_p_value",
                               n_top = 20,
                               adj_p_threshold = 0.05,
                               title = "GSEA Results",
                               colours = c(Enriched = "firebrick", Depleted = "steelblue"),
                               wrap_labels = 30) {

  # Validate inputs
  if (!gene_set_col %in% names(gsea_results)) {
    stop("Column '", gene_set_col, "' not found in gsea_results")
  }
  if (!nes_col %in% names(gsea_results)) {
    stop("Column '", nes_col, "' not found in gsea_results")
  }
  if (!pval_col %in% names(gsea_results)) {
    stop("Column '", pval_col, "' not found in gsea_results")
  }
  if (!adj_pval_col %in% names(gsea_results)) {
    stop("Column '", adj_pval_col, "' not found in gsea_results")
  }

  # Filter by significance
  plot_data <- gsea_results[gsea_results[[adj_pval_col]] < adj_p_threshold, ]

  if (nrow(plot_data) == 0) {
    warning("No gene sets meet the significance threshold (adj_p < ", adj_p_threshold, ")")
    return(NULL)
  }

  # Add direction column
  plot_data$direction <- ifelse(
    plot_data[[nes_col]] > 0,
    "Enriched",
    "Depleted"
  )
  plot_data$direction <- factor(plot_data$direction, levels = c("Depleted", "Enriched"))

  # Select top gene sets (split between enriched and depleted)
  enriched <- plot_data[plot_data$direction == "Enriched", ]
  depleted <- plot_data[plot_data$direction == "Depleted", ]

  n_each <- ceiling(n_top / 2)

  if (nrow(enriched) > 0) {
    enriched <- enriched[order(-abs(enriched[[nes_col]])), ]
    enriched <- head(enriched, n_each)
  }

  if (nrow(depleted) > 0) {
    depleted <- depleted[order(-abs(depleted[[nes_col]])), ]
    depleted <- head(depleted, n_each)
  }

  plot_data <- rbind(enriched, depleted)

  if (nrow(plot_data) == 0) {
    warning("No gene sets to plot after filtering")
    return(NULL)
  }

  # Wrap long labels if requested
  if (!is.null(wrap_labels)) {
    plot_data[[gene_set_col]] <- sapply(
      plot_data[[gene_set_col]],
      function(x) paste(strwrap(x, width = wrap_labels), collapse = "\n")
    )
  }

  # Order gene sets by NES (descending so most positive at top of plot)
  plot_data <- plot_data[order(-plot_data[[nes_col]]), ]
  plot_data[[gene_set_col]] <- factor(
    plot_data[[gene_set_col]],
    levels = rev(unique(plot_data[[gene_set_col]]))
  )

  # Calculate -log10(p-value) for size
  plot_data$neg_log10_p <- -log10(plot_data[[pval_col]])

  # Create plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[nes_col]],
      y = .data[[gene_set_col]],
      size = .data$neg_log10_p,
      colour = .data$direction
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_point() +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_size_continuous(name = expression(-log[10] ~ p - value)) +
    ggplot2::labs(
      x = "Normalised Enrichment Score (NES)",
      y = NULL,
      title = title
    ) +
    theme_jm()

  return(p)
}


#' Plot GSEA results as a network graph
#'
#' Creates a network visualisation of Gene Set Enrichment Analysis results,
#' where nodes are gene sets and edges represent similarity between their
#' leading edge genes. Clusters are detected and can be labelled directly
#' on the plot using convex hulls.
#'
#' @param gsea_results Data frame of GSEA results, typically from fgsea with
#'   cleaned column names (e.g., output from `generate_gsea_boilerplate()`).
#' @param gene_set_col Character string specifying the column containing gene
#'   set names. Default: "gene_set".
#' @param nes_col Character string specifying the column containing Normalised
#'   Enrichment Scores. Default: "nes".
#' @param adj_pval_col Character string specifying the column containing
#'   adjusted p-values for filtering. Default: "adj_p_value".
#' @param leading_edge_col Character string specifying the column containing
#'   leading edge genes (as a list column). Default: "leading_edge".
#' @param n_genes_col Character string specifying the column containing the
#'   number of genes in each gene set. Default: "n_genes".
#' @param adj_p_threshold Adjusted p-value threshold for filtering gene sets.
#'   Default: 0.05.
#' @param similarity_threshold Minimum Jaccard similarity to draw an edge
#'   between two gene sets. Default: 0.3.
#' @param min_cluster_size Minimum number of gene sets in a cluster to draw
#'   a hull and label. Default: 3.
#' @param exclude_singletons Logical. If `TRUE`, singleton clusters (size 1)
#'   are excluded from the plot. Default: `FALSE`.
#' @param show_labels Logical. If `TRUE`, shows cluster labels on the plot.
#'   Default: `TRUE`.
#' @param label_clusters_by Method for selecting the representative label for
#'   each cluster. One of "pagerank" (highest PageRank centrality) or
#'   "n_genes" (largest gene set). Default: "pagerank".
#' @param layout Layout algorithm for the network. Default: "fr"
#'   (Fruchterman-Reingold).
#' @param title Optional plot title. Default: `NULL`.
#' @param node_size_range Numeric vector of length 2 specifying the range of
#'   node sizes. Default: `c(2, 8)`.
#' @param concavity Concavity parameter for hull shapes. Higher values create
#'   more convex hulls. Default: 4.
#' @param expand_hull Amount to expand hulls around points. Default:
#'
#'   `ggplot2::unit(3, "mm")`.
#' @param seed Random seed for layout reproducibility. Default: 42.
#'
#' @return A ggplot2 object. Returns `NULL` with a warning if fewer than 2
#'   gene sets meet the significance threshold.
#'
#'   **Cluster assignments** are stored as an attribute and can be accessed via
#'   `attr(plot, "clusters")`. This returns a data frame with columns:
#'   \describe{
#'     \item{gene_set}{Name of the gene set}
#'     \item{cluster}{Integer cluster ID (from Walktrap algorithm)}
#'     \item{nes}{Normalised Enrichment Score}
#'     \item{cluster_label}{Representative name for clusters meeting
#'       `min_cluster_size` (NA for smaller clusters)}
#'   }
#'
#' @details
#' The network visualisation shows:
#' \itemize{
#'   \item **Nodes**: Gene sets, sized by number of genes, coloured by NES
#'   \item **Edges**: Connections between gene sets with Jaccard similarity
#'     above the threshold, with edge width proportional to similarity
#'   \item **Clusters**: Detected using the Walktrap algorithm, shown as
#'     coloured hulls with labels
#' }
#'
#' Edge weights are calculated as the Jaccard similarity between the leading
#' edge genes of each pair of gene sets. This captures functional similarity
#' based on the genes driving the enrichment signal.
#'
#' Cluster labels are determined by selecting the most "central" gene set
#' within each cluster, either by PageRank centrality or by gene set size.
#'
#' **Note:** Both this function and [cluster_enrichment_results()] use the
#' Walktrap algorithm for clustering, so cluster assignments should be
#' consistent between the two (given the same similarity threshold).
#'
#' @section Required Packages:
#' This function requires the following packages to be installed:
#' \itemize{
#'   \item tidygraph
#'   \item ggraph
#'   \item ggforce (for `geom_mark_hull`)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot GSEA network with default settings
#' plot_gsea_network(gsea_results)
#'
#' # Adjust similarity threshold and cluster size
#' plot_gsea_network(
#'   gsea_results,
#'   similarity_threshold = 0.2,
#'   min_cluster_size = 5
#' )
#'
#' # Without cluster labels
#' plot_gsea_network(gsea_results, show_labels = FALSE)
#'
#' # Access cluster assignments
#' p <- plot_gsea_network(gsea_results)
#' clusters <- attr(p, "clusters")
#' head(clusters)
#' }
#'
#' @seealso [plot_gsea_dots()] for a simpler dot plot visualisation
#'
#' @importFrom ggplot2 aes unit scale_fill_distiller scale_size_continuous labs guides guide_legend theme element_blank element_rect
#' @importFrom rlang .data
plot_gsea_network <- function(gsea_results,
                               gene_set_col = "gene_set",
                               nes_col = "nes",
                               adj_pval_col = "adj_p_value",
                               leading_edge_col = "leading_edge",
                               n_genes_col = "n_genes",
                               adj_p_threshold = 0.05,
                               similarity_threshold = 0.3,
                               min_cluster_size = 3,
                               exclude_singletons = FALSE,
                               show_labels = TRUE,
                               label_clusters_by = c("pagerank", "n_genes"),
                               layout = "fr",
                               title = NULL,
                               node_size_range = c(2, 8),
                               concavity = 4,
                               expand_hull = ggplot2::unit(3, "mm"),
                               seed = 42) {

  # Check for required packages
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("Package 'tidygraph' is required for plot_gsea_network(). ",
         "Please install it with: install.packages('tidygraph')")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Package 'ggraph' is required for plot_gsea_network(). ",
         "Please install it with: install.packages('ggraph')")
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required for plot_gsea_network(). ",
         "Please install it with: install.packages('ggforce')")
  }

  label_clusters_by <- match.arg(label_clusters_by)

  # Validate inputs
  required_cols <- c(gene_set_col, nes_col, adj_pval_col, leading_edge_col, n_genes_col)
  missing_cols <- setdiff(required_cols, names(gsea_results))
  if (length(missing_cols) > 0) {
    stop("Missing columns in gsea_results: ", paste(missing_cols, collapse = ", "))
  }

  # Filter by significance
  plot_data <- gsea_results[gsea_results[[adj_pval_col]] < adj_p_threshold, ]

  if (nrow(plot_data) < 2) {
    warning("Fewer than 2 gene sets meet the significance threshold. ",
            "Cannot create network plot.")
    return(NULL)
  }

  # Create nodes data frame
  nodes <- data.frame(
    name = plot_data[[gene_set_col]],
    nes = plot_data[[nes_col]],
    n_genes = plot_data[[n_genes_col]],
    stringsAsFactors = FALSE
  )

  # Extract leading edge genes
  leading_edges <- plot_data[[leading_edge_col]]
  names(leading_edges) <- plot_data[[gene_set_col]]

  # Calculate pairwise Jaccard similarities
  n_sets <- nrow(nodes)
  edges_list <- list()

  for (i in seq_len(n_sets - 1)) {
    for (j in (i + 1):n_sets) {
      sim <- jaccard_similarity(leading_edges[[i]], leading_edges[[j]])
      if (sim >= similarity_threshold) {
        edges_list[[length(edges_list) + 1]] <- data.frame(
          from = nodes$name[i],
          to = nodes$name[j],
          similarity = sim,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(edges_list) == 0) {
    warning("No edges meet the similarity threshold. ",
            "Try lowering similarity_threshold.")
    # Create graph with no edges
    edges <- data.frame(from = character(0), to = character(0),
                        similarity = numeric(0), stringsAsFactors = FALSE)
  } else {
    edges <- do.call(rbind, edges_list)
  }

  # Create tidygraph object
  graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE,
    node_key = "name"
  )

  # Add cluster membership and centrality
  graph <- graph |>
    tidygraph::activate("nodes") |>
    tidygraph::mutate(
      cluster = as.factor(tidygraph::group_walktrap()),
      pagerank = tidygraph::centrality_pagerank()
    )

  # Extract node data with layout coordinates
  set.seed(seed)
  layout_coords <- ggraph::create_layout(graph, layout = layout)

  # Determine cluster labels
  node_data <- as.data.frame(layout_coords)

  # Count cluster sizes and identify clusters to show
  cluster_sizes <- table(node_data$cluster)
  singleton_clusters <- names(cluster_sizes)[cluster_sizes == 1]
  large_clusters <- names(cluster_sizes)[cluster_sizes >= min_cluster_size]

  # Mark singletons for exclusion from plotting
 node_data$is_singleton <- node_data$cluster %in% singleton_clusters

  # Select representative label for each large cluster
  if (show_labels && length(large_clusters) > 0) {
    cluster_labels <- vapply(large_clusters, function(cl) {
      cluster_nodes <- node_data[node_data$cluster == cl, ]
      if (label_clusters_by == "pagerank") {
        cluster_nodes$name[which.max(cluster_nodes$pagerank)]
      } else {
        cluster_nodes$name[which.max(cluster_nodes$n_genes)]
      }
    }, character(1))

    node_data$cluster_label <- ifelse(
      node_data$cluster %in% large_clusters,
      cluster_labels[as.character(node_data$cluster)],
      NA_character_
    )

    # Only show label for the representative node
    node_data$is_label_node <- node_data$name == node_data$cluster_label
    node_data$display_label <- ifelse(
      node_data$cluster %in% large_clusters,
      cluster_labels[as.character(node_data$cluster)],
      NA_character_
    )
  } else {
    node_data$cluster_label <- NA_character_
    node_data$is_label_node <- FALSE
    node_data$display_label <- NA_character_
    large_clusters <- character(0)
  }

  # Filter node data for hulls (only large clusters) and drop unused factor levels
  hull_data <- node_data[node_data$cluster %in% large_clusters, ]
  hull_data$cluster <- droplevels(hull_data$cluster)

  # Filter node data for plotting (optionally exclude singletons)
  if (exclude_singletons) {
    nodes_to_plot <- node_data[!node_data$is_singleton, ]
  } else {
    nodes_to_plot <- node_data
  }

  # Build the plot
  p <- ggraph::ggraph(layout_coords) +
    # Edges (only between non-singleton nodes)
    ggraph::geom_edge_link(
      ggplot2::aes(width = .data$similarity),
      alpha = 0.3,
      colour = "grey60"
    ) +
    ggraph::scale_edge_width_continuous(
      name = "Jaccard\nsimilarity",
      range = c(0.2, 2),
      breaks = c(0.3, 0.5, 0.7, 0.9)
    )

  # Add hulls for large clusters
  if (nrow(hull_data) > 0 && show_labels) {
    p <- p +
      ggforce::geom_mark_hull(
        data = hull_data,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$cluster,
          label = .data$display_label
        ),
        fill = "grey80",
        concavity = concavity,
        expand = expand_hull,
        alpha = 0.15,
        colour = NA,
        show.legend = FALSE
      )
  }

  # Add nodes (excluding singletons)
  p <- p +
    ggplot2::geom_point(
      data = nodes_to_plot,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        size = .data$n_genes,
        fill = .data$nes
      ),
      shape = 21,
      colour = "grey30",
      stroke = 0.3
    ) +
    ggplot2::scale_fill_distiller(
      palette = "RdYlBu",
      name = "NES",
      direction = -1
    ) +
    ggplot2::scale_size_continuous(
      name = "Gene set\nsize",
      range = node_size_range
    ) +
    ggraph::theme_graph() +
    ggplot2::labs(title = title) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(barwidth = 0.8, barheight = 8)
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )

  # Create cluster assignments data frame
  cluster_df <- data.frame(
    gene_set = node_data$name,
    cluster = as.integer(node_data$cluster),
    nes = node_data$nes,
    stringsAsFactors = FALSE
  )

  # Add cluster labels (representative gene set name for each cluster)
  if (length(large_clusters) > 0) {
    cluster_df$cluster_label <- cluster_labels[as.character(node_data$cluster)]
  } else {
    cluster_df$cluster_label <- NA_character_
  }

  cluster_df <- cluster_df[order(cluster_df$cluster, -abs(cluster_df$nes)), ]

  # Attach as attribute
 attr(p, "clusters") <- cluster_df

  return(p)
}
