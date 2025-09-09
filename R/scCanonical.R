#' Identify Conserved Markers Across All Clusters
#'
#' This function identifies conserved markers for each cluster in a Seurat object
#' across different groups defined by a grouping variable. It uses Seurat's
#' `FindConservedMarkers` function and filters results based on specified criteria.
#'
#' @param obj A Seurat object containing single-cell data.
#' @param grouping.var Character string specifying the metadata column used to
#'   define groups (default: "condition").
#' @param min.pct Numeric value specifying the minimum fraction of cells in which
#'   a gene must be expressed to be considered (default: 0.1).
#' @param logfc.threshold Numeric value specifying the minimum log fold-change
#'   threshold for differential expression (default: 0.25).
#' @param pval_meta_cutoff Numeric value specifying the p-value cutoff for
#'   filtering conserved markers (default: 0.05).
#' @param min.cells.per.group Integer specifying the minimum number of cells per
#'   group required to perform the conserved marker test (default: 3).
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{results}: A data frame containing conserved marker genes for
#'       each cluster, including gene names, cluster identifiers, and relevant
#'       statistics. Returns an empty data frame if no results meet the criteria.
#'     \item \code{skipped}: A data frame listing clusters that were skipped due
#'       to insufficient cells or missing groups, including cluster IDs, group
#'       names, and cell counts. Returns an empty data frame if no clusters were
#'       skipped.
#'   }
#'
#' @details
#' The function sets the identities of the Seurat object to "seurat_clusters" and
#' iterates through each cluster to identify conserved markers using
#' `FindConservedMarkers`. Clusters with insufficient cells (fewer than
#' `min.cells.per.group`) in any group or missing groups are skipped, and a message
#' is printed. Results are filtered based on the specified p-value cutoff
#' (`pval_meta_cutoff`) using the appropriate p-value column (`max_pval`,
#' `p_val_adj`, or `p_val`) depending on the Seurat version. The function uses
#' `dplyr::bind_rows` to combine results and skipped cluster information.
#'
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object
#' result <- get_conserved_for_all(
#'   obj = seurat_obj,
#'   grouping.var = "condition",
#'   min.pct = 0.1,
#'   logfc.threshold = 0.25,
#'   pval_meta_cutoff = 0.05,
#'   min.cells.per.group = 3
#' )
#' print(result$results)
#' print(result$skipped)
#' }
#'
#' @importFrom Seurat Idents Idents<- WhichCells FindConservedMarkers
#' @importFrom dplyr bind_rows
#' @export
get_conserved_for_all <- function(obj,
                                  grouping.var = "condition",
                                  min.pct = 0.1,
                                  logfc.threshold = 0.25,
                                  pval_meta_cutoff = 0.05,
                                  min.cells.per.group = 3) {
  Idents(obj) <- "seurat_clusters"
  clusters <- levels(Idents(obj))
  groups <- unique(obj[[grouping.var]][, 1])

  res_list <- list()
  skipped <- list()

  for (cl in clusters) {
    # Check that this cluster exists with enough cells in every group
    cl_cells <- WhichCells(obj, idents = cl)
    cl_groups <- obj[[grouping.var]][cl_cells, 1]
    tab <- table(cl_groups)

    if (length(tab) < length(groups) || any(tab < min.cells.per.group)) {
      skipped[[cl]] <- data.frame(
        cluster = cl,
        group = names(tab),
        n_cells = as.integer(tab)
      )
      message("Skipping conserved test for cluster ", cl,
              " (some groups missing or < ", min.cells.per.group, " cells)")
      next
    }

    cm <- tryCatch({
      FindConservedMarkers(
        object = obj,
        ident.1 = cl,
        grouping.var = grouping.var,
        only.pos = TRUE,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
      )
    }, error = function(e) {
      message("FindConservedMarkers failed for cluster ", cl, ": ", e$message)
      NULL
    })
    if (is.null(cm)) next

    # Normalize to data.frame & add identifiers
    cm <- as.data.frame(cm)
    cm$gene <- rownames(cm)
    cm$cluster <- cl

    # Filter by a meta p column if present (Seurat versions differ)
    if ("max_pval" %in% colnames(cm)) {
      cm <- subset(cm, max_pval < pval_meta_cutoff)
    } else if ("p_val_adj" %in% colnames(cm)) {
      cm <- subset(cm, p_val_adj < pval_meta_cutoff)
    } else if ("p_val" %in% colnames(cm)) {
      cm <- subset(cm, p_val < pval_meta_cutoff)
    } # else keep as is

    if (nrow(cm) > 0) {
      res_list[[cl]] <- cm
    }
  }

  list(
    results = if (length(res_list)) dplyr::bind_rows(res_list) else data.frame(),
    skipped = if (length(skipped)) dplyr::bind_rows(skipped) else data.frame()
  )
}

#' Calculate Specificity Score for Differential Expression Results
#'
#' This function calculates a specificity score for each gene in a differential
#' expression data frame, based on the log2 fold-change and the difference in
#' expression percentages between two groups.
#'
#' @param df A data frame containing differential expression results, typically
#'   from Seurat's `FindConservedMarkers` or similar functions. It must include
#'   columns `avg_log2FC`, `pct.1`, and `pct.2` for the specificity score to be
#'   calculated; otherwise, the score will be set to `NA`.
#'
#' @return A data frame identical to the input `df`, with an additional column
#'   `spec_score` containing the specificity score for each gene. The score is
#'   calculated as `avg_log2FC * max(0, pct.1 - pct.2)`. If the required columns
#'   (`avg_log2FC`, `pct.1`, `pct.2`) are not present, `spec_score` is set to
#'   `NA_real_`.
#'
#' @details
#' The specificity score quantifies the degree of differential expression by
#' combining the log2 fold-change (`avg_log2FC`) with the difference in the
#' percentage of cells expressing the gene in the target group (`pct.1`) versus
#' the reference group (`pct.2`). The `pmax(0, pct.1 - pct.2)` ensures that only
#' positive differences contribute to the score, emphasizing genes more specific
#' to the target group.
#'
#' @examples
#' \dontrun{
#' # Example data frame with differential expression results
#' de_results <- data.frame(
#'   gene = c("Gene1", "Gene2"),
#'   avg_log2FC = c(1.5, 0.8),
#'   pct.1 = c(0.9, 0.4),
#'   pct.2 = c(0.2, 0.5)
#' )
#' result <- add_specificity(de_results)
#' print(result)
#' }
#'
#' @export
add_specificity <- function(df) {
  if (all(c("avg_log2FC", "pct.1", "pct.2") %in% colnames(df))) {
    df$spec_score <- df$avg_log2FC * pmax(0, df$pct.1 - df$pct.2)
  } else {
    df$spec_score <- NA_real_
  }
  df
}

#' Plot Canonical Markers for Clusters
#'
#' This function generates a faceted scatter plot of canonical markers for each
#' cluster, displaying the difference in expression percentage (`delta_pct`) on
#' the x-axis and specificity score (`spec_score`) on the y-axis. Canonical
#' markers are highlighted with distinct colors and labeled, while non-canonical
#' markers are shown in the background.
#'
#' @param cons_joined A data frame containing differential expression results for
#'   all genes, with columns `delta_pct`, `spec_score`, and `cluster`.
#' @param canonical A data frame containing canonical marker genes, with columns
#'   `delta_pct`, `spec_score`, `gene`, and `cluster`.
#' @param custom_colors A named vector of colors for each cluster, used for
#'   coloring points and labels in the plot.
#'
#' @return A `ggplot` object representing the faceted scatter plot of canonical
#'   markers, which is also printed to the active graphics device.
#'
#' @details
#' The function creates a scatter plot using `ggplot2`, with background points for
#' all genes in `cons_joined` and highlighted points for canonical markers in
#' `canonical`. Points are colored by cluster, and gene labels are added for
#' canonical markers using `geom_text_repel` to avoid overlap. The plot is faceted
#' by cluster with free y-axis scales. Custom y-axis breaks are calculated for each
#' cluster, and a "-" marker is added at the left border of each facet panel to
#' indicate these breaks. The plot uses a minimal theme with customized aesthetics
#' for clarity and readability.
#'
#' @examples
#' \dontrun{
#' # Example data
#' cons_joined <- data.frame(
#'   delta_pct = runif(100, -0.5, 0.5),
#'   spec_score = runif(100, 0, 2),
#'   cluster = rep(1:4, each = 25)
#' )
#' canonical <- data.frame(
#'   delta_pct = runif(10, 0, 0.5),
#'   spec_score = runif(10, 1, 2),
#'   gene = paste0("Gene", 1:10),
#'   cluster = rep(1:4, length.out = 10)
#' )
#' custom_colors <- c("1" = "red", "2" = "blue", "3" = "green", "4" = "purple")
#'
#' # Generate plot
#' plot_canonicals(cons_joined, canonical, custom_colors)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_manual labs theme_minimal theme element_blank element_rect element_text
#' @importFrom dplyr group_by do
#' @importFrom ggrepel geom_text_repel
#' @export
plot_canonicals <- function(cons_joined, canonical, custom_colors) {
  # Create base plot
  p <- ggplot2::ggplot(cons_joined, ggplot2::aes(x = delta_pct, y = spec_score)) +
    ggplot2::geom_point(color = "gray70", size = 1.5, alpha = 0.6) +
    ggplot2::geom_point(data = canonical, ggplot2::aes(color = factor(cluster)), size = 2.5) +
    ggrepel::geom_text_repel(
      data = canonical,
      ggplot2::aes(label = gene, color = factor(cluster)),
      size = 3, box.padding = 0.4, max.overlaps = Inf
    ) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y") +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::labs(
      x = expression(Delta~"Percentage Difference (pct.1 - pct.2)"),
      y = "Specificity Score",
      color = "Cluster"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.8),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "none"
    )

  # Calculate y-axis breaks and add "-" markers
  y_breaks_df <- cons_joined %>%
    dplyr::group_by(cluster) %>%
    dplyr::do({
      yb <- pretty(.$spec_score, 4)
      data.frame(spec_score = yb)
    })

  p <- p + ggplot2::geom_text(
    data = y_breaks_df,
    ggplot2::aes(x = -Inf, y = spec_score, label = "-"),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5, size = 4, color = "black"
  )

  # Print and return the plot
  print(p)
  return(p)
}

#' Plot Canonical Markers for Clusters in a Single Row
#'
#' This function generates a faceted scatter plot of canonical markers for each
#' cluster, with all facets arranged in a single row. The plot displays the
#' difference in expression percentage (`delta_pct`) on the x-axis and specificity
#' score (`spec_score`) on the y-axis. Canonical markers are highlighted with
#' distinct colors and labeled, while non-canonical markers are shown in the
#' background.
#'
#' @param cons_joined A data frame containing differential expression results for
#'   all genes, with columns `delta_pct`, `spec_score`, and `cluster`.
#' @param canonical A data frame containing canonical marker genes, with columns
#'   `delta_pct`, `spec_score`, `gene`, and `cluster`.
#' @param custom_colors A named vector of colors for each cluster, used for
#'   coloring points and labels in the plot.
#'
#' @return A `ggplot` object representing the faceted scatter plot of canonical
#'   markers, with all facets in a single row and fixed x-axis labels. The plot is
#'   also printed to the active graphics device.
#'
#' @details
#' The function creates a scatter plot using `ggplot2`, with background points for
#' all genes in `cons_joined` and highlighted points for canonical markers in
#' `canonical`. Points are colored by cluster, and gene labels are added for
#' canonical markers using `geom_text_repel` to avoid overlap. The plot is faceted
#' by cluster in a single row with free y-axis scales. The x-axis is fixed with
#' breaks at 0, 0.5, and 1, and limits from 0 to 1. Custom y-axis breaks are
#' calculated for each cluster, and a "-" marker is added at the left border of
#' each facet panel to indicate these breaks. The plot uses a minimal theme with
#' customized aesthetics for clarity and readability.
#'
#' @examples
#' \dontrun{
#' # Example data
#' cons_joined <- data.frame(
#'   delta_pct = runif(100, 0, 1),
#'   spec_score = runif(100, 0, 2),
#'   cluster = rep(1:4, each = 25)
#' )
#' canonical <- data.frame(
#'   delta_pct = runif(10, 0, 1),
#'   spec_score = runif(10, 1, 2),
#'   gene = paste0("Gene", 1:10),
#'   cluster = rep(1:4, length.out = 10)
#' )
#' custom_colors <- c("1" = "red", "2" = "blue", "3" = "green", "4" = "purple")
#'
#' # Generate plot
#' plot_canonicals_inline(cons_joined, canonical, custom_colors)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_manual scale_x_continuous labs theme_minimal theme element_blank element_rect element_text
#' @importFrom dplyr group_by do
#' @importFrom ggrepel geom_text_repel
#' @export
plot_canonicals_inline <- function(cons_joined, canonical, custom_colors) {
  # Calculate y-axis breaks
  y_breaks_df <- cons_joined %>%
    dplyr::group_by(cluster) %>%
    dplyr::do({
      yb <- pretty(.$spec_score, 4)
      data.frame(spec_score = yb)
    })

  # Create base plot
  p <- ggplot2::ggplot(cons_joined, ggplot2::aes(x = delta_pct, y = spec_score)) +
    ggplot2::geom_point(color = "gray70", size = 1.5, alpha = 0.6) +
    ggplot2::geom_point(data = canonical, ggplot2::aes(color = factor(cluster)), size = 2.5) +
    ggrepel::geom_text_repel(
      data = canonical,
      ggplot2::aes(label = gene, color = factor(cluster)),
      size = 3, box.padding = 0.4, max.overlaps = Inf
    ) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y", nrow = 1) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    ggplot2::labs(
      x = expression(Delta~"Percentage Difference (pct.1 - pct.2)"),
      y = "Specificity Score",
      color = "Cluster"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.8),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::geom_text(
      data = y_breaks_df,
      ggplot2::aes(x = -Inf, y = spec_score, label = "-"),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5, size = 4, color = "black"
    )

  # Print and return the plot
  print(p)
  return(p)
}
