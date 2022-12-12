#' Perform clustering
#'
#' Cluster similar cells based on their expression profiles.
#'
#' @param x A \linkS4class{SpatialExperiment} object.
#' @param method Clustering method (only `leiden` is implemented currently).
#' @param resolution Value of the parameter controlling the coarseness of the clusters (higher resolution yields more clusters).
#' @param n_neighbors The maximum number of nearest neighbors to compute.
#' @param out_dir Name of the output analysis directory. If specified, the object and the corresponding plots will be saved there.
#'
#' @return A SpatialExperiment object.
#'
#' @import stringr
#' @importFrom glue glue
#' @importFrom igraph cluster_leiden
#' @importFrom janitor tabyl
#' @importFrom readr write_csv
#' @importFrom scran buildSNNGraph
#' @importFrom SummarizedExperiment assay assayNames
#'
#' @export
cluster <- function(x, method = c("leiden"), resolution = 1, n_neighbors = 50, out_dir = NULL) {
  # check if the input is valid
  method <- match.arg(method)
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.numeric(resolution)) {
    stop("`resolution` is not a number")
  }
  if (!is.numeric(n_neighbors)) {
    stop("`n_neighbors` is not a number")
  }

  if (!"exprs" %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `exprs` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
    clusters_dir <- glue("{out_dir}/clusters")
    dir.create(clusters_dir)
  }

  # get the expression matrix
  exprs_mat <- assay(x, "exprs")

  # clustering using igraph::cluster_leiden() (added in v1.2.7)
  if (method == "leiden") {
    # build a nearest-neighbor graph
    # matrix is transposed if rows are cells
    g <- scran::buildSNNGraph(exprs_mat, transposed = FALSE, k = n_neighbors)

    # repeat clustering for every resolution value
    n_clust_prev <- 0
    for (res_num in resolution) {
      message(glue("clustering using resolution of {res_num}"))

      # perform clustering
      set.seed(99)
      clusters <- igraph::cluster_leiden(g, objective_function = "modularity", resolution_parameter = res_num, n_iterations = 10)
      clusters <- clusters$membership

      # generate a label for the cluster set
      res_str <- format(as.numeric(res_num), nsmall = 1)
      res_str <- stringr::str_pad(res_str, width = 3, side = "left", pad = "0")
      res_str <- stringr::str_c("res", res_str)
      n_clust <- length(unique(clusters))
      clusters_label <- glue("cluster_{method}_{res_str}_clust{n_clust}")

      # pad with "C" to avoid downstream numeric conversions
      clusters <- as.character(clusters)
      if (n_clust < 100) {
        clusters <- stringr::str_pad(clusters, width = 2, side = "left", pad = "0")
      } else {
        clusters <- stringr::str_pad(clusters, width = nchar(n_clust), side = "left", pad = "0")
      }
      clusters <- stringr::str_c("C", clusters)

      # add clusters to the original object
      # skip if the increased resolution did not yield more clusters
      if (n_clust > n_clust_prev) {
        x[[clusters_label]] <- factor(clusters)
        n_clust_prev <- n_clust
        # generate summary and plots if output directory is specified
        if (!is.null(out_dir)) {
          # summary table
          cluster_summary <- janitor::tabyl(x[[clusters_label]])
          colnames(cluster_summary) <- c("cluster", "cells_num", "cells_freq")
          readr::write_csv(cluster_summary, glue("{clusters_dir}/{clusters_label}-summary.csv"))
          # plots (skip if too many clusters)
          if (n_clust < 150) {
            plot_heatmap(x, group_by = clusters_label, out_dir = clusters_dir)
            plot_spatial(x, color_by = clusters_label, out_dir = clusters_dir)
            plot_dr(x, dr = "UMAP", color_by = clusters_label, out_dir = clusters_dir)
          }
        }
      }
    }
  }

  # clustering using leiden::leiden() (uses python leidenalg via reticulate)
  # if (method == "leiden") {
  #
  #   # find nearest neighbors for each cell
  #   exprs_mat <- assay(x, "exprs")
  #   snn <- RANN::nn2(t(exprs_mat), k = n_neighbors)$nn.idx
  #
  #   # generate an adjacency matrix of nearest neighbors between samples
  #   adj_mat <- matrix(0L, ncol(exprs_mat), ncol(exprs_mat))
  #   rownames(adj_mat) <- colnames(adj_mat) <- colnames(exprs_mat)
  #   for (i in 1:ncol(exprs_mat)) {
  #     adj_mat[i, colnames(exprs_mat)[snn[i, ]]] <- 1L
  #   }
  #
  #   # repeat clustering for every resolution value
  #   for (res_num in resolution) {
  #     # perform clustering
  #     clusters <- leiden::leiden(adj_mat, resolution_parameter = res_num, n_iterations = 10, seed = 99)
  #   }
  # }

  # save object if output directory is specified
  if (!is.null(out_dir)) {
    message("saving object")
    saveRDS(x, paste0(out_dir, "/spe.rds"))
  }

  return(x)
}
