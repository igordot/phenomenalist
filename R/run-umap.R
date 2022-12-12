#' Run UMAP dimensionality reduction
#'
#' @param x A \linkS4class{SpatialExperiment} object.
#' @param assay .
#' @inheritParams uwot::umap
#' @param out_dir Name of the output analysis directory. If specified, the object will be saved there.
#'
#' @return A SpatialExperiment object.
#'
#' @importFrom methods is
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom uwot umap
#'
#' @export
run_umap <- function(x, assay = "exprs", n_neighbors = 50, min_dist = 0.01, n_threads = NULL, out_dir = NULL) {
  # check if the input is valid
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!assay %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `", assay, "` assay")
  }

  # get expression matrix
  e <- assay(x, assay)

  # get expression matrix
  set.seed(99)
  umap_mat <- uwot::umap(t(e), n_neighbors = n_neighbors, min_dist = min_dist, n_threads = n_threads)
  umap_mat <- round(umap_mat, 5)

  reducedDim(x, type = "UMAP") <- umap_mat

  # save object if output directory is specified
  if (!is.null(out_dir)) {
    message("saving object")
    saveRDS(x, paste0(out_dir, "/spe.rds"))
  }

  return(x)
}
