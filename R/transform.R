#' Transform/normalize expression data
#'
#' Transform/normalize expression data.
#' The options are currently:
#' \itemize{
#'   \item \code{log}: log10 + 1
#'   \item \code{z}: Z normalization of each marker (recommended by Hickey et al.)
#' }
#'
#' @param x A \linkS4class{SpatialExperiment} object.
#' @param method A character string indicating which transformation method should be used.
#' @param out_dir Name of the output directory. If specified, the density plots will be saved there.
#'
#' @return A \linkS4class{SpatialExperiment} object with transformed values in the \code{exprs} assay.
#'
#' @import SpatialExperiment
#' @importFrom methods is
#' @importFrom SingleCellExperiment counts logcounts logcounts<-
#' @importFrom SummarizedExperiment assay assay<- assayNames
#'
#' @export
transform <- function(x, method = c("log", "z"), out_dir = NULL) {
  # check if the input is valid
  method <- match.arg(method)
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!"counts" %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `counts` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }

  # generate log-transformed values to use for visualizations
  if (!"logcounts" %in% assayNames(x)) {
    # plot expression density
    if (!is.null(out_dir)) {
      plot_distribution(x, assay = "counts", out_dir = out_dir)
    }
    logcounts(x) <- log10(counts(x) + 1)
    if (!is.null(out_dir)) {
      plot_distribution(x, assay = "logcounts", out_dir = out_dir)
    }
  }

  if (method == "log") {
    # log-transformed values
    assay(x, "exprs") <- logcounts(x)
  } else if (method == "z") {
    # z-norm of each marker based on raw "counts"
    # recommended by Hickey 2021 and used in Akoya Seurat demo
    # z <- logcounts(x)
    z <- counts(x)
    # scale() centers/scales the columns of a numeric matrix
    z <- t(scale(t(z)))
    # limit outliers
    z[z > 5] <- 5
    z[z < -5] <- -5
    assay(x, "exprs", FALSE) <- z
    # plot expression density
    if (!is.null(out_dir)) {
      plot_distribution(x, assay = "exprs", out_dir = out_dir)
    }
  }

  return(x)
}
