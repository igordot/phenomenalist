#' Plot dimensionality reduction
#'
#' Plot cell-level reduced dimension results stored in a SpatialExperiment object.
#'
#' @param x A \linkS4class{SpatialExperiment} object.
#' @param dr .
#' @param assay A character string indicating which expression values should be used.
#' @inheritParams plot_scatter
#' @param out_dir Name of the output analysis directory. If specified, the plots will be saved there.
#'
#' @return .
#'
#' @import ggplot2
#' @importFrom methods is
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment assay assayNames
#'
#' @export
plot_dr <- function(x, dr, color_by, assay = "logcounts", smooth = FALSE, range = c(0.01, 0.99), out_dir = NULL) {
  # check if the input is valid
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(dr)) {
    stop("`dr` is not a character string")
  }
  if (!dr %in% reducedDimNames(x)) {
    stop("SpatialExperiment object does not have a dimensionality reduction `", dr, "`")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }

  # check if plotting cell metadata or expression and get the relevant values
  if (all(color_by %in% names(colData(x)))) {
    vals <- as.data.frame(colData(x), stringsAsFactors = FALSE)[color_by]
  } else if (all(color_by %in% rownames(x))) {
    if (!is.character(assay)) {
      stop("`assay` is not a character string")
    }
    if (!assay %in% assayNames(x)) {
      stop("input SpatialExperiment object does not have a `", assay, "` assay")
    }
    vals <- assay(x, i = assay)
    vals <- as.data.frame(t(vals), stringsAsFactors = FALSE)[color_by]
  } else {
    stop("not all `color_by` values are present in the object")
  }

  # get coordinates
  coords <- reducedDim(x, type = dr)
  coords <- coords[, 1:2]
  colnames(coords) <- paste0(dr, 1:2)
  coords <- as.data.frame(coords)

  # add values to the coordinates table
  coords <- cbind(coords, vals[rownames(coords), , drop = FALSE])

  # generate the plot
  for (val in colnames(vals)) {
    p <- plot_scatter(data = coords, x = names(coords)[1], y = names(coords)[2], color_by = val, smooth = smooth, range = range, title = val)
    if (is.null(out_dir)) {
      return(p)
    } else {
      out_base <- glue("{out_dir}/{val}-{dr}")
      if (smooth) {
        out_base <- glue("{out_base}-smooth")
      }
      message(glue("generating {dr} plot for {val}"))
      ggsave(filename = glue("{out_base}.png"), plot = p, width = 8, height = 5)
    }
  }

  return(p)
}
