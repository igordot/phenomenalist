#' Plot cells on their spatial coordinates
#'
#' Plot cell-level physical coordinates stored in a SpatialExperiment object.
#'
#' @param x A \linkS4class{SpatialExperiment} object.
#' @param assay A character string indicating which expression values should be used.
#' @inheritParams plot_scatter
#' @param out_dir Name of the output analysis directory. If specified, the plots will be saved there.
#'
#' @return .
#'
#' @import ggplot2 SpatialExperiment
#' @importFrom glue glue
#' @importFrom SummarizedExperiment assay assayNames colData
#'
#' @export
plot_spatial <- function(x, color_by, assay = "logcounts", smooth = FALSE, range = c(0.01, 0.99), out_dir = NULL) {

  # check if the input is valid
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }

  # check if plotting cell metadata or expression and get the relevant values
  if (all(color_by %in% names(colData(x)))) {
    vals <- as.data.frame(colData(x)[color_by], stringsAsFactors = FALSE)
  } else if (all(color_by %in% rownames(x))) {
    if (!is.character(assay)) {
      stop("`assay` is not a character string")
    }
    if (!assay %in% assayNames(x)) {
      stop("input SpatialExperiment object does not have a `", assay, "` assay")
    }
    vals <- assay(x, i = assay)
    vals <- t(vals)[, color_by, drop = FALSE]
    vals <- as.data.frame(vals, stringsAsFactors = FALSE)
  } else {
    stop("not all `color_by` values are present in the object")
  }

  # create a sub-directory for the output
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE)
  }

  # get spatial coordinates
  coords <- spatialCoords(x)
  coords <- as.data.frame(coords)

  # reverse `y` so 0 is at the top of the plot
  coords$y <- coords$y * -1

  # add values to the coordinates table
  coords <- cbind(coords, vals[rownames(coords), , drop = FALSE])

  # get aspect ratio
  ratio <- max(coords$y) / max(coords$x)

  # generate the plot
  for (val in colnames(vals)) {
    p <- plot_scatter(data = coords, x = "x", y = "y", color_by = val, smooth = smooth, range = range, title = val, aspect_ratio = ratio)
    if (is.null(out_dir)) {
      return(p)
    } else {
      out_base <- glue("{out_dir}/{val}-spatial")
      if (smooth) {
        out_base <- glue("{out_base}-smooth")
      }
      message(glue("generating spatial plot for {val}"))
      ggsave(filename = glue("{out_base}.png"), plot = p, width = 8, height = 5)
    }
  }
}
