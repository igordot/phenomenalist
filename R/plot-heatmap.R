#' Generate a heatmap
#'
#' @param x A \linkS4class{SpatialExperiment} object.
#' @param group_by Column metadata field(s) to group by.
#' @param assay A character string indicating which expression values should be used.
#' @param out_dir Name of the output analysis directory. If specified, the plots will be saved there.
#'
#' @return A \linkS4class{Heatmap} object.
#'
#' @importFrom ComplexHeatmap draw Heatmap
#' @importFrom ggsci pal_igv
#' @importFrom glue glue
#' @importFrom grDevices dev.off png
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment assay assayNames colData
#'
#' @export
plot_heatmap <- function(x, group_by, assay = "logcounts", out_dir = NULL) {

  # check if the input is valid
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(group_by)) {
    stop("`group_by` is not a character string")
  }
  if (!all(group_by %in% names(colData(x)))) {
    stop("not all `group_by` values are present in the object")
  }
  if (!is.character(assay)) {
    stop("`assay` is not a character string")
  }
  if (!assay %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `", assay, "` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }

  # color scheme
  gradient_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  # discrete_colors <- ggsci::pal_igv("default")(51)

  # generate the plot
  for (g in group_by) {

    # get expression matrix (clusters are columns)
    e <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]], assay.type = assay, statistics = "median")
    e <- SummarizedExperiment::assay(e, i = "median")

    # rotate and scale each marker
    # clusters are now rows to allow for wider plots (usually there are more markers than clusters)
    e <- scale(t(e))

    # limit outliers
    e[e > 5] <- 5
    e[e < -5] <- -5

    # generate the plot
    hm <-
      ComplexHeatmap::Heatmap(
        e,
        name = "Expression", row_title = g, col = gradient_colors,
        cluster_rows = TRUE, cluster_columns = TRUE
      )

    # return the plot or save if output is specified
    if (is.null(out_dir)) {
      return(hm)
    } else {
      out_base <- glue("{out_dir}/{g}-heatmap")
      png(filename = glue("{out_base}.png"), width = 10, height = 5, units = "in", res = 300)
      ComplexHeatmap::draw(hm)
      dev.off()
    }
  }
}
