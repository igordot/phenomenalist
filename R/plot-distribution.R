#' Plot expression distribution
#'
#' Plot distribution of expression intensities, split by marker.
#'
#' @param x A [SpatialExperiment-class] object.
#' @param assay A character string indicating which values should be used.
#' @param out_dir Name of the output analysis directory. If specified, the plots will be saved there.
#'
#' @return .
#'
#' @import ggplot2
#' @importFrom glue glue
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#'
#' @export
plot_distribution <- function(x, assay = "exprs", out_dir = NULL) {
  # check if the input is valid
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
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
    # create a sub-directory for the output
    out_dir <- paste0(out_dir, "/expression-distribution")
    dir.create(out_dir, showWarnings = FALSE)
    out_base <- glue("{out_dir}/expression-distribution-{assay}")
  }

  # adjust plot title based on assay
  if (assay == "counts") {
    plot_title <- "Original Input Expression Values"
  } else if (assay == "logcounts") {
    plot_title <- "Log-Transformed Expression Values"
  } else if (assay == "exprs") {
    plot_title <- "Transformed Expression Values"
  } else {
    plot_title <- "Expression Values"
  }

  # get expression matrix (cells are columns)
  e <- assay(x, i = assay)

  # convert expression matrix into a tidy tibble
  e <- as.data.frame(t(e), stringsAsFactors = FALSE)
  e <- tidyr::pivot_longer(e, tidyselect::everything())

  # generate the plot and return it
  p <-
    ggplot(e, aes(x = .data$value)) +
    # geom_density(adjust = 1 / 2, fill = "gray40") +
    geom_histogram(bins = 50, color = "gray50", fill = "gray50") +
    geom_freqpoly(bins = 50, color = "black") +
    facet_wrap(vars(.data$name), scales = "free") +
    # labs(title = plot_title, x = "", y = "Density") +
    labs(title = plot_title, x = "", y = "Count") +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  if (is.null(out_dir)) {
    return(p)
  } else {
    ggsave(filename = glue("{out_base}.png"), plot = p, width = 12, height = 9)
  }
}
