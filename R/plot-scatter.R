#' Generate a generic scatter plot
#'
#' Create a scatter plot to display the relationship between two continuous variables.
#' This is essentially a wrapper for \code{\link[ggplot2:geom_point]{ggplot2::geom_point()}}.
#'
#' @param data A data frame.
#' @param x .
#' @param y .
#' @param color_by Column metadata field(s) or feature(s) to color by.
#' @param smooth A logical scalar. Smooth values. Helps to visualize expression patterns in a plot with many overlapping points.
#' @param range A vector of 2 values indicating the minimum and maximum percentiles for the color range. Helps to visualize expression patterns when extreme outliers are present. For example, \code{c(0, 0.99)} will not expand the color scale above 99th percentile.
#' @param title Plot title.
#' @param aspect_ratio Aspect ratio of the panel.
#'
#' @return A ggplot object.
#'
#' @import cowplot ggplot2
#' @importFrom ggsci pal_igv
#' @importFrom FNN knn.reg
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scattermore geom_scattermore
#' @importFrom stats quantile
plot_scatter <- function(data, x, y, color_by, smooth = FALSE, range = c(0.01, 0.99), title = "", aspect_ratio = 1) {
  # check if the input is valid
  if (!is.data.frame(data)) {
    stop("input is not a data frame")
  }
  if (!x %in% names(data)) {
    stop("input data frame does not have a `", x, "` column")
  }
  if (!y %in% names(data)) {
    stop("input data frame does not have a `", y, "` column")
  }
  if (!color_by %in% names(data)) {
    stop("input data frame does not have a `", color_by, "` column")
  }
  if (!is.numeric(range)) {
    stop("`range` is not numeric")
  }
  if (length(range) != 2) {
    stop("`range` is not a vector with 2 elements")
  }
  if (range[2] < range[1]) {
    stop("first `range` value is not smaller than the second")
  }
  if (range[2] > 1) {
    stop("`range` values are greater than 1")
  }
  if (!is.character(title)) {
    stop("`title` is not a character string")
  }
  if (!is.numeric(aspect_ratio)) {
    stop("`aspect_ratio` is not numeric")
  }

  # color scheme
  gradient_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  discrete_colors <- c(ggsci::pal_igv("default")(51), ggsci::pal_igv("default", alpha = 0.6)(51), ggsci::pal_igv("default", alpha = 0.4)(51))

  # determine if plotting continuous numeric/expression values
  continuous <- FALSE
  if (is.numeric(data[[color_by]])) {
    if (length(unique(data[[color_by]])) > 100) {
      continuous <- TRUE

      # rescale values to 0-1 (original distribution)
      # scale(x, center = min(x), scale = diff(range(x)))
      # set values to percentiles (flat/even distribution)
      # range(ecdf(x)(x))

      # adjust range of values (remove outliers)
      if (range[2] - range[1] < 1) {
        range <- quantile(data[[color_by]], range)
        data[[color_by]][data[[color_by]] < range[1]] <- range[1]
        data[[color_by]][data[[color_by]] > range[2]] <- range[2]
      }
    } else {
      data[[color_by]] <- as.factor(data[[color_by]])
    }
  }

  # smooth values
  if (continuous && smooth) {
    # smooth values using kernel density estimation (MASS::kde2d)
    # returns x/y coordinates and a matrix of the estimated density
    # dens <- MASS::kde2d(data[[x]], data[[y]], n = 100)
    # find the index of one vector (data$x) in another (dens$x)
    # intervals <- cbind(findInterval(data[[x]], dens[[x]]), findInterval(data[[y]], dens[[y]]))
    # data[[color_by]] <- dens$z[intervals]

    # k-nearest neighbor regression (ks::kde)
    # weights need to sum to sample size
    # w <- data[[color_by]] / sum(data[[color_by]]) * length(data[[color_by]])
    # dens <- ks::kde(data[, c(x, y)], w = w)
    # intervals <- cbind(findInterval(data[[x]], dens$eval.points[[1]]), findInterval(data[[y]], dens$eval.points[[2]]))
    # data[[color_by]] <- dens$estimate[intervals]

    # k-nearest neighbor regression
    knn <- FNN::knn.reg(as.matrix(data[, c(x, y)]), y = data[[color_by]], k = 10)
    data[[color_by]] <- knn$pred
  }

  # make dots larger for smaller datasets
  # due to rasterization properties it is often beneficial to try non-integer point sizes
  point_size <- 3.2
  if (nrow(data) < 100000) point_size <- 4.2
  if (nrow(data) < 10000) point_size <- 7.2
  if (nrow(data) < 1000) point_size <- 10.2

  # randomize cell order
  set.seed(99)
  data <- data[sample(1:nrow(data)), ]

  # generate the plot
  # rasterize with scattermore (simpler dependencies and more popular than ggrastr)
  p <-
    ggplot(data) +
    scattermore::geom_scattermore(
      aes(x = .data[[x]], y = .data[[y]], color = .data[[color_by]]),
      pointsize = point_size, pixels = c(1800, 1800)
    ) +
    labs(title = title) +
    scale_x_continuous(expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = 0.03)) +
    theme_cowplot(rel_large = 1) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  # adjust plot depending on the type of values
  if (continuous) {
    p <- p +
      scale_color_gradientn(name = "", colors = gradient_colors)
  } else {
    p <- p +
      theme(legend.title = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      scale_color_manual(values = discrete_colors)
  }

  return(p)
}
