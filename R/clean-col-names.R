#' Clean up data frame column names
#'
#' @param x A data frame.
#'
#' @return A data frame.
#'
#' @import tidyselect stringr
#' @importFrom dplyr select
#' @importFrom janitor clean_names
#' @importFrom methods is
clean_col_names <- function(x) {
  # check if the input is valid
  if (!is.data.frame(x)) {
    stop("input is not a data frame")
  }

  # adjust MAV column names
  names(x) <- str_remove(names(x), ":Cyc_.*")
  names(x) <- str_remove(names(x), ".*:")

  # adjust Visiopharm column names
  names(x) <- str_remove(names(x), "MP.*\\(")

  # adjust HALO column names
  names(x) <- str_remove(names(x), " Cell Intensity")

  # clean column names
  x <- janitor::clean_names(x, case = "none")

  # adjust common column names
  if (!"x" %in% names(x)) {
    names(x)[names(x) == "X"] <- "x"
    names(x)[names(x) == "X_X"] <- "x"
    names(x)[names(x) == "CtrX"] <- "x"
    if (all(c("XMin", "XMax") %in% names(x))) {
      x$x <- (x$XMin + x$XMax) / 2
    }
  }
  if (!"y" %in% names(x)) {
    names(x)[names(x) == "Y"] <- "y"
    names(x)[names(x) == "Y_Y"] <- "y"
    names(x)[names(x) == "CtrY"] <- "y"
    if (all(c("YMin", "YMax") %in% names(x))) {
      x$y <- (x$YMin + x$YMax) / 2
    }
  }
  if (!"z" %in% names(x)) {
    names(x)[names(x) == "Z"] <- "z"
    names(x)[names(x) == "Z_Z"] <- "z"
  }
  if (!"reg" %in% names(x)) {
    names(x)[names(x) == "Region"] <- "reg"
    names(x)[names(x) == "region"] <- "reg"
  }
  if (!"tile_num" %in% names(x)) {
    names(x)[names(x) == "tile_nr"] <- "tile_num"
    names(x)[names(x) == "tile_number"] <- "tile_num"
  }
  if (!"size" %in% names(x)) {
    names(x)[names(x) == "Cell_Area_mm2"] <- "size"
  }

  return(x)
}
