#' Detect marker/antibody column names
#'
#' @param x A data frame.
#'
#' @return A vector of column names.
#'
#' @import stringr tidyselect
#' @importFrom dplyr between
#' @importFrom methods is
detect_exprs_cols <- function(x) {
  # check if the input is valid
  if (!is.data.frame(x)) {
    stop("input is not a data frame")
  }

  # exclude any columns with NA values
  x <- x[colSums(is.na(x)) == 0]

  # keep just numeric columns
  # x <- select_if(x, is.numeric)
  # x <- select(x, where(is.numeric))
  x <- x[sapply(x, is.numeric)]

  # get expression intensity columns based on common labels
  exprs_cols <- str_subset(names(x), "Cyc|cyc|MP\\s|Cell.Intensity")

  # use column contents if matching based on name did not work
  if (length(exprs_cols) < 10) {
    exprs_cols <- names(x[sapply(x, function(y) {
      # if labels or coordinates, not likely to have small non-integer or negative values
      (any(between(y, 0.05, 0.45)) || min(y) < -100) &&
        # if coordinates, values are not likely to start close to 0
        min(y) < 1 &&
        # if percentiles, max is 100
        (max(y) < 99.5 || max(y) > 100.5)
    })])
  }

  # check if markers were detected
  if (length(exprs_cols) < 10) {
    message("input columns: ", toString(colnames(x)), "\n")
    if (length(exprs_cols) == 0) {
      stop("no markers detected")
    }
    stop("too few markers detected: ", toString(exprs_cols))
  }

  return(exprs_cols)
}
