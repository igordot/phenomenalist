#' Detect marker/antibody column names
#'
#' The input is generally a table of both expression values and various metadata columns.
#' This function attempts to determine which of the columns contain expression values.
#'
#' @param x A data frame.
#'
#' @return A vector of column names.
#'
#' @import stringr tidyselect
#' @importFrom methods is
#' @importFrom stats quantile
detect_exprs_cols <- function(x) {
  # check if the input is valid
  if (!is.data.frame(x)) {
    stop("input is not a data frame")
  }

  # exclude any columns with NA values
  x <- x[colSums(is.na(x)) == 0]
  if (!ncol(x)) {
    stop("all columns have NAs (expression values should not include NAs)")
  }

  # keep just numeric double-precision columns
  x <- x[sapply(x, is.double)]
  if (!ncol(x)) {
    stop("there are no numeric non-integer (double-precision) columns")
  }

  # remove columns with few unique values
  x <- x[sapply(x, function(y) length(unique(y)) > 100)]

  # get expression intensity columns based on potential labels
  exprs_cols <- str_subset(names(x), "Cyc|cyc|MP\\s|Cell.Intensity")

  # if matching based on name worked, use those columns
  if (length(exprs_cols) >= 10) {
    return(exprs_cols)
  }

  # use column values if matching based on name did not work

  # remove columns with small values and percentile columns where max is 100
  x <- x[sapply(x, function(y) {
    (max(y) > 1.5) && (max(y) < 99.5 || max(y) > 100.5)
  })]

  # get columns that are more likely to be expression values (CD/DAPI)
  sub_cols <- str_subset(names(x), "^C[Dd]\\d|DAPI")
  if (length(sub_cols) < 10) {
    sub_cols <- names(x)
  }

  # get the range of possible values based on likely expression columns
  sub_mat <- as.matrix(x[sub_cols])
  sub_mins <- apply(sub_mat, MARGIN = 2, FUN = min)
  sub_medians <- apply(sub_mat, MARGIN = 2, FUN = quantile, 0.5)

  # expression values are positive
  if (quantile(sub_mins, 0.2) >= 0) {
    x <- x[sapply(x, function(y) min(y) >= 0)]
    # expression values start at 0
    if (quantile(sub_mins, 0.9) == 0) {
      x <- x[sapply(x, function(y) min(y) == 0)]
    }
  }

  # max values should be higher than some of the medians
  x <- x[sapply(x, function(y) max(y) > quantile(sub_medians, 0.2))]

  # assuming the the remaining columns are expression columns
  exprs_cols <- names(x)

  # check if enough columns remain after filtering
  if (length(exprs_cols) < 10) {
    message("input columns: ", toString(colnames(x)), "\n")
    if (length(exprs_cols) == 0) {
      stop("no markers detected")
    }
    stop("too few markers detected: ", toString(exprs_cols))
  }

  return(exprs_cols)
}
