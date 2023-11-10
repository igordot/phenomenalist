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

  # message("cols input: ", ncol(x))

  # exclude any columns with NA values
  x <- x[colSums(is.na(x)) == 0]
  if (!ncol(x)) {
    stop("all columns have NAs (expression values should not include NAs)")
  }
  # message("cols non-na: ", ncol(x))

  # keep just numeric double-precision columns
  x <- x[sapply(x, is.double)]
  if (!ncol(x)) {
    stop("there are no numeric non-integer (double-precision) columns")
  }
  # message("cols double: ", ncol(x))

  # remove columns with few unique values
  x <- x[sapply(x, function(y) length(unique(y)) > max(100, length(y) / 100))]
  # message("cols unique: ", ncol(x))

  # get expression intensity columns based on potential labels
  exprs_cols <- str_subset(names(x), "Cyc|cyc|MP\\s|Cell.Intensity")

  # if matching based on name worked, use those columns
  if (length(exprs_cols) >= 10) {
    return(exprs_cols)
  }

  # use column values if matching based on name did not work

  # remove columns with small values and percentile columns where max is 100
  x <- x[sapply(x, function(y) {
    (max(y) > 1.5) && (max(y) < 99.9 || max(y) > 100.1)
  })]
  # message("cols non-percentile: ", ncol(x))

  # get columns that are more likely to be expression values (start with CD)
  sub_cols <- str_subset(names(x), "^C[Dd]\\d")
  if (length(sub_cols) < 10) {
    sub_cols <- names(x)
  }
  # message(" sub cols: ", length(sub_cols))

  # get the range of possible values based on likely expression columns
  sub_mat <- as.matrix(x[sub_cols])
  sub_mins <- apply(sub_mat, MARGIN = 2, FUN = min)
  sub_medians <- apply(sub_mat, MARGIN = 2, FUN = quantile, 0.5)
  sub_maxs <- apply(sub_mat, MARGIN = 2, FUN = quantile, 0.999)

  # check if expression values should be positive
  if (quantile(sub_mins, 0.2) >= 0) {
    x <- x[sapply(x, function(y) min(y) >= 0)]
    # message("cols positive: ", ncol(x))
  }

  # min values should be lower than most of the max values
  # x <- x[sapply(x, function(y) min(y) < quantile(sub_maxs, 0.8))]
  # message("cols min<max: ", ncol(x))

  # q1 values should be lower than most of the max values
  x <- x[sapply(x, function(y) quantile(y, 0.25) < quantile(sub_maxs, 0.8))]
  # message(" maxs: ", quantile(sub_maxs, 0.8))
  # message("cols med<max: ", ncol(x))

  # max values should be higher than some of the medians
  x <- x[sapply(x, function(y) max(y) > quantile(sub_medians, 0.2))]
  # message(" meds: ", quantile(sub_medians, 0.2))
  # message("cols max>med: ", ncol(x))

  # assuming that the remaining columns are expression columns
  exprs_cols <- names(x)

  # check if enough columns remain after filtering
  if (length(exprs_cols) < 10) {
    message("input columns: ", toString(colnames(x)))
    message("expression columns detected: ", toString(exprs_cols))
    stop("too few expression columns detected")
  }

  return(exprs_cols)
}
