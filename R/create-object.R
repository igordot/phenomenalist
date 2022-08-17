#' Create a SpatialExperiment object
#'
#' Create a \linkS4class{SpatialExperiment} object (used in the same manner as \linkS4class{SingleCellExperiment}) that will store all of the project data.
#'
#' @param x A data frame or a path to a file (can be gzipped).
#' @param expression_cols Column names of markers/antibody expression as a vector or a grep pattern. Auto-detected if not specified.
#' @param metadata_cols Column names of cell metadata (not markers/antibodies) as a vector or a grep pattern. Auto-detected if not specified.
#' @param skip_cols Column names to ignore (such as bad antibodies) as a vector or a grep pattern.
#' @param clean_names A logical scalar. Clean the data frame column names to remove problematic characters and make them unique.
#' @param transformation A character string indicating which transformation method should be used. See \code{\link[=transform]{transform()}}.
#' @param out_dir Name of the output analysis directory. If specified, the object and the corresponding plots will be saved there.
#'
#' @return A \linkS4class{SpatialExperiment} object.
#'
#' @import SpatialExperiment stringr
#' @importFrom data.table fread
#' @importFrom dplyr select
#' @importFrom janitor clean_names
#' @importFrom methods is
#' @importFrom tibble column_to_rownames
#'
#' @export
#'
#' @examples
#' tonsil_csv <- system.file("extdata", "tonsil-akoya-2018-500.csv", package = "phenomenalist")
#' tonsil_spe <- create_object(tonsil_csv, skip_cols = "DAPI|Blank", transformation = "z")
create_object <- function(x, expression_cols = NULL, metadata_cols = NULL, skip_cols = NULL, clean_names = TRUE, transformation = NULL, out_dir = NULL) {

  # check if the input is valid
  if (is.character(x)) {
    x <- data.table::fread(x, stringsAsFactors = FALSE, data.table = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("input is not a file or a data frame")
  }
  if (!is.null(out_dir)) {
    if (dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` already exists")
    }
  }

  # check if the data frame dimensions make sense
  if (nrow(x) < 500) {
    stop("data frame has too few rows/cells")
  }
  if (ncol(x) < 10) {
    stop("data frame has too few columns")
  }

  message("number of input table rows: ", nrow(x))
  message("number of input table columns: ", ncol(x))
  message("")

  # remove columns that should be ignored from the table
  if (!is.null(skip_cols)) {
    # treating as a pattern to get matching column names if length of 1
    if (length(skip_cols) == 1) {
      skip_cols <- str_subset(names(x), pattern = skip_cols)
    }
    x <- x[, setdiff(names(x), skip_cols)]
    message("skipped columns: ", toString(skip_cols), "\n")
  }

  # create the expression data frame
  if (is.null(expression_cols)) {
    expression_cols <- detect_exprs_cols(x)
  }
  # treating as a pattern to get matching column names if length of 1
  if (length(expression_cols) == 1) {
    expression_cols <- str_subset(names(x), pattern = expression_cols)
  }
  # check that specified column names are valid
  if (length(setdiff(expression_cols, names(x))) > 0) {
    stop("missing markers: ", toString(setdiff(expression_cols, names(x))))
  }
  exprs <- x[, expression_cols]

  # use the non-expression columns for the metadata data frame if not specified
  if (is.null(metadata_cols)) {
    metadata_cols <- setdiff(names(x), expression_cols)
  }
  # treating as a pattern to get matching column names if length of 1
  if (length(metadata_cols) == 1) {
    metadata_cols <- str_subset(names(x), pattern = metadata_cols)
  }
  # check that specified column names are valid
  if (length(setdiff(metadata_cols, names(x))) > 0) {
    stop("missing metadata columns: ", toString(setdiff(metadata_cols, names(x))))
  }
  x <- x[, metadata_cols]

  # clean column names
  if (clean_names) {

    # run generic column name cleanup
    exprs <- clean_col_names(exprs)
    x <- clean_col_names(x)

    # identify cell IDs
    if (!"cell_id" %in% names(x)) {
      names(x)[names(x) == "CellID"] <- "cell_id"
      names(x)[names(x) == "label"] <- "cell_id"
      names(x)[names(x) == "Object_Id"] <- "cell_id"
    }
    if (!"cell_id" %in% names(x)) {
      x$cell_id <- rownames(x)
    }

    # force cell IDs to be strings to avoid any downstream issues
    if (is.numeric(x$cell_id)) {
      x$cell_id <- str_pad(as.character(x$cell_id), width = 7, pad = "0")
      x$cell_id <- str_c("C", x$cell_id)
    }

    # make sure cell IDs are unique
    x$cell_id <- make.names(x$cell_id, unique = TRUE)
  }

  # confirm that all the necessary metadata is present
  if (!"cell_id" %in% names(x)) {
    stop("data frame must contain `cell_id` column")
  }
  if (!"x" %in% names(x)) {
    stop("data frame must contain `x` column")
  }
  if (!"y" %in% names(x)) {
    stop("data frame must contain `y` column")
  }

  # set cell IDs as rownames
  rownames(x) <- x$cell_id

  # create the expression matrix
  exprs <- as.matrix(exprs)
  exprs <- exprs[, sort(colnames(exprs))]

  # set cell IDs as rownames
  rownames(exprs) <- rownames(x)

  # remove rows/cells without values (assuming negative values are not possible)
  # exprs <- exprs[rowSums(exprs) > 0, ]
  # x <- x[rownames(exprs), ]

  message("number of expression columns: ", ncol(exprs))
  message("number of metadata columns: ", ncol(x))
  message("")

  message("expression columns: ", toString(colnames(exprs)), "\n")
  message("metadata columns: ", toString(colnames(x)), "\n")

  # create a SpatialExperiment object
  # in v1.5.2 (1/2022) spatialData slot was deprecated (move the contents to colData)
  s <-
    SpatialExperiment::SpatialExperiment(
      assay = list(counts = t(exprs)),
      colData = x,
      spatialCoordsNames = c("x", "y")
    )

  # create output directory if specified
  if (!is.null(out_dir)) {
    dir.create(out_dir)
  }

  # apply transformation if specified
  if (!is.null(transformation)) {
    s <- transform(s, method = transformation, out_dir = out_dir)
    s <- run_umap(s, n_threads = 4)
  }

  # save object if output directory is specified
  if (!is.null(out_dir)) {
    message("saving object")
    saveRDS(s, paste0(out_dir, "/spe.rds"))
  }

  return(s)
}
