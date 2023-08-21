test_that("import tonsil data default", {
  tonsil_csv <- system.file("extdata", "tonsil-akoya-2018-500.csv", package = "phenomenalist")
  tonsil_spe <- suppressMessages(create_object(tonsil_csv))
  expect_s4_class(tonsil_spe, "SpatialExperiment")
  expect_equal(ncol(tonsil_spe), 500)
})

test_that("import tonsil data with skip_cols and transformation", {
  tonsil_csv <- system.file("extdata", "tonsil-akoya-2018-500.csv", package = "phenomenalist")
  tonsil_spe <- suppressMessages(create_object(tonsil_csv, skip_cols = "DAPI|Blank", transformation = "z"))
  expect_s4_class(tonsil_spe, "SpatialExperiment")
  expect_equal(nrow(tonsil_spe), 18)
  expect_equal(ncol(tonsil_spe), 500)
  expect_equal(length(assayNames(tonsil_spe)), 3)
})

test_that("import spleen data default", {
  spleen_csv <- system.file("extdata", "spleen-goltsev-2018-500.csv", package = "phenomenalist")
  spleen_spe <- suppressMessages(create_object(spleen_csv))
  expect_s4_class(spleen_spe, "SpatialExperiment")
  expect_equal(nrow(spleen_spe), 32)
  expect_equal(ncol(spleen_spe), 500)
  expect_equal(length(assayNames(spleen_spe)), 1)
})

test_that("import spleen data with skip_cols", {
  spleen_csv <- system.file("extdata", "spleen-goltsev-2018-500.csv", package = "phenomenalist")
  spleen_spe <- suppressMessages(create_object(spleen_csv, skip_cols = "blank"))
  expect_s4_class(spleen_spe, "SpatialExperiment")
  expect_equal(nrow(spleen_spe), 30)
  expect_equal(ncol(spleen_spe), 500)
  expect_equal(length(assayNames(spleen_spe)), 1)
})
