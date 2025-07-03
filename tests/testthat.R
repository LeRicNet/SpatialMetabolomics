# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

# Unit tests for SpatialMetabolic package

library(testthat)
library(SpatialMetabolics)
library(SpatialExperiment)
library(SingleCellExperiment)

# Note: create_test_data is now defined in helper-test-data.R

# Test class creation and validity
test_that("SpatialMetabolic class creation works", {
  spm <- create_test_data()

  expect_s4_class(spm, "SpatialMetabolic")
  expect_true(inherits(spm, "SpatialExperiment"))
  expect_equal(ncol(spm), 100)
  expect_equal(nrow(spm), 100)
})
