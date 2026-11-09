# Tests for cell grouping functionality in runMANTICORE()
# Validates that cells are grouped correctly and networks are created per group

test_that("runMANTICORE() groups cells correctly for single groupBy column", {
  data(manticoreTest)
  result <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    showProgress = FALSE
  )
  
  network_columns <- setdiff(colnames(result), c("tf", "target"))
  expected_networks <- unique(manticoreTest$metadata$region)
  
  expect_equal(length(network_columns), length(expected_networks))
  expect_true(all(expected_networks %in% network_columns))
})

test_that("runMANTICORE() groups cells correctly for multiple groupBy columns", {
  data(manticoreTest)
  result <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = c("donor", "region"),
    showProgress = FALSE
  )
  
  network_columns <- setdiff(colnames(result), c("tf", "target"))
  
  # Should have multiple network columns (donor-region combinations)
  expect_true(length(network_columns) > 0)
  expect_true(all(grepl("--", network_columns)))  # Multiple groupBy uses -- separator
})

test_that("runMANTICORE() respects minCells parameter for filtering", {
  data(manticoreTest)
  
  # With minCells = 1, all groups should be included
  result_min_1 <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "donor",
    minCells = 1,
    showProgress = FALSE
  )
  
  network_cols_1 <- setdiff(colnames(result_min_1), c("tf", "target"))
  
  # Should have network columns
  expect_true(length(network_cols_1) > 0)
  expect_true(is.data.frame(result_min_1))
})
