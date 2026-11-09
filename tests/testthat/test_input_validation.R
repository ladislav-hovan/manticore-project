# Tests for input validation and error handling
# Validates that functions properly validate inputs and provide meaningful errors

test_that("scorpion() alphaValue parameter validation - negative value", {
  data(manticoreTest)
  
  # alphaValue < 0 should error
  expect_error(
    scorpion(
      tfMotifs = manticoreTest$tf,
      gexMatrix = manticoreTest$gex[, 1:50],
      ppiNet = manticoreTest$ppi,
      alphaValue = -0.1,
      showProgress = FALSE
    )
  )
})

test_that("scorpion() alphaValue parameter validation - value > 1", {
  data(manticoreTest)
  
  # alphaValue > 1 should error
  expect_error(
    scorpion(
      tfMotifs = manticoreTest$tf,
      gexMatrix = manticoreTest$gex[, 1:50],
      ppiNet = manticoreTest$ppi,
      alphaValue = 1.5,
      showProgress = FALSE
    )
  )
})

test_that("runMANTICORE() input validation checks gexMatrix dimensions", {
  data(manticoreTest)
  
  # Create mismatched metadata
  bad_metadata <- manticoreTest$metadata[1:100, ]
  
  # Should error because gexMatrix has more columns than metadata has rows
  expect_error(
    runMANTICORE(
      gexMatrix = manticoreTest$gex,
      tfMotifs = manticoreTest$tf,
      ppiNet = manticoreTest$ppi,
      cellsMetadata = bad_metadata,
      groupBy = "region",
      showProgress = FALSE
    )
  )
})

test_that("runMANTICORE() validates groupBy column existence", {
  data(manticoreTest)
  
  # Should error with invalid groupBy column
  expect_error(
    runMANTICORE(
      gexMatrix = manticoreTest$gex,
      tfMotifs = manticoreTest$tf,
      ppiNet = manticoreTest$ppi,
      cellsMetadata = manticoreTest$metadata,
      groupBy = "nonexistent_column",
      showProgress = FALSE
    )
  )
})

test_that("runMANTICORE() batch effect correction requires batch parameter", {
  data(manticoreTest)
  
  # Should error when removeBatchEffect=TRUE but no batch parameter
  expect_error(
    runMANTICORE(
      gexMatrix = manticoreTest$gex,
      tfMotifs = manticoreTest$tf,
      ppiNet = manticoreTest$ppi,
      cellsMetadata = manticoreTest$metadata,
      groupBy = "region",
      removeBatchEffect = TRUE,
      batch = NULL,
      showProgress = FALSE
    )
  )
})
