# Tests for data preprocessing and parameter effects in runMANTICORE()
# Validates normalization, batch correction, and parameter passing

test_that("runMANTICORE() parameter passing to scorpion() works with alphaValue", {
  data(manticoreTest)
  result_alpha_low <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    alphaValue = 0.1,
    showProgress = FALSE
  )
  
  result_alpha_high <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    alphaValue = 0.8,
    showProgress = FALSE
  )
  
  # Results should differ with different alphaValue
  expect_false(identical(result_alpha_low, result_alpha_high))
  
  # But should have same structure
  expect_equal(colnames(result_alpha_low), colnames(result_alpha_high))
  expect_equal(nrow(result_alpha_low), nrow(result_alpha_high))
})

test_that("runMANTICORE() data normalization affects results", {
  data(manticoreTest)
  result_normalized <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    normalizeData = TRUE,
    showProgress = FALSE
  )
  
  result_unnormalized <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    normalizeData = FALSE,
    showProgress = FALSE
  )
  
  # Results should differ with different normalization
  expect_false(identical(result_normalized, result_unnormalized))
  
  # But should have same structure
  expect_equal(colnames(result_normalized), colnames(result_unnormalized))
})

test_that("runMANTICORE() batch effect correction works with valid batch", {
  data(manticoreTest)
  
  result <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    removeBatchEffect = TRUE,
    batch = manticoreTest$metadata$donor,
    showProgress = FALSE
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("tf", "target") %in% colnames(result)))
})

test_that("runMANTICORE() skips batch correction when batch has a single level", {
  data(manticoreTest)

  args <- list(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    showProgress = FALSE
  )

  set.seed(1)
  uncorrected <- do.call(runMANTICORE, args)
  set.seed(1)
  expect_message(
    single_batch <- do.call(runMANTICORE, c(args, list(
      removeBatchEffect = TRUE,
      batch = rep("b1", ncol(manticoreTest$gex))
    ))),
    "fewer than two levels"
  )

  expect_equal(single_batch, uncorrected)
})
