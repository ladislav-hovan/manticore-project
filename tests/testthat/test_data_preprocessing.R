# Tests for data preprocessing and parameter effects in runMANTICORE()
# Validates normalization, batch correction, and parameter passing

test_that("runMANTICORE() parameter passing to scorpion() works with alphaValue", {
  data(scorpionTest)
  result_alpha_low <- runMANTICORE(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    alphaValue = 0.1,
    showProgress = FALSE
  )
  
  result_alpha_high <- runMANTICORE(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
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
  data(scorpionTest)
  result_normalized <- runMANTICORE(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    normalizeData = TRUE,
    showProgress = FALSE
  )
  
  result_unnormalized <- runMANTICORE(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
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
  data(scorpionTest)
  
  result <- runMANTICORE(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    removeBatchEffect = TRUE,
    batch = scorpionTest$metadata$donor,
    showProgress = FALSE
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("tf", "target") %in% colnames(result)))
})

test_that("runMANTICORE() skips batch correction when batch has a single level", {
  data(scorpionTest)

  args <- list(
    gexMatrix = scorpionTest$gex,
    tfMotifs = scorpionTest$tf,
    ppiNet = scorpionTest$ppi,
    cellsMetadata = scorpionTest$metadata,
    groupBy = "region",
    showProgress = FALSE
  )

  set.seed(1)
  uncorrected <- do.call(runMANTICORE, args)
  set.seed(1)
  expect_message(
    single_batch <- do.call(runMANTICORE, c(args, list(
      removeBatchEffect = TRUE,
      batch = rep("b1", ncol(scorpionTest$gex))
    ))),
    "fewer than two levels"
  )

  expect_equal(single_batch, uncorrected)
})
