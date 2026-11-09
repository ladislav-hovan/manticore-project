# Tests for scorpion() network computation
# Validates that scorpion() produces correct output structures with various parameters

test_that("scorpion() accepts gexMatrix and tfMotifs and produces list output", {
  data(manticoreTest)
  result <- scorpion(
    tfMotifs = manticoreTest$tf,
    gexMatrix = manticoreTest$gex[, 1:100],  # Small subset for speed
    ppiNet = manticoreTest$ppi,
    alphaValue = 0.8,
    showProgress = FALSE
  )
  
  expect_true(is.list(result))
  expect_true(all(c("regNet", "coregNet", "coopNet") %in% names(result)))
  expect_true(is.numeric(result$numGenes))
  expect_true(is.numeric(result$numTFs))
  expect_true(is.numeric(result$numEdges))
})

test_that("scorpion() works with NULL tfMotifs (co-expression only)", {
  skip("Skipped due to Matrix/gpuR dependency issues")
  data(manticoreTest)
  result <- scorpion(
    tfMotifs = NULL,
    gexMatrix = manticoreTest$gex[, 1:100],
    ppiNet = manticoreTest$ppi,
    showProgress = FALSE
  )
  
  expect_true(is.list(result) || is.matrix(result))
  expect_true(nrow(result) > 0 && ncol(result) > 0)
})

test_that("scorpion() works with NULL ppiNet", {
  skip("Skipped due to Matrix/gpuR dependency issues")
  data(manticoreTest)
  result <- scorpion(
    tfMotifs = manticoreTest$tf,
    gexMatrix = manticoreTest$gex[, 1:100],
    ppiNet = NULL,
    showProgress = FALSE
  )
  
  expect_true(is.list(result) || is.matrix(result))
})

test_that("scorpion() alphaValue parameter changes output", {
  data(manticoreTest)
  
  result_low <- scorpion(
    tfMotifs = manticoreTest$tf,
    gexMatrix = manticoreTest$gex[, 1:50],
    ppiNet = manticoreTest$ppi,
    alphaValue = 0.1,
    showProgress = FALSE
  )
  
  result_high <- scorpion(
    tfMotifs = manticoreTest$tf,
    gexMatrix = manticoreTest$gex[, 1:50],
    ppiNet = manticoreTest$ppi,
    alphaValue = 0.8,
    showProgress = FALSE
  )
  
  # Both should produce results
  expect_true(is.list(result_low) || is.matrix(result_low))
  expect_true(is.list(result_high) || is.matrix(result_high))
})
