# Tests for parameter passing between lioness_single() and runMANTICORE()
# Validates that parameters are correctly passed from wrapper to core function

test_that("runMANTICORE() correctly passes computingEngine parameter to lioness_single()", {
  data(manticoreTest)

  result_cpu <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    computingEngine = "cpu",
    showProgress = FALSE
  )

  expect_true(is.data.frame(result_cpu))
  expect_true(all(c("tf", "target") %in% colnames(result_cpu)))
})

test_that("runMANTICORE() correctly passes nPC parameter to lioness_single()", {
  data(manticoreTest)

  result <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    nPC = 20,
    showProgress = FALSE
  )

  expect_true(is.data.frame(result))
})

test_that("runMANTICORE() correctly passes assocMethod parameter to lioness_single()", {
  data(manticoreTest)

  result <- runMANTICORE(
    gexMatrix = manticoreTest$gex,
    tfMotifs = manticoreTest$tf,
    ppiNet = manticoreTest$ppi,
    cellsMetadata = manticoreTest$metadata,
    groupBy = "region",
    assocMethod = "pearson",
    showProgress = FALSE
  )

  expect_true(is.data.frame(result))
})
