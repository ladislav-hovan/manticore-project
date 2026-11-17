# Tests for parameter consistency between lioness_single() and runMANTICORE()
# Ensures both functions share all common parameters with the same names

test_that("Both functions have consistent common parameters", {
  lioness_sig <- formals(lioness_single)
  runmanticore_sig <- formals(runMANTICORE)

  # List of parameters that should exist in both functions
  common_params <- c("gexMatrix", "tfMotifs", "ppiNet", "computingEngine",
                     "nCores", "assocMethod",
                     "alphaValue", "hammingValue", "nIter", "outNet",
                     "zScaling", "showProgress", "randomizationMethod",
                     "scaleByPresent")

  for (param in common_params) {
    expect_true(param %in% names(lioness_sig),
                info = paste("lioness_single() missing parameter:", param))
    expect_true(param %in% names(runmanticore_sig),
                info = paste("runMANTICORE() missing parameter:", param))
  }
})
