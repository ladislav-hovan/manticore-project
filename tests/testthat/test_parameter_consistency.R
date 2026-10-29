# Tests for parameter consistency between scorpion() and runMANTICORE()
# Ensures both functions share all common parameters with the same names

test_that("Both functions have consistent common parameters", {
  scorpion_sig <- formals(scorpion)
  runscoprion_sig <- formals(runMANTICORE)
  
  # List of parameters that should exist in both functions
  common_params <- c("gexMatrix", "tfMotifs", "ppiNet", "computingEngine",
                     "nCores", "gammaValue", "nPC", "assocMethod",
                     "alphaValue", "hammingValue", "nIter", "outNet",
                     "zScaling", "showProgress", "randomizationMethod",
                     "scaleByPresent", "filterExpr")
  
  for (param in common_params) {
    expect_true(param %in% names(scorpion_sig),
                info = paste("scorpion() missing parameter:", param))
    expect_true(param %in% names(runscoprion_sig),
                info = paste("runMANTICORE() missing parameter:", param))
  }
})
