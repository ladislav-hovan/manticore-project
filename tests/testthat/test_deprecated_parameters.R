# Tests for deprecated parameter names
# Ensures old parameter names have been removed from both functions

test_that("scorpion() does NOT have countMatrix parameter", {
  sig <- formals(scorpion)
  expect_false("countMatrix" %in% names(sig))
})

test_that("runMANTICORE() does NOT have countMatrix parameter", {
  sig <- formals(runMANTICORE)
  expect_false("countMatrix" %in% names(sig))
})

test_that("runMANTICORE() does NOT have tfTargets parameter", {
  sig <- formals(runMANTICORE)
  expect_false("tfTargets" %in% names(sig))
})

test_that("Old parameter names do not exist in scorpion()", {
  sig <- formals(scorpion)
  old_names <- c("countMatrix", "gexData")
  
  for (name in old_names) {
    expect_false(name %in% names(sig),
                 info = paste("Old parameter name found:", name))
  }
})

test_that("Old parameter names do not exist in runMANTICORE()", {
  sig <- formals(runMANTICORE)
  old_names <- c("countMatrix", "tfTargets", "tfData")
  
  for (name in old_names) {
    expect_false(name %in% names(sig),
                 info = paste("Old parameter name found:", name))
  }
})
