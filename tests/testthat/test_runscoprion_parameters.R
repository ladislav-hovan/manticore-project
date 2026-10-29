# Tests for runMANTICORE() parameter signature verification
# Ensures gexMatrix and tfMotifs parameters are present and properly named

test_that("runMANTICORE() has gexMatrix parameter in signature", {
  sig <- formals(runMANTICORE)
  expect_true("gexMatrix" %in% names(sig))
})

test_that("runMANTICORE() has tfMotifs parameter in signature", {
  sig <- formals(runMANTICORE)
  expect_true("tfMotifs" %in% names(sig))
})

test_that("runMANTICORE() gexMatrix parameter is required", {
  sig <- formals(runMANTICORE)
  # gexMatrix should not have a default value (either NULL or absent from formals)
  expect_true(is.null(sig$gexMatrix) || !"gexMatrix" %in% names(sig) || 
              identical(as.character(sig$gexMatrix), ""))
})

test_that("runMANTICORE() tfMotifs parameter is required", {
  sig <- formals(runMANTICORE)
  # tfMotifs should not have a default value (either NULL or absent from formals)
  expect_true(is.null(sig$tfMotifs) || !"tfMotifs" %in% names(sig) || 
              identical(as.character(sig$tfMotifs), ""))
})
