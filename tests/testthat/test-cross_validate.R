test_that("Basic cross-validation calls", {
  data_file <- system.file(
    "extdata", "amp_demo_grammar_full.csv", package = "maxent.ot"
  )
  tableaux_df <- read.csv(data_file)

  # Test call with scalar values
  cv_results <- cross_validate(tableaux_df, 10, 0, 1)
  expect_equal(nrow(cv_results), 1)

  # Test call with vectors of scalars
  cv_results <- cross_validate(tableaux_df, 10, c(0, 1, 2), c(0.01, 0.1, 1))
  expect_equal(nrow(cv_results), 3)

  # Test call with grid search
  cv_results <- cross_validate(
    tableaux_df, 10, c(0, 1, 2), c(0.01, 0.1, 1), grid_search=TRUE
  )
  expect_equal(nrow(cv_results), 9)

  # Test call with list of vectors
  cv_results <- cross_validate(
    tableaux_df, 10, list(c(0, 1, 2, 3), c(3, 4, 5, 6)),
    list(c(0.01, 0.1, 1, 2), c(2, 1, 0.1, 0.01))
  )
  expect_equal(nrow(cv_results), 2)
})
