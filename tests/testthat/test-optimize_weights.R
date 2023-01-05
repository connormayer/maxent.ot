test_that("Basic optimize_weights call", {
  data_file <- system.file(
    "extdata", "sample_data_frame.csv", package = "maxent.ot"
  )
  tableaux_df <- read.csv(data_file)
  my_model <- optimize_weights(tableaux_df)
  expect_equal(length(my_model$weights), 2)

  expect_equal(my_model$k, 2)

  expect_equal(my_model$n, 3)

  expect_equal(my_model$bias_params, NA)
})

