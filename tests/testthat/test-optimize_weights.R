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

test_that("Displays harmony and maxent values", {
  data_file <- system.file(
    "extdata", "sample_data_frame.csv", package = "maxent.ot"
  )
  df <- read.csv(data_file)
  model <- optimize_weights(df)
  prediction <- predict_probabilities(
    df, model$weights,
    include_harmony = TRUE,
    include_maxent_values = TRUE
  )
  expect_true("harmony" %in% names(prediction))
  expect_true("maxent_values" %in% names(prediction))
})


