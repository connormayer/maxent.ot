test_that("Basic predict_probabilities call", {
  data_file <- system.file(
    "extdata", "sample_data_frame.csv", package = "maxent.ot"
  )
  tableaux_df <- read.csv(data_file)
  my_model <- optimize_weights(tableaux_df)
  predictions <- predict_probabilities(tableaux_df, my_model$weights)

  expect_equal(
    predictions$predictions$Predicted,
    c(0.5, 0.5, 1, 0),
    tolerance=1e-5
  )
})
