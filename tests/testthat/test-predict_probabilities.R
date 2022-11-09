test_that("Basic predict_probabilities call", {
  data_file <- system.file(
    "extdata", "sample_data_file.txt", package = "maxent.ot"
  )
  my_model <- optimize_weights(data_file)
  predictions <- predict_probabilities(data_file, my_model$weights)

  expect_equal(
    predictions$predictions$`Predicted Probability`,
    c(0.5, 0.5, 1, 0),
    tolerance=1e-5
  )
})
