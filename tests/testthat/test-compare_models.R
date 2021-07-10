test_that("Basic compare_models call", {
  data_file_1 <- system.file(
    "extdata", "sample_data_file.txt", package = "maxent.ot"
  )
  data_file_2 <- system.file(
    "extdata", "sample_data_file_2.txt", package = "maxent.ot"
  )
  model1 <- optimize_weights(data_file_1)
  model2 <- optimize_weights(data_file_2)

  result <- compare_models(model1, model2)

  expect_equal(result$k_delta, 1)
})

