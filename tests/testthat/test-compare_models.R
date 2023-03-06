test_that("Basic compare_models call", {
  data_file_1 <- system.file(
    "extdata", "sample_data_frame.csv", package = "maxent.ot"
  )
  data_file_2 <- system.file(
    "extdata", "sample_data_frame_large.csv", package = "maxent.ot"
  )
  df1 <- read.csv(data_file_1)
  df2 <- read.csv(data_file_2)
  model1 <- optimize_weights(df1)
  model2 <- optimize_weights(df2)

  result <- compare_models(model1, model2, method='bic')

  expect_equal(nrow(result), 2)
})

