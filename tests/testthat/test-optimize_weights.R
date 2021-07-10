test_that("Basic optimize_weights call", {
  data_file <- system.file(
    "extdata", "sample_data_file.txt", package = "maxent.ot"
  )
  my_model <- optimize_weights(data_file)
  expect_equal(length(my_model$weights), 2)
  expect_equal(my_model$weights[['Constraint1']], 14.20336, tolerance=1e-6)
  expect_equal(my_model$weights[['Constraint2']], 14.20341, tolerance=1e-6)

  expect_equal(my_model$loglik, -1.386295, tolerance=1e-6)

  expect_equal(my_model$k, 2)

  expect_equal(my_model$n, 3)

  expect_equal(my_model$bias_params, NA)
})

