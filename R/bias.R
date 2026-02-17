optimization_function <- function(bias, training_data, test_data) {

  model <- optimize_weights(training_data, mu = bias, sigma = bias)

  train_ll <- model$loglik
  weights  <- model$weights

  predictions <- predict_probabilities(test_data, weights)
  test_ll <- predictions$loglik

  out <- list(bias = bias, train_ll = train_ll, test_ll = test_ll, weights = weights)

  return(out)
}

