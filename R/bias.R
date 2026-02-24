# validates mu and sigma according to passed data
is_valid <- function(mu, sigma, train_data, test_data) {

  # potentially unnecessary if assumed to only be used with valid splits
  if (length(train_data) != length(test_data)){
    stop("train and test data must be the same length")
  }
  if (!is.numeric(mu) || !is.numeric(sigma)) {
    stop("mu and sigma must be numeric.")
  }

  if (length(mu) != length(sigma)) {
    stop("mu and sigma must be the same length.")
  }

  constraint_cols <- sapply(train_data, is.numeric)
  n_constraints <- sum(constraint_cols) - 1

  if (length(mu) > 1 && length(mu) != n_constraints) {
    stop("If vectors, mu and sigma must match number of constraints.")
  }

  return(TRUE)
}


#' Helper function that takes a bias vector, potential mu/sigma vector
#' (if mapping to each value), train/test splits,
#' outputs the test loglik to be used in optim later
optimization_function <- function(mu, sigma, train_data, test_data) {

  is_valid(mu, sigma, train_data)

  model <- optimize_weights(train_data, mu = mu, sigma = sigma)

  predictions <- predict_probabilities(test_data, model$weights)

  return(predictions$loglik)
}


