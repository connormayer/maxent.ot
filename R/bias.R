#' @export
optimize_bias <- function(data, k, mu = 0,
                          sigma = 0.5, method = "L-BFGS-B") { # check default method, default sigma

  result <- stats::optim(
    par = c(mu, sigma),
    fn = optimization_function,
    data = data,
    k = k,
    method = method
  )

  return(list(
    mu = result$par[1],
    sigma = result$par[2],
    mean_ll = -result$value
  ))
}

optimization_function <- function(par, data, k) {

  partitions <- maxent.ot:::partition_data(data, k)

  log_liks <- c()

  for (curr in seq_len(k)) {

    training_data <- data
    training_data$Frequency <- 0
    test_data <- training_data

    training_part <- partitions[partitions$partition != curr, ]
    test_part     <- partitions[partitions$partition == curr, ]

    training_tableau <- populate_tableau(training_data, training_part)
    test_tableau     <- populate_tableau(test_data, test_part)

    fold_ll <- optimize_fold(
      par,
      training_tableau,
      test_tableau
    )

    log_liks <- c(log_liks, fold_ll)
  }

  return(-mean(log_liks))
}


#' Helper function that takes a bias vector, potential mu/sigma vector
#' (if mapping to each value), train/test splits,
#' outputs the test loglik to be used in optim later
optimize_fold <- function(par, train_data, test_data) { # par contains vectorized mu, sigma values

  mu    <- par[1]
  sigma <- par[2]

  is_valid(mu, sigma, train_data, test_data)

  model <- optimize_weights(train_data, mu = mu, sigma = sigma)

  predictions <- predict_probabilities(test_data, model$weights)

  return(predictions$loglik)
}


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




