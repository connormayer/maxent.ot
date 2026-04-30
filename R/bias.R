#' @export
optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                          k = NULL, sigma = 0.5, mu = 0,
                          method = "L-BFGS-B") {

  if (!is.null(train_data) && !is.null(test_data)) {
    fn <- optimization_function_fixed
    fn_args <- list(train_data = train_data, test_data = test_data, mu = mu)
  } else {
    fn <- optimization_function
    fn_args <- list(data = data, k = k, mu = mu)
  }

  result <- do.call(
    stats::optim,
    c(
      list(
        par = c(log_sigma = log(sigma)),
        fn = fn,
        method = method,
        lower = c(-Inf)
      ),
      fn_args
    )
  )

  sigma_hat <- exp(result$par[1])

  return(list(
    mu = mu,
    sigma = sigma_hat,
    mean_ll = -result$value,
    convergence = result$convergence
  ))
}


optimization_function_fixed <- function(par, train_data, test_data, mu) {

  val <- tryCatch({

    sigma <- exp(par[1])

    training_data <- train_data
    training_data$Frequency <- 0
    testing_data <- test_data

    train_tableau <- populate_tableau(training_data, train_data)
    test_tableau  <- populate_tableau(testing_data, test_data)

    ll <- optimize_fold(mu, sigma, train_tableau, test_tableau)

    if (!is.finite(ll)) return(1e12)

    -ll

  }, error = function(e) {
    1e12
  })

  return(val)
}

optimization_function <- function(par, data, k, mu) {

  val <- tryCatch({

    sigma <- exp(par[1])

    partitions <- maxent.ot:::partition_data(data, k)
    log_liks <- numeric(k)

    for (curr in seq_len(k)) {

      training_data <- data
      training_data$Frequency <- 0
      test_data <- training_data

      training_part <- partitions[partitions$partition != curr, ]
      test_part     <- partitions[partitions$partition == curr, ]

      training_tableau <- populate_tableau(training_data, training_part)
      test_tableau     <- populate_tableau(test_data, test_part)

      fold_ll <- optimize_fold(
        mu,
        sigma,
        training_tableau,
        test_tableau
      )

      if (!is.finite(fold_ll)) return(1e12)

      log_liks[curr] <- fold_ll
    }

    -mean(log_liks)

  }, error = function(e) {
    1e12
  })

  return(val)
}


optimize_fold <- function(mu, sigma, train_data, test_data) {

  model <- optimize_weights(train_data, mu = mu, sigma = sigma)

  predictions <- predict_probabilities(test_data, model$weights)

  ll <- predictions$loglik

  if (!is.finite(ll)) return(1e12)

  return(ll)
}
