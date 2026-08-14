# Runs optimize_bias n times, returning the sigma that most benefitted the
# log-likelihood
#' @export
optimize_bias_max <- function(train_data, test_data,
                             mu = 0,
                             method = "L-BFGS-B",
                             upper_bound = 100,
                             control_params = NA) {

  n <- length(train_data)

  results <- vector("list", n)

  for (i in seq_len(n)) {
    results[[i]] <- optimize_bias(
      train_data = train_data[[i]],
      test_data = test_data[[i]],
      mu = mu,
      method = method,
      upper_bound = upper_bound,
      control_params = control_params
    )
  }

  sigmas <- sapply(results, function(x) x$sigma)
  mean_lls <- sapply(results, function(x) x$mean_ll)

  best_sigma <- median(sigmas)

  list(
    sigma = best_sigma,
    mu = mu,
    sigmas = sigmas,
    mean_ll = mean(mean_lls),
    results = results
  )
}


#' @export
optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                          k = NULL, sigma = 1, mu = 0,
                          method = "L-BFGS-B", upper_bound = 100,
                          control_params = NA) {

  if (!is.null(train_data) && !is.null(test_data)) {
    fn <- optimization_function
  }

  if (!is.na(control_params)) {
    control_params$fnscale <- -1
  } else {
    control_params <- list(fnscale = -1)
  }

  # Environment used to record every evaluation made by optim()
  history <- new.env()
  history$df <- data.frame(
    sigma = numeric(),
    ll = numeric()
  )

  sigma <- log(sigma)

  if (length(sigma) > 1) {

    result <- stats::optim(
      par = c(sigma),
      fn = fn,
      train_data = train_data,
      test_data = test_data,
      mu = mu,
      history = history,
      control = control_params
      # lower = c(0),
      # upper = c(upper_bound)
    )

  } else {

    result <- stats::optim(
      par = c(sigma),
      fn = fn,
      train_data = train_data,
      test_data = test_data,
      mu = mu,
      history = history,
      control = control_params,
      lower = c(-2),
      upper = c(2),
      method = "Brent"
    )

  }

  sigma_hat <- result$par[1]

  return(list(
    mu = mu,
    sigma = exp(sigma_hat),
    mean_ll = -result$value,
    convergence = result$convergence,
    history = history$df
  ))
}


optimization_function <- function(par, train_data, test_data,mu, history) {

  sigma <- exp(par[1])

  print(sigma)

  model <- optimize_weights(
    train_data,
    mu = mu,
    sigma = sigma
  )

  predictions <- predict_probabilities(
    test_data,
    model$weights
  )

  ll <- predictions$loglik

  print(ll)

  history$df <- rbind(
    history$df,
    data.frame(
      sigma = sigma,
      ll = ll
    )
  )

  if (!is.finite(ll))
    return(1e12)

  return(ll)
}
