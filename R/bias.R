# Runs optimizt_bias n times, returning the sigma that most benefitted the
# log-likelihood
#' @export
optimize_bias_max <- function(data = NULL, train_data = NULL, test_data = NULL,
                                     k = NULL, sigma = 1, mu = 0,
                                     method = "L-BFGS-B", upper_bound = 100,
                                     control_params = NA,
                                     n = 10) {

  results <- vector("list", n)

  for (i in 1:n) {
    results[[i]] <- optimize_bias(
      data = data,
      train_data = train_data,
      test_data = test_data,
      k = k,
      sigma = sigma,
      mu = mu,
      method = method,
      upper_bound = upper_bound,
      control_params = control_params
    )
  }

  best_index <- which.max(sapply(results, function(x) x$mean_ll))
  best_result <- results[[best_index]]

  best_result$n_runs <- n
  best_result$best_run <- best_index
  best_result$all_results <- results

  return(best_result)
}

#' @export
optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                          k = NULL, sigma = 1, mu = 0,
                          method = "L-BFGS-B", upper_bound = 100,
                          control_params=NA) {

  if (!is.null(train_data) && !is.null(test_data)) {
    fn <- optimization_function
    fn_args <- list(train_data = train_data, test_data = test_data, mu = mu)
  }

  if (!is.na(control_params)) {
    control_params$fnscale <- -1
  }
  else {
    control_params <- list(fnscale = -1)
  }

  sigma <- log(sigma)
  if (length(sigma) > 1) {
    result <- stats::optim(
      par = c(sigma),
      fn=fn,
      train_data = train_data,
      test_data = test_data,
      mu=mu,
      control=control_params,
      # lower = c(0),
      # upper = c(upper_bound),
    )
  }

  else{
    result <- stats::optim(
      par = c(sigma),
      fn=fn,
      train_data = train_data,
      test_data = test_data,
      mu=mu,
      control=control_params,
      lower = c(-2),
      upper = c(2),
      method="Brent"
    )
  }


  sigma_hat <- result$par[1]

  return(list(
    mu = mu,
    sigma = exp(sigma_hat),
    mean_ll = -result$value,
    convergence = result$convergence
  ))
}


optimization_function <- function(par, train_data, test_data, mu) {

  sigma <- exp(par[1])
  # sigma <- par[1]
  print(sigma)

  # print("Optimizing model") DEBUG
  model <- optimize_weights(train_data, mu = mu, sigma = sigma)
  # print("Getting probabilities") DEBUG
  predictions <- predict_probabilities(test_data, model$weights)
  # print("Got predictions") DEBUG
  ll <- predictions$loglik
  print(ll)

  if (!is.finite(ll))
    return(1e12)

  return(ll)
}


# Joint Optimization approach, ONLY FOR TESTING
join_optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                               k = NULL, sigma = 1, mu = 0, init_weights = NULL,
                               method = "L-BFGS-B", upper_bound = 100,
                               control_params = NA) {
  if (!is.null(train_data) && !is.null(test_data)) {
    fn <- joint_optimization_function
  }

  if (!is.na(control_params)) {
    control_params$fnscale <- -1
  }
  else {
    control_params <- list(fnscale = -1)
  }

  if (is.null(init_weights)) {
    init_model <- optimize_weights(train_data, mu = mu, sigma = sigma)
    init_weights <- init_model$weights
  }

  par_init <- c(log(sigma), init_weights)

  result <- stats::optim(
    par = par_init,
    fn = fn,
    train_data = train_data,
    test_data = test_data,
    mu = mu,
    control = control_params,
    method = method
  )

  sigma_hat <- exp(result$par[1])
  weights_hat <- result$par[-1]

  return(list(
    mu = mu,
    sigma = sigma_hat,
    weights = weights_hat,
    mean_ll = result$value,
    convergence = result$convergence
  ))
}

joint_optimization_function <- function(par, train_data, test_data, mu) {
  sigma <- exp(par[1])
  weights <- par[-1]

  predictions <- predict_probabilities(test_data, weights)
  ll <- predictions$loglik

  if (!is.finite(ll)) return(-1e12)
  return(ll)
}

