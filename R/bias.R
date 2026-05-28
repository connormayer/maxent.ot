#' @export
optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                          k = NULL, sigma = 0.5, mu = 0,
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


  result <- stats::optim(
    par = c(sigma),
    fn=fn,
    train_data = train_data,
    test_data = test_data,
    mu=mu,
    control=control_params,
    # lower = c(0),
    # upper = c(upper_bound),
    # method="L-BFGS-B"
  )

  sigma_hat <- result$par[1]

  return(list(
    mu = mu,
    sigma = sigma_hat,
    mean_ll = -result$value,
    convergence = result$convergence
  ))
}


optimization_function <- function(par, train_data, test_data, mu) {

  sigma <- exp(par[1])
  # sigma <- par[1]
  # print(sigma) DEBUG

  # print("Optimizing model") DEBUG
  model <- optimize_weights(train_data, mu = mu, sigma = sigma)
  # model <- optimize_weights(train_data, mu = mu, sigma = sigma)
  # print("Getting probabilities") DEBUG
  predictions <- predict_probabilities(test_data, model$weights)
  # print("Got predictions") DEBUG
  ll <- predictions$loglik
  # print(ll) DEBUG

  if (!is.finite(ll))
    return(1e12)

  return(ll)
}
