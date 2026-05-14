#' @export
optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                          k = NULL, sigma = 0.5, mu = 0,
                          method = "L-BFGS-B", upper_bound = 100) {

  if (!is.null(train_data) && !is.null(test_data)) {
    fn <- optimization_function
    fn_args <- list(train_data = train_data, test_data = test_data, mu = mu)
  }

  result <- do.call(
    stats::optim,
    c(
      list(
        par = c(sigma),
        fn = fn,
        method = method,
        lower = c(0),
        upper = c(upper_bound)
      ),
      fn_args
    )
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

  val <- tryCatch(
    {

      sigma <- exp(par[1])

      ll <- optimize_fold(mu, sigma, train_data, test_data)

      if (!is.finite(ll))
        return(1e12)

      return(-ll)
    }, error = function(e) {
        1e12
      }
    )

  return(val)
}



optimize_fold <- function(mu, sigma, train_data, test_data) {

  model <- optimize_weights(train_data, mu = mu, sigma = sigma)

  predictions <- predict_probabilities(test_data, model$weights)

  ll <- predictions$loglik

  if (!is.finite(ll))
    return(1e12)

  return(ll)
}
