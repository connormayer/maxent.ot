#' Finds optimal sigma value across multiple train/test sets
#'
#' Runs `optimize_bias` once for each train/test set and returns
#' the mean sigma across those runs.
#'
#' Valid input train/test data is created through `optimize_bias_helper`.
#'
#' @param train_data A vector of n tableaux that are treated as training data
#' across each internal call to `optimize_bias`
#'
#' @param test_data A vector of n tableaux that are treated as testing data
#' across each internal call to `optimize_bias`
#'
#' @param mu The mu value used for the center of the bias argument variation,
#' defaulting to 0
#'
#' @param sigma The starting sigma value which `optimize_bias` optimizes,
#' defaulting to 1, representing the standard deviation of the bias argument
#'
#' @export
optimize_bias_max <- function(train_data, test_data, mu = 0, sigma = 1) {

  n <- length(train_data)

  results <- vector("list", n)

  for (i in seq_len(n)) {

    results[[i]] <- optimize_bias(
      train_data = train_data[[i]],
      test_data = test_data[[i]],
      mu = mu,
      sigma = sigma
    )
  }

  sigmas <- sapply(results, function(x) x$sigma)
  mean_lls <- sapply(results, function(x) x$mean_ll)

  best_sigma <- mean(sigmas)

  list(
    sigma = best_sigma,
    mu = mu,
    sigmas = sigmas,
    mean_ll = mean(mean_lls),
    results = results
  )
}


#' Finds the optimal sigma value for a single train/test set
#'
#' Optimizes sigma by maximizing the log-likelihood of the predictions
#' made on the test data.
#'
#' Sigma is optimized on the log scale so that sigma is always positive.
#'
#' @param data (optional) Path to dataset, not necessary if being called through
#' `optimize_bias_max`
#'
#' @param train_data Training tableau used to optimize the model weights
#'
#' @param test_data Testing tableau used to evaluate the resulting model
#'
#' @param k (optional) Number of partitions, not necessary if being called
#' through `optimize_bias_max`
#'
#' @param sigma Initial sigma value, defaulting to 1
#'
#' @param mu The mu value used for the center of the bias argument variation,
#' defaulting to 0
#'
#' @export
optimize_bias <- function(data = NULL, train_data = NULL, test_data = NULL,
                          k = NULL, sigma = 1, mu = 0) {

  if (!is.null(train_data) && !is.null(test_data)) {
    fn <- optimization_function
  }

  history <- new.env()
  history$df <- data.frame(
    sigma = numeric(),
    ll = numeric()
  )

  sigma <- log(sigma)

  result <- stats::optim(
    par = c(sigma),
    fn = fn,
    train_data = train_data,
    test_data = test_data,
    mu = mu,
    history = history,
    control = list(fnscale = -1),
    lower = c(-2),
    upper = c(2),
    method = "Brent"
  )

  sigma_hat <- result$par[1]

  return(list(
    mu = mu,
    sigma = exp(sigma_hat),
    mean_ll = -result$value,
    convergence = result$convergence,
    history = history$df
  ))
}


#' Objective function used to optimize sigma
#'
#' For a given sigma:
#'  Optimizes the model weights on the training data
#'  Uses those weights to make predictions on the test data
#'  Calculates the test log-likelihood
#'
#' The resulting log-likelihood is maximized by `optimize_bias`.

optimization_function <- function(par, train_data, test_data, mu, history) {

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
    return(-1e12)

  return(ll)
}


#' Creates train/test sets and finds the optimal sigma
#'
#' Reads the input data, partitions it into n sets, creates
#' corresponding train/test tableaux, and passes them to
#' `optimize_bias_max`.
#'
#' @param path Path to the input CSV file
#'
#' @param n Number of partitions to create
#'
#' @export
optimize_bias_helper <- function(path, n, mu = 0, sigma = 1) {

  input <- read_csv(path)
  processed_input <- maxent.ot:::load_input(input)
  data <- processed_input$data

  partitions <- maxent.ot:::partition_data(data, k = n)

  train_data_list <- vector("list", n)
  test_data_list <- vector("list", n)

  training_data <- data %>% mutate(Frequency = 0)
  test_data <- training_data

  for (hold_out in 1:n) {

    training_part <- partitions %>%
      filter(partition != hold_out)

    test_part <- partitions %>%
      filter(partition == hold_out)

    train_data_list[[hold_out]] <-
      maxent.ot:::populate_tableau(
        training_data,
        training_part
      )

    test_data_list[[hold_out]] <-
      maxent.ot:::populate_tableau(
        test_data,
        test_part
      )
  }

  result <- optimize_bias_max(
    train_data = train_data_list,
    test_data = test_data_list,
    mu = mu,
    sigma = sigma
  )

  cat("Best sigma:", result$sigma, "\n")
  cat("Best log-likelihood:", result$mean_ll, "\n")

  cat("\nSigma for each split:\n")
  print(result$sigmas)
}
