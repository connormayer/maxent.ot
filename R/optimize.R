#' @export
optimize_weights <- function(input_file=NA, starting_weights=NA,
                             bias_file=NA,
                             mu_scalar=NA, mu_vector=NA,
                             sigma_scalar=NA, sigma_vector=NA) {

  # Organize our inputs
  input <- load_data_otsoft(input_file)
  long_names <- input[[1]]
  short_names <- input[[2]]
  data <- input[[3]]
  num_constraints <- length(long_names)
  bias_params <- process_bias_arguments(
    bias_file, mu_scalar, mu_vector, sigma_scalar, sigma_vector,
    num_constraints
  )

  # If starting weights aren't provided, initialize all weights to 1
  if (is.na(starting_weights)) {
    constraint_weights <- rep(1, length(long_names))
  } else {
    constraint_weights <- starting_weights
  }

  # Perform optimization
  best <- optim(
    constraint_weights,
    calculate_log_likelihood,
    data=data,
    bias_params=bias_params,
    control=list(fnscale = -1),
    lower=rep(0, length(constraint_weights)),
    # The default upper bound is Inf, but the function we're optimizing
    # can't be evaluated at Inf. This results in the optimizer finding
    # non-optimal weights, so we pass in a large finite value instead.
    upper=rep(2^64, length(constraint_weights)),
    method="L-BFGS-B"
  )
  out_weights <- best[[1]]
  names(out_weights) <- long_names
  return(out_weights)
}

calculate_log_likelihood <- function(constraint_weights, data, bias_params=NA) {
  ll = 0
  # Sum of likelihoods of each datum
  for (datum in data) {
    datum_likelihood <- calculate_datum_likelihood(constraint_weights, datum)
    ll <- ll + datum_likelihood
  }
  # Minus the bias term
  if (any_not_na(bias_params)) {
    bias_term <- calculate_bias(bias_params, constraint_weights)
    ll <- ll - bias_term
  }
  return(ll)
}

calculate_bias <- function(bias_params, constraint_weights) {
  top <- (constraint_weights - bias_params$mus)^2
  bottom <- 2 * bias_params$sigmas^2
  bias <- sum(top / bottom)
  return(bias)
}

calculate_datum_likelihood <- function(constraint_weights, datum) {
  freqs <- datum[,3]
  violations <- datum[,4:ncol(datum)]
  violations_mat <- data.matrix(violations)
  violations_mat[is.na(violations_mat)] <- 0
  harmony <- violations_mat %*% constraint_weights
  e_harmonies <- exp(1)^-harmony
  z <- sum(e_harmonies)
  candidate_probs <- e_harmonies / z
  log_candidate_probs <- log(e_harmonies / z)
  log_prob <- sum(freqs * log_candidate_probs)
  return(log_prob)
}

any_not_na <- function(...) {
  x <- list(...)
  return(any(!is.na(x)))
}

process_bias_arguments <- function(bias_file=NA, mu_scalar=NA, mu_vector=NA,
                                   sigma_scalar=NA, sigma_vector=NA,
                                   num_constraints=NA) {
  if (any_not_na(bias_file)) {
    # Read bias parameters from provided file location
    if (any_not_na(mu_scalar, mu_vector, sigma_vector, sigma_vector)) {
      warning(
        "Both a bias file and bias scalars/vectors were provided\n",
        "Ignoring scalars/vectors and using parameters from file"
      )
    }
    bias_params <- load_bias_file_otsoft(bias_file)
  } else if (any_not_na(mu_scalar, mu_vector) &
      any_not_na(sigma_scalar,sigma_vector)) {
    # Set bias parameters from provided scalars/vectors
    bias_params <- data.table::data.table()

    # Set our mus
    if (any_not_na(mu_vector)) {
      if (any_not_na(mu_scalar)) {
        warning(
          "Both a vector and a scalar for mu were provided\n",
          "Ignoring scalar and using vector parameters"
        )
      }
      if (length(mu_vector) != num_constraints) {
        stop(sprintf(
          paste(
            "Number of constraint mus (%d) and number of constraints in",
            "data (%d) do not match"
          ),
          length(mu_vector), num_constraints)
        )
      }
      bias_params$mus <- mu_vector
    }
    else {
      bias_params$mus <- rep(mu_scalar, num_constraints)
    }
    # Set our sigmas
    if (any_not_na(sigma_vector)) {
      if (any_not_na(sigma_scalar)) {
        warning(
          "Both a vector and a scalar for sigma were provided\n",
          "Ignoring scalar and using vector parameters"
        )
      }
      if (length(sigma_vector) != num_constraints) {
        stop(sprintf(
          paste(
            "Number of constraint sigmas (%d) and number of constraints in",
            "data (%d) do not match"
          ),
          length(sigma_vector), num_constraints)
        )
      }
      bias_params$sigmas <- sigma_vector
    }
    else {
      bias_params$sigmas <- rep(sigma_scalar, num_constraints)
    }
  } else if (any_not_na(mu_scalar, mu_vector, sigma_scalar, sigma_vector)) {
    stop("You must specify both constraint mus and sigmas, or neither.")
  } else {
    # Don't use a bias parameter
    bias_params <- NA
  }

  # Check that bias params fall within acceptable ranges
  if (any_not_na(bias_params)) {
    if (!all(bias_params$mu >= 0)) {
      stop("All constraint mus must be >= 0")
    }
    if (!all(bias_params$sigma > 0)) {
      stop("All constraint sigmas must be > 0")
    }
  }

  return(bias_params)
}
