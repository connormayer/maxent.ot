#' Optimizes constraint weights given a data set and optional biases. If no
#' bias arguments are provided, the bias term will not be included in the
#' optimization.
#'
#' TODO: Add description of how optimization functions.
#'
#' @param input_file The path to the input data file. This file contains one or
#'   more OT tableaux consisting of mappings between underlying and surface
#'   forms with observed frequency and violation profiles. Constraint
#'   violations must be numeric.
#' @param bias_file (optional) The path to the file containing mus and sigma
#'   for constraint biases. If this argument is provided, the scalar and vector
#'   mu and sigma arguments will be ignored.
#' @param mu_scalar (optional) A single scalar value that will serve as the mu
#'   for each constraint in the bias term. Constraint weights will also be
#'   initialized to this value. This value will not be used if either bias_file
#'   or mu_vector are provided.
#' @param mu_vector (optional) A vector of mus for each constraint in the bias
#'   term. The length of this vector must equal the number of constraints in
#'   the input file. If bias_file is provided, this argument will be ignored.
#'   If this argument is provided, mu_scalar will be ignored.
#' @param sigma_scalar (optional) A single scalar value that will serve as the
#'   sigma for each constraint in the bias term. This value will not be used if
#'   either bias_file or sigma_vector are provided.
#' @param sigma_vector (optional) A vector of sigmas for each constraint in the
#'   bias term. The length of this vector must equal the number of constraints
#'   in the input file. If bias_file is provided, this argument will be ignored.
#'   If this argument is provided, sigma_scalar will be ignored.
#' @param input_format (optional) A string specifying the format of the input
#'   files. Currently only otsoft-style formatting is supported.
#'
#' @return A named vector containing optimized constraint weights.
#'
#' @examples
#'   optimize_weights('my_tableaux.csv')
#'   optimize_weights('my_tableaux.csv', 'my_biases.csv')
#'   optimize_weights('my_tableaux.csv', mu_vector = c(1,2), sigma_vector = c(100, 200))
#'   optimize_weights('my_tableaux.csv', mu_scalar = 0, sigma_scalar = 1000)
#'   optimize_weights('my_tableaux.csv', mu_vector = c(1, 2), sigma_scalar = 1000)
#'
#' @export
optimize_weights <- function(input_file, bias_file=NA,
                             mu_scalar=NA, mu_vector=NA,
                             sigma_scalar=NA, sigma_vector=NA,
                             input_format='otsoft') {

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

  # If mus aren't provided, initialize all weights to 1
  if (any_not_na(bias_params)) {
    constraint_weights <- bias_params$mus
  } else {
    constraint_weights <- rep(1, length(long_names))
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

# Calculate the log likelihood of the data given the current constraint weights
# and bias parameters. This is the function that is optimized.
calculate_log_likelihood <- function(constraint_weights,
                                     data, bias_params=NA) {
  ll = 0
  # Sum of likelihoods of each datum
  for (tableau in data) {
    tableau_likelihood <- calculate_tableau_likelihood(constraint_weights, tableau)
    ll <- ll + tableau_likelihood
  }
  # Minus the bias term
  if (any_not_na(bias_params)) {
    bias_term <- calculate_bias(bias_params, constraint_weights)
    ll <- ll - bias_term
  }
  return(ll)
}

# Calculates the bias term in the optimized function.
calculate_bias <- function(bias_params, constraint_weights) {
  top <- (constraint_weights - bias_params$mus)^2
  bottom <- 2 * bias_params$sigmas^2
  bias <- sum(top / bottom)

  return(bias)
}

calculate_tableau_likelihood <- function(constraint_weights, tableau) {
  freqs <- tableau[,3]
  violations <- tableau[,4:ncol(tableau)]

  log_candidate_probs <- calculate_tableau_probabilities(
    constraint_weights, violations
  )
  log_prob <- sum(freqs * log_candidate_probs)
  return (log_prob)
}

# Calculates the likelihood of a single tableau.
calculate_tableau_probabilities <- function(constraint_weights, tableau) {
  tableau <- data.matrix(tableau)
  # Replace empty constraint violations with 0
  tableau[is.na(tableau)] <- 0

  # Calculate log probability
  harmony <- tableau %*% constraint_weights
  e_harmonies <- exp(1)^-harmony
  z <- sum(e_harmonies)
  candidate_probs <- e_harmonies / z
  log_candidate_probs <- log(e_harmonies / z)

  return(log_candidate_probs)
}

# Helper function that checks whether any arguments are NA
any_not_na <- function(...) {
  x <- list(...)
  return(any(!is.na(x)))
}

# Function to load and validate bias parameters.
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
      bias_params <- cbind(bias_params, mu_vector)
    }
    else {
      bias_params <- cbind(bias_params, rep(mu_scalar, num_constraints))
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
      bias_params <- cbind(bias_params, sigma_vector)
    }
    else {
      bias_params <- cbind(bias_params, rep(sigma_scalar, num_constraints))
    }
    names(bias_params) <- c("mus", "sigmas")
  } else if (any_not_na(mu_scalar, mu_vector, sigma_scalar, sigma_vector)) {
    stop("You must specify both constraint mus and sigmas, or neither.")
  } else {
    # Don't use a bias parameter
    bias_params <- NA
  }

  # Check that bias params fall within acceptable ranges
  if (any_not_na(bias_params)) {
    if (!all(bias_params$mus >= 0)) {
      stop("All constraint mus must be >= 0")
    }
    if (!all(bias_params$sigmas > 0)) {
      stop("All constraint sigmas must be > 0")
    }
  }
  return(bias_params)
}

