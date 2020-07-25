# A few constants
DEFAULT_WEIGHT <- 1

#' Optimize maxent OT constraint weights
#'
#' Optimizes constraint weights given a data set and optional biases. If no
#' bias arguments are provided, the bias term will not be included in the
#' optimization.
#'
#' By default, this function optimizes the log likelihood of the training data
#' by changing the values of \eqn{w}, the vector of constraint weights. The log
#' likelihood of the training data is
#'
#' \deqn{LL_w(D) = \sum_{i=1}^{n}{\ln L_w(x)}
#' - \sum_{k=1}^{m}{\frac{(w_k - \mu_k)^2}{2\sigma_k^2}}}
#'
#' where \eqn{w_k} is the weight of constraint \eqn{k}, and \eqn{\mu_k} and
#' \eqn{\sigma_k} parameterize a normal distribution that serves as a
#' prior bias for \eqn{w_k}.
#'
#' \eqn{L_(x)} is the likelihood of a single tableau whose input is the
#' underlying form \eqn{x}. \eqn{L_(x)} is defined as
#'
#' \deqn{L_w(x) = \prod_{y\in \mathcal{Y}(x)}{\nu(y)P(y|x; w)}}
#'
#' where \eqn{\mathcal{Y}(x)} is the set of observed surface realizations of
#' \eqn{x}, \eqn{\nu(y)} is the number of observed tokens of surface form
#' \eqn{y}, and \eqn{P(y|x; w)} is the probability of realizing underlying
#' \eqn{x} as surface \eqn{y} given the constraint weighting \eqn{w}. This
#' probability is defined as
#'
#' \deqn{P(y|x; w) = \frac{1}{Z_w(x)}\exp(\sum_{k=1}^{m}{w_k f_k(y, x)})}
#'
#' where \eqn{f_k(y, x)} is the number of violatons of constraint \eqn{k}
#' incurred by mapping underlying \eqn{x} to surface \eqn{y}, and \eqn{Z(x)}
#' is a normalization term defined as
#'
#' \deqn{Z(x) = \sum_{y\in \mathcal{Y}(x)}{\exp(\sum_{k=1}^{m}{w_k f_k(y, x)})}}
#'
#' Optimization is done using the `optim` function from the R-core statistics
#' library. By default it uses `L-BFGS-B` optimization, which is a quasi-Newton
#' method that allows upper and lower bounds on variables. Constraint weights
#' are restricted to finite, non-negative values.
#'
#' If no bias parameters are specified (either the `bias_file` argument or some
#' combination of the scalar/vector mu/sigma parameters), optimization will be
#' done without the bias term.
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
#'   initialized to this value. This value will not be used if either
#'   `bias_file` or `mu_vector` are provided.
#' @param mu_vector (optional) A vector of mus for each constraint in the bias
#'   term. The length of this vector must equal the number of constraints in
#'   the input file. If `bias_file` is provided, this argument will be
#'   ignored. If this argument is provided, `mu_scalar` will be ignored.
#' @param sigma_scalar (optional) A single scalar value that will serve as the
#'   sigma for each constraint in the bias term. This value will not be used if
#'   either `bias_file` or `sigma_vector` are provided.
#' @param sigma_vector (optional) A vector of sigmas for each constraint in the
#'   bias term. The length of this vector must equal the number of constraints
#'   in the input file. If `bias_file` is provided, this argument will be ignored.
#'   If this argument is provided, `sigma_scalar` will be ignored.
#' @param input_format (optional) A string specifying the format of the input
#'   files. Currently only OTSoft-style formatting is supported.
#' @param in_sep (optional) The delimiter used in the input files. Defaults to
#'   tabs.
#' @param control_params (optional) A named list of control parameters that
#'   will be passed to the `optim` function. See documentation of that function
#'   for details. Note that some parameter settings may interfere with
#'   optimization. The parameter `fnscale` will be overwritten to `-1` if
#'   specified, since this must be treated as a maximization problem.
#' @param upper_bound (optional) The maximum value for constraint weights.
#' @param model_name (optional) A name for the model. If not provided, the file
#'   name will be used.
#' @return An object with the following named attributes:
#'         * `weights`: A named list of the optimal constraint weights
#'         * `log_lik`: the log likelihood of the data under the discovered
#'           weights
#'         * `k`: the number of constraints
#'         * `n`: the number of data points in the training set
#'
#' @examples
#'   optimize_weights('my_tableaux.csv')
#'   optimize_weights('my_tableaux.csv', 'my_biases.csv')
#'   optimize_weights('my_tableaux.csv', mu_vector = c(1, 2), sigma_vector = c(100, 200))
#'   optimize_weights('my_tableaux.csv', mu_scalar = 0, sigma_scalar = 1000)
#'   optimize_weights('my_tableaux.csv', mu_vector = c(1, 2), sigma_scalar = 1000)
#'   optimize_weights('my_tableau.csv, control_params = list(maxit = 500))
#'
#' @export
optimize_weights <- function(input_file, bias_file = NA,
                             mu_scalar = NA, mu_vector = NA,
                             sigma_scalar = NA, sigma_vector = NA,
                             input_format = 'otsoft', in_sep = '\t',
                             control_params = NA, upper_bound = 100,
                             model_name = NA) {

  # Organize our inputs
  input <- load_data_otsoft(input_file, sep = in_sep)
  long_names <- input$full_names
  short_names <- input$abbr_names
  data <- input$data
  n <- input$n
  num_constraints <- length(long_names)
  bias_params <- process_bias_arguments(
    bias_file, mu_scalar, mu_vector, sigma_scalar, sigma_vector,
    num_constraints
  )

  # If no model name provided, use filename sans extension
  if (is.na(model_name)) {
    model_name <- tools::file_path_sans_ext(basename(input_file))
  }

  # If mus aren't provided, initialize all weights to 1
  # TODO: Does initializing contraints to the mus make sense?
  if (any_not_na(bias_params)) {
    constraint_weights <- bias_params$mus
  } else {
    constraint_weights <- rep(DEFAULT_WEIGHT, length(long_names))
  }

  # Pass through control parameters specified by the user, but fnscale needs to
  # be -1 to turn this into a maximization problem
  if (!is.na(control_params)) {
    control_params$fnscale <- -1
  }
  else {
    control_params <- list(fnscale = -1)
  }

  # Build ourselves a matrix for efficient computation
  # Pre-alocate space
  data_matrix <- matrix(0L, nrow = nrow(data), ncol = ncol(data) + 2)
  # Map URs to integers
  data_matrix[,1] <- as.integer(as.factor(data[,1]))
  # Set the violation profiles
  data_matrix[,2:(ncol(data_matrix) - 3)] <- data.matrix(data[,3:ncol(data)])
  # Replace empty cells with 0
  data_matrix[is.na(data_matrix)] <- 0

  # Perform optimization
  best <- tryCatch({
    optim(
      constraint_weights,
      calculate_log_likelihood_helper,
      data=data_matrix,
      bias_params=bias_params,
      control=control_params,
      lower=rep(0, length(constraint_weights)),
      # The default upper bound is Inf, but the function we're optimizing
      # can't be evaluated at Inf. This results in the optimizer finding
      # non-optimal weights, so we pass in a large finite value instead.
      upper=rep(upper_bound, length(constraint_weights)),
      method="L-BFGS-B"
    )
  },
  error=function(cond){
      if (cond$message == "L-BFGS-B needs finite values of 'fn'") {
        message("\nThis error indicates that the likelihood function has ",
                "returned a non-finite value because the constraint ",
                "weights have become too large. You can resolve this by ",
                "passing in a lower upper bound on maximum weights using ",
                "the `upper_bound' argument, or by introducing a stronger ",
                "bias towards lower weights by manipulating the mu and sigma ",
                "parameters.")
      }
      stop(cond)
    }
  )
  out_weights <- best[[1]]
  names(out_weights) <- long_names
  out_object <- list(
    name = model_name,
    weights = out_weights,
    loglik = best$value,
    k = length(out_weights),
    n = n,
    bias_params = bias_params
  )
  return(out_object)
}

# Calculate the log likelihood of the data given the current constraint weights
# and bias parameters. This is the function that is optimized.
calculate_log_likelihood_helper <- function(constraint_weights,
                                     data, bias_params=NA) {
  # Set a few column indexes
  freq_ix <- 2
  log_prob_ix <- ncol(data) - 1
  lik_ix <- log_prob_ix + 1

  # Caluclate log likelihood of data
  data <- calculate_probabilities(constraint_weights, data)
  ll <- sum(data[,freq_ix] * data[, log_prob_ix])

  # Minus the bias term
  if (any_not_na(bias_params)) {
    bias_term <- calculate_bias(bias_params, constraint_weights)
    ll <- ll - bias_term
  }
  return(ll)
}

# Calculate probabilities for all candidates based on current constraint
# weights
calculate_probabilities <- function(constraint_weights, data) {
  freq_ix <- 2
  harm_ix <- ncol(data) - 2
  log_prob_ix <- harm_ix + 1

  data[, harm_ix] <- data[, (freq_ix + 1):(harm_ix - 1)] %*% constraint_weights
  data[, harm_ix] <- exp(1)^-data[, harm_ix]
  data[, log_prob_ix] <- log(apply(data, 1, normalize_row, data, harm_ix))
  return(data)
}

# Helper function that applies Z normalization
normalize_row <- function(row, m, col_num) {
  return(row[col_num] / sum(m[m[,1] == row[1],][,col_num]))
}

# Calculates the bias term in the optimized function.
calculate_bias <- function(bias_params, constraint_weights) {
  top <- (constraint_weights - bias_params$mus)^2
  bottom <- 2 * bias_params$sigmas^2
  bias <- sum(top / bottom)

  return(bias)
}

# Helper function that checks whether any arguments are NA
any_not_na <- function(...) {
  x <- list(...)
  return(any(!is.na(x)))
}

# Function to load and validate bias parameters.
process_bias_arguments <- function(bias_file = NA,
                                   mu_scalar = NA, mu_vector = NA,
                                   sigma_scalar = NA, sigma_vector=NA,
                                   num_constraints = NA, sep = '\t') {
  if (any_not_na(bias_file)) {
    # Read bias parameters from provided file location
    if (any_not_na(mu_scalar, mu_vector, sigma_vector, sigma_vector)) {
      warning(
        "Both a bias file and bias scalars/vectors were provided\n",
        "Ignoring scalars/vectors and using parameters from file"
      )
    }
    bias_params <- load_bias_file_otsoft(bias_file, sep = sep)
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
