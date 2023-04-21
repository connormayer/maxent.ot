# Constants
DEFAULT_WEIGHT <- 0
DEFAULT_UPPER_BOUND <- 100

#' Optimize MaxEnt OT constraint weights
#'
#' Optimizes constraint weights given a data set and optional biases. If no
#' bias arguments are provided, the bias term(s) will not be included in the
#' optimization.
#'
#' The objective function \eqn{J(w)} that is optimized is defined as
#'
#' \deqn{J(w) = \sum_{i=1}^{n}{\ln P(y_i|x_i; w)}
#' - \sum_{k=1}^{m}{\frac{(w_k - \mu_k)^2}{2\sigma_k^2}}}
#'
#' The first term in this equation calculates the natural logarithm of the
#' conditional likelihood of the training data under the weights \eqn{w}. \eqn{n}
#' is the number of data points (i.e., the sample size or the sum of the frequency
#' column in the input),\eqn{x_i} is the input form of the \eqn{i}th data
#' point, and \eqn{y_i} is the observed surface form corresponding to
#' \eqn{x_i}.\eqn{P(y_i|x_i; w)} represents the probability of realizing
#' underlying \eqn{x_i} as surface \eqn{y_i} given weights \eqn{w}. This
#' probability is defined as
#'
#' \deqn{P(y_i|x_i; w) = \frac{1}{Z_w(x_i)}\exp(-\sum_{k=1}^{m}{w_k f_k(y_i, x_i)})}
#'
#' where \eqn{f_k(y_i, x_i)} is the number of violations of constraint \eqn{k}
#' incurred by mapping underlying \eqn{x_i} to surface \eqn{y_i}. \eqn{Z_w(x_i)}
#' is a normalization term defined as
#'
#' \deqn{Z(x_i) = \sum_{y\in\mathcal{Y}(x_i)}{\exp(-\sum_{k=1}^{m}{w_k f_k(y, x_i)})}}
#'
#' where \eqn{\mathcal{Y}(x_i)} is the set of observed surface realizations of
#' input \eqn{x_i}.
#'
#' The second term of the equation for calculating the objective function is
#' the optional bias term, where \eqn{w_k} is the weight of constraint \eqn{k}, and
#' \eqn{\mu_k} and \eqn{\sigma_k} parameterize a normal distribution that
#' serves as a prior for the value of \eqn{w_k}. \eqn{\mu_k} specifies the mean
#' of this distribution (the expected weight of constraint \eqn{k} before
#' seeing any data) and \eqn{sigma_k} reflects certainty in this value: lower
#' values of \eqn{\sigma_k} penalize deviations from \eqn{\mu_k} more severely,
#' and thus require greater amounts of data to move \eqn{w_k} away from
#' \eqn{mu_k}. While increasing \eqn{\sigma_k} will improve the fit to the
#' training data, it may result in overfitting, particularly for small data
#' sets.
#'
#' A general bias with \eqn{\mu_k = 0} for all \eqn{k} is commonly used as a
#' form of simple regularization to prevent overfitting (see, e.g., Goldwater
#' and Johnson 2003). Bias terms have also been used to model proposed
#' phonological learning biases; see for example Wilson (2006), White (2013),
#' and Mayer (2021, Ch. 4). The choice of \eqn{\sigma} depends on the sample
#' size. As the number of data points increases, \eqn{\sigma} must decrease in
#' order for the effect of the bias to remain constant: specifically,
#' \eqn{n\sigma^2} must be held constant, where \eqn{n} is the number of tokens.
#'
#' Optimization is done using the \link[stats]{optim} function from the R-core
#' statistics library. By default it uses `L-BFGS-B` optimization, which is a
#' quasi-Newtonian method that allows upper and lower bounds on variables.
#' Constraint weights are restricted to finite, non-negative values.
#'
#' If no bias parameters are specified (either the `bias_file` argument or some
#' combination of the scalar/vector mu/sigma parameters), optimization will be
#' done without the bias term.
#'
#' @param input The input data frame/data table/tibble. This should contain one
#'   or more OT tableaux consisting of mappings between underlying and surface
#'   forms with observed frequency and violation profiles. Constraint violations
#'   must be numeric.
#'
#'   For an example of the data frame format, see inst/extdata/sample_data_frame.csv.
#'   You can read this file into a data frame using read.csv or into a tibble
#'   using dplyr::read_csv.
#'
#'   This function also supports the legacy OTSoft file format. You can use this
#'   format by passing in a file path string to the OTSoft file rather than a
#'   data frame.
#'
#'   For examples of OTSoft format, see inst/extdata/sample_data_file.txt.
#' @param bias_input (optional)
#'   A data frame/data table/tibble containing the bias mus and sigmas. Each row
#'   corresponds to an individual constraint, and consists of three columns:
#'   `Constraint`, which contains the constraint name, `Mu`, which contains the
#'   mu, and `Sigma`, which contains the sigma. If this argument is provided,
#'   the scalar and vector mu and sigma arguments will be ignored.

#'   Like the `input` argument, this function also supports the legacy OTSoft
#'   file format for this argument. In this case, `bias_input` should be a path
#'   to the bias parameters in OTSoft format.
#'
#'   For examples of OTSoft bias format, see inst/extdata/sample_bias_file_otsoft.txt.
#'   Each row in this file should be the name of the constraint, followed by the
#'   mu, followed by the sigma (separated by tabs).
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
#' @param control_params (optional) A named list of control parameters that
#'   will be passed to the \link[stats]{optim} function. See the documentation
#'   of that function for details. Note that some parameter settings may
#'   interfere with optimization. The parameter `fnscale` will be overwritten
#'   with `-1` if specified, since this must be treated as a maximization
#'   problem.
#' @param upper_bound (optional) The maximum value for constraint weights.
#'   Defaults to 100.
#' @param encoding (optional) The character encoding of the input file. Defaults
#'  to "unknown".
#' @param model_name (optional) A name for the model. If not provided, the name
#'   of the variable will be used if the input is a data frame. If the input
#'   is a path to an OTSoft file, the filename will be used.
#' @return An object with the following named attributes:
#' \itemize{
#'         \item `weights`: A named list of the optimal constraint weights
#'         \item `log_lik`: the log likelihood of the data under the discovered
#'           weights
#'         \item `k`: the number of constraints
#'         \item `n`: the number of data points in the training set
#' }
#' @examples
#'   # Get paths to toy data and bias files.
#'   df_file <- system.file(
#'       "extdata", "sample_data_frame.csv", package = "maxent.ot"
#'   )
#'   bias_file <- system.file(
#'        "extdata", "sample_bias_data_frame.csv", package = "maxent.ot"
#'   )

#'   # Fit weights to data with no biases
#'   tableaux_df <- read.csv(df_file)
#'   optimize_weights(tableaux_df)
#'
#'   # Fit weights with biases specified in file
#'   bias_df <- read.csv(bias_file)
#'   optimize_weights(tableaux_df, bias_df)
#'
#'   # Fit weights with biases specified in vector form
#'   optimize_weights(
#'       tableaux_df, mu_vector = c(1, 2), sigma_vector = c(100, 200)
#'   )
#'
#'   # Fit weights with biases specified as scalars
#'   optimize_weights(tableaux_df, mu_scalar = 0, sigma_scalar = 1000)
#'
#'   # Fit weights with mix of scalar and vector biases
#'   optimize_weights(tableaux_df, mu_vector = c(1, 2), sigma_scalar = 1000)
#'
#'   # Pass additional arguments to optim function
#'   optimize_weights(tableaux_df, control_params = list(maxit = 500))
#'
#' @export
optimize_weights <- function(input, bias_input = NA,
                             mu_scalar = NA, mu_vector = NA,
                             sigma_scalar = NA, sigma_vector = NA,
                             control_params = NA,
                             upper_bound = DEFAULT_UPPER_BOUND,
                             encoding = 'unknown', model_name = NA) {
  if (is.data.frame(input)) {
    if (is.na(model_name)) {
      # If no provided model name, use name of input variable
      model_name <- toString(substitute(input))
    }
  }
  processed_input <- load_input(
    input, encoding = encoding, model_name = model_name
  )
  long_names <- processed_input$long_names
  data <- processed_input$data
  n <- processed_input$n
  model_name <- processed_input$model_name

  num_constraints <- length(long_names)
  bias_params <- process_bias_arguments(
    long_names, bias_input, mu_scalar, mu_vector, sigma_scalar, sigma_vector,
    num_constraints
  )

  # If mus aren't provided, initialize all weights to 0
  if (any_not_na(bias_params)) {
    constraint_weights <- bias_params$Mu
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

  # Convert strings in violation profiles to numbers
  data[, 3:ncol(data)] <- lapply(data[, 3:ncol(data)], as.numeric)

  # Build ourselves a matrix for efficient computation
  # Pre-allocate space
  data_matrix <- matrix(0L, nrow = nrow(data), ncol = ncol(data) + 2)
  # Map URs to integers
  data_matrix[, 1] <- as.integer(as.factor(data[, 1]))
  # Set the violation profiles
  data_matrix[, 2:(ncol(data_matrix) - 3)] <- data.matrix(data[, 3:ncol(data)])
  # Replace empty cells with 0
  data_matrix[is.na(data_matrix)] <- 0

  # Perform optimization
  best <- tryCatch({
    stats::optim(
      constraint_weights,
      objective_func,
      gr=calculate_gradient,
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
    loglik = calculate_log_likelihood(out_weights, data_matrix),
    k = length(out_weights),
    n = n,
    bias_params = bias_params
  )
  return(out_object)
}

calculate_log_likelihood <- function(constraint_weights, data) {
  # Set a few column indexes
  freq_ix <- 2
  log_prob_ix <- ncol(data) - 1
  lik_ix <- log_prob_ix + 1

  # Calculate log likelihood of data
  data <- calculate_probabilities(constraint_weights, data)
  ll <- sum(data[, freq_ix] * data[, log_prob_ix])

  return(ll)
}

# This is the function that is optimized.
objective_func <- function(constraint_weights, data, bias_params=NA) {
  ll <- calculate_log_likelihood(constraint_weights, data)

  # Minus the bias term
  if (any_not_na(bias_params)) {
    bias_term <- calculate_bias(bias_params, constraint_weights)
    ll <- ll - bias_term
  }

  return(ll)
}

# Calculate probabilities for all candidates based on current constraint
# weights
calculate_probabilities <- function(constraint_weights, data,
                                    temperature = DEFAULT_TEMPERATURE) {
  freq_ix <- 2
  harm_ix <- ncol(data) - 2
  log_prob_ix <- harm_ix + 1

  data[, harm_ix] <- data[, (freq_ix + 1):(harm_ix - 1)] %*% matrix(constraint_weights)
  data[, harm_ix] <- exp(-data[, harm_ix] / temperature)
  data[, log_prob_ix] <- log(apply(data, 1, normalize_row, data, harm_ix))
  return(data)
}

# Calculate gradients
calculate_gradient <- function(constraint_weights,
                               data,
                               temperature=DEFAULT_TEMPERATURE,
                               bias_params=NA){

  # Gradients for log likelihood wrt each constraint
  grad_ll <- calculate_grad_ll(constraint_weights, data, temperature)

  # If there's a bias term
  if (any_not_na(bias_params)) {
    # Calculate gradients for bias term wrt each constraint
    grad_bias <- calculate_grad_bias(bias_params, constraint_weights)
    # Include bias contribution to gradient
    grad <- grad_ll - grad_bias
  } else {
    # Else gradients consist solely of d(ll)/dCon
    grad <- grad_ll
  }

  return(grad)
}

# Calculate gradients for log likelihood
calculate_grad_ll <- function(constraint_weights,
                              data,
                              temperature){

  # Get expected "feature activation" matrix
  expt_mat <- expectation_mat(constraint_weights, data, temperature)

  # Get feature matrix
  feat_mat <- data[, 3:(ncol(data)-3)]

  # Get frequency matrix
  freq_mat <- data[, 2]

  # Get gradients of log likelihood wrt constraints
  # Produces a (1,n) matrix, where n is number of constraints
  grad_ll <- t(freq_mat) %*% (expt_mat - feat_mat)
  # Convert matrix to vector
  grad_ll <- c(grad_ll)

  return(grad_ll)
}

# Make expectation matrix
expectation_mat <- function(constraint_weights,
                            data,
                            temperature){

  # Get current log probabilities for all candidates
  # This returns a matrix with:
  # col1: UR's identified with indices (each unique UR is assigned a unique ix)
  # col2: freq
  # ncol(data_mat) - 1: log prob
  # feature matrix: (ur_ix + 2):(log_prob_ix - 2)
  data_mat <- calculate_probabilities(constraint_weights,
                                      data,
                                      temperature)

  # Set a few column indices
  ur_ix <- 1
  log_prob_ix <- ncol(data_mat) - 1

  # Initialize empty expected "feature activation" matrix
  exp_mat <- matrix(0L, nrow=nrow(data_mat), ncol=ncol(data_mat)-5)

  # For each row
  for (i in 1:nrow(data_mat)){

    # Get current UR's id
    ur_id <- data_mat[i, 1]

    # Get feature matrix for all URs that share the same UR id
    # UR ids are in column 1
    feat_mat <- data_mat[data_mat[, 1] == ur_id, (ur_ix + 2):(log_prob_ix - 2),
                         drop=FALSE]

    # Get probability matrix for candidates that share the same UR
    prob_mat <- data_mat[data_mat[, 1] == ur_id, log_prob_ix,
                         drop=FALSE]
    # Convert from log probs to probs
    prob_mat <- exp(prob_mat)

    # Calculate row-matrix for expected "feature activation" for i-th candidate
    # Store it in i-th row of expected "feature activation" matrix
    exp_mat[i, ] <- t(prob_mat) %*% feat_mat
  }

  # Return matrix of expected "feature activation"
  # Dimensions of matrix: (d, n)
  # Where d=num of candidates & n=num of constraints
  return(exp_mat)
}

# Calculate gradients for bias term
calculate_grad_bias <- function(bias_params, constraint_weights) {
  top <- constraint_weights - bias_params$Mu
  bottom <- bias_params$Sigma^2
  grad_bias <- top / bottom

  # Returns a n-length vector, where n=num of constraints
  return(grad_bias)
}

# Helper function that applies Z normalization
normalize_row <- function(row, m, col_num) {
  return(row[col_num] / sum(m[m[, 1] == row[1],][, col_num]))
}

# Calculates the bias term in the optimized function.
calculate_bias <- function(bias_params, constraint_weights) {
  top <- (constraint_weights - bias_params$Mu)^2
  bottom <- 2 * bias_params$Sigma^2
  bias <- sum(top / bottom)

  return(bias)
}

# Helper function that checks whether any arguments are NA
any_not_na <- function(...) {
  x <- list(...)
  return(any(!is.na(x)))
}

# Function to load and validate bias parameters.
process_bias_arguments <- function(names, bias_input = NA,
                                   mu_scalar = NA, mu_vector = NA,
                                   sigma_scalar = NA, sigma_vector=NA,
                                   num_constraints = NA) {
  if (any_not_na(bias_input)) {
    # Read bias parameters from provided file location
    if (any_not_na(mu_scalar, mu_vector, sigma_vector, sigma_vector)) {
      warning(
        "Both a bias file and bias scalars/vectors were provided\n",
        "Ignoring scalars/vectors and using parameters from file"
      )
    }
    if (is.data.frame(bias_input)) {
      bias_params <- bias_input
    }
    else {
      bias_params <- load_bias_file_otsoft(bias_input)
    }
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
    bias_params <- cbind(names, bias_params)
    names(bias_params) <- c("Constraint", "Mu", "Sigma")

  } else if (any_not_na(mu_scalar, mu_vector, sigma_scalar, sigma_vector)) {
    stop("You must specify both constraint mus and sigmas, or neither.")
  } else {
    # Don't use a bias parameter
    bias_params <- NA
  }

  # Check that bias params fall within acceptable ranges
  if (any_not_na(bias_params)) {
    if (!all(bias_params$Mu >= 0)) {
      stop("All constraint mus must be >= 0")
    }
    if (!all(bias_params$Sigma > 0)) {
      stop("All constraint sigmas must be > 0")
    }
  }
  return(bias_params)
}
