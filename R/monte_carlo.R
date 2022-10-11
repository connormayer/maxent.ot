# Constants
DEFAULT_UPPER_BOUND <- 1000

#' Create simulated data and learn weights for these data
#'
#' Creates a simulated data set by picking an output for each instance of an
#' input.
#' The probability of picking a particular output is guided by its conditional
#' probability given the input.
#' Learns constraint weights for each simulated data set.
#'
#' This function creates multiple simulated data sets, and learns a set of
#' weights that maximizes the likelihood of data for each simulated data set.
#'
#' To create a simulated data set, one output is randomly chosen for each
#' instance of an input.
#' The probability of picking a particular output, \eqn{O_i}, which arises from
#' input \eqn{I_j} depends on \eqn{Pr(O_i|I_j)}.
#'
#' The function `optimize_weights()` is called to find a set of weights that
#' maximize the likelihood of the simulated data.
#' All optional arguments of `optimize_weights()` that were available for the
#' user to specify biases and bounds are likewise available in this function,
#' `monte_carlo_weights()`.
#'
#' The process of simulating a data set and learning weights that optimize the
#' likelihood of the simulated data is repeated as per the number of specified
#' simulations.
#'
#' @section Why use this function?:
#'
#' This function gives us a way to estimate constraint weights via a Monte Carlo
#' process.
#' For example we might be interested in the effect of temperature on polarizing
#' predicted probabilities, and the resulting constraint weights.
#' This function can produce a distribution of constraint weights for the
#' simulated polarized data, as well as a distribution of constraint weights for
#' the simulated non-polarized ones, thereby allowing a comparison of the two.
#'
#' @param pred_prob A data frame with a column for predicted probabilities.
#'   This object should be in the same format as the object returned by the
#'   `predict_probabilities` function.
#' @param num_simul The number of simulations to run.
#' @param bias_file (optional) The path to the file containing mus and sigma
#'   for constraint biases. If this argument is provided, the scalar and vector
#'   mu and sigma arguments will be ignored. Each row in this file should be the
#'   name of the constraint, followed by the mu, followed by the sigma
#'   (separated by whatever the relevant separator is; default is tabs).
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
#' @param penalty_func (optional) ???
#' @param output_path (optional) A string specifying the path to a file to
#'   which the output will be saved. If the file exists it will be overwritten.
#'   If this argument isn't provided, the output will not be written to a file.
#' @param out_sep (optional) The delimiter used in the output files.
#'   Defaults to tabs.
#' @param control_params (optional) A named list of control parameters that
#'   will be passed to the \link[stats]{optim} function. See the documentation
#'   of that function for details. Note that some parameter settings may
#'   interfere with optimization. The parameter `fnscale` will be overwritten
#'   with `-1` if specified, since this must be treated as a maximization
#'   problem.
#' @param upper_bound (optional) The maximum value for constraint weights.
#'   Defaults to 1000.
#' @return A data frame with the following structure:
#' \itemize{
#'         \item rows: As many rows as the number of simulations
#'         \item columns: As many columns as the number of constraints
#' }
#' @examples
#'   # Get paths to toy data file
#'   data_file <- system.file(
#'       "extdata", "sample_data_file.txt", package = "maxent.ot"
#'   )
#'
#'   # Fit weights to data with no biases
#'   fit_model <- optimize_weights(data_file)
#'
#'   # Predict probabilities for the same input with temperature = 2
#'   pred_df <- predict_probabilities(
#'       data_file, fit_model$weights, temperature = 2
#'   )
#'
#'  # Run 5 monte carlo simulations
#'  # based on predicted probabilities when temperature = 2,
#'  # and learn weights for these 5 simulated data sets
#'  monte_carlo_weights(pred_df, 5)
#'
#'  # Save learned weights to a file
#'  tmp_output <- tempfile()
#'  monte_carlo_weights(pred_df, 5, output_path=tmp_output)
#' @export
# Learns constraint weights for multiple randomly generated SR responses
monte_carlo_weights <- function(pred_prob, num_simul,
                                bias_file = NA,
                                mu_scalar = NA, mu_vector = NA,
                                sigma_scalar = NA, sigma_vector = NA,
                                penalty_func = NA,
                                output_path = NA, out_sep = "\t",
                                control_params = NA,
                                upper_bound = DEFAULT_UPPER_BOUND) {

  # Create file that calculates conditional probability over trial
  cdnpred_prob <- cdnProb_trial(pred_prob)

  # Initialize data frame to store learned weights
  num_feats <- ncol(cdnpred_prob)-6
  output <- matrix(nrow = num_simul, ncol = num_feats)

  # Learn weights for each simulation
  for (i in 1:num_simul) {

    # Create a simulated response file
    simul_resp_file <- monte_carlo(cdnpred_prob)

    # Learn weights for simulated response
    curr_model <- optimize_weights(simul_resp_file,
                                   bias_file = NA,
                                   mu_scalar = NA, mu_vector = NA,
                                   sigma_scalar = NA, sigma_vector = NA,
                                   penalty_func = NA,
                                   input_format = 'df',
                                   control_params = NA,
                                   upper_bound = DEFAULT_UPPER_BOUND,
                                   model_name = "curr_model")

    # Record learned weights
    output[i,] <- curr_model$weights
  }

  # Put in names of constraints
  colnames(output) <- colnames(simul_resp_file)[4:ncol(simul_resp_file)]

  # Write to output file if desired
  if (!is.na(output_path)) {
    utils::write.table(output, file = output_path, sep = out_sep,
                       row.names = FALSE)
  }

  return(output)
}

# Function that calculates conditional probability of output given trial
cdnProb_trial <- function (pred_prob) {

  # Build ourselves a matrix for efficient computation
  # Pre-allocate space
  data_matrix <- matrix(0L, nrow = nrow(pred_prob), ncol = ncol(pred_prob))
  # Port over the violation profiles and p(SR|UR)
  data_matrix[, 1:(ncol(data_matrix) - 2)] <- data.matrix(pred_prob[, 1:(ncol(pred_prob) - 2)])
  # Replace empty cells with 0
  data_matrix[is.na(data_matrix)] <- 0

  # Record trial_id in last column of data_matrix
  trial_id <- 0
  curr_ur <- ""
  sr_ls <- list()

  # For each row
  for (i in 1:nrow(pred_prob)) {

    # If still the same trial as previous row
    # Conditions: same UR, SR is is the list of SRs associated with this UR
    if (pred_prob[i, 1] == curr_ur && (pred_prob[i, 2] %in% sr_ls)) {
      # Record current trial_id
      data_matrix[i, ncol(data_matrix)] <- trial_id

      # Else: This row begins a new trial
    } else {
      # Current trial needs a new trial_id
      trial_id = trial_id + 1
      # Record current trial_id in output file
      data_matrix[i, ncol(data_matrix)] <- trial_id

      # Store new UR of current trial in curr_ur
      curr_ur <- pred_prob[i, 1]
      # Start a new list of SRs for current trial in sr_ls
      sr_ls <- list(pred_prob[i, 2])

      # Track following rows that belong to the same trial
      next_row <- i+1
      # While the following rows belong to the same trial as row i
      # Conditions: such rows have the same UR and different SRs
      while (next_row <= nrow(pred_prob) && pred_prob[next_row, 1] == curr_ur && !(pred_prob[next_row, 2] %in% sr_ls)) {
        # Add alternative SR options to SR list
        sr_ls <- append(sr_ls, pred_prob[next_row, 2])
        # Move on to following row
        next_row <- next_row+1
      }
    }
  }

  # Record probability conditioned on trial in "new" 2nd last column of output
  data_matrix[, (ncol(data_matrix)-1)] <- apply(data_matrix, 1,
                                                normalize_trial, data_matrix, (ncol(data_matrix)-2))

  # Change observed frequencies to 0
  data_matrix[, 3] <- 0

  # Re-introduce UR & SR characters
  output <- cbind(pred_prob[, 1:2], data_matrix[, 3:ncol(data_matrix)])

  # Column names
  # Port over column names from pred_prob
  names(output) <- colnames(pred_prob)
  # Change "Predicted Probability" to "Pred p(SR|UR)"
  colnames(output)[colnames(output) == "Predicted Probability"] <- "Pred p(SR|UR)"
  # Change "Observed Probability" to "Pred p(SR|trial)"
  colnames(output)[colnames(output) == "Observed Probability"] <- "Pred p(SR|trial)"
  # Change "Error" to "Trial id"
  colnames(output)[colnames(output) == "Error"] <- "Trial id"

  return(output)
}

# Helper function that applies normalization over trial
normalize_trial <- function (row, m, col_num) {
  return(row[col_num]/sum(m[m[, ncol(m)] == row[ncol(m)], ][, col_num]))
}

# Function that creates 1 simulated response based on probabilities conditioned over trial
# The input_file is the one that has probs over trial
monte_carlo <- function (data_file) {

  # Get total number of trials
  num_trial <- max(data_file[,ncol(data_file)])

  # For each trial, pick 1 SR response
  for (i in 1:num_trial) {
    # Get row indices of all trials with 'Trial id' == i
    id_vec <- which(data_file$'Trial id' == i)

    # Get conditional prob over trial for all trials in id_vec
    prob_vec <- data_file[id_vec, 'Pred p(SR|trial)']

    # Randomly pick 1 SR according to conditional prob
    # Currently only 1 random draw per trial
    # To support multiple draws per trial in future?
    picked_sr_id <- sample(x=id_vec, size=1, prob=prob_vec)

    # Record picked SR
    data_file[picked_sr_id, 'Freq'] <- 1
  }

  # Drop last 3 columns so file has only the columns required for optimize_weights()
  simulated_resp_file <- data_file[, -(ncol(data_file)-2):-(ncol(data_file))]

  return(simulated_resp_file)
}
