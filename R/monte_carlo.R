# Constants


#' Describe function here
#'

# TODO: get optimize_weights() to work with an input file that is a matrix (here)
# in addition to also working with .txt files (only file that is currently supported)
# nb: currently i'm using optimize_weights_fromDF that works with matrix input-type.

cdnProb_trial <- function (prob_file) {

  # Prepare response_file: data_matrix
  # Drop last two columns (observed prob & error)
  data_matrix <- prob_file[, -(ncol(prob_file)-1):-(ncol(prob_file))]

  # Insert a new column for Conditional prob over trial
  data_matrix["Pred p(SR|trial)"] <- 0
  # Insert a new column for Trial id
  data_matrix["Trial id"] <- 0

  # Record trial_id in last column of data_matrix
  trial_id <- 0
  curr_ur <- ""
  sr_ls <- list()

  for (i in 1:nrow(prob_file)) {

    # If still the same trial
    # Conditions: same UR, SR is is the list of SRs associated with this UR
    if (prob_file[i, 1] == curr_ur && (prob_file[i, 2] %in% sr_ls)) {
      # Record current trial_id in output last column of output
      data_matrix[i, ncol(data_matrix)] <- trial_id

      # Else: This row begins a new trial
      # Current trial needs a new trial_id
    } else {
      trial_id = trial_id + 1
      # Record current trial_id in output file
      data_matrix[i, ncol(data_matrix)] <- trial_id

      # Store new UR of current trial in curr_ur
      curr_ur <- prob_file[i, 1]
      # Start a new list of SRs for current trial in sr_ls
      sr_ls <- list(prob_file[i, 2])

      # Track following rows that belong to the same trial
      next_row <- i+1
      # Conditions: such rows have the same UR and different SRs
      while (next_row <= nrow(prob_file) && prob_file[next_row, 1] == curr_ur && !(prob_file[next_row, 2] %in% sr_ls)) {
        # Add alternative SR options to SR list
        sr_ls <- append(sr_ls, prob_file[next_row, 2])
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

  # Change name of "Predicted Probability" to "Pred p(SR|UR)"
  colnames(data_matrix)[colnames(data_matrix) == "Predicted Probability"] <- "Pred p(SR|UR)"

  return(data_matrix)
}

# Helper function that applies normalization over trial
normalize_trial <- function (row, m, col_num) {
  return(as.numeric(row[col_num])/sum(m[m[, ncol(m)] == row[ncol(m)], ][, col_num]))
}

# Function that creates 1 simulated response based on probabilities conditioned over trial
# The input_file is the one that has probs over trial
monte_carlo <- function (data_file) {

  # Get total number of trials
  num_trial <- max(data_file[,ncol(data_file)-2])

  # For each trial, pick 1 SR response
  for (i in 1:num_trial) {
    # Get row indices of all trials with 'Trial id' == i
    id_vec <- which(data_file$'Trial id' == i)

    # Get conditional prob over trial for all trials in id_vec
    prob_vec <- data_file[id_vec, 'Conditional prob over trial']

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

# Learns constraint weights for multiple randomly generated SR responses
# TODO: support writing output to a .txt file?
monte_carlo_weights <- function(prob_file, num_simul, output_path = NA, out_sep = "\t") {

  # Create file that calculates conditional probability over trial
  cdnProb_file <- cdnProb_trial(prob_file)

  # Initialize data frame to store learned weights
  num_feats <- ncol(cdnProb_file)-6
  output <- matrix(nrow = num_simul, ncol = num_feats)

  # Learn weights for each simulation
  for (i in 1:num_simul) {

    # Create a simulated response file
    simul_resp_file <- monte_carlo(cdnProb_file)

    # Learn weights for simulated response
    # TODO: input_format something other than "otsoft"? & model name
    # TODO: learn weights from "optimize_weights" rather than "optimize_weights_fromDF" once input format is fixed
    curr_model <- optimize_weights_fromDF(simul_resp_file, model_name = "curr_model")

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
