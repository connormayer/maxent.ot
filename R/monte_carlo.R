# Constants


#' Describe function here
#'



cdnProb_trial <- function (prob_file, response_file, input_format = "otsoft",
                           in_sep = "\t", output_path = NA, out_sep = "\t") {
  # Prepare response_file
  input <- load_data_otsoft(response_file, sep = in_sep)
  long_names <- input$full_names
  short_names <- input$abbr_names
  data <- input$data
  n <- input$n

  data_matrix <- matrix(0L, nrow = nrow(data), ncol = ncol(data) +
                          2)
  data_matrix[, 1] <- as.integer(as.factor(data[, 1]))
  data_matrix[, 4:(ncol(data_matrix) - 2)] <- apply(as.matrix(data[,
                                                                   4:ncol(data)]), 2, as.numeric)
  data_matrix[is.na(data_matrix)] <- 0


  # Record trial_id in second-last column of output
  trial_id <- 0
  curr_ur <- ""
  sr_ls <- list()

  for (i in 1:nrow(prob_file)) {

    # If still the same trial
    # Record current trial_id in output 2nd last column of output
    # Either new UR
    if (prob_file[i, 1] == curr_ur && (prob_file[i, 2] %in% sr_ls)) {
      data_matrix[i, (ncol(data_matrix)-1)] <- trial_id
      #data_matrix[i, 1] <- trial_id    # Put trial_id in 1st column if using normalize row

      # This row begins a new trial
    } else {
      # Current trial needs a new trial_id
      trial_id = trial_id + 1
      # Record current trial_id in output file
      data_matrix[i, (ncol(data_matrix)-1)] <- trial_id
      #data_matrix[i, 1] <- trial_id    # Put trial_id in 1st column if using normalize row

      # Store new UR of current trial
      curr_ur <- prob_file[i, 1]
      # Start a new list of SRs for current trial
      sr_ls <- list(prob_file[i, 2])

      # Track following rows that belong to the same trial
      next_row <- i+1
      # Such rows have the same UR and different SRs
      while (next_row <= nrow(data) && prob_file[next_row, 1] == curr_ur && !(prob_file[next_row, 2] %in% sr_ls)) {
        # Add alternative SR options to SR list
        sr_ls <- append(sr_ls, prob_file[next_row, 2])
      }
    }
  }

  # Add conditional prob over UR in a new last column
  data_matrix <- cbind(data_matrix, prob_file[, "Predicted Probability"])

  # Record probability conditioned on trial in "new" 2nd last column of output
  data_matrix[, (ncol(data_matrix)-1)] <- apply(data_matrix, 1,
                                                normalize_trial, data_matrix, ncol(data_matrix))

  # Change observed frequencies to 0
  data_matrix[, 3] <- 0

  # Put row and column names back in
  output <- cbind(data[, 1:2], data_matrix[, 3:ncol(data_matrix)])
  names(output) <- c(c(c("UR", "SR", "Freq"), unlist(long_names)),
                     "Trial id", "Conditional prob over trial", "Conditional prob over UR")
  # Write to output file if desired
  if (!is.na(output_path)) {
    utils::write.table(output, file = output_path, sep = out_sep,
                       row.names = FALSE)
  }

  return(output)
}

# Helper function that applies normalization over trial
normalize_trial <- function (row, m, col_num) {
  return(row[col_num]/sum(m[m[, (ncol(m)-2)] == row[(ncol(m)-2)], ][, col_num]))
}

# Function that creates 1 simulated response based on probabilities conditioned over trial
# The input_file is the one that has probs over trial
monte_carlo <- function (data_file) {

  # Get total number of trials
  num_trial <- max(data_file[,ncol(data_file)-2])

  #
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

# TODO: biggest function that calls on smaller helper fns

# Create file that conditions over trial
# prob_per_trial_file <- cdnProb_trial(prob_file, response_file)

# Start looping below
# simulated_resp_file <- monte_carlo()


# Learn weights
# output <- optimize_weights_fromDF <- optimize_weights_fromDF(simulated_resp_file)
