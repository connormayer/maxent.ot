cross_validate <- function(input, k, bias_file = NA,
                           mu_scalar = NA, mu_vector = NA,
                           sigma_scalar = NA, sigma_vector = NA,
                           in_sep = '\t', control_params = NA,
                           upper_bound = DEFAULT_UPPER_BOUND,
                           encoding = 'unknown', model_name = NA) {
  processed_input <- load_input(
    input, sep = in_sep, encoding = encoding, model_name = model_name
  )
  data <- processed_input$data
  partitions <- partition_data(data, k)

  mu_scalar <- 0
  sigma_scalar <- 1

  log_liks <- c()

  for (hold_out in (1:k)) {
    training_data <- data
    training_data$Frequency <- 0
    test_data <- training_data

    training_part <- partitions[partitions$partition != hold_out,]
    test_part <- partitions[partitions$partition == hold_out,]

    training_tableau <- populate_tableau(training_data, training_part)
    test_tableau <- populate_tableau(test_data, test_part)

    m <- optimize_weights(training_tableau, mu_scalar = mu_scalar, sigma_scalar = sigma_scalar)
    predictions <- predict_probabilities(test_tableau, m$weights)
    log_liks <- c(predictions$loglik, log_liks)
  }

  mean_ll <- mean(log_liks)

  print('hi')

}

populate_tableau <- function(tableau, tokens) {
  tokens$count <- seq(nrow(tokens))
  counts <- aggregate(count ~ Input + Output, data = tokens, FUN = length)
  for (i in (1:nrow(counts))) {
    row <- counts[i, ]
    tableau[tableau$Input == row$Input & tableau$Output == row$Output, ]$Frequency <- row$count
  }
  return(tableau)
}

partition_data <- function(data, k) {
  freq <- data$Frequency
  data_vals <- data[, 1:2]

  # Create one row per input token and shuffle them
  randomized_data <- data_vals[sample(rep(seq_len(nrow(data_vals)), freq)),]
  n <- nrow(randomized_data)
  # Assign them partition numbers
  randomized_data$partition <- rep(seq_len(k), ceiling(n / k))[1:n]

  return(randomized_data)
}
