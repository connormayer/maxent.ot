
#' Cross-validate bias parameters for constraint weights.
#'
#' Performs k-fold cross-validation of a data set and a set of input bias
#' parameters. Cross-validation allows the space of bias parameters to be
#' searched to find the settings that best support generalization to unseen data.
#'
#' The cross-validation procedure is as follows:
#'
#' \enumerate{
#'   \item Randomly divide the data into k partitions.
#'   \item Iterate through every combination of mu and sigma specified in the
#'     input arguments (see the documentation for the `grid_search` argument
#'     for details on how this is done).
#'   \item For each combination, for each of the k partitions, train a model
#'     on the other (k-1) partitions using `optimize_weights` and then run
#'     `predict_probabilities` on the remaining partition.
#'   \item Record the mean log likelihood the models apply to the held-out
#'     partitions.
#' }
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
#' @param k The number of folds to use in cross-validation.
#' @param mu_values A vector or list of mu bias parameters to use in
#'   cross-validation. Parameters may either be scalars, in which case the
#'   same mu parameter will be applied to every constraint, or vectors/lists
#'   containing a separate mu bias parameter for each constraint.
#' @param sigma_values A vector or list of sigma bias parameters to use in
#'   cross-validation. Parameters may either be scalars, in which case the
#'   same sigma parameter will be applied to every constraint, or vectors/lists
#'   containing a separate sigma bias parameter for each constraint.
#' @param grid_search (optional) If TRUE, the Cartesian product of the values
#'   in `mu_values` and `sigma_values` will be validated. For example, if
#'   `mu_values = c(0, 1)` and `sigma_values = c(0.1, 1)`, cross-validation will
#'   be done on the mu/sigma pairs `(0, 0.1), (0, 1), (1, 0.1), (1, 1)`. If
#'   FALSE (default), cross-validation will be done on each pair of values at
#'   the same indices in `mu_values` and `sigma_values`. For example, if
#'   `mu_values = c(0, 1)` and `sigma_values = c(0.1, 1)`, cross-validation will
#'   be done on the mu/sigma pairs `(0, 0.1), (1, 1)`.
#' @param output_path (optional) A string specifying the path to a file to
#'   which the cross-validation results will be saved. If the file exists it
#'   will be overwritten. If this argument isn't provided, the output will not
#'   be written to a file.
#' @param out_sep (optional) The delimiter used in the output files.
#'   Defaults to tabs.
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
#' @param model_name (optional) A name for the model. If not provided, the file
#'   name will be used if the input is a file path. If the input is a data frame
#'   the name of the variable will be used.
#' @param allow_negative_weights (optional) Whether the optimizer should allow
#'   negative weights. Defaults to FALSE.
#' @return A data frame with the following columns:
#' \itemize{
#'         \item `model_name`: the name of the model
#'         \item `mu`: the value(s) of mu tested
#'         \item `sigma`: the value(s) of sigma tested
#'         \item `folds`: the number of folds
#'         \item `mean_ll`: the mean log likelihood of k-fold cross-validation
#'                          using these bias parameters
#' }
#' @examples
#'   # Get paths to OTSoft file. Note that you can also pass dataframes into
#'   # this function, as described in the documentation for `optimize`.
#'   data_file <- system.file(
#'       "extdata", "amp_demo_grammar.csv", package = "maxent.ot"
#'   )
#'   tableaux_df <- read.csv(data_file)
#'
#'   # Define mu and sigma parameters to try
#'   mus <- c(0, 1)
#'   sigmas <- c(0.01, 0.1)
#'
#'   # Do 2-fold cross-validation
#'   cross_validate(tableaux_df, 2, mus, sigmas)
#'
#'   # Do 2-fold cross-validation with grid search of parameters
#'   cross_validate(tableaux_df, 2, mus, sigmas, grid_search=TRUE)
#'
#'   # You can also use vectors/lists for some/all of the bias parameters to set
#'   # separate biases for each constraint
#'   mus_v <- list(
#'     c(0, 1),
#'     c(1, 0)
#'   )
#'   sigmas_v <- list(
#'     c(0.01, 0.1),
#'     c(0.1, 0.01)
#'   )
#'
#'   cross_validate(tableaux_df, 2, mus_v, sigmas_v)
#'
#'   # Save cross-validation results to a file
#'   tmp_output <- tempfile()
#'   cross_validate(tableaux_df, 2, mus, sigmas, output_path=tmp_output)
#' @export
cross_validate <- function(input, k, mu_values, sigma_values,
                           grid_search = FALSE, output_path = NA,
                           out_sep = ',', control_params = NA,
                           upper_bound = DEFAULT_UPPER_BOUND,
                           encoding = 'unknown', model_name = NA,
                           allow_negative_weights = FALSE) {
  if (is.data.frame(input)) {
    if (is.na(model_name)) {
      # If no provided model name, use name of input variable
      model_name <- toString(substitute(input))
    }
  }
  processed_input <- load_input(
    input, encoding = encoding, model_name = model_name
  )
  data <- processed_input$data
  model_name <- processed_input$model_name
  partitions <- partition_data(data, k)
  result_df <- data.frame()

  if (grid_search) {
    for (mu in mu_values) {
      for (sigma in sigma_values) {
        result_row <- do_validation(
          k, data, partitions, mu[[1]], sigma[[1]], model_name, control_params,
          upper_bound, allow_negative_weights
        )
        result_df <- rbind(result_df, result_row)
      }
    }
  } else {
    if (length(mu_values) != length(sigma_values)) {
      stop("mu_values and sigma_values must be the same length unless
            gridsearch = TRUE.")
    }

    for (i in (1:length(mu_values))) {
      result_row <- do_validation(
        k, data, partitions, mu_values[[i]], sigma_values[[i]], model_name,
        control_params, upper_bound, allow_negative_weights
      )
      result_df <- rbind(result_df, result_row)
    }
  }

  if (!is.na(output_path)) {
    utils::write.table(
      result_df, file=output_path, sep=out_sep, row.names = FALSE
    )
  }
  return(result_df)
}


do_validation <- function(k, data, partitions, mu, sigma, model_name,
                          control_params, upper_bound, allow_negative_weights) {
  # Performs cross-validation given a particular value of mu and sigma
  log_liks_test <- c()
  log_liks_training <- c()

  for (hold_out in (1:k)) {
    print(sprintf(
      "Training paramers mu:%s, sigma:%s, fold number:%s", toString(mu),
      toString(sigma), hold_out
    ))
    training_data <- data
    training_data$Frequency <- 0
    test_data <- training_data

    training_part <- partitions[partitions$partition != hold_out,]
    test_part <- partitions[partitions$partition == hold_out,]

    training_tableau <- populate_tableau(training_data, training_part)
    test_tableau <- populate_tableau(test_data, test_part)

    m <- optimize_weights(
      training_tableau, mu = mu, sigma = sigma, control_params = control_params,
      upper_bound = upper_bound, allow_negative_weights = allow_negative_weights
    )
    predictions_test <- predict_probabilities(test_tableau, m$weights)
    log_liks_test <- c(predictions_test$loglik, log_liks_test)

    predictions_training <- predict_probabilities(training_tableau, m$weights)
    log_liks_training <-c(predictions_training$loglik, log_liks_training)
  }
  mean_ll_test <- mean(log_liks_test)
  mean_ll_training <- mean(log_liks_training)

  df <- data.frame(
    model_name = model_name, mu = toString(mu), sigma = toString(sigma),
    folds = k, mean_ll_test = mean_ll_test, mean_ll_training = mean_ll_training
  )
  return(df)
}

#A function to turn input-output pairs back into tableaux
populate_tableau <- function(tableau, tokens) {
  # Fills in a tableau with token frequency counts
  tableau$Frequency <- 0 #First, initialize all the frequencies to 0
  tokens$count <- seq(nrow(tokens))
  counts <- stats::aggregate(count ~ Input + Output, data = tokens, FUN = length)
  for (i in (1:nrow(counts))) {
    row <- counts[i, ]
    tableau[tableau$Input == row$Input & tableau$Output == row$Output, ]$Frequency <- row$count
  }
  return(tableau)
}

partition_data <- function(data, k) {
  # Divides data into k roughly equal partitions
  freq <- data$Frequency
  data_vals <- data[, 1:2]
  n <- sum(freq)

  if (n < k) {
    stop(sprintf(
      "Cannot divide %s data points into %s partitions: please choose a
       a smaller value of k.", n, k
    ))
  }

  # Create one row per input token and shuffle them
  randomized_data <- data_vals[sample(rep(seq_len(nrow(data_vals)), freq)),]
  n <- nrow(randomized_data)
  # Assign them partition numbers
  randomized_data$partition <- rep(seq_len(k), ceiling(n / k))[1:n]

  return(randomized_data)
}
