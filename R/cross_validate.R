
#' Cross-validate bias parameters for constraint weights.
#'
#' Performs k-fold cross validation. TODO
#'
#' @param input The path to the input data file or a dataframe/tibble.
#'   This should contain more OT tableaux consisting of
#'   mappings between underlying and surface forms with observed frequency and
#'   violation profiles. Constraint violations must be numeric.
#'
#'   If this is a file path, the file should be in OTSoft format.
#'   For examples of OTSoft format, see inst/extdata/sample_data_file.txt.
#'   For an example of the dataframe format, see inst/extdata/sample_data_frame.txt.
#'   You can read the latter file into a dataframe using read.csv or into a tibble
#'   using dplyr::read_csv.
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
#'   be done on the mu/sigma pairs `c(0, 0.1), c(0, 1), c(1, 0.1), c(1, 1)`. If
#'   FALSE (default), cross-validation will be done on each pair of values at
#'   the same indices in `mu_values` and `sigma_values`. For example, if
#'   `mu_values = c(0, 1)` and `sigma_values = c(0.1, 1)`, cross-validation will
#'   be done on the mu/sigma pairs `c(0, 0.1), c(1, 1)`.
#' @param output_path (optional) A string specifying the path to a file to
#'   which the cross-validation results will be saved. If the file exists it
#'   will be overwritten. If this argument isn't provided, the output will not
#'   be written to a file.
#' @param in_sep (optional) The delimiter used in the input file.
#'   Defaults to tabs.
#' @param out_sep (optional) The delimiter used in the output files.
#'   Defaults to tabs.
#' @param encoding (optional) The character encoding of the input file. Defaults
#'  to "unknown".
#' @param model_name (optional) A name for the model. If not provided, the file
#'   name will be used if the input is a file path. If the input is a data frame
#'   the name of the variable will be used.
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
#'       "extdata", "sample_hu_wug.txt", package = "maxent.ot"
#'   )
#'
#'   # Define mu and sigma parameters to try
#'   mus <- c(0, 1)
#'   sigmas <- c(0.01, 0.1)
#'
#'   # Do 10-fold cross-validation
#'   cross_validate(data_file, 10, mus, sigmas)
#'
#'   # Do 10-fold cross-validation with grid search of parameters
#'   cross_validate(data_file, 10, mus, sigmas, grid_search=TRUE)
#'
#'   # You can also use vectors for some/all of the bias parameters.
#' mus_v <- list(
#'   c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
#'   c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
#' )
#' sigmas_v <- list(
#'   c(0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1),
#'   c(0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01)
#' )
#'
#' cross_validate(data, 10, mus_v, sigmas_v)
#'
#'   # Save cross-validation results to a file
#'   tmp_output <- tempfile()
#'   cross_validate(data_file, 10, mus, sigmas, output_path=tmp_output)
#' @export
cross_validate <- function(input, k, mu_values, sigma_values,
                           grid_search = FALSE, output_path = NA,
                           in_sep = '\t', out_sep = ',', control_params = NA,
                           upper_bound = DEFAULT_UPPER_BOUND,
                           encoding = 'unknown', model_name = NA) {
  if (is.data.frame(input)) {
    if (is.na(model_name)) {
      # If no provided model name, use name of input variable
      model_name <- toString(substitute(input))
    }
  }
  processed_input <- load_input(
    input, sep = in_sep, encoding = encoding, model_name = model_name
  )
  data <- processed_input$data
  model_name <- processed_input$model_name
  partitions <- partition_data(data, k)
  result_df <- data.frame()

  if (grid_search) {
    for (mu in mu_values) {
      for (sigma in sigma_values) {
        result_row <- do_validation(
          k, data, partitions, mu, sigma, model_name, in_sep, control_params,
          upper_bound
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
        k, data, partitions, mu_values[i], sigma_values[i], model_name, in_sep,
        control_params, upper_bound
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

do_validation <- function(k, data, partitions, mu, sigma, model_name, in_sep,
                          control_params, upper_bound) {
  if (length(mu) > 1) {
    mu_scalar <- NA
    mu_vector <- unlist(mu)
  } else {
    mu_scalar <- mu
    mu_vector <- NA
  }

  if (length(sigma) > 1) {
    sigma_scalar <- NA
    sigma_vector <- unlist(sigma)
  } else {
    sigma_scalar <- sigma
    sigma_vector <- NA
  }

  log_liks <- c()

  for (hold_out in (1:k)) {
    sprintf(
      "Training paramers mu:%s, sigma:%s, fold number:%s", mu, sigma, hold_out
    )
    training_data <- data
    training_data$Frequency <- 0
    test_data <- training_data

    training_part <- partitions[partitions$partition != hold_out,]
    test_part <- partitions[partitions$partition == hold_out,]

    training_tableau <- populate_tableau(training_data, training_part)
    test_tableau <- populate_tableau(test_data, test_part)

    m <- optimize_weights(
      training_tableau, mu_scalar = mu_scalar, mu_vector = mu_vector,
      sigma_scalar = sigma_scalar, sigma_vector = sigma_vector,
      in_sep = in_sep, control_params = control_params, upper_bound = upper_bound
    )
    predictions <- predict_probabilities(test_tableau, m$weights)
    log_liks <- c(predictions$loglik, log_liks)
  }
  mean_ll <- mean(log_liks)
  df <- data.frame(
    model_name = model_name, mu = toString(mu), sigma = toString(sigma),
    folds = k, mean_ll = mean_ll
  )
  return(df)
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
