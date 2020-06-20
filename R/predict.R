#' Predict probabilities of candidates based on their violation profiles and
#' constraint weights.
#'
#' @param test_file The path to the input data file. This file contains one or
#'   more OT tableaux consisting of mappings between underlying and surface
#'   forms with observed frequency and violation profiles. Constraint
#'   violations must be numeric.
#' @param constraint_weights The constraint weights to use.
#' @param input_format (optional) A string specifying the format of the input
#'   files. Currently only otsoft-style formatting is supported. Note that this
#'   means a column containing frequencies must be present. However, these
#'   frequencies will be ignored.
#'
#' @return A data table containing all the tableaux, with probabilities
#'   assigned to each candidate.
#'
#' @export
predict_probabilities <- function(test_file, constraint_weights,
                                  input_format = "otsoft") {

  if (input_format == "otsoft") {
    input <- load_data_otsoft(test_file)
    long_names <- input[[1]]
    short_names <- input[[2]]
    tableaux <- input[[3]]
  } else {
    stop(sprintf("Invalid input format %s", input_format))
  }

  if (any(constraint_weights < 0)) {
    stop("Constraint weights must be non-negative")
  }

  output <- data.table::data.table()

  for (tableau in tableaux) {
    violations <- tableau[, 4:ncol(tableau)]
    log_probs <- calculate_tableau_probabilities(constraint_weights, violations)
    probs <- exp(log_probs)
    tableau <- cbind(tableau, probs)
    output <- rbind(output, tableau)
  }

  names(output) <- c(c(c("UR", "SR", "Freq"), unlist(long_names)), "prob")

  return(output)
}
