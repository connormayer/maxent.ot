#' Predict probabilities of OT candidates
#'
#' Predict probabilities of candidates based on their violation profiles and
#' constraint weights. For each input/output pair in the provided, file this
#' will calculate the probability of that output given the input form and the
#' provided weights. This is defined as
#'
#' \deqn{P(y|x; w) = \frac{1}{Z_w(x)}\exp(\sum_{k=1}^{m}{w_k f_k(y, x)})}
#'
#' where \eqn{f_k(y, x)} is the number of violatons of constraint \eqn{k}
#' incurred by mapping underlying \eqn{x} to surface \eqn{y}, \eqn{w_k} is the
#' weight associated with contraint \eqn{k}, and  \eqn{Z(x)} is a normalization
#' term defined as
#'
#' \deqn{Z(x) = \sum_{y\in \mathcal{Y}(x)}{\exp(\sum_{k=1}^{m}{w_k f_k(y, x)})}}
#'
#' The resulting probabilities will be appended to a data frame object
#' representing the input tableaux. This data frame can also be saved to a file
#' if the `output_path` argument is provided.
#'
#' @param test_file The path to the input data file. This file contains one or
#'   more OT tableaux consisting of mappings between underlying and surface
#'   forms with observed frequency and violation profiles. Constraint
#'   violations must be numeric.
#' @param constraint_weights The constraint weights to use.
#' @param output_path (optional) A string specifying the path to a file to
#'   which the output will be saved. If the file exists it will be overwritten.
#'   If this argument isn't provided, the output will not be written to a file.
#' @param input_format (optional) A string specifying the format of the input
#'   files. Currently only otsoft-style formatting is supported. Note that this
#'   means a column containing frequencies must be present. However, these
#'   frequencies will be ignored.
#' @param in_sep (optional) The delimiter used in the input file.
#'   Defaults to tabs.
#' @param out_sep (optional) The delimiter used in the output files.
#'   Defaults to tabs.
#'
#' @return A data table containing all the tableaux, with probabilities
#'   assigned to each candidate.
#'
#' @export
predict_probabilities <- function(test_file, constraint_weights,
                                  output_path=NA,
                                  input_format = "otsoft",
                                  in_sep = "\t", out_sep = '\t') {

  if (input_format == "otsoft") {
    input <- load_data_otsoft(test_file, in_sep)
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

  if (!is.na(output_path)) {
    write.table(output, file=output_path, sep=out_sep, row.names = FALSE)
  }

  return(output)
}
