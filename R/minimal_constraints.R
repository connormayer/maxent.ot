#' Find the optimal minimal subset of weights for a given dataset and model choice.
#'
#' Uses stepwise regression with backwards elimination to find all subsets of
#' weights and compares these subset models to one another, finding the
#' subset that results in the best performance overall
#' @export
minimal_constraints <- function(dataset, method = "aic", approach = "naive") {
  #TODO implement naive approach
  if (approach == "naive") {
    all_models <- step_naive(dataset, method)
  }
  else if (approach == "recursive") {
    all_models <- step(dataset, method)
  }
  comparison <- compare_models(all_models, method = method)
  return(comparison)
}



get_constraint_names <- function(dataset) {
  colnames(dataset)[4:ncol(dataset)]
}



compare_lrt <- function(model, candidate_models, candidate_constraints) {
  best_pval <- 1
  best_subset <- NULL

  for (i in seq_along(candidate_models)) {
    candidate <- candidate_models[[i]]
    res <- compare_models(model, candidate, method = "lrt")
    pval <- res$p_value
    if (pval < best_pval) {
      best_pval <- pval
      best_subset <- candidate_constraints[[i]]
    }
  }
  best_subset
}



compare_aic <- function(model, candidate_models, candidate_constraints, method) {
  res <- do.call(compare_models, c(list(model = model, method = method), candidate_models))
  best_name <- res$model[1]
  for (i in seq_along(candidate_models)) {
    if (candidate_models[[i]]$name == best_name) {
      return(candidate_constraints[[i]])
    }
  }
  return(NULL)
}



#'Recursively generates all subsets of weights and fits optimal models to each
#'subset, ultimately returning all models (to be used in minimal_constraints)
step1 <- function(dataset, method = "aic") {
  all_models <- list()
  constraint_names <- get_constraint_names(dataset)

  generate_subsets <- function(constraints) {
    non_constraint_cols <- colnames(dataset)[!(colnames(dataset) %in% constraint_names)]
    keep_columns <- c(constraints, non_constraint_cols)
    subset_dataset <- dataset[, keep_columns, drop = FALSE]

    model <- optimize_weights(subset_dataset)
    all_models[[length(all_models) + 1]] <<- model

    if (length(constraints) == 1) return()

    candidate_models <- list()
    candidate_constraints <- list()

    for (to_drop in constraints) {
      candidate_subset <- constraints[constraints != to_drop]
      candidate_dataset <- dataset[, c(candidate_subset, non_constraint_cols), drop = FALSE]
      candidate_model <- optimize_weights(candidate_dataset)
      candidate_models[[length(candidate_models) + 1]] <- candidate_model
      candidate_constraints[[length(candidate_constraints) + 1]] <- candidate_subset
    }

    if (method == "lrt") {
      best_subset <- compare_lrt(model, candidate_models, candidate_constraints)
    }
    else if (method %in% c("aic", "aic_c", "bic")) {
      best_subset <- compare_aic(model, candidate_models, candidate_constraints, method)
    }

    if (!is.null(best_subset)) {
      generate_subsets(best_subset)
    }
  }

  generate_subsets(constraint_names)
  return(all_models)
}



# TODO: iterate through constraints, naively dropping the constraint that
# performs worst, output each of these n - 1 subsets
step_naive <- function(dataset, weights){

}







