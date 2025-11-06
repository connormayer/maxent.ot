#' Find the optimal minimal subset of weights for a given dataset and model choice.
#'
#' Uses stepwise regression with backwards elimination to find all subsets of
#' weights and compares these subset models to one another, finding the
#' subset that results in the best performance overall
#' @export
minimal_constraints <- function(dataset, weights, method = "aic", approach = "naive") {
  #TODO implement naive approach
  if (approach == "naive") {
    all_models <- step_naive(dataset, model, weights)
  }
  else if (approach == "recursive") {
    all_models <- step(dataset, weights)
  }
  comparison <- compare_models(all_models, method = method)
  return(comparison)
}


#'Recursively generates all subsets of weights and fits optimal models to each
#'subset, ultimately returning all models (to be used in minimal_constraints)
step1 <- function(dataset, weights) {
  all_models <- list()
  constraint_names <- names(weights)

  #helper that generates optimal subsets from a given superset of constraints
  generate_subsets <- function(constraints) {
    if (length(constraints) == 0){
      return(list())
    }


    non_constraint_cols <- colnames(dataset)[!(colnames(dataset) %in% constraint_names)]
    keep_columns <- c(constraints, non_constraint_cols)
    subset_dataset <- dataset[, keep_columns, drop = FALSE]
    model <- optimize_weights(subset_dataset)
    all_models[[length(all_models) + 1]] <<- model

    if (length(constraints) == 1){
      return(list())
    }

    best_score <- Inf
    best_subset <- NULL
    # iterates through constraints and stores loglik for the subset excluding
    # this constraint to an array, ultimately choosing the best n - 1 subset
    for (to_drop in constraints) {
      candidate_subset <- constraints[constraints != to_drop]
      candidate_dataset <- dataset[, c(candidate_subset, non_constraint_cols), drop = FALSE]
      candidate_model <- optimize_weights(candidate_dataset)
      if (candidate_model$loglik < best_score) {
        best_score <- candidate_model$loglik
        best_subset <- candidate_subset
      }
    }
    # continues the recursion of this optimal subset to find further optimality
    generate_subsets(best_subset)
  }

  # compiles all generated subsets to later use in model comparison
  generate_subsets(constraint_names)
  empty_dataset <- dataset[, !(colnames(dataset) %in% constraint_names), drop = FALSE]
  all_models[[length(all_models) + 1]] <- optimize_weights(empty_dataset)

  return(all_models)
}

# TODO: iterate through constraints, naively dropping the constraint that
# performs worst, output each of these n - 1 subsets
step_naive <- function(dataset, weights){

}







