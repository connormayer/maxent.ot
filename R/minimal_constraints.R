#' Find the optimal minimal subset of weights for a given dataset and model choice.
#'
#' Uses stepwise regression with backwards elimination to find all subsets of
#' weights and compares these subset models to one another, finding the
#' subset that results in the best performance overall
#' @export
minimal_constraints <- function(dataset, model, weights, method = "aic") {
  # Find all possible weight subsets and store them in a vector
  all_models <- step(dataset, model, weights)
  # call compare models, entering in all subsets as the possible models
  comparison <- compare_models(all_models, method = method)
  # outputs a table showcasing the best performing models, of which their
  # weights can be accessed
  return(comparison)
}


#'Recursively generates all subsets of weights and fits optimal models to each
#'subset, ultimately returning all models (to be used in minimal_constraints)
step <- function(dataset, model, weights) {
  all_models <- list()
  constraint_names <- names(weights)
  # Recursively generates all weight subsets
  generate_subsets <- function(constraints) {
    if (length(constraints) == 0) return(list())
    first <- constraints[1]
    rest <- constraints[-1]
    subsets_rest <- generate_subsets(rest)
    subsets <- list(c(first))
    for (subset in subsets_rest) {
      subsets[[length(subsets) + 1]] <- c(first, subset)
    }
    all_subsets <- c(subsets, subsets_rest)
    return(all_subsets)
  }
  all_subsets <- generate_subsets(constraint_names)
  for (subset in all_subsets) {
    # Keep relevant columns
    non_constraint_cols <- colnames(dataset)[!(colnames(dataset) %in% constraint_names)]
    # Subset including only selected constraints
    keep_columns <- c(subset, non_constraint_cols)
    subset_dataset <- dataset[, keep_columns, drop = FALSE]
    # optimize this sub model and append it to the entire list of models
    sub_model <- optimize_weights(subset_dataset)
    all_models[[length(all_models) + 1]] <- sub_model
  }

  return(all_models)
}





