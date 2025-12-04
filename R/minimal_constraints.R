#' Find the optimal minimal subset of weights for a given dataset and model choice.
#'
#' Uses stepwise regression with backwards elimination to find all subsets of
#' weights and compares these subset models to one another.
#' @export
minimal_constraints <- function(dataset, method = "aic", approach = "iterative") {
  if (approach == "iterative") {
    all_models <- step_iterative(dataset, method)
  }
  else if (approach == "recursive") {
    all_models <- step_recursive(dataset, method)
  }
  else {
    stop(sprintf(paste("Approach must be 'iterative' or 'recursive'")))
  }

  comparison <- do.call(compare_models, c(all_models, list(method = method)))
  out_object <- list(
    constraints = lapply(all_models, function(m) {
      w <- unlist(m$weights)
      as.list(w)
    }),
    comparison_table = comparison
  )
  return(out_object)
}


# Returns all constraint names from a given dataset
get_constraint_names <- function(dataset) {
  colnames(dataset)[4:ncol(dataset)]
}


# Compares given models using lrt with compare_models
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


# Compares models using aic, aic_c, or bic
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
step_recursive <- function(dataset, method = "aic") {
  all_models <- list()
  constraint_names <- get_constraint_names(dataset)

  generate_subsets <- function(constraints) {
    non_constraint_cols <- c("Input", "Output", "Frequency")
    keep_columns <- c(non_constraint_cols, constraints)

    subset_dataset <- dataset[, keep_columns, drop = FALSE]

    model <- optimize_weights(subset_dataset)
    all_models[[length(all_models) + 1]] <<- model

    if (length(constraints) == 1) return()

    candidate_models <- list()
    candidate_constraints <- list()

    for (to_drop in constraints) {
      candidate_subset <- constraints[constraints != to_drop]

      keep_cand <- c(non_constraint_cols, candidate_subset)
      candidate_dataset <- dataset[, keep_cand, drop = FALSE]

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
    else {
      stop(sprintf(paste("Chosen method must be one of: 'aic', 'aic_c', 'bic', 'lrt'")))
    }

    if (!is.null(best_subset)) {
      generate_subsets(best_subset)
    }
  }

  generate_subsets(constraint_names)
  return(all_models)
}

# iterative approach
step_iterative <- function(dataset, method = "aic") {
  all_models <- list()
  constraint_names <- get_constraint_names(dataset)

  constraints <- constraint_names
  non_constraint_cols <- c("Input", "Output", "Frequency")

  while (length(constraints) > 1) {
    keep_columns <- c(non_constraint_cols, constraints)
    subset_dataset <- dataset[, keep_columns, drop = FALSE]

    model <- optimize_weights(subset_dataset)
    all_models[[length(all_models) + 1]] <- model

    candidate_models <- list()
    candidate_constraints <- list()

    for (to_drop in constraints) {
      candidate_subset <- constraints[constraints != to_drop]
      keep_cand <- c(non_constraint_cols, candidate_subset)

      candidate_dataset <- dataset[, keep_cand, drop = FALSE]
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
    else {
      stop(sprintf(paste("Chosen method must be 'lrt', 'aic', 'aic_c', or 'bic'")))
    }

    if (is.null(best_subset)) break
    constraints <- best_subset
  }
  if (length(constraints) == 1) {
    keep_columns <- c(non_constraint_cols, constraints)
    final_dataset <- dataset[, keep_columns, drop = FALSE]
    final_model <- optimize_weights(final_dataset)
    all_models[[length(all_models) + 1]] <- final_model
  }

  return(all_models)
}
