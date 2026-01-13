#' Find the optimal minimal subset of weights for a given dataset, method of
#' comparison, and approach.
#'
#' Uses stepwise regression with backwards elimination to find all subsets of
#' weights and compares these subset models to one another.
#'
#' @param dataset The input data frame/data table/tibble. This should contain one
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
#'
#' @param method The chosen method of comparison to be used to determine the
#'  subset of weights. Can be chosen from one of: 'aic', 'bic', 'aic_c', 'lrt'
#'  which specify the usage of:
#'    - Akaike Information Criterion (**aic**)
#'    - Bayesian Information Criterion (**bic**)
#'    - Akaike Information Criterion corrected for small sample sizes (**aic_c**)
#'    - Likelihood Ratio Test (**lrt**)
#'
#'  Refer to the `compare.R` file for further description of these modes of comparison
#'
#' @param approach A string specifying the desired approach, one of:
#'  * `recursive`: Uses pure recursion to explore every possible subset and
#'  drops suboptimal weights according to whichever method is specified. Max
#'  recursion depth affects this approach.
#'  * `recursive_without_pruning`: Uses recursion where suboptimal subsets are
#'  pruned to minimize recursive calls. Max recursion depth affects this approach.
#'  * `iterative`: Uses iteration to drop the current suboptimal weight according
#'  to the specified method. Max recursion depth does not affect this approach.
#' @param max_depth (optional) The maximum possible recursion depth if the
#'  approach is specified as "recursive_without_pruning"
#' @return An object with the following named attributes:
#' \itemize{
#'         \item `constraints`: The constraints used for each model subset.
#'         \item `comparison_table`: A data table containing the results of running
#'         `compare_models` on all created subset models.
#' }
#' @examples
#'   TODO implement
#' @export
minimal_constraints <- function(dataset, method = "aic",
                                approach = "iterative", max_depth = Inf) {

  all_constraint_names <- get_constraint_names(dataset)

  if (approach == "iterative") {
    all_models <- step_iterative(dataset, method)
  }
  else if (approach == "recursive") {
    all_models <- step_recursive(dataset, method)
  }
  else if (approach == "recursive_without_pruning") {
    all_models <- step_without_pruning(dataset, max_depth)
  }
  else {
    stop("Approach must be 'iterative', 'recursive', or 'recursive_without_pruning'")
  }

  if (method %in% c("aic", "aic_c", "bic")) {
    comparison <- do.call(compare_models, c(all_models, list(method = method)))
  }
  else if (method == "lrt") {
    full_model <- all_models[[1]]
    best_index <- 1
    best_ll <- full_model$log_likelihood

    for (i in 2:length(all_models)) {
      ll <- all_models[[i]]$log_likelihood
      if (length(ll) && !is.na(ll) && ll > best_ll) {
        best_ll <- ll
        best_index <- i
      }
    }

    best_subset_model <- all_models[[best_index]]
    full_weights = full_model$weights
    best_weights = best_subset_model$weights
    if (any(full_weights != best_weights)){
      comparison <- compare_models(full_model, best_subset_model, method = "lrt")
    }
    else{
      comparison <- full_model
    }
  }
  else {
    stop("Method must be 'aic','aic_c','bic', or 'lrt'")
  }

  constraints <- lapply(all_models, function(m) {
    w <- m$weights
    full_w <- setNames(rep(NA, length(all_constraint_names)), all_constraint_names)
    for (name in names(w)) {
      full_w[name] <- w[[name]]
    }
    as.list(full_w)
  })

  out_object <- list(
    constraints = constraints,
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

step_without_pruning <- function(dataset, max_depth = Inf) {
  all_models <- list()
  constraint_names <- get_constraint_names(dataset)
  non_constraint_cols <- c("Input", "Output", "Frequency")

  recurse <- function(constraints, depth) {
    if (depth > max_depth) {
      return()
    }

    keep_columns <- c(non_constraint_cols, constraints)
    subset_dataset <- dataset[, keep_columns, drop = FALSE]
    model <- optimize_weights(subset_dataset)
    all_models[[length(all_models) + 1]] <<- model

    if (length(constraints) <= 1){
      return()
    }
    for (to_drop in constraints) {
      new_subset <- constraints[constraints != to_drop]
      recurse(new_subset, depth + 1)
    }
  }
  recurse(constraint_names, depth = 0)
  return(all_models)
}

