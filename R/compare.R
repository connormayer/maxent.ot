#' @export
compare_models_lrt <- function(model1, model2) {
  if (model1$k < model2$k) {
    chi_sq_val <- 2 * (model2$loglik - model1$loglik)
    k <- model2$k - model1$k
  } else {
    stop("Model 1 must have fewer parameters than model two for the ",
         "likelihood ratio test to be applied. These parameters should ",
         "be a subset of Model 2's parameters, but this function does not ",
         "explicitly check for this.")
  }
  result <- list(chiq_sq = chi_sq_val, p_value = pchisq(chi_sq_val, k))
  return(result)
}

#' @export
compare_models_aic <- function(...) {
  models <- list(...)
  if (length(models) < 2) {
    stop("Must provide at least two models to compare.")
  }
  aics <- sapply(models, calculate_aic)
  return(aics)
}

#' @export
compare_models_aic_c <- function(...) {
  models <- list(...)
  if (length(models) < 2) {
    stop("Must provide at least two models to compare.")
  }
  aic_cs <- sapply(models, calculate_aic_c)
  return(aic_cs)
}

#' @export
compare_models_bic <- function(...) {
  models <- list(...)
  if (length(models) < 2) {
    stop("Must provide at least two models to compare.")
  }
  bics <- sapply(models, calculate_bic)
  return(bics)
}

calculate_aic <- function(model) {
  return(2 * model$k - 2 * model$loglik)
}

calculate_aic_c <- function(model) {
  value <- calculate_aic(model) + (2 * model$k^2 + 2 * model$k) / model$n - model$k - 1
  return(value)
}

calculate_bic <- function(model) {
  return(model$k * log(model$n) - 2 * model$loglik)
}
