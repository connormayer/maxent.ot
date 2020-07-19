#' Compare Maxent OT models using a variety of methods
#'
#' Compares two or more fitted models to determine which provides the best fit
#' to the data, using a variety of methods.
#'
#' A few caveats for several of the comparison methods:
#'   * The likelihood ratio test (`lrt`) method applies to exactly two models,
#'     and assumes that the parameters of these models are *nested*: that is,
#'     the constraints in the smaller model are a strict subset of the
#'     constraints in the larger model. This function will verify this to some
#'     extent based on the number and names of constraints.
#'   * The Akaike Information Criterion adjusted for small sample sizes
#'     (`aic_c`) and the Bayesian Information Criterion (`bic`) rely on sample
#'     sizes in their calculations. The same size for a data set is defined as
#'     the sum of the column of surface form frequencies. If you want to apply
#'     these methods, it is important that the values in the column are token
#'     counts, not relative frequencies. Applying these methods to relative
#'     frequencies, which effectively ignore sample size, will produce invalid
#'     results.
#'   * TODO: What about comparing models with different biases?
#'
#' The `aic`, `aicc`, and `bic` comparison methods return raw AIC/AICc/BIC
#' values as well as weights corresponding to these values.These weights
#' are calculated similarly for each model:
#'
#' \deqn{W_i = \frac{\exp(-0.5 \delta_i)}{\sum_{j=1}^{m}{\exp(-0.5 \delta_j)}}}
#'
#' where \eqn{\delta_i} is the difference in score (AIC, AICc, BIC) between model
#' \eqn{i} and the model with the best score and \eqn{m} is the number of models
#' being compared. These weights provide the relative likelihood or conditional
#' probability of this model being the best model (by whatever definition of
#' "best" is assumed by the measurement type) given the data and the selection
#' of models it is being compared to.
#'
#' @param ... Two or more models objects to be compared. These objects should
#'   be in the same format as the objects returned by the `optimize_weights`
#'   function. Note that the likelihood ratio test applies to exactly two
#'   models, which the other comparison methods can be applied to arbitrarily
#'   many models.
#' @param method The method of comparison to use. This currently includes:
#'     * `lrt`: The likelihood ratio test. This method can only be applied to a
#'       maximum of two models, and the parameters of these models (i.e.) their
#'       constraints *must be in a strict subset/superset relationship.* If
#'       your models do not meet these requirements, you should use a different
#'       metod.
#'
#'       The likelihood ratio is calculated as follows:
#'
#'        \deqn{LR = 2(LL_2 - LL_1)}
#'
#'       where \eqn{LL_2} is log likelihood of the model with more parameters.
#'       A p-value is calculated by conducting a chi-squared test with
#'       \eqn{X^2 = LR} and the degrees of freedom set to the difference in
#'       number of parameters between the two models. This p-value tells us
#'       whether the difference in likelihood between the two models is
#'       significant (i.e., whether the extra parameters in the full model
#'       are justified by the increase in model fit).
#'     * `aic`: The Akaike Information Criterion. This is calculated as follows
#'       for each model:
#'
#'       \deqn{AIC = 2k - 2LL}
#'
#'       where \eqn{k} is the number of model parameters and LL is the model's
#'       log likelihood.
#'     * `aic_c`: The Akaike Information Criterion corrected for small sample
#'       This is calculated as follows:
#'
#'       \deqn{AIC_c = 2k - 2LL + \frac{2k^2 + 2k}{n - k - 1}}
#'
#'       where \eqn{n} is the number of samples an the other parameters are
#'       identical to those used in the AIC calculation. As \eqn{n} approaches
#'       infinity, this final term converges to 0. Please see the note earlier
#'       in this functions documentation for information about sample sizes.
#'     * `bic`: The Bayesian Information Criterion.` This is calculated as
#'       follows:
#'
#'       \deqn{BIC = k\ln(n) - 2LL}
#'
#'       As with `aic_c`, this calculation relies on the number of samples.
#'       Please see the discussion on sample sizes earlier in the documentation
#'       of this function before using this method.
#'
#' @return A data frame containing information about the comparison. The
#'   contents and size of this data frame vary depending on the method used.
#'     * `lrt`: A data frame with a single row and the following columns:
#'         + `description`: the names of the two models being compared. The
#'           name  of the model with more parameters will be first.
#'         + `chi_sq`: the chi-squared value calculated during the test,
#'         + `k_delta`: the difference in parameters between the two models
#'           used as degrees of freedom in the chi-squared test
#'         + `p_value` the p-value calculated by the test
#'     * `aic`: A data frame with as many rows as there were models passed
#'       in. The models are sorted in descending order of AIC (i.e., best
#'       first). This data frame has the following columns:
#'         + `model`: The name of the model
#'         + `k`: The number of parameters
#'         + `aic`: The model's AIC value
#'         + `aic.wt`: The model's AIC weight: this reflects the relative
#'           likelihood (or conditional probability) that this model is the
#'           "best" model in the set.
#'         + `cum.wt`: The cumulative sum of AIC weights up to and including
#'           this model.
#'         + `ll`: The log likelihood of this model.
#'     * `aicc`: The data frame returned here is analogous to the structure
#'       of the AIC dataframe, with AICc values replacing AICs and accordingly
#'       modified column names. There is one additional column:
#'         + `n`: The number of samples in the data the model is fit to.
#'     * `bic`: The data frame returned here is analogous to the structure of
#'       the AIC and AICc data frames. Like the AICc data frame, it contains
#'       the `n` column.
#'
#'
#' @export
compare_models <- function(..., method='lrt') {
  models <- list(...)

  if (length(models) < 2) {
    stop("Must provide at least two models to compare.")
  }

  if (method == 'lrt') {
    if (length(models) != 2) {
      stop("Likelihood ratio test requires exactly two models.")
    }

    if (models[[1]]$k < models[[2]]$k) {
      sub_model <- models[[1]]
      full_model <- models[[2]]
    } else if (models[[2]]$k < models[[1]]$k) {
      sub_model <- models[[2]]
      full_model <- models[[1]]
    }
    else {
      stop("Models have the same number of parameters, cannot apply ",
           "likelihood ratio test.")
    }

    if (!all(names(sub_model$weights) %in% names(full_model$weights))) {
      stop("Constraint names in the smaller model are not a subset ",
           "of those in the larger model. Models must be nested to ",
           "apply the likelihood ratio test.")
    }
    chi_sq_val <- 2 * (full_model$loglik - sub_model$loglik)
    k <- full_model$k - sub_model$k

    result <- data.frame(
      description = paste(full_model$name, sub_model$name, sep="~"),
      chi_sq = chi_sq_val,
      k_delta = k,
      p_value = pchisq(chi_sq_val, k)
    )
  } else {
    result <- get_scores(models, method)
  }
  return(result)
}

# Helper function that calculates aics for a set of models, then calculates
# AIC weights and deltas for the models, returning them in a data frame sorted
# by AIC score.
get_scores <- function(models, method='aic') {
  if (method == 'aic') {
    scores <- sapply(models, calculate_aic)
  } else if (method == 'aic_c') {
    scores <- sapply(models, calculate_aic_c)
  } else if (method == 'bic') {
    scores <- sapply(models, calculate_bic)
  } else {
    stop(sprintf("Invalid method %s", method))
  }
  sorted_ix <- sort.list(scores)
  sorted_scores <- scores[sorted_ix]
  sorted_models <- models[sorted_ix]

  for (i in 1:length(sorted_models)) {
    model <- sorted_models[[i]]
    model$score <- sorted_scores[i]
    sorted_models[[i]] <- model
  }

  best_model <- sorted_models[[1]]
  weight_sum <- 0

  for (i in 1:length(sorted_models)) {
    model <- sorted_models[[i]]
    model$delta <- model$score - best_model$score
    model$raw_weight <- exp(-0.5 * model$delta)
    weight_sum <- weight_sum + model$raw_weight
    sorted_models[[i]] <- model
  }

  for (i in 1:length(sorted_models)) {
    sorted_models[[i]]$weight <- sorted_models[[i]]$raw_weight / weight_sum
  }

  if (method == 'aic') {
    new_df <- data.frame(
      model = character(),
      k = integer(),
      aic = double(),
      aic.wt = double(),
      cum.wt = double(),
      ll = double()
    )
  } else if (method == 'aic_c') {
    new_df <- data.frame(
      model = character(),
      k = integer(),
      n = integer(),
      aicc = double(),
      aicc.Wt = double(),
      cum.wt = double(),
      ll = double()
    )
  } else if (method == 'bic') {
    new_df <- data.frame(
      model = character(),
      k = integer(),
      n = integer(),
      bic = double(),
      bic.wt = double(),
      cum.wt = double(),
      ll = double()
    )
  }

  cumul_weight <- 0

  for (model in sorted_models) {
    cumul_weight <- cumul_weight + model$weight
    if (method == 'aic') {
      m_df <- data.frame(
        model = model$name,
        k = model$k,
        aic = model$score,
        aic.wt = model$weight,
        cum.wt = cumul_weight,
        ll = model$loglik
      )
    } else if (method == 'aic_c') {
      m_df <- data.frame(
        model = model$name,
        k = model$k,
        n = model$n,
        aicc = model$score,
        aicc.wt = model$weight,
        cum.wt = cumul_weight,
        ll = model$loglik
      )
    } else if (method == 'bic') {
      m_df <- data.frame(
        model = model$name,
        k = model$k,
        n = model$n,
        bic = model$score,
        bic.wt = model$weight,
        cum.wt = cumul_weight,
        ll = model$loglik
      )
    }
    new_df <- rbind(new_df, m_df)
  }

  return(new_df)
}

# Helper function that calculates AIC for a single model
calculate_aic <- function(model) {
  return(2 * model$k - 2 * model$loglik)
}

# Helper function that calculates AICc for a single model
calculate_aic_c <- function(model) {
  value <- calculate_aic(model) + (2 * model$k^2 + 2 * model$k) / (model$n - model$k - 1)
  return(value)
}

# Helper function that calculates BIC for a single model
calculate_bic <- function(model) {
  return(model$k * log(model$n) - 2 * model$loglik)
}
