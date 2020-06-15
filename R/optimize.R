source("read-data.R")

calculate_log_likelihood <- function(constraint_weights, data, bias_params) {
  # Sum of likelihoods of each datum
  ll = 0
  for (datum in data) {
    ll <- ll + calculate_datum_likelihood(constraint_weights, datum)
  }
  if (!is.na(bias_params)) {
    ll <- ll - calculate_bias(bias_params[,1], bias_params[,2], constraint_weights)
  }
  return(ll)
}

calculate_bias <- function(mus, sigmas, constraint_weights) {
  top <- (constraint_weights - mus)^2
  bottom <- 2 * sigmas^2
  bias <- sum(top / bottom)
  return(bias)
}

calculate_datum_likelihood <- function(constraint_weights, datum) {
  freqs <- datum[,3]
  violations <- datum[,4:ncol(datum)]
  violations_mat <- data.matrix(violations)
  violations_mat[is.na(violations_mat)] <- 0
  harmony <- violations_mat %*% constraint_weights
  e_harmonies <- exp(1)^-harmony
  z <- sum(e_harmonies)
  candidate_probs <- e_harmonies / z
  log_candidate_probs <- log(e_harmonies / z)
  log_prob <- sum(freqs * log_candidate_probs)
  return(log_prob)
}

optimize <- function(input_file, constraint_weights=NA, bias_file=NA) {
  input <- load_data_otsoft(input_file)
  if (!is.na(bias_file)) {
    bias_params <- load_bias_file_otsoft(bias_file)
  } else {
    bias_params <- NA
  }
  long_names <- input[[1]]
  short_names <- input[[2]]
  data <- input[[3]]

  if (is.na(constraint_weights)) {
    constraint_weights <- rep(1, length(long_names))
  }
  best <- optim(
    constraint_weights,
    calculate_log_likelihood,
    data=data,
    bias_params=bias_params,
    control=list(fnscale = -1),
    lower=rep(0, length(constraint_weights)),
    # This needs a very large finite number, otherwise it lands on suboptimal weights
    upper=rep(2^64, length(constraint_weights)),
    method="L-BFGS-B"
  )
  out_weights <- best[[1]]
  names(out_weights) <- long_names
  return(out_weights)
}

#bar <- optimize("test-files/SampleDataFile.txt", bias_file="test-files/SampleConstraintFile.txt" )
