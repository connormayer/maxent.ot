source("read-data.R")

calculate_log_likelihood <- function(constraint_weights, data) {
  # Sum of likelihoods of each datum
  ll = 0
  for (datum in data) {
    ll <- ll + calculate_datum_likelihood(constraint_weights, datum)
  }
  return(ll)
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

optimize <- function(input_file, constraint_weights) {
  input <- load_data_otsoft(input_file)
  long_names <- input[[1]]
  short_names <- input[[2]]
  data <- input[[3]]

  weights <- rep(1, length(long_names))
  best <- optim(
    weights,
    calculate_log_likelihood,
    data=data,
    control=list(fnscale = -1),
    lower=rep(0, length(weights)),
    # This needs a very large finite number, otherwise it lands on suboptimal weights
    upper=rep(2^64, length(weights)),
    method="L-BFGS-B")
  return(best)
}

optimize("test-files/SampleDataFile.txt")
