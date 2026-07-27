library(tidyverse)
library(maxent.ot)

shona_input <- read_csv("inst/extdata/shona.csv")
processed_input <- maxent.ot:::load_input(shona_input)
data <- processed_input$data

partitions <- maxent.ot:::partition_data(data, k = 5)
hold_out <- 1

training_part <- partitions %>% filter(partition != hold_out)
test_part <- partitions %>% filter(partition == hold_out)

training_data <- data %>% mutate(Frequency = 0)
test_data <- training_data

training_tableau <- maxent.ot:::populate_tableau(training_data, training_part)
test_tableau <- maxent.ot:::populate_tableau(test_data, test_part)

result <- optimize_bias_max(
  train_data = training_tableau,
  test_data = test_tableau,
  sigma = 0.5,
  mu = 0,
  method = "L-BFGS-B",
  n = 5
)

cat("Best sigma:", result$sigma, "\n")
cat("Best log-likelihood:", result$mean_ll, "\n")
cat("Best run:", result$best_run, "of", result$n_runs, "\n")
cat("Convergence:", result$convergence, "\n\n")
