library(tidyverse)
library(maxent.ot)

shona_input <- read_csv("inst/extdata/shona.csv")
processed_input <- maxent.ot:::load_input(shona_input)
data <- processed_input$data

partitions <- maxent.ot:::partition_data(data, k = 5)

train_data_list <- vector("list", 5)
test_data_list <- vector("list", 5)

training_data <- data %>% mutate(Frequency = 0)
test_data <- training_data

for (hold_out in 1:5) {

  training_part <- partitions %>%
    filter(partition != hold_out)

  test_part <- partitions %>%
    filter(partition == hold_out)

  train_data_list[[hold_out]] <-
    maxent.ot:::populate_tableau(
      training_data,
      training_part
    )

  test_data_list[[hold_out]] <-
    maxent.ot:::populate_tableau(
      test_data,
      test_part
    )
}

result <- optimize_bias_max(
  train_data = train_data_list,
  test_data = test_data_list,
  mu = 0,
  method = "L-BFGS-B"
)

cat("Best sigma:", result$sigma, "\n")
cat("Best log-likelihood:", result$mean_ll, "\n")
cat("Convergence:", result$convergence, "\n")

cat("\nSigma for each split:\n")
print(result$sigmas)
