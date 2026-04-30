library(maxent.ot)
library(dplyr)

shona_input <- read.csv("inst/extdata/shona.csv")

processed_input <- maxent.ot:::load_input(shona_input)
data <- processed_input$data


k <- 5

inputs <- unique(data$Input)

set.seed(123)
fold_ids <- sample(rep(1:k, length.out = length(inputs)))

input_folds <- data.frame(
  Input = inputs,
  fold = fold_ids
)

data_with_folds <- merge(data, input_folds, by = "Input")

hold_out <- 1

training_part <- data_with_folds %>%
  filter(fold != hold_out)

test_part <- data_with_folds %>%
  filter(fold == hold_out)

training_data <- data %>% mutate(Frequency = 0)
test_data     <- training_data

training_tableau <- maxent.ot:::populate_tableau(training_data, training_part)
test_tableau     <- maxent.ot:::populate_tableau(test_data, test_part)

result <- optimize_bias(
  train_data = training_part,
  test_data  = test_part,
  sigma = 0.5,
  mu = 0
)

print(result)
