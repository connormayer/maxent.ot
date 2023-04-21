library(maxent.ot)
library(tidyverse)
library(xtable)

setwd("C:/Users/conno/git_repos/maxent.ot/misc/pda_tutorial")
simple_input <- read_csv("data/simple_input_df.csv")
# simple_input <- "data/simple_input_otsoft.txt"

result1 <- predict_probabilities(simple_input, c(2,1))
result1$loglik
result1$predictions

best <- optimize_weights(simple_input)

result2 <- predict_probabilities(simple_input, best$weights)

