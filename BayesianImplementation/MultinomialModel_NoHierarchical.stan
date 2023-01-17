data {
  int<lower=1> N;  // number of tableaux -- that is to say, number of datapoints
  int<lower=0> max_num_candidates;
  array[N] int<lower=2, upper=max_num_candidates> J; // vector of num_candidate in each tableau
  int<lower=1> K; //  number of constraints; this is constant across tableaux
  array[N] int<lower=1,upper=max_num_candidates> Y; // index of winning candidate in each tableau
  array[N] matrix[max_num_candidates, K] violations; // violations
  vector[K] constraint_mus;
  vector[K] constraint_sigmas;
}

parameters {
  vector<lower=0>[K] beta;  // attribute effects  // constraint weights; should be just one
}


model {
  //beta ~ normal(0,1);
  beta ~ normal(constraint_mus,constraint_sigmas);
  for (i in 1:N) // look up "block" - just looks at a specific part of a matrix
    //local_violations = ; // this is supposed to get the i'th tableau's candidates and violations from the stack violations (which is num_tableaux * max_num_candidates (padded out) * num_constraints in size), and pull out only the first J[i] rows, which are just those that are nonzero and not padded out so relevant here, and then all columns (the K constraints, which are constant for all tableaux)
    Y[i] ~ categorical_logit(block(violations[i], 1, 1, J[i], K)*beta); // so this gets just one instance, one token, at a time. could we ramp up to do multinomial and not categorical?
}

generated quantities {
  array[N] int sim_winners;
  for (i in 1:N)
    sim_winners[i] = categorical_logit_rng(block(violations[i], 1, 1, J[i], K)*beta);
}
