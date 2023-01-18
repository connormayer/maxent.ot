data {
  int<lower=1> N;  // number of tableaux -- that is to say, number of datapoints
  int<lower=0> max_num_candidates;
  array[N] int<lower=2, upper=max_num_candidates> J; // vector of num_candidate in each tableau
  int<lower=1> K; //  number of constraints; this is constant across tableaux
  array[K] int<lower=0, upper=1> is_constraint_hierarhical_by_items;
  int<lower=1> L; // number of levels of the random effect
  array[N] int<lower=1, upper=L> ll; // an array telling you which level of the random effect a given UR belongs to.
  array[N] int<lower=1,upper=max_num_candidates> Y; // index of winning candidate in each tableau
  array[N] matrix[max_num_candidates, K] violations; // violations
  vector[K] constraint_mus;
  vector[K] constraint_sigmas;
}


transformed data {
  int num_constraints_hierarchical_on_items = sum(is_constraint_hierarhical_by_items);
  array[K] int t_s;
  for (i in 1:K){
    t_s[i] = (i > 1 ? t_s[i-1] : 0) + is_constraint_hierarhical_by_items[i];
  }
}


parameters {
  array[K] real mu; // the vector of constraint weights, coefficients - group level
  array[K] real<lower=0> sigma; // the variance on these guys
  matrix[L,num_constraints_hierarchical_on_items] beta_params; // a vector of random-effect-level coefficients that are drawn from real_mu, one for each leve lof the random effect, for only those constraints that are hierarchical. these are the random-effect-level offsets // means drawn from the group-level real mu
  //matrix[L, K] beta;
}


transformed parameters {
  array[L] vector[K] beta;
{
    // local variable to index which parameter value to use
    //int param_ind = 1;
    for (constraint_index in 1:K) {
      if (is_constraint_hierarhical_by_items[constraint_index] == 1) {
        for (level_of_item_hierarchical_structure in 1:L) {

        beta[level_of_item_hierarchical_structure,constraint_index] = beta_params[level_of_item_hierarchical_structure,t_s[constraint_index]];
        //param_ind += 1;
      }} else {
                for (level_of_item_hierarchical_structure in 1:L) {

        beta[level_of_item_hierarchical_structure,constraint_index] = mu[constraint_index];
                }
      }
    }
  }
}




model {
mu ~ normal(constraint_mus,constraint_sigmas); // prior on mean of each mu
//sigma ~ normal(constraint_sigmas,1); // prior on sigma for that mu
to_vector(beta_params) ~ normal(mu,1); // this is the variance on how much indidivual random offets can be from the mu

for (i in 1:N) { // then for each datapoint
    Y[i] ~ categorical_logit(block(violations[i], 1, 1, J[i], K)* beta[ll[i]]); // take the datapoint n, get the random effect group it belongs to - ll - then index into the matrix of betas to find that row, and multiply it by the constraint violation matrix - this is what needs to be blocked on my end - for that one.
  }
}


generated quantities {
  vector[N] sim_winners;
  vector[N] log_lik;
  for (i in 1:N){
    sim_winners[i] = categorical_logit_rng(block(violations[i], 1, 1, J[i], K)*beta[ll[i]]);
    log_lik[i] = categorical_logit_lpmf(Y[i] | block(violations[i], 1, 1, J[i], K)*beta[ll[i]]);
  }
}
