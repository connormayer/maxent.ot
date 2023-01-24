data {
  int<lower=1> N;  // number of tableaux -- that is to say, number of datapoints
  int<lower=0> max_num_candidates;
  array[N] int<lower=2, upper=max_num_candidates> J; // vector of num_candidate in each tableau
  int<lower=1> K; //  number of constraints; this is constant across tableaux
  array[K] int<lower=0, upper=1> is_constraint_hierarchical_by_items;
  array[K] int<lower=0, upper=1> is_constraint_hierarchical_by_subjects;
  int<lower=1> L_items; // number of levels of the random effect by items
  int<lower=1> L_subjects; // number of levels of the random effect by subjects

  array[N] int<lower=1, upper=L_items> ll_items; // an array telling you which level of the random effect a given UR belongs to.
  array[N] int<lower=1, upper=L_subjects> ll_subjects; // an array telling you which level of the random effect a given UR belongs to.

  array[N] int<lower=1,upper=max_num_candidates> Y; // index of winning candidate in each tableau
  array[N] matrix[max_num_candidates, K] violations; // violations
  vector[K] constraint_mus;
  vector[K] constraint_sigmas;
  int<lower=0,upper=1> run_estimation;
  int<lower=0,upper=1> fit_mus;
  int<lower=0,upper=1> fit_sigmas;


}


transformed data {

  int num_constraints_hierarchical_on_items = sum(is_constraint_hierarchical_by_items);
  int num_constraints_hierarchical_on_subjects = sum(is_constraint_hierarchical_by_subjects);

  array[K] int t_s;
  for (i in 1:K){
    t_s[i] = (i > 1 ? t_s[i-1] : 0) + is_constraint_hierarchical_by_items[i];
  }
  array[K] int t_i;
  for (i in 1:K){
    t_i[i] = (i > 1 ? t_i[i-1] : 0) + is_constraint_hierarchical_by_subjects[i];
  }

}


parameters {

  //if (fit_mus == 1){
  array[K] real<lower=0> distribution_of_likely_mus;
  //}

  //if (fit_sds == 1){
  array[K] real<lower=0> distribution_of_likely_sigmas;
  //}

  array[K] real<lower=0>  mu; // the vector of constraint weights, coefficients - group level
  //array[K] real<lower=0> sigma; // the variance on these guys
  matrix[L_items, num_constraints_hierarchical_on_items] beta_params_items; // a vector of random-effect-level coefficients that are drawn from real_mu, one for each leve lof the random effect, for only those constraints that are hierarchical. these are the random-effect-level offsets // means drawn from the group-level real mu
  matrix[L_subjects, num_constraints_hierarchical_on_subjects] beta_params_subjects; // a vector of random-effect-level coefficients that are drawn from real_mu, one for each leve lof the random effect, for only those constraints that are hierarchical. these are the random-effect-level offsets // means drawn from the group-level real mu

}


transformed parameters {
  array[L_items, L_subjects] vector[K] beta;

  {

  for (constraint_index in 1:K)	{
    if (is_constraint_hierarchical_by_items[constraint_index] == 1) {
      for (level_of_hierarchical_structure_items in 1:L_items) {
        if (is_constraint_hierarchical_by_subjects[constraint_index] == 1) {
          for (level_of_hierarchical_structure_subjects in 1:L_subjects) {
            beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = (beta_params_items[level_of_hierarchical_structure_items, t_i[constraint_index]]+beta_params_subjects[level_of_hierarchical_structure_subjects, t_s[constraint_index]] - mu[constraint_index]);

          }}
          else {
            for (level_of_hierarchical_structure_subjects in 1:L_subjects){
              beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = beta_params_items[level_of_hierarchical_structure_items, t_i[constraint_index]];
            }}}}
            else{
              for (level_of_hierarchical_structure_items in 1:L_items){
                if (is_constraint_hierarchical_by_subjects[constraint_index] == 1){
                  for (level_of_hierarchical_structure_subjects in 1:L_subjects)	{
                    beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = beta_params_subjects[level_of_hierarchical_structure_subjects, t_s[constraint_index]];
                  }}
                  else		{
                    for (level_of_hierarchical_structure_subjects in 1:L_subjects)	{
                      beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = mu[constraint_index];
                    }}}}


  }
}
}

model {
  if (fit_mus == 1){
    distribution_of_likely_mus ~ normal(0,10);

  }
  if (fit_sigmas == 1){
    distribution_of_likely_sigmas ~ exponential(1);
}

  if(fit_mus == 1 &&  fit_sigmas == 1){
      mu ~ normal(distribution_of_likely_mus,distribution_of_likely_sigmas); // prior on mean of each mu

  } else{ if (fit_mus == 1 &&  fit_sigmas == 0){
          mu ~ normal(distribution_of_likely_mus,constraint_sigmas); // prior on mean of each mu

  } else{ if (fit_mus == 0 &&  fit_sigmas == 1){
          mu ~ normal(constraint_mus,distribution_of_likely_sigmas); // prior on mean of each mu

  }else {
      mu ~ normal(constraint_mus,constraint_sigmas); // prior on mean of each mu

  }}}

  for (i in 1:num_constraints_hierarchical_on_items){
    col(beta_params_items,i) ~ normal(mu,1);
  }
  for (i in 1:num_constraints_hierarchical_on_subjects){
    col(beta_params_subjects,i) ~ normal(mu,1);
  }

  if (run_estimation == 1){
  for (i in 1:N) { // then for each datapoint
    Y[i] ~ categorical_logit(block(violations[i], 1, 1, J[i], K)* beta[ll_items[i], ll_subjects[i]]); // take the datapoint n, get the random effect group it belongs to - ll - then index into the matrix of betas to find that row, and multiply it by the constraint violation matrix - this is what needs to be blocked on my end - for that one.
  }}
}


generated quantities {
  vector[N] sim_winners;
  vector[N] log_lik;
  for (i in 1:N){
    sim_winners[i] = categorical_logit_rng(block(violations[i], 1, 1, J[i], K)*beta[ll_items[i], ll_subjects[i]]);
    log_lik[i] = categorical_logit_lpmf(Y[i] | block(violations[i], 1, 1, J[i], K)*beta[ll_items[i],ll_subjects[i]]);
  }
}
