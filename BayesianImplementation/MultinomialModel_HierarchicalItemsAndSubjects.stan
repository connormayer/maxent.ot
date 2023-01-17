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
}


transformed data {
  //  print("in transformed parameters!");
 // vector[(L_items*L_subjects)]transformed_constraint_mus = to_vector(append_col(constraint_mus,constraint_mus));
//  vector[(L_items*L_subjects)]transformed_constraint_sigmas = to_vector(append_col(constraint_sigmas,constraint_sigmas));
  //print(dims(transformed_constraint_mus),"are transformed constraint mus dims");
  //vector[L_items*L_subjects] transformed_constraint_mus;

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
  array[K] real mu; // the vector of constraint weights, coefficients - group level
  array[K] real<lower=0> sigma; // the variance on these guys
  matrix[L_items, num_constraints_hierarchical_on_items] beta_params_items; // a vector of random-effect-level coefficients that are drawn from real_mu, one for each leve lof the random effect, for only those constraints that are hierarchical. these are the random-effect-level offsets // means drawn from the group-level real mu
  matrix[L_subjects, num_constraints_hierarchical_on_subjects] beta_params_subjects; // a vector of random-effect-level coefficients that are drawn from real_mu, one for each leve lof the random effect, for only those constraints that are hierarchical. these are the random-effect-level offsets // means drawn from the group-level real mu

}


transformed parameters {
  array[L_items, L_subjects] vector[K] beta;
  //arrray[2] mu_longer;

  //print("in transformed params, dims beta is ", dims(beta));
  {
  //int param_ind_items = 1;
  //int param_ind_subjects = 1;
  for (constraint_index in 1:K)	{
    if (is_constraint_hierarchical_by_items[constraint_index] == 1) {
      for (level_of_hierarchical_structure_items in 1:L_items) {
        if (is_constraint_hierarchical_by_subjects[constraint_index] == 1) {
          for (level_of_hierarchical_structure_subjects in 1:L_subjects) {
            beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = (beta_params_items[level_of_hierarchical_structure_items, t_i[constraint_index]]+beta_params_subjects[level_of_hierarchical_structure_subjects, t_s[constraint_index]] - mu[constraint_index]);
            //param_ind_items += 1
            //	param_ind_subjects += 1
          }}
          else {
            for (level_of_hierarchical_structure_subjects in 1:L_subjects){
              beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = beta_params_items[level_of_hierarchical_structure_items, t_i[constraint_index]];
              //param_ind_items += 1
            }}}}
            else{
              for (level_of_hierarchical_structure_items in 1:L_items){
                if (is_constraint_hierarchical_by_subjects[constraint_index] == 1){
                  for (level_of_hierarchical_structure_subjects in 1:L_subjects)	{
                    beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = beta_params_subjects[level_of_hierarchical_structure_subjects, t_s[constraint_index]];
                    //param_ind_subjects += 1
                  }}
                  else		{
                    //print("in last else");
                    for (level_of_hierarchical_structure_subjects in 1:L_subjects)	{
                      //print("level of hierarchial structure subjects is", level_of_hierarchical_structure_subjects);
                      //print("mu constraint index is",dims(mu[constraint_index]));
                      beta[level_of_hierarchical_structure_items, level_of_hierarchical_structure_subjects, constraint_index] = mu[constraint_index];
                    }}}}


  }
}
}

model {
  //print("in model!");
  mu ~ normal(constraint_mus,constraint_sigmas); // prior on mean of each mu
  //print("dims of constrain_mus are", dims(constraint_mus),"constraint sigmas are",dims(constraint_sigmas),"mu is ", dims(mu));
  //sigma ~ normal(constraint_sigmas,1); // prior on sigma for that mu
  //  print("dims beta params subjects ",dims((beta_params_subjects)));

  //print("dims beta params subjects to vector",dims(to_vector(beta_params_subjects)));

  //mu_longer = append(mu, mu);
  //print(dims(transformed_constraint_mus),"are dims of transformed mu");

  for (i in 1:num_constraints_hierarchical_on_items){
    col(beta_params_items,i) ~ normal(mu,1);
  }
  for (i in 1:num_constraints_hierarchical_on_subjects){
    col(beta_params_subjects,i) ~ normal(mu,1);
  }
  //to_vector(beta_params_items) ~ normal(mu,1); // this is the variance on how much indidivual random offets can be from the mu
 // to_vector(beta_params_subjects) ~ normal(mu,1); // this is the variance on how much indidivual random offets can be from the mu
 //  // print("dims beta params ites",dims(beta_params_items));

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
