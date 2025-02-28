data {
  int<lower=1> Y;                // Number of years
  array[Y] int<lower=0> N_obs;    // Observed total lake abundance
  array[Y] int<lower=0> A1_obs;   // Observed age-1 migrants in samples
  array[Y] int<lower=0> A_total;  // Total aged fish in outmigration samples
  
  real<lower=0> obs_error_prior;        // obs error on total abundance negative binomial
  real mu_theta_prior;   // Mean prior for age 1 outmigration proportion
  real<lower=0> sigma_theta_prior;// SD prior for age1 outmigration proportion
  real mu_M_prior;       // Mean prior for mortality of age 2s in year 2
  real<lower=0> sigma_M_prior;    // SD prior for mortality
  
  real<lower=0> N2_init_prior;    //prior for mean of N2s in first year
  real<lower=0> N2_init_sigma_prior; //prior for sd of N2s in first year
}

parameters {
  real<lower=0> obs_error;  // Observation error for total abundance negative binomial
  real mu_theta;         // Population-level mean of age 1 outmigration proportion in logit space
  real mu_M;             // Population-level mean of mortality in logit space
  
  real<lower=0> sigma_theta;      // Yearly variation in outmigration
  real<lower=0> sigma_M;          // Yearly variation in mortality
  
  array[Y] real<lower=0, upper=1> theta_raw;  // Yearly outmigration deviation - z score
  array[Y] real<lower=0> M_raw;   // Yearly mortality deviation - z score
  
  real<lower=0> init_N2; //initialization for age 2 fish in year 1
}

transformed parameters {
  array[Y] real<lower=0, upper=1> theta;
  array[Y] real<lower=0> M;
  array[Y] real<lower=0> N1; // Age-1 in lake
  array[Y] real<lower=0> N2; // Age-2 in lake
  array[Y] real<lower=0> N_lake; // Total lake abundance
  array[Y] real<lower=0> O1; // Age-1 outmigrants
  array[Y] real<lower=0> O2; // Age-2 outmigrants
  
  // Hierarchical transformation
  for (y in 1:Y) {
    theta[y] = inv_logit(mu_theta + sigma_theta * theta_raw[y]);
    M[y] = inv_logit(mu_M + sigma_M * M_raw[y]);
  }
  
  // State dynamics
  for (y in 1:Y) {
    if (y == 1) {
      N2[y] = init_N2;  // No age-2s in first year
    } else {
      N2[y] = N1[y-1] * (1 - theta[y-1]) * M[y-1];
    }
    
    N1[y] = N_obs[y] - N2[y];  // Implicit calculation from total
    N_lake[y] = N1[y] + N2[y];
    
    // Outmigrants
    O1[y] = theta[y] * N1[y];
    O2[y] = N2[y];
  }
}

model {
  // Priors
  obs_error ~ exponential(obs_error_prior);
  mu_theta ~ normal(mu_theta_prior, sigma_theta_prior);
  mu_M ~ normal(mu_M_prior, sigma_M_prior);
  sigma_theta ~ normal(0, sigma_theta_prior);
  sigma_M ~ normal(0, sigma_M_prior);
  
  //initialization of age 2s
  init_N2 ~ lognormal(log(N2_init_prior), N2_init_sigma_prior);
  
  // Hierarchical deviations - z scores
  theta_raw ~ normal(0,1);
  M_raw ~ normal(0, 1);
  
  // Observation model: Negative binomial for total lake abundance
  for (y in 1:Y) {
    N_obs[y] ~ neg_binomial_2(N_lake[y], 1/obs_error);
  }
  
  // Observation model: Binomial for age composition in outmigrants
  for (y in 1:Y) {
    A1_obs[y] ~ binomial(A_total[y], O1[y] / (O1[y] + O2[y]));
  }
}
