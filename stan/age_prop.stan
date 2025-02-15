data {
  int<lower=1> N;
  int<lower=1> L; // Number of lakes
  int<lower=1, upper=L> lake[N]; // Lake index
  int<lower=0, upper=1> observed[N];
  int age_1s[N]; // Observed proportion
  int sample_size[N]; // Sample sizes
  real mu_prior[L];
}

parameters {
  real mu[L]; // Lake-level mean proportion
  vector<lower=0>[L] p_sigma;
  real p_z[N];
}
transformed parameters{
  real<lower=0, upper=1> p[N]; // True proportions per observation
  for (i in 1:N) {
    p[i] = inv_logit(mu[lake[i]] +  p_sigma[lake[i]] * p_z[i]); // Hierarchical beta prior
  }
}
model {
  for(i in 1:L){
  mu[i] ~ normal(mu_prior[i], 0.3);
  }
  p_sigma ~ gamma(3, 3);
  p_z ~ normal(0, 1);
  for (i in 1:N) {
    if(observed[i] == 1){
      age_1s[i] ~ binomial(sample_size[i], p[i]); // Proper binomial likelihood
    }
  }
}

