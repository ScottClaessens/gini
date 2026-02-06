data {
  int<lower=0> N;              // total number of sites
  vector[N] gini_obs;          // observed gini estimates
  vector<lower=0>[N] gini_se;  // standard errors for gini estimates
  vector<lower=0>[N] pop_size; // population size estimates
}
parameters {
  real phi;             // location for gini 
  real<lower=0> tau;    // scale for gini
  vector[N] gini_true;  // unknown true gini values
  real mu;              // location for population size
  real<lower=0> sigma;  // scale for population size
}
model {
  // priors
  phi ~ normal(0.5, 0.05);
  tau ~ exponential(10);
  mu ~ normal(10, 2);
  sigma ~ exponential(2);
  // likelihoods
  gini_true ~ normal(phi, tau);
  gini_obs ~ normal(gini_true, gini_se);
  pop_size ~ lognormal(mu, sigma);
}
