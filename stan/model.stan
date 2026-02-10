data {
  int<lower=0> N;                         // total number of sites
  vector<lower=0, upper=1>[N] gini_obs;   // observed gini estimates
  vector<lower=0>[N] gini_se;             // standard errors for gini estimates
  vector<lower=0>[N] pop_size;            // population size estimates
  vector[N] start_year;                   // start year
  vector[N] end_year;                     // end year
}
parameters {
  real<lower=0, upper=1> phi;             // location for gini 
  real<lower=1e-12> tau;                  // scale for gini
  vector<lower=0, upper=1>[N] gini_true;  // unknown true gini values
  real alpha;                             // intercept for population size
  real beta;                              // slope for population size
  real<lower=1e-12, upper=10> sigma;      // scale for population size
  vector<lower=0, upper=1>[N] year_raw;   // unknown true year (scaled 0-1)
}
transformed parameters {
  vector[N] year_true;
  for (n in 1:N) {
    year_true[n] = start_year[n] + year_raw[n] * (end_year[n] - start_year[n]);
  }
}
model {
  // priors
  phi ~ normal(0.5, 0.05);
  tau ~ exponential(10);
  alpha ~ normal(10, 2);
  beta ~ normal(0, 1);
  sigma ~ exponential(2);
  year_raw ~ uniform(0, 1);
  // gini model
  gini_true ~ normal(phi, tau);
  gini_obs ~ normal(gini_true, gini_se);
  // population size model
  pop_size ~ lognormal(alpha + beta * (year_true / 1000), sigma);
}
