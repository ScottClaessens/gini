data {
  int<lower=0> N;                          // number of sites
  vector<lower=0, upper=1>[N] gini_lower;  // lower bound for gini estimate
  vector<lower=0, upper=1>[N] gini_upper;  // upper bound for gini estimate
  vector<lower=0, upper=1>[N] time_lower;  // lower bound for site date
  vector<lower=0, upper=1>[N] time_upper;  // upper bound for site date
}

parameters {
  real alpha;                                          // mean for gini
  real<lower=0> phi;                                   // dispersion for gini
  vector<lower=gini_lower, upper=gini_upper>[N] gini;  // true gini values
  vector<lower=time_lower, upper=time_upper>[N] time;  // true site dates
}

model {
  
  // initialise
  real mu;
  real shape1;
  real shape2;
  
  // priors
  alpha ~ normal(0, 1);
  phi ~ normal(10, 1);
  
  // gini model
  mu = inv_logit(alpha);
  
  shape1 = mu * phi + 1e-06;
  shape2 = (1.0 - mu) * phi + 1e-06;
  
  gini ~ beta(shape1, shape2);
  
}

generated quantities {

  real mu;
  real shape1;
  real shape2;

  vector[N] gini_rep;
  vector[N] log_lik;

  // recompute transformed quantities
  mu = inv_logit(alpha);

  shape1 = mu * phi + 1e-6;
  shape2 = (1.0 - mu) * phi + 1e-6;

  for (n in 1:N) {
    
    // posterior predictive simulation
    gini_rep[n] = beta_rng(shape1, shape2);

    // pointwise log-likelihood
    log_lik[n] = beta_lpdf(gini[n] | shape1, shape2);
    
  }
}
