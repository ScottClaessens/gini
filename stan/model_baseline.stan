data {
  
  int<lower=1> N_sites;                          // number of sites
  vector<lower=0, upper=1>[N_sites] gini_lower;  // lower bound for gini
  vector<lower=0, upper=1>[N_sites] gini_upper;  // upper bound for gini
  
}

parameters {
  
  // beta parameters
  real alpha;
  real<lower=0> phi;
  
  // latent "true" gini values
  vector<lower=gini_lower, upper=gini_upper>[N_sites] gini;
  
}

model {
  
  real mu;
  real shape1;
  real shape2;
  
  // ─────────────────────────────────────────
  // Priors
  // ─────────────────────────────────────────
  
  alpha ~ normal(0, 1);
  phi ~ normal(10, 1);
  
  // ─────────────────────────────────────────
  // Beta likelihood
  // ─────────────────────────────────────────
  
  mu = inv_logit(alpha);
  
  shape1 = mu * phi + 1e-06;
  shape2 = (1.0 - mu) * phi + 1e-06;
  
  gini ~ beta(shape1, shape2);
  
}

generated quantities {

  real mu;
  real shape1;
  real shape2;

  vector[N_sites] gini_rep;
  vector[N_sites] log_lik;

  // ─────────────────────────────────────────
  // Recompute transformed quantities
  // ─────────────────────────────────────────
  
  mu = inv_logit(alpha);

  shape1 = mu * phi + 1e-6;
  shape2 = (1.0 - mu) * phi + 1e-6;

  for (n in 1:N_sites) {
    
    // ─────────────────────────────────────────
    // Posterior predictive simulation
    // ─────────────────────────────────────────
    
    gini_rep[n] = beta_rng(shape1, shape2);

    // ─────────────────────────────────────────
    // Pointwise log-likelihood
    // ─────────────────────────────────────────
    
    log_lik[n] = beta_lpdf(gini[n] | shape1, shape2);
    
  }
}
