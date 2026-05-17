data {
  
  int<lower=1> N_sites;                    // number of sites
  vector<lower=0, upper=1>[N_sites] gini;  // gini estimates
  vector[N_sites] pop_scaled;              // population size (scaled)
  int<lower=0, upper=1> prior_only;        // sample the prior only?
  
}

parameters {
  
  real alpha;          // gini intercept
  real<lower=0> phi;   // gini precision parameter
  real theta;          // population size intercept
  real<lower=0> shape; // population size shape parameter
  
}

model {
  
  // ─────────────────────────────────────────
  // Priors
  // ─────────────────────────────────────────
  
  alpha ~ normal(0, 1);
  phi ~ lognormal(3, 1);
  theta ~ normal(0, 1);
  shape ~ exponential(1);
  
  // ─────────────────────────────────────────
  // Likelihood
  // ─────────────────────────────────────────
  
  if (!prior_only) {
    
    for (n in 1:N_sites) {
    
      // ─────────────────────────────────────────
      // Gini
      // ─────────────────────────────────────────
      
      real mu;
      real shape1;
      real shape2;
    
      mu = inv_logit(alpha);
      shape1 = mu * phi;
      shape2 = (1.0 - mu) * phi;
      
      gini[n] ~ beta(shape1, shape2);
      
      // ─────────────────────────────────────────
      // Log population size
      // ─────────────────────────────────────────
      
      if (pop_scaled[n] != -9999) {
        pop_scaled[n] ~ gamma(shape, shape / exp(theta));
      }
    
    }
    
  }
  
}

generated quantities {

  vector[N_sites] gini_rep;
  vector[N_sites] pop_scaled_rep;
  vector[N_sites] log_lik;

  for (n in 1:N_sites) {
    
    // ─────────────────────────────────────────
    // Posterior predictive simulation
    // ─────────────────────────────────────────
    
    real mu;
    real shape1;
    real shape2;
    
    mu = inv_logit(alpha);
    shape1 = mu * phi;
    shape2 = (1.0 - mu) * phi;
  
    gini_rep[n] = beta_rng(shape1, shape2);
    
    pop_scaled_rep[n] = gamma_rng(shape, shape / exp(theta));

    // ─────────────────────────────────────────
    // Pointwise log-likelihood for Gini
    // ─────────────────────────────────────────
    
    log_lik[n] = beta_lpdf(gini[n] | shape1, shape2);
    
  }
}
