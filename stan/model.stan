functions {
  vector ode(real time,      // time
             vector state,   // states
             vector theta) { // parameters
    // states
    real G    = state[1];  // inv_logit(gini)
    real logP = state[2];  // log(pop_size)
    real logC = state[3];  // log(cropland)
    // parameters
    real rG = exp(theta[1]);          // rate of inequality mean reversion
    real rP = exp(theta[2]);          // rate of population growth 
    real rC = exp(theta[3]);          // rate of cropland production
    real mu = theta[4] +              // equilibrium for inequality
      (theta[5] * log1p(exp(logP))) +    
      (theta[6] * log1p(exp(logC)));
    // differential equations
    real dG = rG * (mu - G);
    real dlogP = rP;                     // dP = rP * P
    real dlogC = rC * exp(logP - logC);  // dC = rC * P
    return to_vector({dG, dlogP, dlogC});
  }
}
data {
  int<lower=0> N;                                // total number of sites
  int<lower=0, upper=N> N_regions;               // number of regions
  int<lower=0, upper=N> N_times;                 // number of unique time points
  array[N] real<lower=0, upper=1> gini;          // gini
  array[N] real<lower=0> pop_size;               // population size
  array[N] real<lower=0> cropland;               // cropland
  array[N] int<lower=1, upper=N_regions> region; // region id
  array[N_times] real<lower=0, upper=1> t;       // unique time points (0-1)
  array[N] int<lower=1, upper=N_times> t_idx;    // index for time points
}
parameters {
  real init_gini;                // initial state for gini (logit scale)
  real init_pop_size;            // initial state for population size (log scale)
  real init_cropland;            // initial state for cropland (log scale)
  vector[6] theta;               // ode parameters
  array[9] real<lower=0> tau;    // region SDs
  array[9] vector[N_regions] z;  // region-specific effects
  real<lower=0> phi;             // beta precision for gini
  real<lower=0> sigma;           // lognormal variance for population size
  real<lower=0> omega;           // lognormal variance for cropland
}
transformed parameters{
  // construct region-specific effects
  vector[N_regions] init_gini_r     = init_gini     + (tau[1] * z[1]);
  vector[N_regions] init_pop_size_r = init_pop_size + (tau[2] * z[2]);
  vector[N_regions] init_cropland_r = init_cropland + (tau[3] * z[3]);
  array[N_regions] vector[6] theta_r;
  for (r in 1:N_regions) {
    for (k in 1:6) {
      theta_r[r][k] = theta[k] + (tau[k + 3] * z[k + 3][r]);
    }
  }
  
  // solve ode
  array[N_regions, N_times] vector[3] latent;
  for (r in 1:N_regions) {
    latent[r, 1, 1] = init_gini_r[r];
    latent[r, 1, 2] = init_pop_size_r[r];
    latent[r, 1, 3] = init_cropland_r[r];
    latent[r, 2:N_times] = ode_rk45(
      ode, to_vector(latent[r, 1, ]), 0, t[2:N_times], theta_r[r]
    );
  }
}
model {
  // priors for global parameters
  init_gini ~ normal(-1, 0.5);
  init_pop_size ~ normal(0, 1);
  init_cropland ~ normal(-5, 1);
  theta ~ normal(0, 1);
  
  // priors for region-specific parameters
  for (i in 1:9) {
    tau[i] ~ exponential(2);
    z[i] ~ normal(0, 1);
  }
  
  // priors for measurement error
  phi ~ exponential(2);
  sigma ~ exponential(2);
  omega ~ exponential(2);
  
  // likelihood
  for (i in 1:N) {
    gini[i] ~ beta_proportion(inv_logit(latent[region[i], t_idx[i], 1]), phi);
    pop_size[i] ~ lognormal(latent[region[i], t_idx[i], 2], sigma);
    cropland[i] ~ lognormal(latent[region[i], t_idx[i], 3], omega);
  }
}
generated quantities {
  array[N] real gini_rep;
  array[N] real pop_size_rep;
  array[N] real cropland_rep;
  array[100] real t_rep = linspaced_array(100, 0, 1);
  array[100] vector[3] global_latent_rep;
  array[N_regions, 100] vector[3] regional_latent_rep;
  
  // posterior predictive checks
  for (i in 1:N) {
    gini_rep[i] = beta_proportion_rng(
      inv_logit(latent[region[i], t_idx[i], 1]), phi
    );
    pop_size_rep[i] = lognormal_rng(latent[region[i], t_idx[i], 2], sigma);
    cropland_rep[i] = lognormal_rng(latent[region[i], t_idx[i], 3], omega);
  }
  
  // global ode prediction over 0-1 range
  global_latent_rep[1, 1] = init_gini;
  global_latent_rep[1, 2] = init_pop_size;
  global_latent_rep[1, 3] = init_cropland;
  global_latent_rep[2:100] = ode_rk45(
    ode, to_vector(global_latent_rep[1, ]), 0, t_rep[2:100], theta
  );
  
  // regional ode prediction over 0-1 range
  for (r in 1:N_regions) {
    regional_latent_rep[r, 1, 1] = init_gini_r[r];
    regional_latent_rep[r, 1, 2] = init_pop_size_r[r];
    regional_latent_rep[r, 1, 3] = init_cropland_r[r];
    regional_latent_rep[r, 2:100] = ode_rk45(
      ode, to_vector(regional_latent_rep[r, 1, ]), 0, t_rep[2:100], theta_r[r]
    );
  }
}
