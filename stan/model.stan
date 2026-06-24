functions {
  vector ode(real time,      // time
             vector state,   // states
             vector theta) { // parameters
    // states
    real G = state[1];  // G = inv_logit(gini)
    real P = state[2];  // P = population size
    real C = state[3];  // C = cropland
    // parameters
    real rG = exp(theta[1]);     // rate of inequality mean reversion
    real rP = exp(theta[2]);     // rate of population increase
    real rC = exp(theta[3]);     // rate of cropland production
    real K  = exp(theta[4]) * C; // population carrying capacity
    real mu = theta[5] +         // equilibrium for inequality
      (theta[6] * log(P + 1)) + 
      (theta[7] * log(C + 1));
    // differential equations
    real dG = rG * (mu - G);
    real dP = rP * P * (1 - (P / K));
    real dC = rC * P;
    return to_vector({dG, dP, dC});
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
  real init_gini;               // initial state for gini (logit scale)
  real init_pop_size;           // initial state for population size (log scale)
  real init_cropland;           // initial state for cropland (log scale)
  vector[7] theta;              // ode parameters
  array[3] real<lower=0> tau;   // region SDs
  array[3] vector[N_regions] z; // region-specific effects
  real<lower=0> phi;            // beta precision for gini
  real<lower=0> sigma;          // lognormal variance for population size
  real<lower=0> omega;          // lognormal variance for cropland
}
transformed parameters{
  // construct region-specific effects
  vector[N_regions] init_gini_r     = init_gini     + (tau[1] * z[1]);
  vector[N_regions] init_pop_size_r = init_pop_size + (tau[2] * z[2]);
  vector[N_regions] init_cropland_r = init_cropland + (tau[3] * z[3]);
  
  // solve ode
  array[N_regions, N_times] vector[3] latent;
  for (r in 1:N_regions) {
    latent[r, 1, 1] = init_gini_r[r];
    latent[r, 1, 2] = exp(init_pop_size_r[r]);
    latent[r, 1, 3] = exp(init_cropland_r[r]);
    latent[r, 2:N_times] = ode_rk45(
      ode, to_vector(latent[r, 1, ]), 0, t[2:N_times], theta
    );
  }
}
model {
  // priors for global parameters
  init_gini ~ normal(-1, 0.5);
  init_pop_size ~ normal(0, 1);
  init_cropland ~ normal(-5, 1);
  theta[{1,2,3,6,7}] ~ normal(0, 1);
  theta[5] ~ normal(0, 0.2);
  theta[4] ~ normal(6.5, 0.5);
  
  // priors for region-specific parameters
  for (i in 1:3) {
    tau[i] ~ exponential(1);
    z[i] ~ normal(0, 1);
  }
  
  // priors for measurement error
  phi ~ exponential(1);
  sigma ~ exponential(1);
  omega ~ exponential(1);
  
  // likelihood
  for (i in 1:N) {
    gini[i] ~ beta_proportion(inv_logit(latent[region[i], t_idx[i], 1]), phi);
    pop_size[i] ~ lognormal(log(latent[region[i], t_idx[i], 2]), sigma);
    cropland[i] ~ lognormal(log(latent[region[i], t_idx[i], 3]), omega);
  }
}
generated quantities {
  array[N] real gini_rep;
  array[N] real pop_size_rep;
  array[N] real cropland_rep;
  array[100] real t_rep = linspaced_array(100, 0, 1);
  array[100] vector[3] latent_rep;
  
  // posterior predictive checks
  for (i in 1:N) {
    gini_rep[i] = beta_proportion_rng(
      inv_logit(latent[region[i], t_idx[i], 1]), phi
    );
    pop_size_rep[i] = lognormal_rng(
      log(latent[region[i], t_idx[i], 2]), sigma
    );
    cropland_rep[i] = lognormal_rng(
      log(latent[region[i], t_idx[i], 3]), omega
    );
  }
  
  // smooth ode prediction over 0-1 range
  latent_rep[1, 1] = init_gini;
  latent_rep[1, 2] = exp(init_pop_size);
  latent_rep[1, 3] = exp(init_cropland);
  latent_rep[2:100] = ode_rk45(
    ode, to_vector(latent_rep[1, ]), 0, t_rep[2:100], theta
  );
}
