functions {
  vector ode(real time,      // time
             vector state,   // states
             vector theta) { // parameters
    // states
    real G = state[1]; // G = inv_logit(gini)
    real P = state[2]; // P = population size
    real C = state[3]; // C = cropland
    // parameters
    real rP = exp(theta[1]);                // rate of population change
    real K  = exp(theta[2] + theta[3] * C); // additional carrying capacity
    real rC = exp(theta[4]);                // rate of cropland production
    real rG = exp(theta[5]);                // rate of inequality mean reversion
    real mu = theta[6] +                    // equilibrium mean for inequality
      (theta[7] * log(P)) + 
      (theta[8] * log(C));
    // differential equations
    real dG = rG * (mu - G);
    real dP = rP * P * (1 - (P / (P + K)));
    real dC = rC * P;
    return to_vector({dG, dP, dC});
  }
}
data {
  int<lower=0> N;                                 // total number of sites
  int<lower=0> N_times;                           // n unique time points
  array[N] real<lower=0, upper=1> gini;           // gini
  array[N] real<lower=0> pop_size;                // population size
  array[N] real<lower=0> cropland;                // cropland
  array[N_times] real<lower=0, upper=1> t;        // unique time points (0-1)
  array[N] int<lower=1, upper=N_times> t_idx;     // index for time points
}
parameters {
  real init_gini;              // initial state for gini
  real<lower=0> init_pop_size; // initial state for population size
  real<lower=0> init_cropland; // initial state for cropland
  vector[8] theta;             // ode parameters
  real<lower=0> phi;           // beta precision for gini
  real<lower=0> sigma;         // lognormal variance for population size
  real<lower=0> omega;         // lognormal variance for cropland
}
transformed parameters{
  array[N_times] vector[3] latent;
  latent[1][1] = init_gini;
  latent[1][2] = init_pop_size;
  latent[1][3] = init_cropland;
  latent[2:N_times] = ode_rk45(
    ode, to_vector(latent[1, ]), 0, t[2:N_times], theta
  );
}
model {
  // priors
  init_gini ~ normal(0, 1);
  init_pop_size ~ lognormal(0, 1);
  init_cropland ~ lognormal(0, 1);
  theta ~ normal(0, 1);
  phi ~ exponential(1);
  sigma ~ exponential(1);
  omega ~ exponential(1);
  
  // likelihood
  for (i in 1:N) {
    
    // gini
    real mu = inv_logit(latent[t_idx[i], 1]);
    gini[i] ~ beta(mu * phi, (1 - mu) * phi);
    
    // population size
    pop_size[i] ~ lognormal(log(latent[t_idx[i], 2]), sigma);
    
    // cropland
    cropland[i] ~ lognormal(log(latent[t_idx[i], 3]), omega);
    
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
    real mu = inv_logit(latent[t_idx[i], 1]);
    gini_rep[i] = beta_rng(mu * phi, (1 - mu) * phi);
    pop_size_rep[i] = lognormal_rng(log(latent[t_idx[i], 2]), sigma);
    cropland_rep[i] = lognormal_rng(log(latent[t_idx[i], 3]), omega);
  }
  
  // smooth ode prediction over 0-1 range
  latent_rep[1][1] = init_gini;
  latent_rep[1][2] = init_pop_size;
  latent_rep[1][3] = init_cropland;
  latent_rep[2:100] = ode_rk45(
    ode, to_vector(latent_rep[1, ]), 0, t_rep[2:100], theta
  );
}
