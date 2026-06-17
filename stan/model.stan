functions {
  vector ode(real time,      // time
             vector state,   // states {y = gini, x = pop size}
             vector theta) { // parameters
    // states
    real y = state[1];
    real x = state[2];
    // parameters
    real k     = theta[1];
    real alpha = theta[2];
    real beta  = theta[3];
    real r     = theta[4];
    // differential equations
    real dy_dt = -exp(k) * (y - (alpha + (beta * log(x))));
    real dx_dt = exp(r) * x;
    return to_vector({dy_dt, dx_dt});
  }
}
data {
  int<lower=0> N;                             // total number of sites
  int<lower=0> N_times;                       // number of unique time points
  array[N] real<lower=0, upper=1> gini;       // gini
  array[N] real<lower=0> pop_size;            // population size
  array[N_times] real<lower=0, upper=1> t;    // unique time points (0-1)
  array[N] int<lower=1, upper=N_times> t_idx; // index for time points
}
parameters {
  real init_gini;      // initial state for gini
  real init_pop_size;  // initial state for population size
  vector[4] theta;     // ode parameters
  real<lower=0> phi;   // beta precision
  real<lower=0> sigma; // lognormal variance
}
transformed parameters{
  array[N_times] vector[2] latent;
  latent[1][1] = init_gini;
  latent[1][2] = init_pop_size;
  latent[2:N_times] = ode_rk45(
    ode, to_vector(latent[1, ]), 0.0, t[2:N_times], theta
  );
}
model {
  // priors
  init_gini ~ normal(0, 1);
  init_pop_size ~ normal(0, 1);
  theta[1] ~ normal(0, 1);
  theta[2] ~ normal(0, 1);
  theta[3] ~ normal(0, 1);
  theta[4] ~ normal(0, 1);
  phi ~ exponential(1);
  sigma ~ exponential(1);
  
  // likelihood
  for (i in 1:N) {
    real mu = inv_logit(latent[t_idx[i], 1]);
    real nu = latent[t_idx[i], 2];
    gini[i] ~ beta(mu * phi, (1 - mu) * phi);
    pop_size[i] ~ lognormal(nu, sigma);
  }
}
generated quantities {
  array[N] real gini_rep;
  array[N] real pop_size_rep;
  for (i in 1:N) {
    real mu = inv_logit(latent[t_idx[i], 1]);
    real nu = latent[t_idx[i], 2];
    gini_rep[i] = beta_rng(mu * phi, (1 - mu) * phi);
    pop_size_rep[i] = lognormal_rng(nu, sigma);
  }
}
