functions {
  vector ode(real time,      // time
             vector state,   // states
             vector theta) { // parameters
    // states
    real logP = state[1]; // log(pop_size)
    // parameters
    real rP = exp(theta[1]); // rate of population growth 
    // differential equation
    real dlogP = rP; // dP = rP * P
    return to_vector({dlogP});
  }
}
data {
  int N;                                         // total number of records
  int<lower=1, upper=N> N_dates;                 // number of unique dates
  int<lower=1, upper=N> N_regions;               // number of unique regions
  int<lower=1, upper=N> N_obs_pop;               // number of observed pop_size
  array[N_dates] real date;                      // dates (in millenia)
  array[N] int<lower=1, upper=N_regions> region; // region ids
  array[N_obs_pop] real<lower=0> pop_size;       // population size
  array[N] int date_idx;                         // link records to dates
  array[N_obs_pop] int pop_idx;                  // link pop_size to records
}
parameters {
  real init_pop_size;            // initial state for population size (log)
  vector[1] theta;               // ode parameters
  array[2] real<lower=0> tau;    // region SDs
  array[2] vector[N_regions] z;  // region-specific effects
  real<lower=0> sigma;           // lognormal variance for population size
}
transformed parameters{
  // construct region-specific effects
  vector[N_regions] init_pop_size_r = init_pop_size + (tau[1] * z[1]);
  array[N_regions] vector[1] theta_r;
  for (r in 1:N_regions) {
    theta_r[r][1] = theta[1] + (tau[2] * z[2][1]);
  }
  
  // solve ode
  array[N_regions, N_dates] vector[1] latent;
  for (r in 1:N_regions) {
    latent[r, 1][1] = init_pop_size_r[r];
    latent[r, 2:N_dates] = ode_rk45(
      ode, to_vector(latent[r, 1]), date[1], date[2:N_dates], theta_r[r]
    );
  }
}
model {
  // priors
  init_pop_size ~ normal(0, 1);
  theta ~ normal(0, 1);
  tau ~ exponential(2);
  for (i in 1:2) {
    z[i] ~ normal(0, 1);
  }
  sigma ~ exponential(2);
  
  // likelihood
  for (i in 1:N_obs_pop) {
    pop_size[i] ~ lognormal(latent[region[pop_idx[i]], date_idx[i]][1], sigma);
  }
}
