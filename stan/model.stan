functions {
  vector ode(real date,      // date
             vector state,   // states
             vector theta) { // parameters
    // states
    real logP = state[1]; // log(pop_size)
    real logC = state[2]; // log(cropland)
    // parameters
    real rP = exp(theta[1]); // rate of population growth
    real rC = exp(theta[2]); // rate of cropland production
    // differential equations
    real dlogP;
    real dlogC;
    dlogP = rP;                    // dP/dt = rP * P
    dlogC = rC * exp(logP - logC); // dC/dt = rC * P
    return to_vector({dlogP, dlogC});
  }
}
data {
  int N;                                         // total number of records
  int<lower=1, upper=N> N_dates;                 // number of unique dates
  int<lower=1, upper=N> N_regions;               // number of unique regions
  int<lower=1, upper=N> N_obs_pop;               // number of observed pop_size
  int<lower=1, upper=N> N_obs_crop;              // number of observed cropland
  array[N_dates] real date;                      // dates (in millenia)
  int i0;                                        // index for 0 CE
  int i1600;                                     // index for 1600 CE
  array[N] int<lower=1, upper=N_regions> region; // region ids
  array[N_obs_pop] real<lower=0> pop_size;       // population size
  array[N_obs_crop] real<lower=0> cropland;      // population size
  array[N] int date_idx;                         // link records to dates
  array[N_obs_pop] int pop_idx;                  // link pop_size to records
  array[N_obs_crop] int crop_idx;                // link cropland to records
}
parameters {
  real init_pop_size;              // initial state for population size (log)
  real init_cropland;              // initial state for cropland (log)
  array[3] vector[2] theta;        // ode parameters over three time periods
  array[8] real<lower=0> tau;   // region SDs
  array[8] vector[N_regions] z; // region-specific effects
  real<lower=0> sigma;             // lognormal variance for population size
  real<lower=0> omega;             // lognormal variance for cropland
}
transformed parameters{
  // construct region-specific effects
  vector[N_regions] init_pop_size_r = init_pop_size + (tau[1] * z[1]);
  vector[N_regions] init_cropland_r = init_cropland + (tau[2] * z[2]);
  array[3, N_regions] vector[2] theta_r;
  for (p in 1:3) {
    for (r in 1:N_regions) {
      for (k in 1:2) {
        int i = 2 + (p - 1) * 2 + k;
        theta_r[p, r][k] = theta[p][k] + (tau[i] * z[i][r]);
      }
    }
  }
  
  // solve ode
  array[N_regions, N_dates] vector[2] latent;
  for (r in 1:N_regions) {
    // initial values
    latent[r, 1][1] = init_pop_size_r[r];
    latent[r, 1][2] = init_cropland_r[r];
    // first period = 10,000 BCE to 0 CE
    latent[r, 2:i0] = ode_rk45(
      ode, to_vector(latent[r, 1]), 
      date[1], date[2:i0], theta_r[1, r]
    );
    // second period = 0 CE to 1700 CE
    latent[r, (i0 + 1):i1600] = ode_rk45(
      ode, to_vector(latent[r, i0]),
      date[i0], date[(i0 + 1):i1600], theta_r[2, r]
    );
    // third period = 1700 CE to 1980 CE
    latent[r, (i1600 + 1):N_dates] = ode_rk45(
      ode, to_vector(latent[r, i1600]),
      date[i1600], date[(i1600 + 1):N_dates], theta_r[3, r]
    );
  }
}
model {
  // priors
  init_pop_size ~ normal(-1, 1);
  init_cropland ~ normal(-1, 1);
  for (p in 1:3) theta[p] ~ normal(-1, 1);
  for (i in 1:8) z[i] ~ normal(0, 1);
  tau ~ exponential(1);
  sigma ~ exponential(1);
  omega ~ exponential(1);
  
  // likelihood for population size
  for (i in 1:N_obs_pop) {
    real mu = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][1];
    pop_size[i] ~ lognormal(mu, sigma);
  }
  
  // likelihood for cropland
  for (i in 1:N_obs_crop) {
    real nu = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][2];
    cropland[i] ~ lognormal(nu, omega);
  }
}
generated quantities {
  array[N_obs_pop] real pop_size_rep;
  array[N_obs_crop] real cropland_rep;
  array[121] real date_rep = linspaced_array(121, -100, 20);
  array[121] vector[2] global_latent_rep;
  array[N_regions, 121] vector[2] regional_latent_rep;
  
  // posterior predictive check for population size
  for (i in 1:N_obs_pop) {
    pop_size_rep[i] = lognormal_rng(
      latent[region[pop_idx[i]], date_idx[pop_idx[i]]][1], sigma
    );
  }
  
  // posterior predictive check for cropland
  for (i in 1:N_obs_crop) {
    cropland_rep[i] = lognormal_rng(
      latent[region[crop_idx[i]], date_idx[crop_idx[i]]][2], omega
    );
  }
  
  // global ode prediction across three periods
  global_latent_rep[1][1] = init_pop_size;
  global_latent_rep[1][2] = init_cropland;
  global_latent_rep[2:101] = ode_rk45(
    ode, to_vector(global_latent_rep[1]), 
    date_rep[1], date_rep[2:101], theta[1]
  );
  global_latent_rep[102:117] = ode_rk45(
    ode, to_vector(global_latent_rep[101]), 
    date_rep[101], date_rep[102:117], theta[2]
  );
  global_latent_rep[118:121] = ode_rk45(
    ode, to_vector(global_latent_rep[117]), 
    date_rep[117], date_rep[118:121], theta[3]
  );
  
  // regional ode prediction across three periods
  for (r in 1:N_regions) {
    regional_latent_rep[r, 1][1] = init_pop_size_r[r];
    regional_latent_rep[r, 1][2] = init_cropland_r[r];
    regional_latent_rep[r, 2:101] = ode_rk45(
      ode, to_vector(regional_latent_rep[r, 1]), 
      date_rep[1], date_rep[2:101], theta_r[1, r]
    );
    regional_latent_rep[r, 102:117] = ode_rk45(
      ode, to_vector(regional_latent_rep[r, 101]), 
      date_rep[101], date_rep[102:117], theta_r[2, r]
    );
    regional_latent_rep[r, 118:121] = ode_rk45(
      ode, to_vector(regional_latent_rep[r, 117]), 
      date_rep[117], date_rep[118:121], theta_r[3, r]
    );
  }
}
