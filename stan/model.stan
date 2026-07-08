functions {
  vector ode(real date,      // date
             vector state,   // states
             vector theta) { // parameters
    // states
    real logitP = state[1]; // logit( Pr(pop_size > 0) )
    real logP   = state[2]; // log(pop_size)
    real logitC = state[3]; // logit( Pr(cropland > 0) )
    real logC   = state[4]; // log(cropland)
    real logitG = state[5]; // logit(gini)
    // parameters
    real bP = exp(theta[1]); // increase in probability of positive population
    real rP = exp(theta[2]); // rate of population growth
    real bC = exp(theta[3]); // increase in probability of positive cropland
    real rC = exp(theta[4]); // rate of cropland production
    // differential equations
    real dlogitP = bP;                    // dlogitP/dt = bP
    real dlogP   = rP;                    // dP/dt      = rP * P
    real dlogitC = bC;                    // dlogitC/dt = bC
    real dlogC   = rC * exp(logP - logC); // dC/dt      = rC * P
    real dlogitG = 0;                     // dlogitG/dt = 0
    return to_vector({dlogitP, dlogP, dlogitC, dlogC, dlogitG});
  }
}
data {
  int N;                                         // total number of records
  int<lower=1, upper=N> N_dates;                 // number of unique dates
  int<lower=1, upper=N> N_regions;               // number of unique regions
  int<lower=1, upper=N> N_obs_pop;               // number of observed pop_size
  int<lower=1, upper=N> N_obs_crop;              // number of observed cropland
  int<lower=1, upper=N> N_obs_gini;              // number of observed gini
  array[N_dates] real date;                      // dates (in millenia)
  int i0;                                        // index for 0 CE
  int i1600;                                     // index for 1600 CE
  array[N] int<lower=1, upper=N_regions> region; // region ids
  array[N_obs_pop] real<lower=0> pop_size;       // population size
  array[N_obs_crop] real<lower=0> cropland;      // cropland
  array[N_obs_gini] real<lower=0, upper=1> gini; // gini
  array[N] int date_idx;                         // link records to dates
  array[N_obs_pop] int pop_idx;                  // link pop_size to records
  array[N_obs_crop] int crop_idx;                // link cropland to records
  array[N_obs_gini] int gini_idx;                // link gini to records
}
parameters {
  real init_logit_pop;           // initial state: logit prob of population > 0
  real init_pop_size;            // initial state: population size (log)
  real init_logit_crop;          // initial state: logit prob of cropland > 0
  real init_cropland;            // initial state: cropland (log)
  real init_gini;                // initial state: gini (logit)
  array[3] vector[4] theta;      // ode parameters over three time periods
  array[17] real<lower=0> tau;   // region SDs
  array[17] vector[N_regions] z; // region-specific effects
  real<lower=0> sigma;           // lognormal variance for population size
  real<lower=0> omega;           // lognormal variance for cropland
  real<lower=0> phi;             // beta precision for gini
}
transformed parameters{
  // construct region-specific effects
  vector[N_regions] init_logit_pop_r  = init_logit_pop  + (tau[1] * z[1]);
  vector[N_regions] init_pop_size_r   = init_pop_size   + (tau[2] * z[2]);
  vector[N_regions] init_logit_crop_r = init_logit_crop + (tau[3] * z[3]);
  vector[N_regions] init_cropland_r   = init_cropland   + (tau[4] * z[4]);
  vector[N_regions] init_gini_r       = init_gini       + (tau[5] * z[5]);
  array[3, N_regions] vector[4] theta_r;
  for (p in 1:3) {
    for (r in 1:N_regions) {
      for (k in 1:4) {
        int i = 5 + (p - 1) * 4 + k;
        theta_r[p, r][k] = theta[p][k] + (tau[i] * z[i][r]);
      }
    }
  }
  
  // solve ode
  array[N_regions, N_dates] vector[5] latent;
  for (r in 1:N_regions) {
    // initial values
    latent[r, 1][1] = init_logit_pop_r[r];
    latent[r, 1][2] = init_pop_size_r[r];
    latent[r, 1][3] = init_logit_crop_r[r];
    latent[r, 1][4] = init_cropland_r[r];
    latent[r, 1][5] = init_gini_r[r];
    // first period = 10,000 BCE to 0 CE
    latent[r, 2:i0] = ode_rk45(
      ode, latent[r, 1], date[1], date[2:i0], theta_r[1, r]
    );
    // second period = 0 CE to 1600 CE
    latent[r, (i0+1):i1600] = ode_rk45(
      ode, latent[r, i0], date[i0], date[(i0+1):i1600], theta_r[2, r]
    );
    // third period = 1600 CE to 1980 CE
    latent[r, (i1600+1):N_dates] = ode_rk45(
      ode, latent[r, i1600], date[i1600], date[(i1600+1):N_dates], theta_r[3, r]
    );
  }
}
model {
  // priors
  init_logit_pop ~ normal(-4, 1);
  init_pop_size ~ normal(-1, 1);
  init_logit_crop ~ normal(-4, 1);
  init_cropland ~ normal(-1, 1);
  init_gini ~ normal(-1, 1);
  for (p in 1:3) theta[p] ~ normal(-2, 0.5);
  for (i in 1:17) z[i] ~ normal(0, 1);
  tau ~ exponential(1);
  sigma ~ exponential(1);
  omega ~ exponential(1);
  phi ~ exponential(1);
  
  // likelihood for population size
  for (i in 1:N_obs_pop) {
    real logit_prob = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][1];
    if (pop_size[i] == 0) {
      target += bernoulli_logit_lpmf(0 | logit_prob);
    } else {
      real mu = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][2];
      target += bernoulli_logit_lpmf(1 | logit_prob);
      target += lognormal_lpdf(pop_size[i] | mu, sigma);
    }
  }
  
  // likelihood for cropland
  for (i in 1:N_obs_crop) {
    real logit_prob = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][3];
    if (cropland[i] == 0) {
      target += bernoulli_logit_lpmf(0 | logit_prob);
    } else {
      real nu = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][4];
      target += bernoulli_logit_lpmf(1 | logit_prob);
      target += lognormal_lpdf(cropland[i] | nu, omega);
    }
  }
  
  // likelihood for gini
  for (i in 1:N_obs_gini) {
    real mu = latent[region[gini_idx[i]], date_idx[gini_idx[i]]][5];
    target += beta_proportion_lpdf(gini[i] | inv_logit(mu), phi);
  }
}
generated quantities {
  array[N_obs_pop] real pop_size_rep;
  array[N_obs_crop] real cropland_rep;
  array[N_obs_gini] real gini_rep;
  array[121] real date_rep = linspaced_array(121, -100, 20);
  array[121] vector[5] global_latent_rep;
  array[N_regions, 121] vector[5] regional_latent_rep;
  
  // posterior predictive check for population size
  for (i in 1:N_obs_pop) {
    real logit_prob = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][1];
    if (bernoulli_logit_rng(logit_prob) == 0) {
      pop_size_rep[i] = 0;
    } else {
      real mu = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][2];
      pop_size_rep[i] = lognormal_rng(mu, sigma);
    }
  }
  
  // posterior predictive check for cropland
  for (i in 1:N_obs_crop) {
    real logit_prob = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][3];
    if (bernoulli_logit_rng(logit_prob) == 0) {
      cropland_rep[i] = 0;
    } else {
      real nu = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][4];
      cropland_rep[i] = lognormal_rng(nu, omega);
    }
  }
  
  // posterior predictive check for gini
  for (i in 1:N_obs_gini) {
    real mu = latent[region[gini_idx[i]], date_idx[gini_idx[i]]][5];
    gini_rep[i] = beta_proportion_rng(inv_logit(mu), phi);
  }
  
  // global ode prediction across three periods
  global_latent_rep[1][1] = init_logit_pop;
  global_latent_rep[1][2] = init_pop_size;
  global_latent_rep[1][3] = init_logit_crop;
  global_latent_rep[1][4] = init_cropland;
  global_latent_rep[1][5] = init_gini;
  global_latent_rep[2:101] = ode_rk45(
    ode, global_latent_rep[1], date_rep[1], date_rep[2:101], theta[1]
  );
  global_latent_rep[102:117] = ode_rk45(
    ode, global_latent_rep[101], date_rep[101], date_rep[102:117], theta[2]
  );
  global_latent_rep[118:121] = ode_rk45(
    ode, global_latent_rep[117], date_rep[117], date_rep[118:121], theta[3]
  );
  
  // regional ode prediction across three periods
  for (r in 1:N_regions) {
    regional_latent_rep[r, 1][1] = init_logit_pop_r[r];
    regional_latent_rep[r, 1][2] = init_pop_size_r[r];
    regional_latent_rep[r, 1][3] = init_logit_crop_r[r];
    regional_latent_rep[r, 1][4] = init_cropland_r[r];
    regional_latent_rep[r, 1][5] = init_gini_r[r];
    regional_latent_rep[r, 2:101] = ode_rk45(
      ode, regional_latent_rep[r, 1], 
      date_rep[1], date_rep[2:101], theta_r[1, r]
    );
    regional_latent_rep[r, 102:117] = ode_rk45(
      ode, regional_latent_rep[r, 101], 
      date_rep[101], date_rep[102:117], theta_r[2, r]
    );
    regional_latent_rep[r, 118:121] = ode_rk45(
      ode, regional_latent_rep[r, 117], 
      date_rep[117], date_rep[118:121], theta_r[3, r]
    );
  }
}
