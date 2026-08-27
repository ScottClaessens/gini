functions {
  // construct cholesky for gaussian process
  matrix construct_gp_cholesky(array[] vector coords,
                               real tau,
                               real lscale) {
    int N = dims(coords)[1];
    matrix[N, N] K = gp_exp_quad_cov(coords, tau, lscale);
    matrix[N, N] L = cholesky_decompose(K + diag_matrix(rep_vector(1e-8, N)));
    return(L);
  }
  
  // ode function
  vector ode(real date,      // date
             vector state,   // states
             vector theta) { // parameters
    // states
    real logitP = state[1];     // logit( Pr(pop_size > 0) )
    real logP   = state[2];     // log(pop_size)
    real logitC = state[3];     // logit( Pr(cropland > 0) )
    real logC   = state[4];     // log(cropland)
    real logitI = state[5];     // logit( Pr(irrigated > 0) )
    real logI   = state[6];     // log(irrigated)
    real logitU = state[7];     // logit( Pr(urban > 0) )
    real logU   = state[8];     // log(urban)
    real logitG = state[9];     // logit(gini)
    // parameters
    real bP = exp(theta[1]);     // increase in probability of population > 0
    real rP = exp(theta[2]);     // rate of population growth
    real bC = exp(theta[3]);     // increase in probability of cropland > 0
    real rC = exp(theta[4]);     // rate of cropland production
    real bI = exp(theta[5]);     // increase in probability of irrigated > 0
    real rI = exp(theta[6]);     // rate of irrigated area growth
    real bU = exp(theta[7]);     // increase in probability of urban > 0
    real rU = exp(theta[8]);     // rate of urban area production
    real alpha = theta[9];       // continuous-time intercept for gini
    real betaP = theta[10];      // effect of population size on gini
    real betaC = theta[11];      // effect of cropland on gini
    real betaI = theta[12];      // effect of irrigated area on gini
    real betaU = theta[13];      // effect of urban area on gini
    real gamma = exp(theta[14]); // damping parameter
    // transformed parameters
    real rG =
      alpha
      + betaP * (inv_logit(logitP) * log1p_exp(logP))
      + betaC * (inv_logit(logitC) * log1p_exp(logC))
      + betaI * (inv_logit(logitI) * log1p_exp(logI))
      + betaU * (inv_logit(logitU) * log1p_exp(logU))
      - gamma * logitG;
    // differential equations
    real dlogitP = bP;                    // dlogitP/dt = bP
    real dlogP   = rP;                    // dP/dt      = rP * P
    real dlogitC = bC;                    // dlogitC/dt = bC
    real dlogC   = rC * exp(logP - logC); // dC/dt      = rC * P
    real dlogitI = bI;                    // dlogitI/dt = bI
    real dlogI   = rI * exp(logP - logI); // dI/dt      = rI * P
    real dlogitU = bU;                    // dlogitU/dt = bU
    real dlogU   = rU * exp(logP - logU); // dU/dt      = rU * P
    real dlogitG = rG;                    // dlogitG/dt = rG
    return to_vector({dlogitP, dlogP, dlogitC, dlogC, dlogitI, dlogI,
                      dlogitU, dlogU, dlogitG});
  }
}
data {
  int N;                                         // total number of records
  int<lower=1, upper=N> N_dates;                 // number of unique dates
  int<lower=1, upper=N> N_regions;               // number of unique regions
  int<lower=1, upper=N> N_obs_pop;               // number of observed pop_size
  int<lower=1, upper=N> N_obs_crop;              // number of observed cropland
  int<lower=1, upper=N> N_obs_irr;               // number of observed irrigated
  int<lower=1, upper=N> N_obs_urb;               // number of observed urban
  int<lower=1, upper=N> N_obs_gini;              // number of observed gini
  array[N_dates] real date;                      // dates (in centuries)
  int i0;                                        // index for 0 CE
  int i1600;                                     // index for 1600 CE
  array[N] int<lower=1, upper=N_regions> region; // region ids
  array[N_obs_pop] real<lower=0> pop_size;       // population size
  array[N_obs_crop] real<lower=0> cropland;      // cropland
  array[N_obs_irr] real<lower=0> irrigated;      // irrigated area
  array[N_obs_urb] real<lower=0> urban;          // urban area
  array[N_obs_gini] real<lower=0, upper=1> gini; // gini
  array[N] int date_idx;                         // link records to dates
  array[N_obs_pop] int pop_idx;                  // link pop_size to records
  array[N_obs_crop] int crop_idx;                // link cropland to records
  array[N_obs_irr] int irr_idx;                  // link irrigated to records
  array[N_obs_urb] int urb_idx;                  // link urban to records
  array[N_obs_gini] int gini_idx;                // link gini to records
  array[N_regions] vector[3] coords;             // x,y,z region coordinates
}
parameters {
  real init_logit_pop;           // initial state: logit prob of population > 0
  real init_pop_size;            // initial state: population size (log)
  real init_logit_crop;          // initial state: logit prob of cropland > 0
  real init_cropland;            // initial state: cropland (log)
  real init_logit_irr;           // initial state: logit prob of irrigated > 0
  real init_irrigated;           // initial state: irrigated (log)
  real init_logit_urb;           // initial state: logit prob of urban > 0
  real init_urban;               // initial state: urban (log)
  real init_gini;                // initial state: gini (logit)
  array[3] vector[14] theta;     // ode parameters over three time periods
  array[38] real<lower=0> tau;   // region SDs
  array[38] vector[N_regions] z; // region-specific effects
  array[6] real<lower=0> lscale; // length-scale for spatial GPs
  array[4] real<lower=0> sigma;  // lognormal variances
  real<lower=0> phi;             // beta precision for gini
}
transformed parameters{
  // construct gp vectors
  vector[N_regions] init_gini_gp; 
  vector[N_regions] alpha_gp; 
  vector[N_regions] betaP_gp; 
  vector[N_regions] betaC_gp; 
  vector[N_regions] betaI_gp; 
  vector[N_regions] betaU_gp; 
  init_gini_gp = construct_gp_cholesky(coords, tau[1], lscale[1]) * z[1];
  alpha_gp     = construct_gp_cholesky(coords, tau[2], lscale[2]) * z[2];
  betaP_gp     = construct_gp_cholesky(coords, tau[3], lscale[3]) * z[3];
  betaC_gp     = construct_gp_cholesky(coords, tau[4], lscale[4]) * z[4];
  betaI_gp     = construct_gp_cholesky(coords, tau[5], lscale[5]) * z[5];
  betaU_gp     = construct_gp_cholesky(coords, tau[6], lscale[6]) * z[6];
  
  // construct region-specific initial values
  vector[N_regions] init_logit_pop_r  = init_logit_pop  + (tau[7] * z[7]);
  vector[N_regions] init_pop_size_r   = init_pop_size   + (tau[8] * z[8]);
  vector[N_regions] init_logit_crop_r = init_logit_crop + (tau[9] * z[9]);
  vector[N_regions] init_cropland_r   = init_cropland   + (tau[10] * z[10]);
  vector[N_regions] init_logit_irr_r  = init_logit_irr  + (tau[11] * z[11]);
  vector[N_regions] init_irrigated_r  = init_irrigated  + (tau[12] * z[12]);
  vector[N_regions] init_logit_urb_r  = init_logit_urb  + (tau[13] * z[13]);
  vector[N_regions] init_urban_r      = init_urban      + (tau[14] * z[14]);
  vector[N_regions] init_gini_r       = init_gini       + init_gini_gp;
  
  // construct region-specific ode parameters
  array[3, N_regions] vector[14] theta_r;
  for (r in 1:N_regions) {
    // period 1
    theta_r[1, r][1]  = theta[1][1]  + (tau[15] * z[15][r]);  // bP
    theta_r[1, r][2]  = theta[1][2]  + (tau[16] * z[16][r]);  // rP
    theta_r[1, r][3]  = theta[1][3]  + (tau[17] * z[17][r]);  // bC
    theta_r[1, r][4]  = theta[1][4]  + (tau[18] * z[18][r]);  // rC
    theta_r[1, r][5]  = theta[1][5]  + (tau[19] * z[19][r]);  // bI
    theta_r[1, r][6]  = theta[1][6]  + (tau[20] * z[20][r]);  // rI
    theta_r[1, r][7]  = theta[1][7]  + (tau[21] * z[21][r]);  // bU
    theta_r[1, r][8]  = theta[1][8]  + (tau[22] * z[22][r]);  // rU
    theta_r[1, r][9]  = theta[1][9]  + alpha_gp[r];           // alpha
    theta_r[1, r][10] = theta[1][10] + betaP_gp[r];           // betaP
    theta_r[1, r][11] = theta[1][11] + betaC_gp[r];           // betaC
    theta_r[1, r][12] = theta[1][12] + betaI_gp[r];           // betaI
    theta_r[1, r][13] = theta[1][13] + betaU_gp[r];           // betaU
    theta_r[1, r][14] = theta[1][14];                         // gamma
    // period 2
    theta_r[2, r][1]  = theta[2][1]  + (tau[23] * z[23][r]);  // bP
    theta_r[2, r][2]  = theta[2][2]  + (tau[24] * z[24][r]);  // rP
    theta_r[2, r][3]  = theta[2][3]  + (tau[25] * z[25][r]);  // bC
    theta_r[2, r][4]  = theta[2][4]  + (tau[26] * z[26][r]);  // rC
    theta_r[2, r][5]  = theta[2][5]  + (tau[27] * z[27][r]);  // bI
    theta_r[2, r][6]  = theta[2][6]  + (tau[28] * z[28][r]);  // rI
    theta_r[2, r][7]  = theta[2][7]  + (tau[29] * z[29][r]);  // bU
    theta_r[2, r][8]  = theta[2][8]  + (tau[30] * z[30][r]);  // rU
    theta_r[2, r][9]  = theta[1][9]  + alpha_gp[r];           // alpha
    theta_r[2, r][10] = theta[1][10] + betaP_gp[r];           // betaP
    theta_r[2, r][11] = theta[1][11] + betaC_gp[r];           // betaC
    theta_r[2, r][12] = theta[1][12] + betaI_gp[r];           // betaI
    theta_r[2, r][13] = theta[1][13] + betaU_gp[r];           // betaU
    theta_r[2, r][14] = theta[1][14];                         // gamma
    // period 3
    theta_r[3, r][1]  = theta[3][1]  + (tau[31] * z[31][r]);  // bP
    theta_r[3, r][2]  = theta[3][2]  + (tau[32] * z[32][r]);  // rP
    theta_r[3, r][3]  = theta[3][3]  + (tau[33] * z[33][r]);  // bC
    theta_r[3, r][4]  = theta[3][4]  + (tau[34] * z[34][r]);  // rC
    theta_r[3, r][5]  = theta[3][5]  + (tau[35] * z[35][r]);  // bI
    theta_r[3, r][6]  = theta[3][6]  + (tau[36] * z[36][r]);  // rI
    theta_r[3, r][7]  = theta[3][7]  + (tau[37] * z[37][r]);  // bU
    theta_r[3, r][8]  = theta[3][8]  + (tau[38] * z[38][r]);  // rU
    theta_r[3, r][9]  = theta[1][9]  + alpha_gp[r];           // alpha
    theta_r[3, r][10] = theta[1][10] + betaP_gp[r];           // betaP
    theta_r[3, r][11] = theta[1][11] + betaC_gp[r];           // betaC
    theta_r[3, r][12] = theta[1][12] + betaI_gp[r];           // betaI
    theta_r[3, r][13] = theta[1][13] + betaU_gp[r];           // betaU
    theta_r[3, r][14] = theta[1][14];                         // gamma
  }
  
  // solve ode
  array[N_regions, N_dates] vector[9] latent;
  for (r in 1:N_regions) {
    // initial values
    latent[r, 1][1] = init_logit_pop_r[r];
    latent[r, 1][2] = init_pop_size_r[r];
    latent[r, 1][3] = init_logit_crop_r[r];
    latent[r, 1][4] = init_cropland_r[r];
    latent[r, 1][5] = init_logit_irr_r[r];
    latent[r, 1][6] = init_irrigated_r[r];
    latent[r, 1][7] = init_logit_urb_r[r];
    latent[r, 1][8] = init_urban_r[r];
    latent[r, 1][9] = init_gini_r[r];
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
  // priors for initial states
  init_logit_pop ~ normal(-4, 0.5);
  init_pop_size ~ normal(-1, 0.5);
  init_logit_crop ~ normal(-4, 0.5);
  init_cropland ~ normal(-1, 0.5);
  init_logit_irr ~ normal(-4, 0.5);
  init_irrigated ~ normal(-1, 0.5);
  init_logit_urb ~ normal(-4, 0.5);
  init_urban ~ normal(-1, 0.5);
  init_gini ~ normal(-1, 0.5);
  
  // priors for ode parameters
  for (p in 1:3) {
    theta[p][1] ~ normal(-2, 0.5);  // bP
    theta[p][2] ~ normal(-2, 0.5);  // rP
    theta[p][3] ~ normal(-2, 0.5);  // bC
    theta[p][4] ~ normal(-2, 0.5);  // rC
    theta[p][5] ~ normal(-2, 0.5);  // bC
    theta[p][6] ~ normal(-2, 0.5);  // rC
    theta[p][7] ~ normal(-2, 0.5);  // bC
    theta[p][8] ~ normal(-2, 0.5);  // rC
    theta[p][9] ~ normal(0, 0.5);   // alpha
    theta[p][10] ~ normal(0, 0.5);  // betaP
    theta[p][11] ~ normal(0, 0.5);  // betaC
    theta[p][12] ~ normal(0, 0.5);  // betaI
    theta[p][13] ~ normal(0, 0.5);  // betaU
    theta[p][14] ~ normal(-2, 0.5); // gamma
  }
  
  // priors for standardised varying effects, GP, and measurement error
  for (i in 1:38) z[i] ~ normal(0, 1);
  tau ~ exponential(2);
  lscale ~ exponential(2);
  sigma ~ exponential(2);
  phi ~ exponential(2);
  
  // likelihood for population size
  for (i in 1:N_obs_pop) {
    real logit_prob = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][1];
    if (pop_size[i] == 0) {
      target += bernoulli_logit_lpmf(0 | logit_prob);
    } else {
      real mu = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][2];
      target += bernoulli_logit_lpmf(1 | logit_prob);
      target += lognormal_lpdf(pop_size[i] | mu, sigma[1]);
    }
  }
  
  // likelihood for cropland
  for (i in 1:N_obs_crop) {
    real logit_prob = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][3];
    if (cropland[i] == 0) {
      target += bernoulli_logit_lpmf(0 | logit_prob);
    } else {
      real mu = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][4];
      target += bernoulli_logit_lpmf(1 | logit_prob);
      target += lognormal_lpdf(cropland[i] | mu, sigma[2]);
    }
  }
  
  // likelihood for irrigated area
  for (i in 1:N_obs_irr) {
    real logit_prob = latent[region[irr_idx[i]], date_idx[irr_idx[i]]][5];
    if (irrigated[i] == 0) {
      target += bernoulli_logit_lpmf(0 | logit_prob);
    } else {
      real mu = latent[region[irr_idx[i]], date_idx[irr_idx[i]]][6];
      target += bernoulli_logit_lpmf(1 | logit_prob);
      target += lognormal_lpdf(irrigated[i] | mu, sigma[3]);
    }
  }
  
  // likelihood for urban area
  for (i in 1:N_obs_urb) {
    real logit_prob = latent[region[urb_idx[i]], date_idx[urb_idx[i]]][7];
    if (urban[i] == 0) {
      target += bernoulli_logit_lpmf(0 | logit_prob);
    } else {
      real mu = latent[region[urb_idx[i]], date_idx[urb_idx[i]]][8];
      target += bernoulli_logit_lpmf(1 | logit_prob);
      target += lognormal_lpdf(urban[i] | mu, sigma[4]);
    }
  }
  
  // likelihood for gini
  for (i in 1:N_obs_gini) {
    real mu = latent[region[gini_idx[i]], date_idx[gini_idx[i]]][9];
    target += beta_proportion_lpdf(gini[i] | inv_logit(mu), phi);
  }
}
generated quantities {
  array[N_obs_pop] real pop_size_rep;
  array[N_obs_crop] real cropland_rep;
  array[N_obs_irr] real irrigated_rep;
  array[N_obs_urb] real urban_rep;
  array[N_obs_gini] real gini_rep;
  array[121] real date_rep = linspaced_array(121, -100, 20);
  array[N_regions, 121] vector[9] regional_latent_rep;
  
  // posterior predictive check for population size
  for (i in 1:N_obs_pop) {
    real logit_prob = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][1];
    if (bernoulli_logit_rng(logit_prob) == 0) {
      pop_size_rep[i] = 0;
    } else {
      real mu = latent[region[pop_idx[i]], date_idx[pop_idx[i]]][2];
      pop_size_rep[i] = lognormal_rng(mu, sigma[1]);
    }
  }
  
  // posterior predictive check for cropland
  for (i in 1:N_obs_crop) {
    real logit_prob = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][3];
    if (bernoulli_logit_rng(logit_prob) == 0) {
      cropland_rep[i] = 0;
    } else {
      real mu = latent[region[crop_idx[i]], date_idx[crop_idx[i]]][4];
      cropland_rep[i] = lognormal_rng(mu, sigma[2]);
    }
  }
  
  // posterior predictive check for irrigated area
  for (i in 1:N_obs_irr) {
    real logit_prob = latent[region[irr_idx[i]], date_idx[irr_idx[i]]][5];
    if (bernoulli_logit_rng(logit_prob) == 0) {
      irrigated_rep[i] = 0;
    } else {
      real mu = latent[region[irr_idx[i]], date_idx[irr_idx[i]]][6];
      irrigated_rep[i] = lognormal_rng(mu, sigma[3]);
    }
  }
  
  // posterior predictive check for urban area
  for (i in 1:N_obs_urb) {
    real logit_prob = latent[region[urb_idx[i]], date_idx[urb_idx[i]]][7];
    if (bernoulli_logit_rng(logit_prob) == 0) {
      urban_rep[i] = 0;
    } else {
      real mu = latent[region[urb_idx[i]], date_idx[urb_idx[i]]][8];
      urban_rep[i] = lognormal_rng(mu, sigma[4]);
    }
  }
  
  // posterior predictive check for gini
  for (i in 1:N_obs_gini) {
    real mu = latent[region[gini_idx[i]], date_idx[gini_idx[i]]][9];
    gini_rep[i] = beta_proportion_rng(inv_logit(mu), phi);
  }
  
  // regional ode prediction across three periods
  for (r in 1:N_regions) {
    regional_latent_rep[r, 1][1] = init_logit_pop_r[r];
    regional_latent_rep[r, 1][2] = init_pop_size_r[r];
    regional_latent_rep[r, 1][3] = init_logit_crop_r[r];
    regional_latent_rep[r, 1][4] = init_cropland_r[r];
    regional_latent_rep[r, 1][5] = init_logit_irr_r[r];
    regional_latent_rep[r, 1][6] = init_irrigated_r[r];
    regional_latent_rep[r, 1][7] = init_logit_urb_r[r];
    regional_latent_rep[r, 1][8] = init_urban_r[r];
    regional_latent_rep[r, 1][9] = init_gini_r[r];
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
