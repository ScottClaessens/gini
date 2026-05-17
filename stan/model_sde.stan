functions {
  
  // ──────────────────────────────────────────
  // Solver for the asymptotic Q matrix
  // ──────────────────────────────────────────
  
  matrix ksolve (matrix A, matrix Q) {
    int d = rows(A);
    int d2 = (d * d - d) %/% 2;
    matrix [d + d2, d + d2] O;
    vector [d + d2] triQ;
    matrix[d,d] AQ;
    int z = 0;         // z is row of output
    for (j in 1:d) {   // for column reference of solution vector
      for (i in 1:j) { // and row reference...
        if (j >= i) {  // if i and j denote a covariance parameter
          int y = 0;   // start new output row
          z += 1;      // shift current output row down
          for (ci in 1:d) {   // for columns and
            for (ri in 1:d) { // rows of solution
              if (ci >= ri) { // when in upper tri (inc diag)
                y += 1;       // move to next column of output
                if (i == j) { // if output row is a diag element
                  if (ri == i) O[z, y] = 2 * A[ri, ci];
                  if (ci == i) O[z, y] = 2 * A[ci, ri];
                }
                if (i != j) { // if output row is not a diag element
                  //if column matches row, sum both A diags
                  if (y == z) O[z, y] = A[ri, ri] + A[ci, ci];
                  if (y != z) { // otherwise...
                    // if solution element is related to output row...
                    if (ci == ri) { // if solution element is variance
                      // if variance of solution corresponds to row
                      if (ci == i) O[z, y] = A[j, ci];
                      // if variance of solution corresponds to col
                      if (ci == j) O[z, y] = A[i, ci];
                    }
                    //if solution element is a related covariance
                    if (ci != ri && (ri == i || ri == j || 
                                     ci == i || ci == j)) {
                      // for row 1,2 / 2,1 of output,
                      // if solution row ri 1 (match)
                      // and column ci 3, we need A[2,3]
                      if (ri == i) O[z, y] = A[j, ci];
                      if (ri == j) O[z, y] = A[i, ci];
                      if (ci == i) O[z, y] = A[j, ri];
                      if (ci == j) O[z, y] = A[i, ri];
                    }
                  }
                }
                if (is_nan(O[z, y])) O[z, y] = 0;
              }
            }
          }
        }
      }
    }
    z = 0; // get upper tri of Q
    for (j in 1:d) {
      for (i in 1:j) {
        z += 1;
        triQ[z] = Q[i, j];
      }
    }
    triQ = -O \ triQ; // get upper tri of asymQ
    z = 0; // put upper tri of asymQ into matrix
    for (j in 1:d) {
      for (i in 1:j) {
        z += 1;
        AQ[i, j] = triQ[z];
        if (i != j) AQ[j, i] = triQ[z];
      }
    }
    return AQ;
  }
  
}

data {
  
  int<lower=1> N_sites;                    // number of sites
  int<lower=1> N_times;                    // number of unique times
  vector<lower=0, upper=1>[N_sites] gini;  // gini estimates
  vector[N_sites] pop_scaled;              // population size (scaled)
  vector<lower=0>[N_times] ts;             // time differences
  array[N_sites] int ts_idx;               // id mapping sites to times
  int<lower=0, upper=1> prior_only;        // sample the prior only?
  
}

parameters {
  
  vector<upper=0>[2] A_diag;             // autoregressive terms of A matrix
  vector[2] A_offdiag;                   // cross-lagged terms of A matrix
  vector<lower=0>[2] Q_sigma;            // standard deviations of Q matrix
  vector[2] b;                           // SDE intercepts
  vector[2] eta_initial;                 // initial states
  real<lower=0> phi;                     // gini precision parameter
  real<lower=0> shape;                   // population size shape parameter
  array[N_times - 1] vector[2] z_drift;  // stochastic drift

}

transformed parameters {

  // ──────────────────────────────────────────
  // SDE process
  // ──────────────────────────────────────────

  array[N_times] vector[2] eta;                  // latent states
  matrix[2, 2] A = diag_matrix(A_diag);          // selection matrix
  matrix[2, 2] Q = diag_matrix(square(Q_sigma)); // drift matrix
  matrix[2, 2] Q_inf;                            // asymptotic covariance matrix
  
  // fill off diagonal of A matrix
  A[1, 2] = A_offdiag[1];
  A[2, 1] = A_offdiag[2];
  
  // calculate asymptotic covariance
  Q_inf = ksolve(A, Q);
  
  // set initial states
  eta[1] = eta_initial;
  
  // loop over unique time points
  for (i in 2:N_times) {
    
    matrix[2, 2] A_delta;
    matrix[2, 2] VCV;
    matrix[2, 2] L_VCV;
    matrix[2, 2] A_solve;
    
    A_delta = matrix_exp(A * ts[i]);
    VCV = Q_inf - quad_form_sym(Q_inf, A_delta');
    L_VCV = cholesky_decompose(VCV);
    A_solve = A \ add_diag(A_delta, -1);
    
    eta[i] = to_vector(
      A_delta * eta[i - 1] + (A_solve * b) + (L_VCV * z_drift[i - 1])
    );
    
  }

}

model {
  
  // ─────────────────────────────────────────
  // Priors
  // ─────────────────────────────────────────
  
  A_diag ~ normal(0, 1);
  A_offdiag ~ normal(0, 1);
  Q_sigma ~ normal(0, 1);
  b ~ normal(0, 1);
  eta_initial ~ normal(0, 1);
  phi ~ lognormal(3, 1);
  shape ~ exponential(1);
  for (i in 1:(N_times - 1)) {
    z_drift[i] ~ normal(0, 1);
  }
  
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
    
      mu = inv_logit(eta[ts_idx[n], 1]);
      shape1 = mu * phi;
      shape2 = (1.0 - mu) * phi;
      
      gini[n] ~ beta(shape1, shape2);
      
      // ─────────────────────────────────────────
      // Log population size
      // ─────────────────────────────────────────
      
      if (pop_scaled[n] != -9999) {
        pop_scaled[n] ~ gamma(shape, shape / exp(eta[ts_idx[n], 2]));
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

    mu = inv_logit(eta[ts_idx[n], 1]);
    shape1 = mu * phi;
    shape2 = (1.0 - mu) * phi;
    
    gini_rep[n] = beta_rng(shape1, shape2);
    
    pop_scaled_rep[n] = gamma_rng(shape, shape / exp(eta[ts_idx[n], 2]));

    // ─────────────────────────────────────────
    // Pointwise log-likelihood for Gini
    // ─────────────────────────────────────────
    
    log_lik[n] = beta_lpdf(gini[n] | shape1, shape2);
    
  }
}
