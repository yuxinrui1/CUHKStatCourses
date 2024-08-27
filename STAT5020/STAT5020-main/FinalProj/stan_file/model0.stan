data {
  int<lower = 1> N;
  real y1[N];
  real y2[N];
  real y3[N];
  real y4[N];
  real y5[N];
  real y6[N];
  real y7[N];
  real mu0[7];
}

transformed data {
  vector[2] zero2 = rep_vector(0.0, 2);
  cov_matrix[2] mat_phi0 = [[0.04, 0.05], [0.05, 0.7]];
}

parameters {
  real mu[7];
  real lam[4];
  real gam[2];
  vector[2] xi[N];
  real eta[N];
  real<lower = 0> psi_e[7];
  real<lower = 0> psi_d;
  cov_matrix[2] mat_phi;
}

model {
  // prior
  mat_phi ~ inv_wishart(5, mat_phi0);
  
  psi_e[1] ~ inv_gamma(9, 4);
  psi_e[2] ~ inv_gamma(600, 50000);
  psi_e[3] ~ inv_gamma(9, 4);
  psi_e[4] ~ inv_gamma(9, 4);
  psi_e[5] ~ inv_gamma(9, 4);
  psi_e[6] ~ inv_gamma(9, 4);
  psi_e[7] ~ inv_gamma(9, 4);
  psi_d ~ inv_gamma(300, 12000);

  for (i in 1:7) {
    mu[i] ~ normal(mu0[i], .5);
  }
  lam[1] ~ normal(3, .5 * sqrt(psi_e[2]));
  lam[2] ~ normal(1, .5 * sqrt(psi_e[4]));
  lam[3] ~ normal(.9, .5 * sqrt(psi_e[5]));
  lam[4] ~ normal(.05, .5 * sqrt(psi_e[7]));
  gam[1] ~ normal(24, .5 * sqrt(psi_d));
  gam[2] ~ normal(-4, .2 * sqrt(psi_d));
  
  // explanatory latent variables
  for (i in 1:N) {
    xi[i] ~ multi_normal(zero2, mat_phi);
  }
  
  // measurement equation
  for (i in 1:N) {
    y1[i] ~ normal(mu[1] + eta[i],            sqrt(psi_e[1]));
    y2[i] ~ normal(mu[2] + eta[i]   * lam[1], sqrt(psi_e[2]));
    y3[i] ~ normal(mu[3] + xi[i, 1],          sqrt(psi_e[3]));
    y4[i] ~ normal(mu[4] + xi[i, 1] * lam[2], sqrt(psi_e[4]));
    y5[i] ~ normal(mu[5] + xi[i, 1] * lam[3], sqrt(psi_e[5]));
    y6[i] ~ normal(mu[6] + xi[i, 2],          sqrt(psi_e[6]));
    y7[i] ~ normal(mu[7] + xi[i, 2] * lam[4], sqrt(psi_e[7]));
  }
  
  // structural equation
  for (i in 1:N) {
    eta[i] ~ normal(gam[1] * xi[i, 1] + gam[2] * xi[i, 2], sqrt(psi_d));
  }
}
