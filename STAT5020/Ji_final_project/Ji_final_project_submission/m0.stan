data {
  int<lower = 1> N;
  real y1[N];
  real y2[N];
  real y3[N];
  real y4[N];
  real y5[N];
  real y6[N];
  real y7[N];
  real y8[N];
  real y9[N];
  real y10[N];
  real y11[N];
  real y12[N];
  real g[N];
  real mu0[12];
}

transformed data {
  vector[2] zero2 = rep_vector(0.0, 2);
  cov_matrix[2] mat_phi0 = [[0.04, 0.05], [0.05, 0.7]];
}

parameters {
  real mu[12];
  real lam[9];
  real gam[2];
  vector[2] xi[N];
  real eta[N];
  real<lower = 0> psi_e[12];
  real<lower = 0> psi_d;
  cov_matrix[2] mat_phi;
}

model {
  // prior
  mat_phi ~ inv_wishart(5, mat_phi0);
  
  psi_e[1] ~ inv_gamma(9, 4);
  psi_e[2] ~ inv_gamma(9, 4);
  psi_e[3] ~ inv_gamma(9, 4);
  psi_e[4] ~ inv_gamma(9, 4);
  psi_e[5] ~ inv_gamma(9, 4);
  psi_e[6] ~ inv_gamma(9, 4);
  psi_e[7] ~ inv_gamma(9, 4);
  psi_e[8] ~ inv_gamma(9, 4);
  psi_e[9] ~ inv_gamma(9, 4);
  psi_e[10] ~ inv_gamma(9, 4);
  psi_e[11] ~ inv_gamma(9, 4);
  psi_e[12] ~ inv_gamma(9, 4);
  psi_d ~ inv_gamma(9, 4);

  for (i in 1:12) {
    mu[i] ~ normal(mu0[i], .5);
  }
  lam[1] ~ normal(0, 5);
  lam[2] ~ normal(0, 5);
  lam[3] ~ normal(0, 5);
  lam[4] ~ normal(0, 5);
  lam[5] ~ normal(0, 5);
  lam[6] ~ normal(0, 5);
  lam[7] ~ normal(0, 5);
  lam[8] ~ normal(0, 5);
  lam[9] ~ normal(0, 5);
  gam[1] ~ normal(0, 5);
  gam[2] ~ normal(0, 5);
  
  // explanatory latent variables
  for (i in 1:N) {
    xi[i] ~ multi_normal(zero2, mat_phi);
  }
  
  // measurement equation
  for (i in 1:N) {
    y1[i] ~ normal(mu[1] + eta[i],            sqrt(psi_e[1]));
    y2[i] ~ normal(mu[2] + eta[i]   * lam[1], sqrt(psi_e[2]));
    y3[i] ~ normal(mu[3] + eta[i]   * lam[2], sqrt(psi_e[3]));
    y4[i] ~ normal(mu[4] + xi[i, 1]         , sqrt(psi_e[4]));
    y5[i] ~ normal(mu[5] + xi[i, 1] * lam[3], sqrt(psi_e[5]));
    y6[i] ~ normal(mu[6] + xi[i, 1] * lam[4], sqrt(psi_e[6]));
    y7[i] ~ normal(mu[7] + xi[i, 1] * lam[5], sqrt(psi_e[7]));
    y8[i] ~ normal(mu[8] + xi[i, 1] * lam[6], sqrt(psi_e[8]));
    y9[i] ~ normal(mu[9] + xi[i, 2],          sqrt(psi_e[9]));
    y10[i] ~ normal(mu[10] + xi[i, 2] * lam[7], sqrt(psi_e[10]));
    y11[i] ~ normal(mu[11] + xi[i, 2] * lam[8], sqrt(psi_e[11]));
    y12[i] ~ normal(mu[12] + xi[i, 2] * lam[9], sqrt(psi_e[12]));
  }
  
  // structural equation
  for (i in 1:N) {
    eta[i] ~ normal(gam[1] * xi[i, 1] + gam[2] * xi[i, 2], sqrt(psi_d));
  }
}
