data {
  int<lower = 1> N;
  int<lower = 1, upper = 2> g[N];
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
  real mu[7, 2];
  real lam[4, 2];
  real gam[2, 2];
  vector[2] xi[N];
  real eta[N];
  real<lower = 0> psi_e[7, 2];
  real<lower = 0> psi_d[2];
  cov_matrix[2] mat_phi[2];
}

model {
  // prior
  for (k in 1:2) {
    mat_phi[k] ~ inv_wishart(5, mat_phi0);
    psi_e[1, k] ~ inv_gamma(9, 4);
    psi_e[2, k] ~ inv_gamma(600, 50000);
    psi_e[3, k] ~ inv_gamma(9, 4);
    psi_e[4, k] ~ inv_gamma(9, 4);
    psi_e[5, k] ~ inv_gamma(9, 4);
    psi_e[6, k] ~ inv_gamma(9, 4);
    psi_e[7, k] ~ inv_gamma(9, 4);
    psi_d[k] ~ inv_gamma(300, 12000);
  }
  
  for (k in 1:2) {
    for (i in 1:7) {
      mu[i, k] ~ normal(mu0[i], .5);
    }
    lam[1, k] ~ normal(3,   .5);
    lam[2, k] ~ normal(1,   .5);
    lam[3, k] ~ normal(.9,  .5);
    lam[4, k] ~ normal(.05, .5);
    gam[1, k] ~ normal(24,   1);
    gam[2, k] ~ normal(-4,  .5);
  }
  
  
  // explanatory latent variables
  for (i in 1:N) {
    xi[i] ~ multi_normal(zero2, mat_phi[g[i]]);
  }
  
  // measurement equation
  for (i in 1:N) {
    y1[i] ~ normal(mu[1, g[i]] + eta[i],                  sqrt(psi_e[1, g[i]]));
    y2[i] ~ normal(mu[2, g[i]] + eta[i]   * lam[1, g[i]], sqrt(psi_e[2, g[i]]));
    y3[i] ~ normal(mu[3, g[i]] + xi[i, 1],                sqrt(psi_e[3, g[i]]));
    y4[i] ~ normal(mu[4, g[i]] + xi[i, 1] * lam[2, g[i]], sqrt(psi_e[4, g[i]]));
    y5[i] ~ normal(mu[5, g[i]] + xi[i, 1] * lam[3, g[i]], sqrt(psi_e[5, g[i]]));
    y6[i] ~ normal(mu[6, g[i]] + xi[i, 2],                sqrt(psi_e[6, g[i]]));
    y7[i] ~ normal(mu[7, g[i]] + xi[i, 2] * lam[4, g[i]], sqrt(psi_e[7, g[i]]));
  }
  
  // structural equation
  for (i in 1:N) {
    eta[i] ~ normal(gam[1, g[i]] * xi[i, 1] + gam[2, g[i]] * xi[i, 2], sqrt(psi_d[g[i]]));
  }
}
