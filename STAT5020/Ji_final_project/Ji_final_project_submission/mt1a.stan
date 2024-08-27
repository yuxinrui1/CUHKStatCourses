data {
  int<lower = 1> N;
  int g[N];
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
  real mu0[12];
  real<lower = 0, upper = 1> t;
}

transformed data {
  vector[2] zero2 = rep_vector(0.0, 2);
  cov_matrix[2] mat_phi0 = [[0.04, 0.05], [0.05, 0.7]];
}

parameters {
  real mu[12, 2];
  real lam[9, 2];
  real gam[2, 2];
  vector[2] xi[N];
  real eta[N];
  real<lower = 0> psi_e[12, 2];
  real<lower = 0> psi_d[2];
  cov_matrix[2] mat_phi[2];
}

model {
  // prior
  for (k in 1:2) {
    mat_phi[k] ~ inv_wishart(5, mat_phi0);
    psi_e[1, k] ~ inv_gamma(9, 4);
    psi_e[2, k] ~ inv_gamma(9, 4);
    psi_e[3, k] ~ inv_gamma(9, 4);
    psi_e[4, k] ~ inv_gamma(9, 4);
    psi_e[5, k] ~ inv_gamma(9, 4);
    psi_e[6, k] ~ inv_gamma(9, 4);
    psi_e[7, k] ~ inv_gamma(9, 4);
    psi_e[8, k] ~ inv_gamma(9, 4);
    psi_e[9, k] ~ inv_gamma(9, 4);
    psi_e[10, k] ~ inv_gamma(9, 4);
    psi_e[11, k] ~ inv_gamma(9, 4);
    psi_e[12, k] ~ inv_gamma(9, 4);
    psi_d[k] ~ inv_gamma(9, 4);
  }
  
  for (k in 1:2) {
    for (i in 1:12) {
      mu[i, k] ~ normal(mu0[i], .5);
    }
    lam[1, k] ~ normal(0, 5);
    lam[2, k] ~ normal(0, 5);
    lam[3, k] ~ normal(0, 5);
    lam[4, k] ~ normal(0, 5);
    lam[5, k] ~ normal(0, 5);
    lam[6, k] ~ normal(0, 5);
    lam[7, k] ~ normal(0, 5);
    lam[8, k] ~ normal(0, 5);
    lam[9, k] ~ normal(0, 5);
    gam[1, k] ~ normal(0, 5);
    gam[2, k] ~ normal(0, 5);
  }
  
  
  // explanatory latent variables
  for (i in 1:N) {
    xi[i] ~ multi_normal(zero2, mat_phi[g[i]]);
  }
  
  // measurement equation
  for (i in 1:N) {
    y1[i] ~ normal(mu[1, g[i]] + eta[i],                  sqrt(psi_e[1, g[i]]));
    y2[i] ~ normal(mu[2, g[i]] + eta[i]   * lam[1, g[i]], sqrt(psi_e[2, g[i]]));
    y3[i] ~ normal(mu[3, g[i]] + eta[i]   * lam[2, g[i]], sqrt(psi_e[3, g[i]]));
    y4[i] ~ normal(mu[4, g[i]] + xi[i, 1]               , sqrt(psi_e[4, g[i]]));
    y5[i] ~ normal(mu[5, g[i]] + xi[i, 1] * lam[3, g[i]], sqrt(psi_e[5, g[i]]));
    y6[i] ~ normal(mu[6, g[i]] + xi[i, 1] * lam[4, g[i]], sqrt(psi_e[6, g[i]]));
    y7[i] ~ normal(mu[7, g[i]] + xi[i, 1] * lam[5, g[i]], sqrt(psi_e[7, g[i]]));
    y8[i] ~ normal(mu[8, g[i]] + xi[i, 1] * lam[6, g[i]], sqrt(psi_e[8, g[i]]));
    y9[i] ~ normal(mu[9, g[i]] + xi[i, 2]               , sqrt(psi_e[9, g[i]]));
    y10[i] ~ normal(mu[10, g[i]] + xi[i, 2] * lam[7, g[i]], sqrt(psi_e[10, g[i]]));
    y11[i] ~ normal(mu[11, g[i]] + xi[i, 2] * lam[8, g[i]], sqrt(psi_e[11, g[i]]));
    y12[i] ~ normal(mu[12, g[i]] + xi[i, 2] * lam[9, g[i]], sqrt(psi_e[12, g[i]]));
  }
  
  // structural equation
  for (i in 1:N) {
    eta[i] ~ normal(t * gam[1, g[i]] * xi[i, 1] + t * gam[2, g[i]] * xi[i, 2], sqrt(psi_d[g[i]]));
  }
}

generated quantities {
  real U = 0;
  for (i in 1:N) {
    U += (eta[i] - t * gam[1, g[i]] * xi[i, 1] - t * gam[2, g[i]] * xi[i, 2]) * 
         (gam[1, g[i]] * xi[i, 1] + gam[2, g[i]] * xi[i, 2]) / psi_d[g[i]];
  }
}
