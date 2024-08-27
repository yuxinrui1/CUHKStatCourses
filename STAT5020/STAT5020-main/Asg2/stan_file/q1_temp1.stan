data {
  int<lower = 1> N;
  int<lower = 0, upper = 1> y1[N];
  int<lower = 0, upper = 1> y2[N];
  int<lower = 0, upper = 1> y3[N];
  vector[N] y4;
  vector[N] y5;
  vector[N] y6;
  vector[N] y7;
  int<lower = 0, upper = 1> y8[N];
  int<lower = 0, upper = 1> y9[N];
  int<lower = 0, upper = 1> d[N];
}
transformed data {
  vector[2] lwbd;
  vector[2] upbd;
  vector[2] zero2;
  cov_matrix[2] mat_phi0;
  lwbd[1] = negative_infinity();
  lwbd[2] = 0;
  upbd[1] = 0;
  upbd[2] = positive_infinity();
  zero2 = rep_vector(0, 2);
  mat_phi0 = diag_matrix(rep_vector(1.0, 2));
}
parameters {
  real lam2;
  real lam3;
  real lam5;
  real lam6;
  real lam7;
  real lam9;
  real b;
  real gam[2];
  vector[N] y1s;
  vector[N] y2s;
  vector[N] y3s;
  vector[N] eta;
  vector[2] xi[N];
  vector[9] mu;
  cov_matrix[2] mat_phi;
  real<lower = 0> psi_d;
  real<lower = 0> psi_e[4];
}
transformed parameters {
  vector[N] vart1;
  vector[N] vart2;
  for (i in 1:N) {
    vart1[i] = mu[8] + xi[i, 2];
    vart2[i] = mu[9] + lam9 * xi[i, 2];
  }
}
model {
  for (j in 1:4)
    psi_e[j] ~ inv_gamma(9, 3);
  psi_d ~ inv_gamma(9, 3);
  
  for (j in 1:3)
    mu[j] ~ normal(0, 0.25);
  for (j in 4:7)
    mu[j] ~ normal(0, 0.25 * sqrt(psi_e[j-3]));
  
  lam2 ~ normal(0.8, 0.25);
  lam3 ~ normal(0.8, 0.25);
  lam5 ~ normal(0.7, 0.25 * sqrt(psi_e[2]));
  lam6 ~ normal(0.9, 0.25 * sqrt(psi_e[3]));
  lam7 ~ normal(0.7, 0.25 * sqrt(psi_e[4]));
  lam9 ~ normal(0.8, 0.25);
  
  gam[1] ~ normal(1, sqrt(psi_d));
  gam[2] ~ normal(1, sqrt(psi_d));
  
  mat_phi ~ inv_wishart(10, mat_phi0);
  
  // measurement equation model
  for (i in 1:N) {
    y1s[i] ~ normal(mu[1] + eta[i], 1) T[lwbd[(y1[i] + 1)], upbd[(y1[i] + 1)]];
    y2s[i] ~ normal(mu[2] + lam2 * eta[i], 1) T[lwbd[(y2[i] + 1)], upbd[(y2[i] + 1)]];
    y3s[i] ~ normal(mu[3] + lam3 * eta[i], 1) T[lwbd[(y3[i] + 1)], upbd[(y3[i] + 1)]];
    y4[i] ~ normal(mu[4] + xi[i, 1], sqrt(psi_e[1]));
    y5[i] ~ normal(mu[5] + lam5 * xi[i, 1], sqrt(psi_e[2]));
    y6[i] ~ normal(mu[6] + lam6 * xi[i, 1], sqrt(psi_e[3]));
    y7[i] ~ normal(mu[7] + lam7 * xi[i, 1], sqrt(psi_e[4]));
    y8[i] ~ bernoulli_logit(vart1[i]);
    y9[i] ~ bernoulli_logit(vart2[i]);
  }
  for (i in 1:N) {
    xi[i] ~ multi_normal(zero2, mat_phi);
  }
  // structural equation model
  for (i in 1:N) {
    eta[i] ~ normal(b * d[i] + gam[1] * xi[i, 1] + gam[2] * xi[i, 2], sqrt(psi_d));
  }
}









