data {
  int<lower = 1> N;
  int<lower = 0, upper = N> n1;
  int<lower = 0, upper = N> n2;
  int<lower = 0, upper = N> n3;
  int<lower = 1, upper = N> sub1[n1];       // loc which elem equal to 1
  int<lower = 1, upper = N> sub1c[N - n1];  // loc which elem equal to 0
  int<lower = 1, upper = N> sub2[n2];
  int<lower = 1, upper = N> sub2c[N - n2];
  int<lower = 1, upper = N> sub3[n3];
  int<lower = 1, upper = N> sub3c[N - n3];
  vector[N] y4;
  vector[N] y5;
  vector[N] y6;
  vector[N] y7;
  int<lower = 0, upper = 1> y8[N];
  int<lower = 0, upper = 1> y9[N];
  int<lower = 0, upper = 1> d[N];
}
transformed data {
  vector[2] zero2;
  cov_matrix[2] mat_phi0;
  zero2 = rep_vector(0, 2);
  mat_phi0 = [[7, 0], [0, 7]];
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
  real<lower = 0> y1p[n1];
  real<upper = 0> y1n[(N - n1)];
  real<lower = 0> y2p[n2];
  real<upper = 0> y2n[(N - n2)];
  real<lower = 0> y3p[n3];
  real<upper = 0> y3n[(N - n3)];
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
  // prior
  for (j in 1:4)
    psi_e[j] ~ inv_gamma(9, 3);
  psi_d ~ inv_gamma(9, 3);
  
  mat_phi ~ inv_wishart(10, mat_phi0);
  
  for (j in 1:3)
    mu[j] ~ normal(0, 0.4);
  for (j in 4:7)
    mu[j] ~ normal(0, 0.4 * sqrt(psi_e[j-3]));
  for (j in 8:9) 
    mu[j] ~ normal(0, 0.4);
  
  lam2 ~ normal(0.8, 0.4);
  lam3 ~ normal(0.8, 0.4);
  lam5 ~ normal(0.7, 0.4 * sqrt(psi_e[2]));
  lam6 ~ normal(0.9, 0.4 * sqrt(psi_e[3]));
  lam7 ~ normal(0.7, 0.4 * sqrt(psi_e[4]));
  lam9 ~ normal(0.8, 0.4);
  
  b ~ normal(0.3, 0.4 * sqrt(psi_d));
  gam[1] ~ normal(0.2, 0.4 * sqrt(psi_d));
  gam[2] ~ normal(0.5, 0.4 * sqrt(psi_d));
  
  for (i in 1:N) {
    xi[i] ~ multi_normal(zero2, mat_phi);
  }
  
  // measurement equation model
  for (i in 1:n1)
    y1p[i] ~ normal(mu[1] + eta[sub1[i]], 1) T[0, ];
  for (i in 1:(N - n1))
    y1n[i] ~ normal(mu[1] + eta[sub1c[i]], 1) T[, 0];
  for (i in 1:n2)
    y2p[i] ~ normal(mu[2] + lam2 * eta[sub2[i]], 1) T[0, ];
  for (i in 1:(N - n2))
    y2n[i] ~ normal(mu[2] + lam2 * eta[sub2c[i]], 1) T[, 0];
  for (i in 1:n3)
    y3p[i] ~ normal(mu[3] + lam3 * eta[sub3[i]], 1) T[0, ];
  for (i in 1:(N - n3))
    y3n[i] ~ normal(mu[3] + lam3 * eta[sub3c[i]], 1) T[, 0];
  for (i in 1:N) {
    y4[i] ~ normal(mu[4] + xi[i, 1], sqrt(psi_e[1]));
    y5[i] ~ normal(mu[5] + lam5 * xi[i, 1], sqrt(psi_e[2]));
    y6[i] ~ normal(mu[6] + lam6 * xi[i, 1], sqrt(psi_e[3]));
    y7[i] ~ normal(mu[7] + lam7 * xi[i, 1], sqrt(psi_e[4]));
    y8[i] ~ bernoulli_logit(vart1[i]);
    y9[i] ~ bernoulli_logit(vart2[i]);
  }
  // structural equation model
  for (i in 1:N) {
    eta[i] ~ normal(b * d[i] + gam[1] * xi[i, 1] + gam[2] * xi[i, 2], sqrt(psi_d));
  }
}





