data {
  int<lower = 1> N; // sample size
  int<lower = 1> p; // dim of `y`
  vector[N] Y[p];  // observed data
  real c[N];       // fixed effect in m.e.
  real d[N];       // fixed effect in s.e.
  real t;          
}
parameters {
  real xi1[N];  // 1st explanatory latent variable
  real xi2[N];  // 2nd explanatory latent variable
  real eta[N];  // outcome latent variable
  real mu[p];   // mean in m.e.
  real a[p];    // fixed effect in m.e.
  real lam_21;      // effect of latent variable in m.e.
  real lam_31;
  real lam_52;
  real lam_62;
  real lam_83;
  real lam_93;
  real<lower = 0> psi_eps[p];  // variance in m.e.
  real<lower = 0> phi[2]; // variance of e.l.v.
  real b;                 // fixed effect in s.e.
  real gam[5];      // effect in s.e.
  real<lower = 0> psi_delta;    // variance in s.e.
}
model {
  psi_delta ~ inv_gamma(9, 4);  // priors
  for (k in 1:p) {
    psi_eps[k] ~ inv_gamma(9, 4);
  }
  for (k in 1:p) {
    mu[k] ~ normal(0, sqrt(psi_eps[k]));
  }
  a[1] ~ normal(0.5, sqrt(psi_eps[2]));
  a[2] ~ normal(0.5, sqrt(psi_eps[2]));
  a[3] ~ normal(0.5, sqrt(psi_eps[2]));
  a[4] ~ normal(-0.5, sqrt(psi_eps[2]));
  a[5] ~ normal(-0.5, sqrt(psi_eps[2]));
  a[6] ~ normal(-0.5, sqrt(psi_eps[2]));
  a[7] ~ normal(0.8, sqrt(psi_eps[2]));
  a[8] ~ normal(0.8, sqrt(psi_eps[2]));
  a[9] ~ normal(0.8, sqrt(psi_eps[2]));
  lam_21 ~ normal(0.8, sqrt(psi_eps[2]));
  lam_31 ~ normal(0.5, sqrt(psi_eps[3]));
  lam_52 ~ normal(0.8, sqrt(psi_eps[5]));
  lam_62 ~ normal(0.5, sqrt(psi_eps[6]));
  lam_83 ~ normal(0.8, sqrt(psi_eps[8]));
  lam_93 ~ normal(0.5, sqrt(psi_eps[9]));
  b ~ normal(-1, psi_delta);
  gam[1] ~ normal(0.5, psi_delta);
  gam[2] ~ normal(0.5, psi_delta);
  gam[3] ~ normal(0.3, psi_delta);
  gam[4] ~ normal(0.3, psi_delta);
  gam[5] ~ normal(0, psi_delta);
  for (i in 1:N) {
    Y[1, i] ~ normal(mu[1] + a[1] * c[i] + eta[i], sqrt(psi_eps[1]));
    Y[2, i] ~ normal(mu[2] + a[2] * c[i] + lam_21 * eta[i], sqrt(psi_eps[2]));
    Y[3, i] ~ normal(mu[3] + a[3] * c[i] + lam_31 * eta[i], sqrt(psi_eps[3]));
    Y[4, i] ~ normal(mu[4] + a[4] * c[i] + xi1[i], sqrt(psi_eps[4]));
    Y[5, i] ~ normal(mu[5] + a[5] * c[i] + lam_52 * xi1[i], sqrt(psi_eps[5]));
    Y[6, i] ~ normal(mu[6] + a[6] * c[i] + lam_62 * xi1[i], sqrt(psi_eps[6]));
    Y[7, i] ~ normal(mu[7] + a[7] * c[i] + xi2[i], sqrt(psi_eps[7]));
    Y[8, i] ~ normal(mu[8] + a[8] * c[i] + lam_83 * xi2[i], sqrt(psi_eps[8]));
    Y[9, i] ~ normal(mu[9] + a[9] * c[i] + lam_93 * xi2[i], sqrt(psi_eps[9]));
  }
  for (i in 1:N) {
    xi1[i] ~ normal(0, sqrt(phi[1]));
    xi2[i] ~ normal(0, sqrt(phi[2]));
  }
  for (i in 1:N) {
    eta[i] ~ normal(b * d[i] + gam[1] * xi1[i] + gam[2] * xi2[i] + 
                    gam[3] * xi1[i] * xi2[i] + gam[4] * xi1[i]^2 + 
                    t * gam[5] * xi2[i]^2, 
                    sqrt(psi_delta));
  }
}
generated quantities {
  real U = 0;
  for (i in 1:N) {
    U = U - (eta[i] - (b * d[i] + gam[1] * xi1[i] + gam[2] * xi2[i] + 
                    gam[3] * xi1[i] * xi2[i] + gam[4] * xi1[i]^2 + 
                    t * gam[5] * xi2[i]^2) ) * (- gam[5] * xi2[i]^2) / psi_delta;
  }
}



