data {
  int<lower=0> N;
  matrix[N, 9] Y;
  vector[N] C;
  vector[N] D;
  vector[N] xi1;
  vector[N] xi2;
  vector[N] eta;
}

parameters {
  vector[9] mu;
  real b;
  vector<lower=0>[9] a;
  vector[4] gamma;
  vector[6] lambda;
  real<lower=0> sigma_eta;
  vector<lower=0>[9] sigma_eps;
}

model {
  mu ~ normal(0, 5);
  a ~ normal(0, 5);
  lambda ~ normal(0, 5);
  b ~ normal(0, 5);
  gamma ~ normal(0, 5);
  sigma_eta ~ cauchy(0, 5);
  sigma_eps ~ cauchy(0, 5);
  eta ~ normal(b * D + gamma[1] * xi1 + gamma[2] * xi2  , sigma_eta);
  
  for (n in 1:N) {
    Y[n, 1] ~ normal(mu[1] + a[1] * C[n] + eta[n], sigma_eps[1]);
    Y[n, 4] ~ normal(mu[4] + a[4] * C[n] + xi1[n], sigma_eps[4]);
    Y[n, 7] ~ normal(mu[7] + a[7] * C[n] + xi2[n], sigma_eps[7]);
    
    for (i in 2:3) {
      Y[n, i] ~ normal(mu[i] + a[i] * C[n] + lambda[i - 1] * eta[n], sigma_eps[i]);
    }
    
    for (i in 5:6) {
      Y[n, i] ~ normal(mu[i] + a[i] * C[n] + lambda[i - 2] * xi1[n], sigma_eps[i]);
    }
    
    for (i in 8:9) {
      Y[n, i] ~ normal(mu[i] + a[i] * C[n] + lambda[i - 3] * xi2[n], sigma_eps[i]);
    }
  }
}