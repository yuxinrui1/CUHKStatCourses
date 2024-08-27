parameters {
  real x;
}

model {
  x ~ normal(0, 1);
}

generated quantities {
  cov_matrix[2] Sigma[2];
  Sigma[1] = [[2, 0], [0, 1]];
  Sigma[2] = [[3, 0], [0, 4]];
}