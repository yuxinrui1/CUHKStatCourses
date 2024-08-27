data {
  int<lower = 1> p;
  real a[p];
}
transformed data {
  real lbd = 0;
  real ubd;
  ubd = positive_infinity();
}
parameters {
  real alpha;
}
model {
  alpha ~ normal(0, 1) T[lbd, ubd];
}
generated quantities {
  matrix[2, 2] x = [[1, 2], [3, 4]];
  vector[2] y;
  y = x[:, 1];
}