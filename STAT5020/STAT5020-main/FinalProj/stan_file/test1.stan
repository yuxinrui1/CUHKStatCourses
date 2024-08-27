functions {
  int count1(data vector x) {
    int n = 0;
    for (xi in x) {
      if (xi == 1) {
        n += 1;
      }
    }
    return n;
  }
}

data {
  int<lower = 0> N;
  vector[N] y;
}
transformed data {
  int n = count1(y);
  int<lower = 1, upper = N> sub[n];
  int<lower = 1, upper = N> subc[N - n];
  int j = 1;
  int jc = 1;
  for (i in 1:N) {
    if (y[i] == 1) {
      sub[j] = i;
      j += 1;
    }
    else {
      subc[jc] = i;
      jc += 1;
    }
  }
}
parameters {
  real z[N];
  real<lower = 0> psi0;
}
// transformed parameters {
//   real<lower = 0> psi[2];
//   psi[1] = psi0;
//   psi[2] = 1;
// }
model {
  psi0 ~ normal(2, .1) T[0, ];
  for (i in sub) 
    z[i] ~ normal(0, psi0) T[0, ];
  for (i in subc) 
    z[i] ~ normal(0, 1) T[, 0];
}