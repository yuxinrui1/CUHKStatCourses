data {
int<lower = 0> N; // number of observations 
int<lower = 1> L; // number of groups 
int<lower = 1, upper = L> ll[N]; // group info
vector[N] score;
vector<lower = 0>[N] iq;
vector<lower = 0>[N] ses;
vector<lower = 0>[L] iqcl;
int<lower = 0, upper = 1> gender[N]; // gender 
int<lower = 1> D1; // number of random effects 
int<lower = 1> D2; // number of non-random effects 
int<lower = 1> D3; // number of covariates for variance
}
transformed data {
  vector[N] iq_cen;
  vector[N] ses_cen;
  iq_cen = iq - mean(iq); // centerized iq 
  ses_cen = ses - mean(ses); // centerized ses
}
parameters {
  vector[D1] m; vector[D2] beta; vector[D3] theta;
}
transformed parameters {
// coeffecient of non-random effects
// variance structure parameters
// random effects for each group
  vector[N] V;
  vector[N] mu;
  V = theta[1] + theta[2] * iq; // variance structure 
  for (i in 1:N)
    mu[i] = m[1] + m[2] * iq_cen[i] +
    beta[1] * ses_cen[i] + beta[2] * gender[i] + beta[3] * iqcl[ll[i]];
}
model {
  for (i in 1:N)
    score[i] ~ normal(mu[i], sqrt(V[i])); // main regression model
}