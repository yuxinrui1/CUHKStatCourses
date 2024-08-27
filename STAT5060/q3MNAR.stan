data {
int<lower = 0> N; // number of observations 
int<lower = 0> N_obs; // number of obsevered data points 
int<lower = 0> N_mis; // number of missing data points
int<lower = 0> L;
int<lower = 0> D1;
vector[N_obs] y_obs;
vector[D1] x[N];
vector[N] z;
int<lower = 1, upper = L> ll[N]; 
int<lower = 0, upper = 1> R[N]; 
int<lower = 1, upper = N> ii_obs[N_obs]; 
int<lower = 1, upper = N> ii_mis[N_mis];
}
parameters {
real mu;
vector[D1] beta; real<lower = 0> sigma; vector[L] u;
real c;
real alpha;
vector[D1] gamma; vector[N_mis] y_mis;
}
transformed parameters {
vector[N] y; vector[N] mu_y; vector[N] pi; y[ii_obs] = y_obs; y[ii_mis] = y_mis; 
for (n in 1:N) {

mu_y[n] = mu + beta' * x[n] + u[ll[n]] * z[n]; }
for (n in 1:N) {
pi[n] = c + alpha * y[n] + gamma' * x[n];
} }
model {
for (n in 1:N) {
y[n] ~ normal(mu_y[n], sigma); }
for (n in 1:N) {
R[n] ~ bernoulli_logit(pi[n]);
} }