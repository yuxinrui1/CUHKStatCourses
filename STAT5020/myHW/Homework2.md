# Homework2

## Question1

### a. Set true values for the model parameters. Generate data from the model and conduct Bayesian analysis on the basis of 10 replication.

$\begin{gathered}{\left[\begin{array}{l}y_{i 1} \\ y_{i 2} \\ y_{i 3} \\ y_{i 4} \\ y_{i 5} \\ y_{i 6} \\ y_{i 7} \\ y_{i 8} \\ y_{i 9}\end{array}\right]=\left[\begin{array}{ll}\mu_1 & a_1 \\ \mu_2 & a_2 \\ \mu_3 & a_3 \\ \mu_4 & a_4 \\ \mu_5 & a_5 \\ \mu_6 & a_6 \\ \mu_7 & a_7 \\ \mu_8 & a_8 \\ \mu_9 & a_9\end{array}\right]\left[\begin{array}{c}1 \\ c_i\end{array}\right]+\left[\begin{array}{ccc}1 & 0 & 0 \\ \lambda_{21} & 0 & 0 \\ \lambda_{31} & 0 & 0 \\ 0 & 1 & 0 \\ 0 & \lambda_{52} & 0 \\ 0 & \lambda_{62} & 0 \\ 0 & 0 & 1 \\ 0 & 0 & \lambda_{83} \\ 0 & 0 & \lambda_{93}\end{array}\right]\left[\begin{array}{c}\eta_i \\ \xi_{i 1} \\ \xi_{i 2}\end{array}\right]+\left[\begin{array}{c}\varepsilon_{i 1} \\ \varepsilon_{i 2} \\ \varepsilon_{i 3} \\ \varepsilon_{i 4} \\ \varepsilon_{i 5} \\ \varepsilon_{i 6} \\ \varepsilon_{i 7} \\ \varepsilon_{i 8} \\ \varepsilon_{i 9}\end{array}\right]} \\ \eta_i=b d_i+\left[\begin{array}{lll}\gamma_1 & \gamma_2 & \gamma_3 &\gamma_4\end{array}\right]\left[\begin{array}{c}\xi_{i 1} \\ \xi_{i 2} \\ \xi_{i 1} \xi_{i 2} \\ \xi_{i 2}^2\end{array}\right]+\delta_i\end{gathered}$

The true values of parameters set for this question are listed as follow, and 10 data sets are generated based on the true parameters. 

```
mu <- c(3.0, 1.5, 2.0, 1.0, 2.5, 1.8, 3.2, 2.3, 2.8)
a <- c(0.8, 0.6, 0.7, 0.5, 0.9, 0.8, 1.0, 0.9, 0.7)
lambda <- c(0.6, 0.7, 0.4, 0.5, 0.8, 0.6) ##lambda_21, 31, 52, 62, 83, 93
b <- 1.2
gamma <- c(0.4, 0.6, 0.2, 0.3)
```

### b. Demonstrate how to check convergence of the model.

1. Check the Rhat of the 10 replications. If Rhat is close to 1, then the model converges well, otherwise it does not converge. In the following 10 replications, all estimations converge with Rhat close to 1.
   ![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-17-23-03-image.png)
2. Check the estimation process plot of chains. If different chains meet together as the iteration number grows, the model converges.
   
   Below are two figures depict the $\mu$ and $\gamma$ estimation process of 4 chains. We can see that the model converges.

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-17-25-07-image.png)

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-17-25-19-image.png)

### c. Use Bias and RMSE to summarize the estimation results.

The mean of posterior means, estimation bias and RMSE are listed below.

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-17-48-09-image.png)

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-17-47-32-image.png)



### d. Show your prior inputs and check whether the Bayesian analysis is sensitive to the inputs.

My prior inputs are as follows:

```
  mu ~ normal(0, 5);
  a ~ normal(0, 5);
  lambda ~ normal(0, 5);
  b ~ normal(0, 5);
  gamma ~ normal(0, 5);
  sigma_eta ~ cauchy(0, 5);
  sigma_eps ~ cauchy(0, 5);
```

To check whether the Bayesian analysis is sensitive to the inputs, we modify the prior inputs to:

```
  mu ~ normal(10, 5);
  a ~ normal(10, 5);
  lambda ~ normal(10, 5);
  b ~ normal(10, 5);
  gamma ~ normal(10, 5);
  sigma_eta ~ cauchy(10, 5);
  sigma_eps ~ cauchy(10, 5);
```

The model with new priori is in `HW3m2.stan` file. The metrics are listed below:

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-11-06-image.png)

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-10-31-image.png)

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-11-23-image.png)

From the Rhat, bias, and RMSE of the estimation, we can conclude that the convergence and estimation precision are not affected by the different prior setting.

## Question2

### a. Compare the non-linear SEM in Q1 with its linear SEM counterpart.

The new non-linear model is in `HW3m3.stan` file. After we construct the linear SEM model, the metric is listed below:
![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-13-16-image.png)

#### Compare Bayes factor and DIC

The average DIC of the non-linear model is 7756.662, and the average DIC of the linear model is 7443.317. The average Bayes factor of non-linear vs. linear is 44.16. 

This indicates that the Bayes factor and DIC are in favor of the original nonlinear model since the Bayes factor is pretty large. 

### b.Compare the non-linear SEM in Q1 with this new model.

The new non-linear model is in `HW3m4.stan` file. The Bayes factor is 10.2, and the average DIC of the new model is 7749.231. This indicates that the Bayes factor and the DIC are in favor of the original nonlinear model since the Bayes factor is large. 



## Question3

The true parameters are:

```
b <- 0.3
vec_lambda <- c(1, 0.8, 0.8, 1, 0.7, 0.9, 0.7, 1, 0.8)
vec_gamma <- c(0.2, 0.5)
mat_phi <- matrix(data = c(1, 0.2, 0.2, 0.81), ncol = 2)
vec_psi_eps <- c(1, 1, 1, 0.3, 0.3, 0.4, 0.4)
psi_delta <- 0.36
```

The posterior mean of parameters are:

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-28-28-image.png) 

![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-28-43-image.png)

The bias and RMSE of parameter estimations are:
![](/Users/jiqi/Library/Application%20Support/marktext/images/2024-04-10-18-29-17-image.png)

From the bias and RMSE of the estimation, we can conclude that the model estimation is robust to different simulated data.
