
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/ethan-alt/hdbayes/workflows/R-CMD-check/badge.svg)](https://github.com/ethan-alt/hdbayes/actions)
<!-- badges: end -->

## hdbayes

The goal of `hdbayes` is to make it easier for users to conduct Bayesian
analysis methods that leverage historical data.

## Setting up MCMC parameters

We begin by setting up parameters for the MCMC. We detect the number of
cores on our machine as well as set up the warmup and total number of
desired samples.

``` r
## obtain number of cores
ncores = max(1, parallel::detectCores() - 1)
warmup  = 1000          ## warmup for MCMC sampling
total.samples = 10000   ## number of samples post warmup
samples = ceiling(warmup + total.samples / ncores)  ## outputs approx total.samples samples
```

## Simulating Poisson regression data

We now simulate some Poisson data.

``` r
## simulate logistic regression data
set.seed(391)
n  = 200
n0 = 100
N  = n + n0

X    = cbind(1, 'z' = rbinom(N, 1, 0.5), 'x' = rnorm(N, mean = 1, sd = 1) )
beta = c(1, 0.5, -1)
mean = poisson('log')$linkinv(X %*% beta)
y    = rpois(n = N, lambda = mean)

## create current and historical data sets
data = data.frame('y' = y, 'z' = X[, 'z'], 'x' = X[, 'x'])
histdata = data[1:n0, ]
data     = data[-(1:n0), ]
```

## Frequentist analysis of the data sets

We use the `stats::glm` function to conduct a frequentist GLM of the
data sets

``` r
formula = y ~ z + x
family  = poisson('log')
fit.mle.cur  = glm(formula, family, data)
fit.mle.hist = glm(formula, family, histdata)

summary(fit.mle.cur)
#> 
#> Call:
#> glm(formula = formula, family = family, data = data)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -2.1527  -1.0211  -0.3206   0.5750   2.0951  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.95408    0.08015  11.904  < 2e-16 ***
#> z            0.54353    0.11197   4.854 1.21e-06 ***
#> x           -1.04626    0.05938 -17.621  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 512.24  on 199  degrees of freedom
#> Residual deviance: 204.72  on 197  degrees of freedom
#> AIC: 560.36
#> 
#> Number of Fisher Scoring iterations: 5
```

We will compare these with the Bayesian analysis results later

## Bayesian analysis methods

We now utilize the functions in this package.

### Bayesian hierarchical model

The Bayesian hierarchical model (BHM) is the following model:

![
\\begin{align\*}
  y\_i \| x\_i, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta\_0 &\\sim \\text{Poisson}\\left( \\exp(x\_{0i}'\\beta\_0) \\right) \\\\
  \\beta, \\beta\_0 &\\sim N\_p(\\mu, \\Sigma) \\\\
  \\mu &\\sim N\_p(\\mu\_0, \\Sigma\_0) \\\\
  \\Sigma &\\sim \\text{IW}\_p(\\nu\_0, \\Psi\_0)
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta%2C%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%5C%5C%0A%20%20%5CSigma%20%26%5Csim%20%5Ctext%7BIW%7D_p%28%5Cnu_0%2C%20%5CPsi_0%29%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Poisson}\left( \exp(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Poisson}\left( \exp(x_{0i}'\beta_0) \right) \\
  \beta, \beta_0 &\sim N_p(\mu, \Sigma) \\
  \mu &\sim N_p(\mu_0, \Sigma_0) \\
  \Sigma &\sim \text{IW}_p(\nu_0, \Psi_0)
\end{align*}
")

where
![\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
is the vector of regression coefficients of the current data set,
![\\beta\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0")
is the vector of regression coefficients for the historical data set,
![\\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
is the common prior mean of
![\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
and
![\\beta\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0"),
which is treated as random with a normal hyperprior having mean
![\\mu\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_0 "\mu_0"),
and covariance
![\\Sigma\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_0 "\Sigma_0"),
and
![\\Sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma "\Sigma")
is also treated as random, having an inverse-Wishart hyperprior with
![\\nu\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_0 "\nu_0")
degrees of freedom and scale matrix
![\\Psi\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPsi_0 "\Psi_0").

The defaults in `hdbayes` are

-   ![\\mu\_0 = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_0%20%3D%200 "\mu_0 = 0")
-   ![\\Sigma\_0 = I\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_0%20%3D%20I_p "\Sigma_0 = I_p")
-   ![\\nu\_0 = p + 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_0%20%3D%20p%20%2B%2010 "\nu_0 = p + 10")
-   ![\\Psi\_0 = I\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPsi_0%20%3D%20I_p "\Psi_0 = I_p")
    where
    ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
    is the number of predictors (including the intercept if applicable).

We fit this model as follows

    #>                      mean se_mean    sd     2.5%      25%      50%      75%
    #> (Intercept)         0.957   0.001 0.077    0.803    0.904    0.959    1.009
    #> z                   0.529   0.001 0.108    0.317    0.456    0.528    0.602
    #> x                  -1.042   0.001 0.059   -1.157   -1.081   -1.042   -1.003
    #> (Intercept)_hist    1.042   0.001 0.098    0.844    0.978    1.043    1.108
    #> z_hist              0.407   0.002 0.132    0.150    0.318    0.407    0.494
    #> x_hist             -1.018   0.001 0.080   -1.174   -1.073   -1.019   -0.964
    #> beta_mean[1]        0.953   0.003 0.224    0.507    0.809    0.953    1.101
    #> beta_mean[2]        0.446   0.003 0.235   -0.030    0.293    0.449    0.600
    #> beta_mean[3]       -0.982   0.002 0.222   -1.413   -1.129   -0.986   -0.840
    #> beta_cov[1,1]       0.102   0.001 0.050    0.043    0.068    0.090    0.121
    #> beta_cov[1,2]      -0.001   0.000 0.035   -0.070   -0.018    0.000    0.017
    #> beta_cov[1,3]       0.000   0.000 0.034   -0.068   -0.017    0.000    0.017
    #> beta_cov[2,1]      -0.001   0.000 0.035   -0.070   -0.018    0.000    0.017
    #> beta_cov[2,2]       0.101   0.001 0.049    0.043    0.068    0.089    0.119
    #> beta_cov[2,3]       0.000   0.000 0.033   -0.064   -0.017    0.000    0.017
    #> beta_cov[3,1]       0.000   0.000 0.034   -0.068   -0.017    0.000    0.017
    #> beta_cov[3,2]       0.000   0.000 0.033   -0.064   -0.017    0.000    0.017
    #> beta_cov[3,3]       0.100   0.001 0.049    0.043    0.068    0.088    0.120
    #> lp__             -384.264   0.048 3.037 -391.246 -386.052 -383.913 -382.054
    #>                     97.5%     n_eff  Rhat
    #> (Intercept)         1.103  6479.859 1.000
    #> z                   0.741  6331.749 1.000
    #> x                  -0.930 10270.598 1.000
    #> (Intercept)_hist    1.232  7801.067 1.000
    #> z_hist              0.666  7516.881 1.000
    #> x_hist             -0.864 10579.902 1.000
    #> beta_mean[1]        1.391  7318.793 1.000
    #> beta_mean[2]        0.903  7561.702 1.001
    #> beta_mean[3]       -0.530  7919.533 1.000
    #> beta_cov[1,1]       0.232  6439.783 1.000
    #> beta_cov[1,2]       0.068  5420.873 1.000
    #> beta_cov[1,3]       0.071  6392.890 1.001
    #> beta_cov[2,1]       0.068  5420.873 1.000
    #> beta_cov[2,2]       0.230  5519.679 1.001
    #> beta_cov[2,3]       0.066  6135.939 1.001
    #> beta_cov[3,1]       0.071  6392.890 1.001
    #> beta_cov[3,2]       0.066  6135.939 1.001
    #> beta_cov[3,3]       0.225  5450.707 1.002
    #> lp__             -379.396  3958.805 1.002

### Commensurate prior

The commensurate prior assumes the following hierarchical model

![
\\begin{align\*}
  y\_i \| x\_i, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta\_0 &\\sim \\text{Poisson}\\left( \\exp(x\_{0i}'\\beta\_0) \\right) \\\\
  \\beta\_0 &\\sim N\_p(\\mu\_0, \\Sigma\_0) \\\\
  \\beta\_j &\\sim N\_1\\left( \\beta\_{0j}, \\tau\_j^{-1} \\right), j = 1, \\ldots, p
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%5C%5C%0A%20%20%5Cbeta_j%20%26%5Csim%20N_1%5Cleft%28%20%5Cbeta_%7B0j%7D%2C%20%5Ctau_j%5E%7B-1%7D%20%5Cright%29%2C%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Poisson}\left( \exp(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Poisson}\left( \exp(x_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu_0, \Sigma_0) \\
  \beta_j &\sim N_1\left( \beta_{0j}, \tau_j^{-1} \right), j = 1, \ldots, p
\end{align*}
")

where the
![\\tau\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctau_j "\tau_j")’s
are elicited by the user. The defaults in `hdbayes` are

-   ![\\mu\_0 = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_0%20%3D%200 "\mu_0 = 0")
-   ![\\Sigma\_0 = 100 \\times I\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_0%20%3D%20100%20%5Ctimes%20I_p "\Sigma_0 = 100 \times I_p")

This method can be fit as follows

``` r
fit.commensurate = glm.commensurate(
  formula, family, data, histdata, tau = rep(5, 3),
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.commensurate)$summary, 3 )
#>                      mean se_mean    sd     2.5%      25%      50%      75%
#> (Intercept)         0.957   0.001 0.078    0.800    0.905    0.959    1.011
#> z                   0.534   0.001 0.110    0.317    0.458    0.534    0.608
#> x                  -1.045   0.001 0.058   -1.159   -1.084   -1.044   -1.005
#> (Intercept)_hist    1.050   0.001 0.099    0.850    0.984    1.051    1.117
#> z_hist              0.400   0.002 0.134    0.140    0.310    0.400    0.489
#> x_hist             -1.023   0.001 0.080   -1.182   -1.077   -1.022   -0.968
#> lp__             -422.388   0.025 1.727 -426.642 -423.302 -422.057 -421.138
#>                     97.5%     n_eff  Rhat
#> (Intercept)         1.106  7010.399 1.000
#> z                   0.749  6984.113 1.000
#> x                  -0.932 10515.411 1.000
#> (Intercept)_hist    1.241  7917.160 1.000
#> z_hist              0.665  7615.888 1.001
#> x_hist             -0.869 10695.038 1.000
#> lp__             -420.011  4667.489 1.002
```

### Robust Meta-Analytic Predictive (MAP) Prior

The Robust MAP prior is a generalization of the Bayesian Hierarchical
Model (BHM), and takes the form

![
\\begin{align\*}
  y\_i \| x\_i, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta\_0 &\\sim \\text{Poisson}\\left( \\exp(x\_{0i}'\\beta\_0) \\right) \\\\
  \\beta\_0 &\\sim N\_p(\\mu, \\Sigma) \\\\
  \\beta   &\\sim w \\times N\_p(\\mu, \\Sigma) + (1 - w) N\_p(\\mu\_v, \\Sigma\_v) \\\\
  \\mu &\\sim N\_p(\\mu\_0, \\Sigma\_0) \\\\
  \\Sigma &\\sim \\text{IW}\_p(\\nu\_0, \\Psi\_0)
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%20%5C%5C%0A%20%20%5Cbeta%20%20%20%26%5Csim%20w%20%5Ctimes%20N_p%28%5Cmu%2C%20%5CSigma%29%20%2B%20%281%20-%20w%29%20N_p%28%5Cmu_v%2C%20%5CSigma_v%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%5C%5C%0A%20%20%5CSigma%20%26%5Csim%20%5Ctext%7BIW%7D_p%28%5Cnu_0%2C%20%5CPsi_0%29%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Poisson}\left( \exp(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Poisson}\left( \exp(x_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu, \Sigma) \\
  \beta   &\sim w \times N_p(\mu, \Sigma) + (1 - w) N_p(\mu_v, \Sigma_v) \\
  \mu &\sim N_p(\mu_0, \Sigma_0) \\
  \Sigma &\sim \text{IW}_p(\nu_0, \Psi_0)
\end{align*}
")

where
![w \\in (0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%5Cin%20%280%2C1%29 "w \in (0,1)")
controls for the level of borrowing of the historical data. Note that
when
![w = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%3D%201 "w = 1"),
the robust MAP prior effectively becomes the BHM. The defaults are the
same as in the BHM except the default value for w is 0.1.

``` r
fit.robustmap = glm.robustmap(
  formula, family, data, histdata,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.robustmap)$summary, 3 )
#>                      mean se_mean    sd     2.5%      25%      50%      75%
#> (Intercept)         0.958   0.001 0.078    0.806    0.905    0.958    1.010
#> z                   0.529   0.001 0.109    0.318    0.455    0.529    0.604
#> x                  -1.042   0.001 0.059   -1.158   -1.083   -1.043   -1.002
#> (Intercept)_hist    1.044   0.001 0.099    0.845    0.979    1.047    1.111
#> z_hist              0.404   0.002 0.134    0.145    0.314    0.404    0.494
#> x_hist             -1.019   0.001 0.080   -1.178   -1.073   -1.019   -0.965
#> coef_mean[1]        0.957   0.003 0.228    0.495    0.813    0.963    1.104
#> coef_mean[2]        0.442   0.003 0.238   -0.036    0.288    0.445    0.598
#> coef_mean[3]       -0.981   0.002 0.221   -1.400   -1.125   -0.986   -0.841
#> coef_cov[1,1]       0.101   0.001 0.048    0.043    0.068    0.089    0.119
#> coef_cov[1,2]      -0.001   0.000 0.033   -0.069   -0.018   -0.001    0.017
#> coef_cov[1,3]       0.000   0.000 0.034   -0.070   -0.017   -0.001    0.017
#> coef_cov[2,1]      -0.001   0.000 0.033   -0.069   -0.018   -0.001    0.017
#> coef_cov[2,2]       0.102   0.001 0.050    0.044    0.069    0.090    0.120
#> coef_cov[2,3]      -0.001   0.000 0.034   -0.071   -0.018   -0.001    0.016
#> coef_cov[3,1]       0.000   0.000 0.034   -0.070   -0.017   -0.001    0.017
#> coef_cov[3,2]      -0.001   0.000 0.034   -0.071   -0.018   -0.001    0.016
#> coef_cov[3,3]       0.101   0.001 0.049    0.044    0.068    0.089    0.120
#> lp__             -392.070   0.051 3.045 -399.194 -393.829 -391.700 -389.863
#>                     97.5%     n_eff  Rhat
#> (Intercept)         1.109  6552.829 1.001
#> z                   0.742  6834.968 1.001
#> x                  -0.924 10514.495 1.000
#> (Intercept)_hist    1.231  7371.996 1.001
#> z_hist              0.666  7296.883 1.001
#> x_hist             -0.862 10471.413 0.999
#> coef_mean[1]        1.392  7206.707 1.000
#> coef_mean[2]        0.899  7574.338 1.000
#> coef_mean[3]       -0.530  8022.096 1.001
#> coef_cov[1,1]       0.227  6390.709 1.001
#> coef_cov[1,2]       0.066  5923.615 1.001
#> coef_cov[1,3]       0.069  6003.206 1.000
#> coef_cov[2,1]       0.066  5923.615 1.001
#> coef_cov[2,2]       0.231  5976.786 1.000
#> coef_cov[2,3]       0.067  5944.263 1.001
#> coef_cov[3,1]       0.069  6003.206 1.000
#> coef_cov[3,2]       0.067  5944.263 1.001
#> coef_cov[3,3]       0.229  4656.111 1.000
#> lp__             -387.273  3609.124 1.001
```

### Power prior

The Power Prior takes the form

![
\\begin{align\*}
  y\_i \| x\_i, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_{0i}'\\beta) \\right) \\\\
  \\pi(\\beta \| a\_0) &\\propto L(\\beta \| y\_0)^{a\_0} \\pi\_0(\\beta)
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cpi%28%5Cbeta%20%7C%20a_0%29%20%26%5Cpropto%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Poisson}\left( \exp(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta &\sim \text{Poisson}\left( \exp(x_{0i}'\beta) \right) \\
  \pi(\beta | a_0) &\propto L(\beta | y_0)^{a_0} \pi_0(\beta)
\end{align*}
")

where
![L(\\beta \| y\_0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;L%28%5Cbeta%20%7C%20y_0%29 "L(\beta | y_0)")
is the likelihood of the GLM based on the historical data,
![a\_0 \\in (0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0%20%5Cin%20%280%2C1%29 "a_0 \in (0,1)")
is a fixed hyperaparameter controlling the effective sample size
contributed by the data (e.g.,
![a\_0 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0%20%3D%201 "a_0 = 1")
borrows the whole sample size), and
![\\pi\_0(\\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpi_0%28%5Cbeta%29 "\pi_0(\beta)")
is an “initial prior” on
![\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta").

The default in `hdbayes` is a (noninformative) normal prior on
![\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta"):

![
\\beta \\sim N\_p(0, 100 \\times I\_p)
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbeta%20%5Csim%20N_p%280%2C%20100%20%5Ctimes%20I_p%29%0A "
\beta \sim N_p(0, 100 \times I_p)
")

The power prior (with
![a\_0 = 0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0%20%3D%200.5 "a_0 = 0.5"))
may be fit as follows:

``` r
fit.pp = glm.pp(
  formula, family, data, histdata, a0 = 0.5,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.pp)$summary, 3 )
#>                 mean se_mean    sd     2.5%      25%      50%      75%    97.5%
#> (Intercept)    0.976   0.001 0.070    0.838    0.930    0.977    1.023    1.111
#> z              0.505   0.001 0.097    0.318    0.440    0.506    0.572    0.696
#> x             -1.038   0.001 0.052   -1.140   -1.073   -1.037   -1.003   -0.936
#> lp__        -350.029   0.019 1.236 -353.283 -350.578 -349.707 -349.136 -348.641
#>                n_eff  Rhat
#> (Intercept) 5324.335 1.002
#> z           5121.288 1.001
#> x           5979.019 1.000
#> lp__        4166.155 1.002
```

### Normalized power prior (NPP)

The NPP treats the hyperparameter
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0")
as random, allowing the data to decide what is the best value. For
non-Gaussian models, this requires estimating the normalizing constant
![Z(a\_0) = \\int L(\\beta \| y\_0)^{a\_0} \\pi\_0(\\beta) d\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z%28a_0%29%20%3D%20%5Cint%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%20d%5Cbeta "Z(a_0) = \int L(\beta | y_0)^{a_0} \pi_0(\beta) d\beta").

In `hdbayes`, there is one function to estimate the normalizing constant
across a grid of values for
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0")
and another to obtain posterior samples of the normalized power prior.

The NPP may be summarized as

![
\\begin{align\*}
  y\_i \| x\_i, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta &\\sim \\text{Poisson}\\left( \\exp(x\_{0i}'\\beta) \\right) \\\\
  \\pi(\\beta \| a\_0) &\\propto \\frac{1}{Z(a\_0)} L(\\beta \| y\_0)^{a\_0} \\pi\_0(\\beta) \\\\
  \\pi(a\_0)         &\\propto a\_0^{\\alpha\_0 - 1} (1 - a\_0)^{\\gamma\_0 - 1}
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BPoisson%7D%5Cleft%28%20%5Cexp%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cpi%28%5Cbeta%20%7C%20a_0%29%20%26%5Cpropto%20%5Cfrac%7B1%7D%7BZ%28a_0%29%7D%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%20%5C%5C%0A%20%20%5Cpi%28a_0%29%20%20%20%20%20%20%20%20%20%26%5Cpropto%20a_0%5E%7B%5Calpha_0%20-%201%7D%20%281%20-%20a_0%29%5E%7B%5Cgamma_0%20-%201%7D%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Poisson}\left( \exp(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta &\sim \text{Poisson}\left( \exp(x_{0i}'\beta) \right) \\
  \pi(\beta | a_0) &\propto \frac{1}{Z(a_0)} L(\beta | y_0)^{a_0} \pi_0(\beta) \\
  \pi(a_0)         &\propto a_0^{\alpha_0 - 1} (1 - a_0)^{\gamma_0 - 1}
\end{align*}
")

The defaults in `hdbayes` are

-   ![\\pi\_0(\\beta) \\propto N(\\beta \| 0, 100 \\times I\_p)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpi_0%28%5Cbeta%29%20%5Cpropto%20N%28%5Cbeta%20%7C%200%2C%20100%20%5Ctimes%20I_p%29 "\pi_0(\beta) \propto N(\beta | 0, 100 \times I_p)")
-   ![\\alpha\_0 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_0%20%3D%201 "\alpha_0 = 1")
-   ![\\gamma\_0 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_0%20%3D%201 "\gamma_0 = 1")

when
![\\alpha\_0 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_0%20%3D%201 "\alpha_0 = 1")
and
![\\gamma\_0 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_0%20%3D%201 "\gamma_0 = 1"),
the prior on
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0")
is a
![U(0,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;U%280%2C1%29 "U(0,1)")
prior.

### Estimating the normalizing constant

We begin by estimating the normalizing constant over a fine grid of
values. On machines with multiple cores, this may be accomplished more
efficiently using parallel computing.

``` r
## parallelize estimation of log normalizing constant
library(parallel)
a0     = seq(0, 1, length.out = ncores * 5)

## wrapper to obtain log normalizing constant in parallel package
logncfun = function(a0, ...)
  hdbayes::glm.npp.lognc(
    formula = formula, family = family, histdata = histdata, a0 = a0, ...
  )

cl = makeCluster(ncores)
  clusterSetRNGStream(cl, 123)
  clusterExport(cl, varlist = c('formula', 'family', 'histdata'))

  ## call created function
  a0.lognc = parLapply(
    cl = cl, X = a0, fun = logncfun, iter = 5000, warmup = warmup, refresh = 0
  )
stopCluster(cl)

a0.lognc = data.frame( do.call(rbind, a0.lognc) )
head(a0.lognc)
#>           a0       lognc min_n_eff max_Rhat
#> 1 0.00000000   9.6642941  8391.864 1.000460
#> 2 0.01851852  -0.9177977  5370.667 1.000581
#> 3 0.03703704  -4.6320434  5741.243 1.000498
#> 4 0.05555556  -7.8908335  6342.905 1.000519
#> 5 0.07407407 -10.9560923  6385.807 1.000492
#> 6 0.09259259 -13.9331685  6761.353 1.000032
```

The provided function `glm.npp.lognc` estimates the logarithm of the
normalizing constant,
![\\log Z(a\_0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clog%20Z%28a_0%29 "\log Z(a_0)"),
for one specific value of
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0").
We created the function `logncfun` so that the first argument would be
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0"),
allowing us to use the `parLapply` function in the `parallel` package.

The `hdbayes` function `glm.npp.lognc` outputs
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0"),
![Z(a\_0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z%28a_0%29 "Z(a_0)"),
and the minimum effective sample size and maximum R-hat value of the
MCMC sampling of the power prior. It is a good idea to check that the
minimum effective sample size is at least 1,000 and the maximum R-hat
value is less than 1.10

``` r
min(a0.lognc$min_n_eff) ## lowest effective sample size
#> [1] 5370.667
max(a0.lognc$max_Rhat)  ## highest R-hat value
#> [1] 1.001475
```

We can then plot the logarithm of the normalizing constant

``` r
plot(a0.lognc$a0, a0.lognc$lognc)
```

![](README_files/figure-gfm/npp_lognc_plot-1.png)<!-- -->

### Sampling the posterior distribution

We can now sample from the posterior distribution. The function
`glm.npp` takes, as input, values of
![a\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_0 "a_0")
and the estimated logarithm of the normalizing constant. Linear
interpolation is used to estimate
![Z(a\_0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z%28a_0%29 "Z(a_0)")
for values not in the fine grid. Thus, it may be a good idea to conduct
smoothing of the function such as using LOESS, but we ignore that here.

``` r
fit.npp = glm.npp(
  formula, family, data, histdata, a0 = a0.lognc$a0, lognc = a0.lognc$lognc,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.npp)$summary, 3 )
#>                 mean se_mean    sd     2.5%      25%      50%      75%    97.5%
#> (Intercept)    0.982   0.001 0.069    0.841    0.937    0.985    1.029    1.112
#> z              0.498   0.001 0.096    0.310    0.433    0.495    0.561    0.687
#> x             -1.037   0.001 0.051   -1.137   -1.070   -1.037   -1.002   -0.937
#> a0             0.667   0.002 0.233    0.180    0.498    0.705    0.865    0.987
#> lp__        -277.415   0.024 1.477 -281.140 -278.126 -277.067 -276.338 -275.564
#>                n_eff  Rhat
#> (Intercept) 6122.080 1.000
#> z           5530.929 1.001
#> x           8599.024 1.000
#> a0          9312.138 1.000
#> lp__        3849.028 1.001
```

### Normalized asymptotic power prior (NAPP)

NAPP uses a large sample theory argument to formulate a normal
approximation to the power prior, i.e., the prior is given by

![
\\beta \| a\_0 \\sim N(\\hat{\\beta}\_0, a\_0^{-1} \[I\_n(\\beta)\]^{-1}),
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbeta%20%7C%20a_0%20%5Csim%20N%28%5Chat%7B%5Cbeta%7D_0%2C%20a_0%5E%7B-1%7D%20%5BI_n%28%5Cbeta%29%5D%5E%7B-1%7D%29%2C%0A "
\beta | a_0 \sim N(\hat{\beta}_0, a_0^{-1} [I_n(\beta)]^{-1}),
")

where
![\\hat{\\beta}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Cbeta%7D_0 "\hat{\beta}_0")
is the maximum likelihood estimate (MLE) of
![\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
based on the historical data and
![I\_n(\\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I_n%28%5Cbeta%29 "I_n(\beta)")
is the associated information matrix (negative Hessian).

In this case, the normalizing constant is known, so we do not need to
estimate the normalizing constant before sampling.

The following code analyzes the data set using the NAPP:

``` r
fit.napp = glm.napp(
  formula, family, data, histdata,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.napp)$summary, 3 )
#>                 mean se_mean    sd     2.5%      25%      50%      75%    97.5%
#> (Intercept)    0.983   0.001 0.068    0.847    0.938    0.984    1.030    1.114
#> z              0.497   0.001 0.095    0.311    0.433    0.497    0.563    0.680
#> x             -1.036   0.001 0.051   -1.135   -1.071   -1.036   -1.002   -0.936
#> a0             0.659   0.002 0.231    0.177    0.486    0.693    0.857    0.984
#> lp__        -277.400   0.023 1.455 -281.071 -278.080 -277.087 -276.336 -275.582
#>                n_eff  Rhat
#> (Intercept) 5967.763 1.001
#> z           5886.828 1.001
#> x           8855.402 1.000
#> a0          8552.108 1.001
#> lp__        3918.570 1.001
```

## Comparison of methods

We now can compare the point estimate (MLE / posterior mean) and
uncertainty (SE / posterior standard deviation) of the methods.

``` r
##
## COMPARE METHODS
##
fit.list = list('bhm' = fit.bhm, 'commensurate' = fit.commensurate,
                'robustmap' = fit.robustmap, 'napp' = fit.napp,
                'npp' = fit.npp, 'pp' = fit.pp)

post.mean = sapply(
  fit.list, function(x) summary(x)$summary[names(coef(fit.mle.cur)), 'mean']
)
post.sd = sapply(
  fit.list, function(x) summary(x)$summary[names(coef(fit.mle.cur)), 'sd']
)

post.mean = cbind(
  'truth' = beta,
  'mle.cur' = coef(fit.mle.cur),
  'mle.hist' = coef(fit.mle.hist),
  post.mean
)

post.sd = cbind(
  'mle.cur'  = summary(fit.mle.cur)$coefficients[, 'Std. Error'],
  'mle.hist' = summary(fit.mle.hist)$coefficients[, 'Std. Error'],
  post.sd
)

## posterior means
round( post.mean, 3 )
#>             truth mle.cur mle.hist    bhm commensurate robustmap   napp    npp
#> (Intercept)   1.0   0.954    1.067  0.957        0.957     0.958  0.983  0.982
#> z             0.5   0.544    0.384  0.529        0.534     0.529  0.497  0.498
#> x            -1.0  -1.046   -1.021 -1.042       -1.045    -1.042 -1.036 -1.037
#>                 pp
#> (Intercept)  0.976
#> z            0.505
#> x           -1.038

## posterior std dev.
round( post.sd, 3 )
#>             mle.cur mle.hist   bhm commensurate robustmap  napp   npp    pp
#> (Intercept)   0.080    0.103 0.077        0.078     0.078 0.068 0.069 0.070
#> z             0.112    0.142 0.108        0.110     0.109 0.095 0.096 0.097
#> x             0.059    0.082 0.059        0.058     0.059 0.051 0.051 0.052
```
