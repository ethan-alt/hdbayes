
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

## Simulating logistic regression data

We now simulate some logistic regression data.

``` r
## simulate logistic regression data
set.seed(391)
n  = 200
n0 = 100
N  = n + n0

X    = cbind(1, 'z' = rbinom(N, 1, 0.5), 'x' = rnorm(N, mean = 1, sd = 1) )
beta = c(1, 0.5, -1)
mean = binomial('logit')$linkinv(X %*% beta)
y    = rbinom(N, 1, mean)

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
family  = binomial('logit')
fit.mle.cur  = glm(formula, family, data)
fit.mle.hist = glm(formula, family, histdata)

summary(fit.mle.cur)
#> 
#> Call:
#> glm(formula = formula, family = family, data = data)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.7907  -1.1360   0.6741   1.0005   1.9255  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   0.7600     0.2744   2.769  0.00562 ** 
#> z             0.6769     0.3078   2.199  0.02787 *  
#> x            -0.7500     0.1814  -4.134 3.56e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 273.33  on 199  degrees of freedom
#> Residual deviance: 250.60  on 197  degrees of freedom
#> AIC: 256.6
#> 
#> Number of Fisher Scoring iterations: 4
```

We will compare these with the Bayesian analysis results later

## Bayesian analysis methods

We now utilize the functions in this package.

### Bayesian hierarchical model

The Bayesian hierarchical model (BHM) is the following model:

![
\\begin{align\*}
  y_i \| x_i, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta_0 &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x\_{0i}'\\beta_0) \\right) \\\\
  \\beta, \\beta_0 &\\sim N_p(\\mu, \\Sigma) \\\\
  \\mu &\\sim N_p(\\mu_0, \\Sigma_0) \\\\
  \\Sigma &\\sim \\text{IW}\_p(\\nu_0, \\Psi_0)
\\end{align\*}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta%2C%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%5C%5C%0A%20%20%5CSigma%20%26%5Csim%20%5Ctext%7BIW%7D_p%28%5Cnu_0%2C%20%5CPsi_0%29%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_0) \right) \\
  \beta, \beta_0 &\sim N_p(\mu, \Sigma) \\
  \mu &\sim N_p(\mu_0, \Sigma_0) \\
  \Sigma &\sim \text{IW}_p(\nu_0, \Psi_0)
\end{align*}
")

where ![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") is
the vector of regression coefficients of the current data set,
![\\beta_0](https://latex.codecogs.com/png.latex?%5Cbeta_0 "\beta_0") is
the vector of regression coefficients for the historical data set,
![\\mu](https://latex.codecogs.com/png.latex?%5Cmu "\mu") is the common
prior mean of
![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") and
![\\beta_0](https://latex.codecogs.com/png.latex?%5Cbeta_0 "\beta_0"),
which is treated as random with a normal hyperprior having mean
![\\mu_0](https://latex.codecogs.com/png.latex?%5Cmu_0 "\mu_0"), and
covariance
![\\Sigma_0](https://latex.codecogs.com/png.latex?%5CSigma_0 "\Sigma_0"),
and ![\\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma")
is also treated as random, having an inverse-Wishart hyperprior with
![\\nu_0](https://latex.codecogs.com/png.latex?%5Cnu_0 "\nu_0") degrees
of freedom and scale matrix
![\\Psi_0](https://latex.codecogs.com/png.latex?%5CPsi_0 "\Psi_0").

The defaults in `hdbayes` are

-   ![\\mu_0 = 0](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%200 "\mu_0 = 0")
-   ![\\Sigma_0 = I_p](https://latex.codecogs.com/png.latex?%5CSigma_0%20%3D%20I_p "\Sigma_0 = I_p")
-   ![\\nu_0 = p + 10](https://latex.codecogs.com/png.latex?%5Cnu_0%20%3D%20p%20%2B%2010 "\nu_0 = p + 10")
-   ![\\Psi_0 = I_p](https://latex.codecogs.com/png.latex?%5CPsi_0%20%3D%20I_p "\Psi_0 = I_p")
    where ![p](https://latex.codecogs.com/png.latex?p "p") is the number
    of predictors (including the intercept if applicable).

We fit this model as follows

    #>                      mean se_mean    sd     2.5%      25%      50%      75%
    #> (Intercept)         0.805   0.003 0.240    0.337    0.647    0.805    0.965
    #> z                   0.590   0.003 0.269    0.065    0.408    0.588    0.772
    #> x                  -0.751   0.002 0.162   -1.073   -0.859   -0.750   -0.640
    #> (Intercept)_hist    0.877   0.004 0.280    0.339    0.689    0.874    1.062
    #> z_hist              0.453   0.004 0.324   -0.182    0.233    0.449    0.671
    #> x_hist             -0.800   0.002 0.209   -1.221   -0.936   -0.801   -0.660
    #> beta_mean[1]        0.803   0.004 0.303    0.205    0.603    0.803    1.007
    #> beta_mean[2]        0.497   0.004 0.326   -0.138    0.279    0.498    0.715
    #> beta_mean[3]       -0.735   0.003 0.261   -1.244   -0.909   -0.738   -0.561
    #> beta_cov[1,1]       0.105   0.001 0.052    0.043    0.070    0.092    0.125
    #> beta_cov[1,2]      -0.002   0.000 0.036   -0.075   -0.020   -0.001    0.017
    #> beta_cov[1,3]      -0.002   0.000 0.035   -0.075   -0.019   -0.002    0.016
    #> beta_cov[2,1]      -0.002   0.000 0.036   -0.075   -0.020   -0.001    0.017
    #> beta_cov[2,2]       0.106   0.001 0.052    0.045    0.071    0.094    0.126
    #> beta_cov[2,3]      -0.001   0.000 0.035   -0.071   -0.018    0.000    0.017
    #> beta_cov[3,1]      -0.002   0.000 0.035   -0.075   -0.019   -0.002    0.016
    #> beta_cov[3,2]      -0.001   0.000 0.035   -0.071   -0.018    0.000    0.017
    #> beta_cov[3,3]       0.102   0.001 0.050    0.044    0.069    0.090    0.121
    #> lp__             -152.074   0.052 3.092 -159.202 -153.910 -151.692 -149.849
    #>                     97.5%    n_eff  Rhat
    #> (Intercept)         1.282 6020.666 1.000
    #> z                   1.122 7888.703 1.000
    #> x                  -0.436 7039.384 1.001
    #> (Intercept)_hist    1.442 6059.024 1.001
    #> z_hist              1.088 7564.887 1.000
    #> x_hist             -0.394 7900.572 1.000
    #> beta_mean[1]        1.396 6299.741 1.000
    #> beta_mean[2]        1.125 6916.607 1.001
    #> beta_mean[3]       -0.221 7677.676 1.000
    #> beta_cov[1,1]       0.239 7016.053 1.000
    #> beta_cov[1,2]       0.067 6294.136 1.001
    #> beta_cov[1,3]       0.068 6232.746 1.000
    #> beta_cov[2,1]       0.067 6294.136 1.001
    #> beta_cov[2,2]       0.240 5600.817 1.001
    #> beta_cov[2,3]       0.070 5421.308 1.000
    #> beta_cov[3,1]       0.068 6232.746 1.000
    #> beta_cov[3,2]       0.070 5421.308 1.000
    #> beta_cov[3,3]       0.231 4875.843 1.002
    #> lp__             -147.047 3581.083 1.001

### Commensurate prior

The commensurate prior assumes the following hierarchical model

![
\\begin{align\*}
  y_i \| x_i, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta_0 &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x\_{0i}'\\beta_0) \\right) \\\\
  \\beta_0 &\\sim N_p(\\mu_0, \\Sigma_0) \\\\
  \\beta_j &\\sim N_1\\left( \\beta\_{0j}, \\tau_j^{-1} \\right), j = 1, \\ldots, p
\\end{align\*}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%5C%5C%0A%20%20%5Cbeta_j%20%26%5Csim%20N_1%5Cleft%28%20%5Cbeta_%7B0j%7D%2C%20%5Ctau_j%5E%7B-1%7D%20%5Cright%29%2C%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu_0, \Sigma_0) \\
  \beta_j &\sim N_1\left( \beta_{0j}, \tau_j^{-1} \right), j = 1, \ldots, p
\end{align*}
")

where the
![\\tau_j](https://latex.codecogs.com/png.latex?%5Ctau_j "\tau_j")’s are
elicited by the user. The defaults in `hdbayes` are

-   ![\\mu_0 = 0](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%200 "\mu_0 = 0")
-   ![\\Sigma_0 = 100 \\times I_p](https://latex.codecogs.com/png.latex?%5CSigma_0%20%3D%20100%20%5Ctimes%20I_p "\Sigma_0 = 100 \times I_p")

This method can be fit as follows

``` r
fit.commensurate = glm.commensurate(
  formula, family, data, histdata, tau = rep(5, 3),
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.commensurate)$summary, 3 )
#>                      mean se_mean    sd     2.5%      25%      50%      75%
#> (Intercept)         0.845   0.003 0.247    0.371    0.675    0.842    1.010
#> z                   0.611   0.003 0.279    0.064    0.423    0.611    0.801
#> x                  -0.784   0.002 0.170   -1.124   -0.898   -0.781   -0.670
#> (Intercept)_hist    0.933   0.004 0.300    0.341    0.727    0.931    1.133
#> z_hist              0.478   0.004 0.341   -0.190    0.248    0.475    0.705
#> x_hist             -0.847   0.003 0.223   -1.297   -0.998   -0.844   -0.692
#> lp__             -190.434   0.026 1.747 -194.752 -191.335 -190.117 -189.167
#>                     97.5%    n_eff  Rhat
#> (Intercept)         1.331 6681.533 1.000
#> z                   1.162 7977.236 1.000
#> x                  -0.457 7013.581 1.001
#> (Intercept)_hist    1.530 6543.495 1.001
#> z_hist              1.154 7824.225 1.000
#> x_hist             -0.423 6468.872 1.000
#> lp__             -188.006 4349.878 1.000
```

### Robust Meta-Analytic Predictive (MAP) Prior

The Robust MAP prior is a generalization of the Bayesian Hierarchical
Model (BHM), and takes the form

![
\\begin{align\*}
  y_i \| x_i, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta_0 &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x\_{0i}'\\beta_0) \\right) \\\\
  \\beta_0 &\\sim N_p(\\mu, \\Sigma) \\\\
  \\beta   &\\sim w \\times N_p(\\mu, \\Sigma) + (1 - w) N_p(\\mu_v, \\Sigma_v) \\\\
  \\mu &\\sim N_p(\\mu_0, \\Sigma_0) \\\\
  \\Sigma &\\sim \\text{IW}\_p(\\nu_0, \\Psi_0)
\\end{align\*}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%20%5C%5C%0A%20%20%5Cbeta%20%20%20%26%5Csim%20w%20%5Ctimes%20N_p%28%5Cmu%2C%20%5CSigma%29%20%2B%20%281%20-%20w%29%20N_p%28%5Cmu_v%2C%20%5CSigma_v%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%5C%5C%0A%20%20%5CSigma%20%26%5Csim%20%5Ctext%7BIW%7D_p%28%5Cnu_0%2C%20%5CPsi_0%29%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu, \Sigma) \\
  \beta   &\sim w \times N_p(\mu, \Sigma) + (1 - w) N_p(\mu_v, \Sigma_v) \\
  \mu &\sim N_p(\mu_0, \Sigma_0) \\
  \Sigma &\sim \text{IW}_p(\nu_0, \Psi_0)
\end{align*}
")

where
![w \\in (0,1)](https://latex.codecogs.com/png.latex?w%20%5Cin%20%280%2C1%29 "w \in (0,1)")
controls for the level of borrowing of the historical data. Note that
when ![w = 1](https://latex.codecogs.com/png.latex?w%20%3D%201 "w = 1"),
the robust MAP prior effectively becomes the BHM. The defaults are the
same as in the BHM except the default value for w is 0.1.

``` r
fit.robustmap = glm.robustmap(
  formula, family, data, histdata,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.robustmap)$summary, 3 )
#>                      mean se_mean    sd     2.5%      25%      50%      75%
#> (Intercept)         0.810   0.003 0.238    0.352    0.651    0.805    0.970
#> z                   0.590   0.003 0.265    0.070    0.410    0.589    0.767
#> x                  -0.755   0.002 0.165   -1.085   -0.865   -0.752   -0.642
#> (Intercept)_hist    0.882   0.003 0.280    0.335    0.695    0.880    1.066
#> z_hist              0.454   0.004 0.324   -0.182    0.236    0.454    0.671
#> x_hist             -0.805   0.002 0.211   -1.230   -0.947   -0.799   -0.662
#> coef_mean[1]        0.803   0.004 0.303    0.192    0.605    0.805    1.005
#> coef_mean[2]        0.497   0.004 0.324   -0.137    0.278    0.499    0.717
#> coef_mean[3]       -0.741   0.003 0.261   -1.255   -0.913   -0.742   -0.571
#> coef_cov[1,1]       0.105   0.001 0.053    0.045    0.070    0.092    0.125
#> coef_cov[1,2]      -0.002   0.000 0.036   -0.076   -0.019   -0.001    0.017
#> coef_cov[1,3]      -0.003   0.000 0.036   -0.076   -0.019   -0.001    0.015
#> coef_cov[2,1]      -0.002   0.000 0.036   -0.076   -0.019   -0.001    0.017
#> coef_cov[2,2]       0.106   0.001 0.054    0.044    0.070    0.093    0.125
#> coef_cov[2,3]       0.000   0.001 0.036   -0.072   -0.018    0.000    0.017
#> coef_cov[3,1]      -0.003   0.000 0.036   -0.076   -0.019   -0.001    0.015
#> coef_cov[3,2]       0.000   0.001 0.036   -0.072   -0.018    0.000    0.017
#> coef_cov[3,3]       0.102   0.001 0.054    0.044    0.069    0.090    0.120
#> lp__             -157.079   0.054 3.124 -164.220 -158.924 -156.714 -154.833
#>                     97.5%    n_eff  Rhat
#> (Intercept)         1.279 6412.166 1.000
#> z                   1.109 7707.498 1.000
#> x                  -0.440 6751.157 1.000
#> (Intercept)_hist    1.435 6719.929 1.000
#> z_hist              1.088 7932.033 1.000
#> x_hist             -0.405 7509.958 1.001
#> coef_mean[1]        1.399 6541.830 1.000
#> coef_mean[2]        1.135 6859.063 1.000
#> coef_mean[3]       -0.231 7489.061 1.001
#> coef_cov[1,1]       0.237 6235.431 1.000
#> coef_cov[1,2]       0.066 5709.664 1.000
#> coef_cov[1,3]       0.065 5721.097 1.000
#> coef_cov[2,1]       0.066 5709.664 1.000
#> coef_cov[2,2]       0.246 5130.365 1.000
#> coef_cov[2,3]       0.069 5156.491 1.000
#> coef_cov[3,1]       0.065 5721.097 1.000
#> coef_cov[3,2]       0.069 5156.491 1.000
#> coef_cov[3,3]       0.231 3735.297 1.002
#> lp__             -152.095 3323.737 1.002
```

### Power prior

The Power Prior takes the form

![
\\begin{align\*}
  y_i \| x_i, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x\_{0i}'\\beta) \\right) \\\\
  \\pi(\\beta \| a_0) &\\propto L(\\beta \| y_0)^{a_0} \\pi_0(\\beta)
\\end{align\*}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cpi%28%5Cbeta%20%7C%20a_0%29%20%26%5Cpropto%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta) \right) \\
  \pi(\beta | a_0) &\propto L(\beta | y_0)^{a_0} \pi_0(\beta)
\end{align*}
")

where
![L(\\beta \| y_0)](https://latex.codecogs.com/png.latex?L%28%5Cbeta%20%7C%20y_0%29 "L(\beta | y_0)")
is the likelihood of the GLM based on the historical data,
![a_0 \\in (0,1)](https://latex.codecogs.com/png.latex?a_0%20%5Cin%20%280%2C1%29 "a_0 \in (0,1)")
is a fixed hyperaparameter controlling the effective sample size
contributed by the data (e.g.,
![a_0 = 1](https://latex.codecogs.com/png.latex?a_0%20%3D%201 "a_0 = 1")
borrows the whole sample size), and
![\\pi_0(\\beta)](https://latex.codecogs.com/png.latex?%5Cpi_0%28%5Cbeta%29 "\pi_0(\beta)")
is an “initial prior” on
![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta").

The default in `hdbayes` is a (noninformative) normal prior on
![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta"):

![
\\beta \\sim N_p(0, 100 \\times I_p)
](https://latex.codecogs.com/png.latex?%0A%5Cbeta%20%5Csim%20N_p%280%2C%20100%20%5Ctimes%20I_p%29%0A "
\beta \sim N_p(0, 100 \times I_p)
")

The power prior (with
![a_0 = 0.5](https://latex.codecogs.com/png.latex?a_0%20%3D%200.5 "a_0 = 0.5"))
may be fit as follows:

``` r
fit.pp = glm.pp(
  formula, family, data, histdata, a0 = 0.5,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.pp)$summary, 3 )
#>                 mean se_mean    sd     2.5%      25%      50%      75%    97.5%
#> (Intercept)    0.834   0.004 0.245    0.362    0.668    0.829    0.996    1.325
#> z              0.611   0.004 0.285    0.068    0.420    0.602    0.798    1.202
#> x             -0.788   0.002 0.163   -1.119   -0.895   -0.785   -0.673   -0.484
#> lp__        -157.980   0.020 1.255 -161.252 -158.548 -157.653 -157.062 -156.568
#>                n_eff  Rhat
#> (Intercept) 4496.063 1.001
#> z           5362.722 1.000
#> x           4734.142 1.001
#> lp__        3788.551 1.001
```

### Normalized power prior (NPP)

The NPP treats the hyperparameter
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") as random,
allowing the data to decide what is the best value. For non-Gaussian
models, this requires estimating the normalizing constant
![Z(a_0) = \\int L(\\beta \| y_0)^{a_0} \\pi_0(\\beta) d\\beta](https://latex.codecogs.com/png.latex?Z%28a_0%29%20%3D%20%5Cint%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%20d%5Cbeta "Z(a_0) = \int L(\beta | y_0)^{a_0} \pi_0(\beta) d\beta").

In `hdbayes`, there is one function to estimate the normalizing constant
across a grid of values for
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") and another to
obtain posterior samples of the normalized power prior.

The NPP may be summarized as

![
\\begin{align\*}
  y_i \| x_i, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x_i'\\beta) \\right) \\\\
  y\_{0i} \| x\_{0i}, \\beta &\\sim \\text{Bernoulli}\\left( \\text{logit}^{-1}(x\_{0i}'\\beta) \\right) \\\\
  \\pi(\\beta \| a_0) &\\propto \\frac{1}{Z(a_0)} L(\\beta \| y_0)^{a_0} \\pi_0(\\beta) \\\\
  \\pi(a_0)         &\\propto a_0^{\\alpha_0 - 1} (1 - a_0)^{\\gamma_0 - 1}
\\end{align\*}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cpi%28%5Cbeta%20%7C%20a_0%29%20%26%5Cpropto%20%5Cfrac%7B1%7D%7BZ%28a_0%29%7D%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%20%5C%5C%0A%20%20%5Cpi%28a_0%29%20%20%20%20%20%20%20%20%20%26%5Cpropto%20a_0%5E%7B%5Calpha_0%20-%201%7D%20%281%20-%20a_0%29%5E%7B%5Cgamma_0%20-%201%7D%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta) \right) \\
  \pi(\beta | a_0) &\propto \frac{1}{Z(a_0)} L(\beta | y_0)^{a_0} \pi_0(\beta) \\
  \pi(a_0)         &\propto a_0^{\alpha_0 - 1} (1 - a_0)^{\gamma_0 - 1}
\end{align*}
")

The defaults in `hdbayes` are

-   ![\\pi_0(\\beta) \\propto N(\\beta \| 0, 100 \\times I_p)](https://latex.codecogs.com/png.latex?%5Cpi_0%28%5Cbeta%29%20%5Cpropto%20N%28%5Cbeta%20%7C%200%2C%20100%20%5Ctimes%20I_p%29 "\pi_0(\beta) \propto N(\beta | 0, 100 \times I_p)")
-   ![\\alpha_0 = 1](https://latex.codecogs.com/png.latex?%5Calpha_0%20%3D%201 "\alpha_0 = 1")
-   ![\\gamma_0 = 1](https://latex.codecogs.com/png.latex?%5Cgamma_0%20%3D%201 "\gamma_0 = 1")

when
![\\alpha_0 = 1](https://latex.codecogs.com/png.latex?%5Calpha_0%20%3D%201 "\alpha_0 = 1")
and
![\\gamma_0 = 1](https://latex.codecogs.com/png.latex?%5Cgamma_0%20%3D%201 "\gamma_0 = 1"),
the prior on ![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") is a
![U(0,1)](https://latex.codecogs.com/png.latex?U%280%2C1%29 "U(0,1)")
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
#>           a0        lognc min_n_eff max_Rhat
#> 1 0.00000000  9.451078691  4237.621 1.001168
#> 2 0.02941176  3.076353738  4720.327 1.000376
#> 3 0.05882353  0.003102988  4318.647 1.000745
#> 4 0.08823529 -2.518622980  5116.200 1.000483
#> 5 0.11764706 -4.829003994  5109.575 1.000962
#> 6 0.14705882 -7.016724008  5372.618 1.000905
```

The provided function `glm.npp.lognc` estimates the logarithm of the
normalizing constant,
![\\log Z(a_0)](https://latex.codecogs.com/png.latex?%5Clog%20Z%28a_0%29 "\log Z(a_0)"),
for one specific value of
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0"). We created the
function `logncfun` so that the first argument would be
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0"), allowing us to
use the `parLapply` function in the `parallel` package.

The `hdbayes` function `glm.npp.lognc` outputs
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0"),
![Z(a_0)](https://latex.codecogs.com/png.latex?Z%28a_0%29 "Z(a_0)"), and
the minimum effective sample size and maximum R-hat value of the MCMC
sampling of the power prior. It is a good idea to check that the minimum
effective sample size is at least 1,000 and the maximum R-hat value is
less than 1.10

``` r
min(a0.lognc$min_n_eff) ## lowest effective sample size
#> [1] 4237.621
max(a0.lognc$max_Rhat)  ## highest R-hat value
#> [1] 1.002398
```

We can then plot the logarithm of the normalizing constant

``` r
plot(a0.lognc$a0, a0.lognc$lognc)
```

![](README_files/figure-gfm/npp_lognc_plot-1.png)<!-- -->

### Sampling the posterior distribution

We can now sample from the posterior distribution. The function
`glm.npp` takes, as input, values of
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") and the estimated
logarithm of the normalizing constant. Linear interpolation is used to
estimate
![Z(a_0)](https://latex.codecogs.com/png.latex?Z%28a_0%29 "Z(a_0)") for
values not in the fine grid. Thus, it may be a good idea to conduct
smoothing of the function such as using LOESS, but we ignore that here.

``` r
fit.npp = glm.npp(
  formula, family, data, histdata, a0 = a0.lognc$a0, lognc = a0.lognc$lognc,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup, refresh = 0
)
round( summary(fit.npp)$summary, 3 )
#>                 mean se_mean    sd     2.5%      25%      50%      75%    97.5%
#> (Intercept)    0.850   0.003 0.237    0.390    0.692    0.852    1.008    1.320
#> z              0.596   0.003 0.266    0.083    0.416    0.594    0.773    1.115
#> x             -0.795   0.002 0.161   -1.121   -0.900   -0.792   -0.686   -0.485
#> a0             0.679   0.002 0.220    0.206    0.524    0.713    0.864    0.987
#> lp__        -128.846   0.023 1.480 -132.562 -129.597 -128.518 -127.749 -126.999
#>                n_eff  Rhat
#> (Intercept) 5777.536 1.000
#> z           6196.569 1.000
#> x           6028.055 1.000
#> a0          7989.241 1.001
#> lp__        4155.382 1.002
```

### Normalized asymptotic power prior (NAPP)

NAPP uses a large sample theory argument to formulate a normal
approximation to the power prior, i.e., the prior is given by

![
\\beta \| a_0 \\sim N(\\hat{\\beta}\_0, a_0^{-1} \[I_n(\\beta)\]^{-1}),
](https://latex.codecogs.com/png.latex?%0A%5Cbeta%20%7C%20a_0%20%5Csim%20N%28%5Chat%7B%5Cbeta%7D_0%2C%20a_0%5E%7B-1%7D%20%5BI_n%28%5Cbeta%29%5D%5E%7B-1%7D%29%2C%0A "
\beta | a_0 \sim N(\hat{\beta}_0, a_0^{-1} [I_n(\beta)]^{-1}),
")

where
![\\hat{\\beta}\_0](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cbeta%7D_0 "\hat{\beta}_0")
is the maximum likelihood estimate (MLE) of
![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") based on
the historical data and
![I_n(\\beta)](https://latex.codecogs.com/png.latex?I_n%28%5Cbeta%29 "I_n(\beta)")
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
#> (Intercept)    0.883   0.003 0.209    0.472    0.745    0.883    1.023    1.293
#> z              0.532   0.003 0.235    0.079    0.373    0.526    0.687    1.010
#> x             -0.801   0.002 0.142   -1.088   -0.893   -0.801   -0.706   -0.522
#> a0             0.760   0.002 0.186    0.324    0.647    0.802    0.911    0.992
#> lp__        -125.736   0.024 1.494 -129.480 -126.493 -125.406 -124.620 -123.861
#>                n_eff  Rhat
#> (Intercept) 5137.779 1.000
#> z           6460.579 0.999
#> x           5160.494 1.000
#> a0          7319.763 1.000
#> lp__        3869.660 1.002
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
#> (Intercept)   1.0   0.760    1.036  0.805        0.845     0.810  0.883  0.850
#> z             0.5   0.677    0.313  0.590        0.611     0.590  0.532  0.596
#> x            -1.0  -0.750   -0.856 -0.751       -0.784    -0.755 -0.801 -0.795
#>                 pp
#> (Intercept)  0.834
#> z            0.611
#> x           -0.788

## posterior std dev.
round( post.sd, 3 )
#>             mle.cur mle.hist   bhm commensurate robustmap  napp   npp    pp
#> (Intercept)   0.274    0.377 0.240        0.247     0.238 0.209 0.237 0.245
#> z             0.308    0.442 0.269        0.279     0.265 0.235 0.266 0.285
#> x             0.181    0.266 0.162        0.170     0.165 0.142 0.161 0.163
```
