
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/hdbayes)](https://CRAN.R-project.org/package=hdbayes)
[![R-CMD-check](https://github.com/ethan-alt/hdbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ethan-alt/hdbayes/actions)
<!-- badges: end -->

# hdbayes: Bayesian Analysis of Generalized Linear Models with Historical Data

The goal of `hdbayes` is to make it easier for users to conduct Bayesian
analysis methods that leverage historical data.

# Installation

Using `hdbayes` package requires installing the R package
[`cmdstanr`](https://mc-stan.org/cmdstanr/) (not available on CRAN) and
the command-line interface to Stan:
[`CmdStan`](https://mc-stan.org/users/interfaces/cmdstan.html). You may
follow the instructions in [Getting started with
CmdStanR](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) to
install both.

You can install the latest version of `hdbayes` package on CRAN using
the following R command:

`install.packages("hdbayes", type = "source")`

Please note that it is important to set `type = "source"`. Otherwise,
the ‘CmdStan’ models in the package may not be compiled during
installation.

Alternatively, you can install the development version of `hdbayes` from
GitHub via

`remotes::install_github("ethan-alt/hdbayes")`

# Example

We now demonstrate how to use `hdbayes` via a simple simulation study.

## Setting up MCMC parameters

We begin by setting up parameters for the MCMC. We detect the number of
cores on our machine as well as set up the warmup and total number of
desired samples.

``` r
## obtain number of cores
ncores        = max(1, parallel::detectCores() - 1)
chains        = 4      ## number of Markov chains to run
iter_warmup   = 1000   ## warmup per chain for MCMC sampling
iter_sampling = 2000   ## number of samples post warmup per chain
```

## Simulating logistic regression data

Then we simulate some logistic regression data. For simplicity, we
generate one current and one historical data set. Note that the package
allows for using multiple historical data sets.

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
data      = data.frame('y' = y, 'z' = X[, 'z'], 'x' = X[, 'x'])
histdata  = data[1:n0, ]
data      = data[-(1:n0), ]
data.list = list(data, histdata)
```

## Frequentist analysis of the data sets

We use the `stats::glm` function to conduct a frequentist GLM of the
data sets

``` r
formula = y ~ z + x
family  = binomial('logit')
fit.mle.cur  = glm(formula, family, data)
fit.mle.hist = glm(formula, family, histdata)
pars         = names(coefficients(fit.mle.cur))
summary(fit.mle.cur)
#> 
#> Call:
#> glm(formula = formula, family = family, data = data)
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

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0i} \| x\_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta_0) \right) \\
  \beta, \beta_0 &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \mu &\sim N_p(\mu_0, \Sigma_0), \text{ where }\Sigma_0 \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu\_{0,j}, \psi\_{0,j}^2),\\ j = 1, \ldots, p
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta%2C%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma%20%3D%20%5Ctext%7Bdiag%7D%28%5Csigma_1%5E2%2C%20%5Cldots%2C%20%5Csigma_p%5E2%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma_0%20%5Ctext%7B%20is%20a%20diagonal%20matrix%7D%20%5C%5C%0A%20%20%5Csigma_j%20%26%5Csim%20N%5E%7B%2B%7D%28%5Cnu_%7B0%2Cj%7D%2C%20%5Cpsi_%7B0%2Cj%7D%5E2%29%2C%5C%3A%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_0) \right) \\
  \beta, \beta_0 &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \mu &\sim N_p(\mu_0, \Sigma_0), \text{ where }\Sigma_0 \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu_{0,j}, \psi_{0,j}^2),\: j = 1, \ldots, p
\end{align*}")

where ![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") is
the vector of regression coefficients of the current data set,
![\beta_0](https://latex.codecogs.com/png.latex?%5Cbeta_0 "\beta_0") is
the vector of regression coefficients for the historical data set,
![\mu](https://latex.codecogs.com/png.latex?%5Cmu "\mu") is the common
prior mean of
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") and
![\beta_0](https://latex.codecogs.com/png.latex?%5Cbeta_0 "\beta_0"),
which is treated as random with a normal hyperprior having mean
![\mu_0](https://latex.codecogs.com/png.latex?%5Cmu_0 "\mu_0"), and a
diagonal covariance matrix
![\Sigma_0](https://latex.codecogs.com/png.latex?%5CSigma_0 "\Sigma_0").
![\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma") is
also treated as random and assumed to be a diagonal matrix. Half-normal
hyperpriors are elicited on the diagonal entries of
![\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma").

The defaults in `hdbayes` are

- ![\mu_0 = 0](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%200 "\mu_0 = 0")
- ![\Sigma_0 = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_0%20%3D%20100%20%5Ctimes%20I_p "\Sigma_0 = 100 \times I_p")
- ![\nu\_{0,j} = 0](https://latex.codecogs.com/png.latex?%5Cnu_%7B0%2Cj%7D%20%3D%200 "\nu_{0,j} = 0")
  for
  ![j = 1, \ldots, p](https://latex.codecogs.com/png.latex?j%20%3D%201%2C%20%5Cldots%2C%20p "j = 1, \ldots, p")
- ![\psi\_{0,j} = 1](https://latex.codecogs.com/png.latex?%5Cpsi_%7B0%2Cj%7D%20%3D%201 "\psi_{0,j} = 1")
  for
  ![j = 1, \ldots, p](https://latex.codecogs.com/png.latex?j%20%3D%201%2C%20%5Cldots%2C%20p "j = 1, \ldots, p")
  where ![p](https://latex.codecogs.com/png.latex?p "p") is the number
  of predictors (including the intercept if applicable).

We fit this model as follows

``` r
## Bayesian Hierarchical Model
fit.bhm = glm.bhm(
  formula = formula, family = family, data.list = data.list,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, adapt_delta = 0.98,
  parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 5.4 seconds.
#> Chain 3 finished in 6.2 seconds.
#> Chain 4 finished in 6.6 seconds.
#> Chain 2 finished in 7.3 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 6.4 seconds.
#> Total execution time: 7.4 seconds.

suppressWarnings(
  fit.bhm[, 2:7] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
)
#> # A tibble: 6 × 10
#>   variable         mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)     0.846  0.844 0.247 0.241  0.437  1.25   1       9511.    6818.
#> 2 z               0.616  0.616 0.286 0.28   0.153  1.10   1       8318.    6887.
#> 3 x              -0.788 -0.786 0.169 0.168 -1.06  -0.513  1       9259.    7419.
#> 4 (Intercept)_h…  0.938  0.931 0.3   0.288  0.461  1.44   1       8168.    7152.
#> 5 z_hist_1        0.453  0.46  0.36  0.351 -0.138  1.03   1       8292.    7244.
#> 6 x_hist_1       -0.842 -0.834 0.217 0.208 -1.21  -0.497  1.00    8681.    7295.
```

### Commensurate prior

The commensurate prior assumes the following hierarchical model

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0i} \| x\_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu_0, \Sigma_0) , \text{ where }\Sigma_0 \text{ is a diagonal matrix} \\
  \beta_j &\sim N_1\left( \beta\_{0j}, \tau_j^{-1} \right), j = 1, \ldots, p
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%20%2C%20%5Ctext%7B%20where%20%7D%5CSigma_0%20%5Ctext%7B%20is%20a%20diagonal%20matrix%7D%20%5C%5C%0A%20%20%5Cbeta_j%20%26%5Csim%20N_1%5Cleft%28%20%5Cbeta_%7B0j%7D%2C%20%5Ctau_j%5E%7B-1%7D%20%5Cright%29%2C%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu_0, \Sigma_0) , \text{ where }\Sigma_0 \text{ is a diagonal matrix} \\
  \beta_j &\sim N_1\left( \beta_{0j}, \tau_j^{-1} \right), j = 1, \ldots, p
\end{align*}")

where the
![\tau_j](https://latex.codecogs.com/png.latex?%5Ctau_j "\tau_j")’s are
elicited by the user. The defaults in `hdbayes` are

- ![\mu_0 = 0](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%200 "\mu_0 = 0")
- ![\Sigma_0 = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_0%20%3D%20100%20%5Ctimes%20I_p "\Sigma_0 = 100 \times I_p")

This method can be fit as follows

``` r
fit.commensurate = glm.commensurate(
  formula = formula, family = family, data.list = data.list,
  tau = rep(5, 3),
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.7 seconds.
#> Chain 2 finished in 0.7 seconds.
#> Chain 3 finished in 0.7 seconds.
#> Chain 4 finished in 0.7 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.7 seconds.
#> Total execution time: 0.9 seconds.

fit.commensurate[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 6 × 10
#>   variable         mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)     0.847  0.841 0.251 0.245  0.441  1.27   1.00    4970.    4803.
#> 2 z               0.607  0.603 0.278 0.272  0.156  1.08   1       5846.    5657.
#> 3 x              -0.786 -0.782 0.169 0.17  -1.07  -0.513  1       5863.    5660.
#> 4 (Intercept)_h…  0.936  0.931 0.296 0.295  0.46   1.43   1       5207.    5165.
#> 5 z_hist          0.479  0.475 0.339 0.34  -0.077  1.04   1.00    5903.    5601.
#> 6 x_hist         -0.854 -0.852 0.221 0.22  -1.22  -0.501  1       5515.    5819.
```

### Robust meta-analytic predictive (MAP) prior

The Robust MAP prior is a generalization of the Bayesian Hierarchical
Model (BHM), and takes the form

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0i} \| x\_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \beta   &\sim w \times N_p(\mu, \Sigma) + (1 - w) N_p(\mu_v, \Sigma_v) \\
  \mu &\sim N_p(\mu_0, \Sigma_0), \text{ where }\Sigma_0 \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu\_{0,j}, \psi\_{0,j}^2),\\ j = 1, \ldots, p
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_0%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_0%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_0%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma%20%3D%20%5Ctext%7Bdiag%7D%28%5Csigma_1%5E2%2C%20%5Cldots%2C%20%5Csigma_p%5E2%29%20%5C%5C%0A%20%20%5Cbeta%20%20%20%26%5Csim%20w%20%5Ctimes%20N_p%28%5Cmu%2C%20%5CSigma%29%20%2B%20%281%20-%20w%29%20N_p%28%5Cmu_v%2C%20%5CSigma_v%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma_0%20%5Ctext%7B%20is%20a%20diagonal%20matrix%7D%20%5C%5C%0A%20%20%5Csigma_j%20%26%5Csim%20N%5E%7B%2B%7D%28%5Cnu_%7B0%2Cj%7D%2C%20%5Cpsi_%7B0%2Cj%7D%5E2%29%2C%5C%3A%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_0 &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_0) \right) \\
  \beta_0 &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \beta   &\sim w \times N_p(\mu, \Sigma) + (1 - w) N_p(\mu_v, \Sigma_v) \\
  \mu &\sim N_p(\mu_0, \Sigma_0), \text{ where }\Sigma_0 \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu_{0,j}, \psi_{0,j}^2),\: j = 1, \ldots, p
\end{align*}")

where
![w \in (0,1)](https://latex.codecogs.com/png.latex?w%20%5Cin%20%280%2C1%29 "w \in (0,1)")
controls for the level of borrowing of the historical data. Note that
when ![w = 1](https://latex.codecogs.com/png.latex?w%20%3D%201 "w = 1"),
the robust MAP prior effectively becomes the BHM. The defaults are the
same as in the BHM except the default value for w is 0.1, and the
default vague prior
![N_p(\mu_v, \Sigma_v)](https://latex.codecogs.com/png.latex?N_p%28%5Cmu_v%2C%20%5CSigma_v%29 "N_p(\mu_v, \Sigma_v)")
is specified as:

- ![\mu_v = 0](https://latex.codecogs.com/png.latex?%5Cmu_v%20%3D%200 "\mu_v = 0")
- ![\Sigma_v = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_v%20%3D%20100%20%5Ctimes%20I_p "\Sigma_v = 100 \times I_p")

The Robust MAP prior is implemented via a three-step procedure in
hdbayes as follows:

1.  sampling from the prior induced by the BHM (equivalent to sampling
    from the MAP prior) using the `glm.rmap.bhm()` function;
2.  approximating the distribution of the samples from step 1. by a
    mixture of multivariate normal distributions using the
    `glm.rmap.bhm.approx()` function;
3.  sampling from the posterior distribution of Robust MAP prior using
    the output from step 2. via the `glm.rmap()` function.

``` r
## step 1.
fit.hist.bhm = glm.rmap.bhm(
   formula, family,
   hist.data.list = list(histdata),
   iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
   chains = chains, adapt_delta = 0.98,
   parallel_chains = ncores,
   refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 2.4 seconds.
#> Chain 4 finished in 2.8 seconds.
#> Chain 3 finished in 3.0 seconds.
#> Chain 2 finished in 3.1 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 2.8 seconds.
#> Total execution time: 3.2 seconds.
## fit.hist.bhm$hist_bhm can be used for assessing MCMC convergence
samples_bhm = fit.hist.bhm$beta_pred

## step 2.
res_approx = glm.rmap.bhm.approx(
  samples.bhm = samples_bhm,
  G = 1:9, verbose = FALSE
)

## step 3.
fit.rmap = glm.rmap(
  formula, family,
  curr.data = data,
  probs = res_approx$probs,
  means = res_approx$means,
  covs = res_approx$covs,
  w = 0.1,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.6 seconds.
#> Chain 2 finished in 0.5 seconds.
#> Chain 3 finished in 0.6 seconds.
#> Chain 4 finished in 0.6 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.6 seconds.
#> Total execution time: 0.8 seconds.
fit.rmap[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 3 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.841  0.841 0.252 0.249  0.421  1.25   1       3908.    3840.
#> 2 z            0.625  0.621 0.283 0.286  0.156  1.09   1.00    4365.    4116.
#> 3 x           -0.79  -0.79  0.167 0.171 -1.07  -0.519  1.00    3938.    3942.
```

### Power prior

The Power Prior takes the form

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0i} \| x\_{0i}, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta) \right) \\
  \pi(\beta \| a_0) &\propto L(\beta \| y_0)^{a_0} \pi_0(\beta)
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cpi%28%5Cbeta%20%7C%20a_0%29%20%26%5Cpropto%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta) \right) \\
  \pi(\beta | a_0) &\propto L(\beta | y_0)^{a_0} \pi_0(\beta)
\end{align*}")

where
![L(\beta \| y_0)](https://latex.codecogs.com/png.latex?L%28%5Cbeta%20%7C%20y_0%29 "L(\beta | y_0)")
is the likelihood of the GLM based on the historical data,
![a_0 \in (0,1)](https://latex.codecogs.com/png.latex?a_0%20%5Cin%20%280%2C1%29 "a_0 \in (0,1)")
is a fixed hyperaparameter controlling the effective sample size
contributed by the data (e.g.,
![a_0 = 1](https://latex.codecogs.com/png.latex?a_0%20%3D%201 "a_0 = 1")
borrows the whole sample size), and
![\pi_0(\beta)](https://latex.codecogs.com/png.latex?%5Cpi_0%28%5Cbeta%29 "\pi_0(\beta)")
is an “initial prior” on
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta").

The default in `hdbayes` is a (noninformative) normal prior on
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta"):

![\beta \sim N_p(0, 100 \times I_p)](https://latex.codecogs.com/png.latex?%5Cbeta%20%5Csim%20N_p%280%2C%20100%20%5Ctimes%20I_p%29 "\beta \sim N_p(0, 100 \times I_p)")

The power prior (with
![a_0 = 0.5](https://latex.codecogs.com/png.latex?a_0%20%3D%200.5 "a_0 = 0.5"))
may be fit as follows:

``` r
fit.pp = glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = 0.5,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.5 seconds.
#> Chain 3 finished in 0.5 seconds.
#> Chain 2 finished in 0.6 seconds.
#> Chain 4 finished in 0.6 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.6 seconds.
#> Total execution time: 0.8 seconds.

fit.pp[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 3 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.834  0.83  0.246 0.247  0.437  1.24   1.00    3560.    3834.
#> 2 z            0.616  0.615 0.276 0.272  0.157  1.06   1       4675.    4336.
#> 3 x           -0.788 -0.784 0.166 0.165 -1.06  -0.518  1       3765.    4041.
```

### Normalized power prior (NPP)

The NPP treats the hyperparameter
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") as random,
allowing the data to decide what is the best value. For non-Gaussian
models, this requires estimating the normalizing constant
![Z(a_0) = \int L(\beta \| y_0)^{a_0} \pi_0(\beta) d\beta](https://latex.codecogs.com/png.latex?Z%28a_0%29%20%3D%20%5Cint%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%20d%5Cbeta "Z(a_0) = \int L(\beta | y_0)^{a_0} \pi_0(\beta) d\beta").

In `hdbayes`, there is one function to estimate the normalizing constant
across a grid of values for
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") and another to
obtain posterior samples of the normalized power prior.

The NPP may be summarized as

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0i} \| x\_{0i}, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta) \right) \\
  \pi(\beta \| a_0) &\propto \frac{1}{Z(a_0)} L(\beta \| y_0)^{a_0} \pi_0(\beta) \\
  \pi(a_0)         &\propto a_0^{\alpha_0 - 1} (1 - a_0)^{\gamma_0 - 1}
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cpi%28%5Cbeta%20%7C%20a_0%29%20%26%5Cpropto%20%5Cfrac%7B1%7D%7BZ%28a_0%29%7D%20L%28%5Cbeta%20%7C%20y_0%29%5E%7Ba_0%7D%20%5Cpi_0%28%5Cbeta%29%20%5C%5C%0A%20%20%5Cpi%28a_0%29%20%20%20%20%20%20%20%20%20%26%5Cpropto%20a_0%5E%7B%5Calpha_0%20-%201%7D%20%281%20-%20a_0%29%5E%7B%5Cgamma_0%20-%201%7D%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta) \right) \\
  \pi(\beta | a_0) &\propto \frac{1}{Z(a_0)} L(\beta | y_0)^{a_0} \pi_0(\beta) \\
  \pi(a_0)         &\propto a_0^{\alpha_0 - 1} (1 - a_0)^{\gamma_0 - 1}
\end{align*}")

The defaults in `hdbayes` are

- ![\pi_0(\beta) \propto N(\beta \| 0, 100 \times I_p)](https://latex.codecogs.com/png.latex?%5Cpi_0%28%5Cbeta%29%20%5Cpropto%20N%28%5Cbeta%20%7C%200%2C%20100%20%5Ctimes%20I_p%29 "\pi_0(\beta) \propto N(\beta | 0, 100 \times I_p)")
- ![\alpha_0 = 1](https://latex.codecogs.com/png.latex?%5Calpha_0%20%3D%201 "\alpha_0 = 1")
- ![\gamma_0 = 1](https://latex.codecogs.com/png.latex?%5Cgamma_0%20%3D%201 "\gamma_0 = 1")

when
![\alpha_0 = 1](https://latex.codecogs.com/png.latex?%5Calpha_0%20%3D%201 "\alpha_0 = 1")
and
![\gamma_0 = 1](https://latex.codecogs.com/png.latex?%5Cgamma_0%20%3D%201 "\gamma_0 = 1"),
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
logncfun = function(a0, ...){
  hdbayes::glm.npp.lognc(
    formula = formula, family = family, histdata = histdata, a0 = a0, ...
  )
}

cl = makeCluster(ncores)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'histdata'))

## call created function
a0.lognc = parLapply(
  cl = cl, X = a0, fun = logncfun, iter_warmup = iter_warmup, iter_sampling = 5000, 
  chains = chains, refresh = 0
)
stopCluster(cl)

a0.lognc = data.frame( do.call(rbind, a0.lognc) )
head(
  a0.lognc %>% 
    mutate(across(where(is.numeric), round, 3))
)
#>      a0   lognc min_ess_bulk max_Rhat
#> 1 0.000   0.000     9607.838    1.000
#> 2 0.014  -4.229     6501.686    1.001
#> 3 0.027  -6.290     6768.113    1.001
#> 4 0.041  -7.859     6756.573    1.000
#> 5 0.054  -9.212     6579.915    1.000
#> 6 0.068 -10.444     6727.823    1.000
```

The provided function `glm.npp.lognc` estimates the logarithm of the
normalizing constant,
![\log Z(a_0)](https://latex.codecogs.com/png.latex?%5Clog%20Z%28a_0%29 "\log Z(a_0)"),
for one specific value of
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0") and one
historical data set. We created the function `logncfun` so that the
first argument would be
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0"), allowing us to
use the `parLapply` function in the `parallel` package.

The `hdbayes` function `glm.npp.lognc` outputs
![a_0](https://latex.codecogs.com/png.latex?a_0 "a_0"),
![Z(a_0)](https://latex.codecogs.com/png.latex?Z%28a_0%29 "Z(a_0)"), and
the minimum bulk effective sample size and maximum R-hat value of the
MCMC sampling of the power prior. It is a good idea to check that the
minimum bulk effective sample size is at least 1,000 and the maximum
R-hat value is less than 1.10

``` r
min(a0.lognc$min_ess_bulk) ## lowest bulk effective sample size
#> [1] 6418.662
max(a0.lognc$max_Rhat)  ## highest R-hat value
#> [1] 1.001215
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
  formula = formula, family = family, data.list = data.list,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.8 seconds.
#> Chain 2 finished in 0.7 seconds.
#> Chain 3 finished in 0.7 seconds.
#> Chain 4 finished in 0.8 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.8 seconds.
#> Total execution time: 1.0 seconds.
fit.npp[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 4 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.845  0.844 0.244 0.243  0.454  1.25   1.00    4453.    4320.
#> 2 z            0.598  0.593 0.273 0.274  0.155  1.05   1.00    5437.    5044.
#> 3 x           -0.793 -0.789 0.164 0.167 -1.07  -0.532  1.00    4785.    4421.
#> 4 a0_hist_1    0.676  0.713 0.222 0.249  0.266  0.972  1.00    5552.    3683.
```

### Normalized asymptotic power prior (NAPP)

NAPP uses a large sample theory argument to formulate a normal
approximation to the power prior, i.e., the prior is given by

![\beta \| a_0 \sim N(\hat{\beta}\_0, a_0^{-1} \[I_n(\beta)\]^{-1}),](https://latex.codecogs.com/png.latex?%5Cbeta%20%7C%20a_0%20%5Csim%20N%28%5Chat%7B%5Cbeta%7D_0%2C%20a_0%5E%7B-1%7D%20%5BI_n%28%5Cbeta%29%5D%5E%7B-1%7D%29%2C "\beta | a_0 \sim N(\hat{\beta}_0, a_0^{-1} [I_n(\beta)]^{-1}),")

where
![\hat{\beta}\_0](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cbeta%7D_0 "\hat{\beta}_0")
is the maximum likelihood estimate (MLE) of
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") based on
the historical data and
![I_n(\beta)](https://latex.codecogs.com/png.latex?I_n%28%5Cbeta%29 "I_n(\beta)")
is the associated information matrix (negative Hessian).

In this case, the normalizing constant is known, so we do not need to
estimate the normalizing constant before sampling.

The following code analyzes the data set using the NAPP:

``` r
fit.napp = glm.napp(
  formula = formula, family = family, data.list = data.list,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.5 seconds.
#> Chain 2 finished in 0.6 seconds.
#> Chain 3 finished in 0.6 seconds.
#> Chain 4 finished in 0.6 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.6 seconds.
#> Total execution time: 0.8 seconds.
fit.napp[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 4 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.842  0.837 0.241 0.24   0.453  1.24   1.00    4511.    4506.
#> 2 z            0.592  0.587 0.275 0.27   0.14   1.06   1.00    4990.    4505.
#> 3 x           -0.787 -0.785 0.159 0.159 -1.05  -0.528  1       4638.    4442.
#> 4 a0_hist_1    0.671  0.704 0.225 0.253  0.258  0.97   1.00    5304.    3853.
```

## Comparison of methods

We now can compare the point estimate (MLE / posterior mean) and
uncertainty (SE / posterior standard deviation) of the methods.

``` r
##
## COMPARE METHODS
##
fit.list = list('bhm' = fit.bhm, 'commensurate' = fit.commensurate,
                'robustmap' = fit.rmap, 'napp' = fit.napp,
                'npp' = fit.npp, 'pp' = fit.pp)

post.mean = suppressWarnings(
  sapply(
    fit.list, function(x){
      d = x %>% 
        select(names(coef(fit.mle.cur))) %>% 
        summarise_draws()
      d$mean
    } 
  )
)
post.sd = suppressWarnings(
  sapply(
    fit.list, function(x){
      d = x %>% 
        select(names(coef(fit.mle.cur))) %>% 
        summarise_draws()
      d$sd
    } 
  )
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
#> (Intercept)   1.0   0.760    1.036  0.846        0.847     0.841  0.842  0.845
#> z             0.5   0.677    0.313  0.616        0.607     0.625  0.592  0.598
#> x            -1.0  -0.750   -0.856 -0.788       -0.786    -0.790 -0.787 -0.793
#>                 pp
#> (Intercept)  0.834
#> z            0.616
#> x           -0.788

## posterior std dev.
round( post.sd, 3 )
#>             mle.cur mle.hist   bhm commensurate robustmap  napp   npp    pp
#> (Intercept)   0.274    0.377 0.247        0.251     0.252 0.241 0.244 0.246
#> z             0.308    0.442 0.286        0.278     0.283 0.275 0.273 0.276
#> x             0.181    0.266 0.169        0.169     0.167 0.159 0.164 0.166
```
