
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

You can install the released version of `hdbayes` from CRAN with:

``` r
install.packages("hdbayes", type = "source")
```

Please note that it is important to set `type = "source"`. Otherwise,
the ‘CmdStan’ models in the package may not be compiled during
installation.

Alternatively, you can install the development version from GitHub with

``` r
remotes::install_github("ethan-alt/hdbayes")
```

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
allows for using multiple historical data sets. Let the current data set
be denoted by
![D = \\ (y_i, x_i), i = 1, \ldots, n \\](https://latex.codecogs.com/png.latex?D%20%3D%20%5C%7B%20%28y_i%2C%20x_i%29%2C%20i%20%3D%201%2C%20%5Cldots%2C%20n%20%5C%7D "D = \{ (y_i, x_i), i = 1, \ldots, n \}"),
where ![n](https://latex.codecogs.com/png.latex?n "n") denotes the
sample size of the current data. Suppose we have
![H](https://latex.codecogs.com/png.latex?H "H") historical data sets,
and let the
![h^{th}](https://latex.codecogs.com/png.latex?h%5E%7Bth%7D "h^{th}")
historical data set,
![h \in \\1, \ldots, H\\](https://latex.codecogs.com/png.latex?h%20%5Cin%20%5C%7B1%2C%20%5Cldots%2C%20H%5C%7D "h \in \{1, \ldots, H\}"),
be denoted by
![D\_{0h} = \\ (y\_{0hi}, x\_{0hi}), i = 1, \ldots, n\_{0h} \\](https://latex.codecogs.com/png.latex?D_%7B0h%7D%20%3D%20%5C%7B%20%28y_%7B0hi%7D%2C%20x_%7B0hi%7D%29%2C%20i%20%3D%201%2C%20%5Cldots%2C%20n_%7B0h%7D%20%5C%7D "D_{0h} = \{ (y_{0hi}, x_{0hi}), i = 1, \ldots, n_{0h} \}"),
where
![n\_{0h}](https://latex.codecogs.com/png.latex?n_%7B0h%7D "n_{0h}")
denotes the sample size of the
![h^{th}](https://latex.codecogs.com/png.latex?h%5E%7Bth%7D "h^{th}")
historical data. Let
![D_0 = \\D\_{01}, \ldots, D\_{0H}\\](https://latex.codecogs.com/png.latex?D_0%20%3D%20%5C%7BD_%7B01%7D%2C%20%5Cldots%2C%20D_%7B0H%7D%5C%7D "D_0 = \{D_{01}, \ldots, D_{0H}\}")
denote all the historical data.

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

<!-- We will compare these with the Bayesian analysis results later -->

## Bayesian analysis methods

We now utilize the functions in this package.

### Bayesian hierarchical model

The Bayesian hierarchical model (BHM) may be expressed hierarhically as

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right), \\\\i = 1, \ldots, n \\
  y\_{0hi} \| x\_{0hi}, \beta\_{0h} &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0hi}'\beta\_{0h}) \right), \\\\h = 1, \ldots, H, \\\\i = 1, \ldots, n\_{0h} \\
  \beta, \beta\_{01}, \ldots, \beta\_{0H} &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \mu &\sim N_p(\mu_0, \Sigma\_{0}), \text{ where } \Sigma\_{0} \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu\_{0,j}, \psi\_{0,j}^2),\\ j = 1, \ldots, p
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n%20%5C%5C%0A%20%20y_%7B0hi%7D%20%7C%20x_%7B0hi%7D%2C%20%5Cbeta_%7B0h%7D%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0hi%7D%27%5Cbeta_%7B0h%7D%29%20%5Cright%29%2C%20%5C%20%5C%20h%20%3D%201%2C%20%5Cldots%2C%20H%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n_%7B0h%7D%20%5C%5C%0A%20%20%5Cbeta%2C%20%5Cbeta_%7B01%7D%2C%20%5Cldots%2C%20%5Cbeta_%7B0H%7D%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma%20%3D%20%5Ctext%7Bdiag%7D%28%5Csigma_1%5E2%2C%20%5Cldots%2C%20%5Csigma_p%5E2%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_%7B0%7D%29%2C%20%5Ctext%7B%20where%20%7D%20%5CSigma_%7B0%7D%20%5Ctext%7B%20is%20a%20diagonal%20matrix%7D%20%5C%5C%0A%20%20%5Csigma_j%20%26%5Csim%20N%5E%7B%2B%7D%28%5Cnu_%7B0%2Cj%7D%2C%20%5Cpsi_%7B0%2Cj%7D%5E2%29%2C%5C%3A%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right), \ \ i = 1, \ldots, n \\
  y_{0hi} | x_{0hi}, \beta_{0h} &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0hi}'\beta_{0h}) \right), \ \ h = 1, \ldots, H, \ \ i = 1, \ldots, n_{0h} \\
  \beta, \beta_{01}, \ldots, \beta_{0H} &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \mu &\sim N_p(\mu_0, \Sigma_{0}), \text{ where } \Sigma_{0} \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu_{0,j}, \psi_{0,j}^2),\: j = 1, \ldots, p
\end{align*}")

where ![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") is
the vector of regression coefficients of the current data set,
![\beta\_{0h}](https://latex.codecogs.com/png.latex?%5Cbeta_%7B0h%7D "\beta_{0h}")
is the vector of regression coefficients for historical data set
![D\_{0h}](https://latex.codecogs.com/png.latex?D_%7B0h%7D "D_{0h}"),
![\mu](https://latex.codecogs.com/png.latex?%5Cmu "\mu") is the common
prior mean of
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") and
![\beta\_{0h}](https://latex.codecogs.com/png.latex?%5Cbeta_%7B0h%7D "\beta_{0h}"),
which is treated as random with a normal hyperprior having mean
![\mu_0](https://latex.codecogs.com/png.latex?%5Cmu_0 "\mu_0"), and a
diagonal covariance matrix
![\Sigma\_{0}](https://latex.codecogs.com/png.latex?%5CSigma_%7B0%7D "\Sigma_{0}").
![\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma") is
also treated as random and assumed to be a diagonal matrix. Half-normal
hyperpriors are elicited on the diagonal entries of
![\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma").

The defaults in `hdbayes` are

- ![\mu_0 = \textbf{0}\_p](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%20%5Ctextbf%7B0%7D_p "\mu_0 = \textbf{0}_p"),
  where
  ![\textbf{0}\_p](https://latex.codecogs.com/png.latex?%5Ctextbf%7B0%7D_p "\textbf{0}_p")
  denotes a ![p](https://latex.codecogs.com/png.latex?p "p")-dimensional
  vector of 0s
- ![\Sigma\_{0} = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_%7B0%7D%20%3D%20100%20%5Ctimes%20I_p "\Sigma_{0} = 100 \times I_p")
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
#> Chain 2 finished in 8.0 seconds.
#> Chain 1 finished in 8.8 seconds.
#> Chain 3 finished in 9.3 seconds.
#> Chain 4 finished in 9.7 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 9.0 seconds.
#> Total execution time: 9.9 seconds.

suppressWarnings(
  fit.bhm[, 2:7] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
)
#> # A tibble: 6 × 10
#>   variable         mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)     0.842  0.84  0.246 0.246  0.438  1.25      1    9912.    7113.
#> 2 z               0.619  0.615 0.286 0.283  0.151  1.1       1    8977.    6988.
#> 3 x              -0.788 -0.786 0.166 0.159 -1.07  -0.519     1    9548.    6723.
#> 4 (Intercept)_h…  0.937  0.926 0.302 0.293  0.458  1.44      1    8523.    6766.
#> 5 z_hist_1        0.46   0.461 0.368 0.362 -0.159  1.05      1    8395.    6735.
#> 6 x_hist_1       -0.842 -0.836 0.215 0.207 -1.21  -0.497     1    8746.    7471.
```

### Commensurate prior

The commensurate prior assumes the following hierarchical model

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0hi} \| x\_{0hi}, \beta\_{0} &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0hi}'\beta\_{0}) \right) \\
  \beta\_{0} &\sim N_p(\mu_0, \Sigma\_{0}) , \text{ where }\Sigma\_{0} \text{ is a diagonal matrix} \\
  \beta_j &\sim N_1\left( \beta\_{0j}, \tau_j^{-1} \right), j = 1, \ldots, p \\
  \tau_j &\overset{\text{i.i.d.}}{\sim} p\_{\text{spike}} \cdot N^{+}(\mu\_{\text{spike}}, \sigma\_{\text{spike}}^2) + (1 - p\_{\text{spike}}) \cdot N^{+}(\mu\_{\text{slab}}, \sigma\_{\text{slab}}^2),\\ j = 1, \ldots, p 
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0hi%7D%20%7C%20x_%7B0hi%7D%2C%20%5Cbeta_%7B0%7D%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0hi%7D%27%5Cbeta_%7B0%7D%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_%7B0%7D%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_%7B0%7D%29%20%2C%20%5Ctext%7B%20where%20%7D%5CSigma_%7B0%7D%20%5Ctext%7B%20is%20a%20diagonal%20matrix%7D%20%5C%5C%0A%20%20%5Cbeta_j%20%26%5Csim%20N_1%5Cleft%28%20%5Cbeta_%7B0j%7D%2C%20%5Ctau_j%5E%7B-1%7D%20%5Cright%29%2C%20j%20%3D%201%2C%20%5Cldots%2C%20p%20%5C%5C%0A%20%20%5Ctau_j%20%26%5Coverset%7B%5Ctext%7Bi.i.d.%7D%7D%7B%5Csim%7D%20p_%7B%5Ctext%7Bspike%7D%7D%20%5Ccdot%20N%5E%7B%2B%7D%28%5Cmu_%7B%5Ctext%7Bspike%7D%7D%2C%20%5Csigma_%7B%5Ctext%7Bspike%7D%7D%5E2%29%20%2B%20%281%20-%20p_%7B%5Ctext%7Bspike%7D%7D%29%20%5Ccdot%20N%5E%7B%2B%7D%28%5Cmu_%7B%5Ctext%7Bslab%7D%7D%2C%20%5Csigma_%7B%5Ctext%7Bslab%7D%7D%5E2%29%2C%5C%3A%20j%20%3D%201%2C%20%5Cldots%2C%20p%20%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0hi} | x_{0hi}, \beta_{0} &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0hi}'\beta_{0}) \right) \\
  \beta_{0} &\sim N_p(\mu_0, \Sigma_{0}) , \text{ where }\Sigma_{0} \text{ is a diagonal matrix} \\
  \beta_j &\sim N_1\left( \beta_{0j}, \tau_j^{-1} \right), j = 1, \ldots, p \\
  \tau_j &\overset{\text{i.i.d.}}{\sim} p_{\text{spike}} \cdot N^{+}(\mu_{\text{spike}}, \sigma_{\text{spike}}^2) + (1 - p_{\text{spike}}) \cdot N^{+}(\mu_{\text{slab}}, \sigma_{\text{slab}}^2),\: j = 1, \ldots, p 
\end{align*}")

where
![\beta = (\beta_1, \ldots, \beta_p)'](https://latex.codecogs.com/png.latex?%5Cbeta%20%3D%20%28%5Cbeta_1%2C%20%5Cldots%2C%20%5Cbeta_p%29%27 "\beta = (\beta_1, \ldots, \beta_p)'"),
![\beta\_{0} = (\beta\_{01}, \ldots, \beta\_{0p})'](https://latex.codecogs.com/png.latex?%5Cbeta_%7B0%7D%20%3D%20%28%5Cbeta_%7B01%7D%2C%20%5Cldots%2C%20%5Cbeta_%7B0p%7D%29%27 "\beta_{0} = (\beta_{01}, \ldots, \beta_{0p})'").
The commensurability parameters (i.e.,
![\tau_j](https://latex.codecogs.com/png.latex?%5Ctau_j "\tau_j")’s) are
treated as random with a spike-and-slab prior, which is specified as a
mixture of two half-normal priors. The defaults in `hdbayes` are

- ![\mu_0 = \textbf{0}\_p](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%20%5Ctextbf%7B0%7D_p "\mu_0 = \textbf{0}_p")
- ![\Sigma\_{0} = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_%7B0%7D%20%3D%20100%20%5Ctimes%20I_p "\Sigma_{0} = 100 \times I_p")
- ![p\_{\text{spike}}  = 0.1](https://latex.codecogs.com/png.latex?p_%7B%5Ctext%7Bspike%7D%7D%20%20%3D%200.1 "p_{\text{spike}}  = 0.1")
- ![\mu\_{\text{spike}} = 200](https://latex.codecogs.com/png.latex?%5Cmu_%7B%5Ctext%7Bspike%7D%7D%20%3D%20200 "\mu_{\text{spike}} = 200")
- ![\sigma\_{\text{spike}} = 0.1](https://latex.codecogs.com/png.latex?%5Csigma_%7B%5Ctext%7Bspike%7D%7D%20%3D%200.1 "\sigma_{\text{spike}} = 0.1")
- ![\mu\_{\text{slab}} = 0](https://latex.codecogs.com/png.latex?%5Cmu_%7B%5Ctext%7Bslab%7D%7D%20%3D%200 "\mu_{\text{slab}} = 0")
- ![\sigma\_{\text{slab}} = 5](https://latex.codecogs.com/png.latex?%5Csigma_%7B%5Ctext%7Bslab%7D%7D%20%3D%205 "\sigma_{\text{slab}} = 5")

This method can be fit as follows

``` r
fit.commensurate = glm.commensurate(
  formula = formula, family = family, data.list = data.list,
  p.spike = 0.1,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 1.3 seconds.
#> Chain 2 finished in 1.2 seconds.
#> Chain 3 finished in 1.3 seconds.
#> Chain 4 finished in 1.2 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 1.3 seconds.
#> Total execution time: 1.4 seconds.

suppressWarnings(
  fit.commensurate[, 2:10] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
)
#> # A tibble: 9 × 10
#>   variable         mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)     0.838  0.837 0.251 0.251  0.435  1.26   1       5711.    5176.
#> 2 z               0.625  0.624 0.28  0.285  0.168  1.08   1       6389.    5850.
#> 3 x              -0.786 -0.784 0.17  0.167 -1.07  -0.509  1.00    6249.    5381.
#> 4 (Intercept)_h…  0.942  0.933 0.306 0.301  0.443  1.46   1       5755.    5623.
#> 5 z_hist          0.468  0.471 0.357 0.351 -0.124  1.04   1       6387.    5579.
#> 6 x_hist         -0.854 -0.849 0.229 0.226 -1.24  -0.486  1       6247.    5721.
#> 7 comm_prec[1]    4.73   4.16  3.06  3.03   0.724 10.5    1       5532.    3831.
#> 8 comm_prec[2]    4.62   4.11  2.98  2.94   0.736 10.2    1       5644.    3474.
#> 9 comm_prec[3]    4.93   4.42  3.08  3.08   0.834 10.7    1       5757.    3410.
```

### Robust meta-analytic predictive (MAP) prior

The Robust MAP prior is a generalization of the Bayesian Hierarchical
Model (BHM), and takes the form

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y\_{0i} \| x\_{0i}, \beta\_{0h} &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta\_{0h}) \right) \\
  \beta\_{0h} &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \beta   &\sim w \times N_p(\mu, \Sigma) + (1 - w) N_p(\mu_v, \Sigma_v) \\
  \mu &\sim N_p(\mu_0, \Sigma\_{0}), \text{ where }\Sigma\_{0} \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu\_{0,j}, \psi\_{0,j}^2),\\ j = 1, \ldots, p
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%20%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta_%7B0h%7D%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_%7B0h%7D%29%20%5Cright%29%20%5C%5C%0A%20%20%5Cbeta_%7B0h%7D%20%26%5Csim%20N_p%28%5Cmu%2C%20%5CSigma%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma%20%3D%20%5Ctext%7Bdiag%7D%28%5Csigma_1%5E2%2C%20%5Cldots%2C%20%5Csigma_p%5E2%29%20%5C%5C%0A%20%20%5Cbeta%20%20%20%26%5Csim%20w%20%5Ctimes%20N_p%28%5Cmu%2C%20%5CSigma%29%20%2B%20%281%20-%20w%29%20N_p%28%5Cmu_v%2C%20%5CSigma_v%29%20%5C%5C%0A%20%20%5Cmu%20%26%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_%7B0%7D%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma_%7B0%7D%20%5Ctext%7B%20is%20a%20diagonal%20matrix%7D%20%5C%5C%0A%20%20%5Csigma_j%20%26%5Csim%20N%5E%7B%2B%7D%28%5Cnu_%7B0%2Cj%7D%2C%20%5Cpsi_%7B0%2Cj%7D%5E2%29%2C%5C%3A%20j%20%3D%201%2C%20%5Cldots%2C%20p%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right) \\
  y_{0i} | x_{0i}, \beta_{0h} &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_{0h}) \right) \\
  \beta_{0h} &\sim N_p(\mu, \Sigma), \text{ where }\Sigma = \text{diag}(\sigma_1^2, \ldots, \sigma_p^2) \\
  \beta   &\sim w \times N_p(\mu, \Sigma) + (1 - w) N_p(\mu_v, \Sigma_v) \\
  \mu &\sim N_p(\mu_0, \Sigma_{0}), \text{ where }\Sigma_{0} \text{ is a diagonal matrix} \\
  \sigma_j &\sim N^{+}(\nu_{0,j}, \psi_{0,j}^2),\: j = 1, \ldots, p
\end{align*}")

where
![w \in (0,1)](https://latex.codecogs.com/png.latex?w%20%5Cin%20%280%2C1%29 "w \in (0,1)")
controls for the level of borrowing of the historical data. Note that
when ![w = 1](https://latex.codecogs.com/png.latex?w%20%3D%201 "w = 1"),
the robust MAP prior effectively becomes the BHM. The defaults are the
same as in the BHM except the default value for w is 0.1, and the
default vague/non-informative prior
![N_p(\mu_v, \Sigma_v)](https://latex.codecogs.com/png.latex?N_p%28%5Cmu_v%2C%20%5CSigma_v%29 "N_p(\mu_v, \Sigma_v)")
is specified as:

- ![\mu_v = \textbf{0}\_p](https://latex.codecogs.com/png.latex?%5Cmu_v%20%3D%20%5Ctextbf%7B0%7D_p "\mu_v = \textbf{0}_p")
- ![\Sigma_v = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_v%20%3D%20100%20%5Ctimes%20I_p "\Sigma_v = 100 \times I_p")

The posterior samples under the Robust MAP prior is obtained by using
the marginal likelihoods of the vague and meta-analytic predictive (MAP)
priors. Specifically, note that the posterior density under the Robust
MAP prior can be written as

![\begin{align\*}
    p\_{\text{RMAP}}(\beta \| D, D_0, w) &= \frac{L(\beta \| D) \left\[ w \\\pi_I(\beta \| D_0) + (1 - w) \\\pi_V(\beta) \right\]}{\int L(\beta^\* \| D) \left\[ w \\\pi_I( \beta^\* \| D_0) + (1 - w) \\\pi_V( \beta^\* ) \right\] d\beta^\*}, \\
    &= \tilde{w} \\p_I(\beta \| D, D_0) + (1 - \tilde{w}) \\p_V(\beta \| D),
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20%20%20p_%7B%5Ctext%7BRMAP%7D%7D%28%5Cbeta%20%7C%20D%2C%20D_0%2C%20w%29%20%26%3D%20%5Cfrac%7BL%28%5Cbeta%20%7C%20D%29%20%5Cleft%5B%20w%20%5C%3A%5Cpi_I%28%5Cbeta%20%7C%20D_0%29%20%2B%20%281%20-%20w%29%20%5C%3A%5Cpi_V%28%5Cbeta%29%20%5Cright%5D%7D%7B%5Cint%20L%28%5Cbeta%5E%2A%20%7C%20D%29%20%5Cleft%5B%20w%20%5C%3A%5Cpi_I%28%20%5Cbeta%5E%2A%20%7C%20D_0%29%20%2B%20%281%20-%20w%29%20%5C%3A%5Cpi_V%28%20%5Cbeta%5E%2A%20%29%20%5Cright%5D%20d%5Cbeta%5E%2A%7D%2C%20%5C%5C%0A%20%20%20%20%26%3D%20%5Ctilde%7Bw%7D%20%5C%3Ap_I%28%5Cbeta%20%7C%20D%2C%20D_0%29%20%2B%20%281%20-%20%5Ctilde%7Bw%7D%29%20%5C%3Ap_V%28%5Cbeta%20%7C%20D%29%2C%0A%5Cend%7Balign%2A%7D "\begin{align*}
    p_{\text{RMAP}}(\beta | D, D_0, w) &= \frac{L(\beta | D) \left[ w \:\pi_I(\beta | D_0) + (1 - w) \:\pi_V(\beta) \right]}{\int L(\beta^* | D) \left[ w \:\pi_I( \beta^* | D_0) + (1 - w) \:\pi_V( \beta^* ) \right] d\beta^*}, \\
    &= \tilde{w} \:p_I(\beta | D, D_0) + (1 - \tilde{w}) \:p_V(\beta | D),
\end{align*}")

where
![p_I(\beta \| D, D_0) = L(\beta \| D) \pi_I(\beta \| D_0) / Z_I(D, D_0)](https://latex.codecogs.com/png.latex?p_I%28%5Cbeta%20%7C%20D%2C%20D_0%29%20%3D%20L%28%5Cbeta%20%7C%20D%29%20%5Cpi_I%28%5Cbeta%20%7C%20D_0%29%20%2F%20Z_I%28D%2C%20D_0%29 "p_I(\beta | D, D_0) = L(\beta | D) \pi_I(\beta | D_0) / Z_I(D, D_0)")
is the posterior density under the prior induced by the BHM (referred to
as the meta-analytic predictive (MAP) prior),
![p_V(\beta \| D) = L(\beta \| D) \pi_V(\beta) / Z_V(D)](https://latex.codecogs.com/png.latex?p_V%28%5Cbeta%20%7C%20D%29%20%3D%20L%28%5Cbeta%20%7C%20D%29%20%5Cpi_V%28%5Cbeta%29%20%2F%20Z_V%28D%29 "p_V(\beta | D) = L(\beta | D) \pi_V(\beta) / Z_V(D)")
is the posterior density under the vague prior, and

![\tilde{w} = \frac{w \\Z_I(D, D_0)}{w \\Z_I(D, D_0) + (1 - w) \\Z_V(D)}](https://latex.codecogs.com/png.latex?%5Ctilde%7Bw%7D%20%3D%20%5Cfrac%7Bw%20%5C%3AZ_I%28D%2C%20D_0%29%7D%7Bw%20%5C%3AZ_I%28D%2C%20D_0%29%20%2B%20%281%20-%20w%29%20%5C%3AZ_V%28D%29%7D "\tilde{w} = \frac{w \:Z_I(D, D_0)}{w \:Z_I(D, D_0) + (1 - w) \:Z_V(D)}")

is the updated mixture weight. The normalizing constants
![Z_I(D, D_0)](https://latex.codecogs.com/png.latex?Z_I%28D%2C%20D_0%29 "Z_I(D, D_0)")
and ![Z_V(D)](https://latex.codecogs.com/png.latex?Z_V%28D%29 "Z_V(D)")
are estimated via the R package
[`bridgesampling`](https://cran.r-project.org/web/packages/bridgesampling/)
in `hdbayes`.

``` r
res.rmap = glm.rmap(
  formula = formula, family = family, data.list = data.list,
  w = 0.1,
  bridge.args = list(silent = T),
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores, adapt_delta = 0.98,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 2 finished in 7.5 seconds.
#> Chain 3 finished in 9.0 seconds.
#> Chain 4 finished in 9.0 seconds.
#> Chain 1 finished in 11.4 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 9.2 seconds.
#> Total execution time: 11.5 seconds.
#> Warning: 3 of 8000 (0.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 3 finished in 3.6 seconds.
#> Chain 2 finished in 4.2 seconds.
#> Chain 4 finished in 4.2 seconds.
#> Chain 1 finished in 4.5 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 4.1 seconds.
#> Total execution time: 4.6 seconds.
#> 
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.8 seconds.
#> Chain 2 finished in 0.8 seconds.
#> Chain 3 finished in 0.8 seconds.
#> Chain 4 finished in 0.9 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.9 seconds.
#> Total execution time: 1.1 seconds.
fit.rmap = res.rmap$post.samples
fit.rmap[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 3 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.849  0.85  0.253 0.251  0.438  1.26      1    4111.    3289.
#> 2 z            0.617  0.616 0.28  0.28   0.161  1.08      1    3921.    3705.
#> 3 x           -0.791 -0.787 0.169 0.17  -1.08  -0.521     1    4315.    3679.
```

### Power prior

With ![H](https://latex.codecogs.com/png.latex?H "H") historical data
sets, the power prior takes the form

![\begin{align\*}
  &y_i \| x_i, \beta \sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \\\\i = 1, \ldots, n\\
  &\pi\_{\text{PP}}(\beta \| D_0, a\_{0h}) \propto 
    \pi_0(\beta) \prod\_{h=1}^H L(\beta \| D\_{0h})^{a\_{0h}}
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20%26y_i%20%7C%20x_i%2C%20%5Cbeta%20%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%0A%20%20%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n%5C%5C%0A%20%20%26%5Cpi_%7B%5Ctext%7BPP%7D%7D%28%5Cbeta%20%7C%20D_0%2C%20a_%7B0h%7D%29%20%5Cpropto%20%0A%20%20%20%20%5Cpi_0%28%5Cbeta%29%20%5Cprod_%7Bh%3D1%7D%5EH%20L%28%5Cbeta%20%7C%20D_%7B0h%7D%29%5E%7Ba_%7B0h%7D%7D%0A%5Cend%7Balign%2A%7D "\begin{align*}
  &y_i | x_i, \beta \sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \ \ i = 1, \ldots, n\\
  &\pi_{\text{PP}}(\beta | D_0, a_{0h}) \propto 
    \pi_0(\beta) \prod_{h=1}^H L(\beta | D_{0h})^{a_{0h}}
\end{align*}")

where
![L(\beta \| D)](https://latex.codecogs.com/png.latex?L%28%5Cbeta%20%7C%20D%29 "L(\beta | D)")
is the likelihood of the GLM based on data set
![D](https://latex.codecogs.com/png.latex?D "D"),
![a\_{0h} \in (0,1)](https://latex.codecogs.com/png.latex?a_%7B0h%7D%20%5Cin%20%280%2C1%29 "a_{0h} \in (0,1)")
is a fixed hyperparameter controlling the effective sample size
contributed by historical data set
![D\_{0h}](https://latex.codecogs.com/png.latex?D_%7B0h%7D "D_{0h}") ,
and
![\pi_0(\beta)](https://latex.codecogs.com/png.latex?%5Cpi_0%28%5Cbeta%29 "\pi_0(\beta)")
is an “initial prior” for
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta").

The default in `hdbayes` is a (non-informative) normal initial prior on
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta"):

![\beta \sim N_p(\textbf{0}\_p, 100 \times I_p)](https://latex.codecogs.com/png.latex?%5Cbeta%20%5Csim%20N_p%28%5Ctextbf%7B0%7D_p%2C%20100%20%5Ctimes%20I_p%29 "\beta \sim N_p(\textbf{0}_p, 100 \times I_p)")

With ![H = 1](https://latex.codecogs.com/png.latex?H%20%3D%201 "H = 1"),
the power prior (with
![a\_{01} = 0.5](https://latex.codecogs.com/png.latex?a_%7B01%7D%20%3D%200.5 "a_{01} = 0.5"))
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
#> Chain 1 finished in 0.7 seconds.
#> Chain 2 finished in 0.7 seconds.
#> Chain 3 finished in 0.8 seconds.
#> Chain 4 finished in 0.8 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.8 seconds.
#> Total execution time: 1.0 seconds.

fit.pp[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 3 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.834  0.828 0.251 0.249  0.43   1.26   1       3788.    3746.
#> 2 z            0.613  0.616 0.28  0.28   0.155  1.07   1.00    4106.    4060.
#> 3 x           -0.788 -0.785 0.165 0.165 -1.06  -0.526  1.00    4306.    4689.
```

### Normalized power prior (NPP)

The NPP treats the hyperparameter
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}") as
random, allowing the data to decide what is the best value. For
non-Gaussian models, this requires estimating the normalizing constant
![Z(a\_{0h} \| D\_{0h}) = \int L(\beta \| D\_{0h})^{a\_{0h}} \pi_0(\beta)^{1/H} d\beta](https://latex.codecogs.com/png.latex?Z%28a_%7B0h%7D%20%7C%20D_%7B0h%7D%29%20%3D%20%5Cint%20L%28%5Cbeta%20%7C%20D_%7B0h%7D%29%5E%7Ba_%7B0h%7D%7D%20%5Cpi_0%28%5Cbeta%29%5E%7B1%2FH%7D%20d%5Cbeta "Z(a_{0h} | D_{0h}) = \int L(\beta | D_{0h})^{a_{0h}} \pi_0(\beta)^{1/H} d\beta").

In `hdbayes`, there is one function to estimate the normalizing constant
across a grid of values for
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}") and
another to obtain posterior samples of the normalized power prior.

The posterior under the NPP may be summarized as

![\begin{align\*}
  &y_i \| x_i, \beta \sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \\\\i = 1, \ldots, n\\
  &\pi\_{\text{NPP}}(\beta, a_0 \| D_0) = \pi_0(\beta) \prod\_{h=1}^H \frac{ L(\beta \| D\_{0h})^{a\_{0h}} }{Z(a\_{0h} \| D\_{0h})} \pi(a\_{0h})
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20%26y_i%20%7C%20x_i%2C%20%5Cbeta%20%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%0A%20%20%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n%5C%5C%0A%20%20%26%5Cpi_%7B%5Ctext%7BNPP%7D%7D%28%5Cbeta%2C%20a_0%20%7C%20D_0%29%20%3D%20%5Cpi_0%28%5Cbeta%29%20%5Cprod_%7Bh%3D1%7D%5EH%20%5Cfrac%7B%20L%28%5Cbeta%20%7C%20D_%7B0h%7D%29%5E%7Ba_%7B0h%7D%7D%20%7D%7BZ%28a_%7B0h%7D%20%7C%20D_%7B0h%7D%29%7D%20%5Cpi%28a_%7B0h%7D%29%0A%5Cend%7Balign%2A%7D "\begin{align*}
  &y_i | x_i, \beta \sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \ \ i = 1, \ldots, n\\
  &\pi_{\text{NPP}}(\beta, a_0 | D_0) = \pi_0(\beta) \prod_{h=1}^H \frac{ L(\beta | D_{0h})^{a_{0h}} }{Z(a_{0h} | D_{0h})} \pi(a_{0h})
\end{align*}")

The default initial prior on
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") for NPP
is the same as that for PP:

![\beta \sim N_p(\textbf{0}\_p, 100 \times I_p)](https://latex.codecogs.com/png.latex?%5Cbeta%20%5Csim%20N_p%28%5Ctextbf%7B0%7D_p%2C%20100%20%5Ctimes%20I_p%29 "\beta \sim N_p(\textbf{0}_p, 100 \times I_p)")

The default prior on each
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}") in
`hdbayes` is:

![\pi_0(a\_{0h}) \propto a\_{0h}^{\alpha_0 - 1} (1 - a\_{0h})^{\gamma_0 - 1}](https://latex.codecogs.com/png.latex?%5Cpi_0%28a_%7B0h%7D%29%20%5Cpropto%20a_%7B0h%7D%5E%7B%5Calpha_0%20-%201%7D%20%281%20-%20a_%7B0h%7D%29%5E%7B%5Cgamma_0%20-%201%7D "\pi_0(a_{0h}) \propto a_{0h}^{\alpha_0 - 1} (1 - a_{0h})^{\gamma_0 - 1}")

with
![\alpha\_{0} = \gamma\_{0} = 1](https://latex.codecogs.com/png.latex?%5Calpha_%7B0%7D%20%3D%20%5Cgamma_%7B0%7D%20%3D%201 "\alpha_{0} = \gamma_{0} = 1")

When
![\alpha\_{0} = 1](https://latex.codecogs.com/png.latex?%5Calpha_%7B0%7D%20%3D%201 "\alpha_{0} = 1")
and
![\gamma\_{0} = 1](https://latex.codecogs.com/png.latex?%5Cgamma_%7B0%7D%20%3D%201 "\gamma_{0} = 1"),
the prior on
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}") is
a ![U(0,1)](https://latex.codecogs.com/png.latex?U%280%2C1%29 "U(0,1)")
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
#> 1 0.000   0.000     9947.309    1.000
#> 2 0.014  -4.230     6370.440    1.001
#> 3 0.027  -6.285     6449.416    1.001
#> 4 0.041  -7.859     6204.290    1.001
#> 5 0.054  -9.212     6420.296    1.000
#> 6 0.068 -10.439     6630.803    1.001
```

The provided function `glm.npp.lognc` estimates the logarithm of the
normalizing constant,
![\log Z(a\_{0h} \| D\_{0h})](https://latex.codecogs.com/png.latex?%5Clog%20Z%28a_%7B0h%7D%20%7C%20D_%7B0h%7D%29 "\log Z(a_{0h} | D_{0h})"),
for one specific value of
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}") and
one historical data set
![D\_{0h}](https://latex.codecogs.com/png.latex?D_%7B0h%7D "D_{0h}"). We
created the function `logncfun` so that the first argument would be
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}"),
allowing us to use the `parLapply` function in the `parallel` package.

The `hdbayes` function `glm.npp.lognc` outputs
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}"),
![Z(a\_{0h} \| D\_{0h})](https://latex.codecogs.com/png.latex?Z%28a_%7B0h%7D%20%7C%20D_%7B0h%7D%29 "Z(a_{0h} | D_{0h})"),
and the minimum bulk effective sample size and maximum R-hat value of
the MCMC sampling of the power prior. It is a good idea to check that
the minimum bulk effective sample size is at least 1,000 and the maximum
R-hat value is less than 1.10

``` r
min(a0.lognc$min_ess_bulk) ## lowest bulk effective sample size
#> [1] 6059.067
max(a0.lognc$max_Rhat)  ## highest R-hat value
#> [1] 1.001153
```

We can then plot the logarithm of the normalizing constant

``` r
plot(a0.lognc$a0, a0.lognc$lognc)
```

![](README_files/figure-gfm/npp_lognc_plot-1.png)<!-- -->

### Sampling from the posterior distribution

We can now sample from the posterior distribution. The function
`glm.npp` takes, as input, values of
![a\_{0h}](https://latex.codecogs.com/png.latex?a_%7B0h%7D "a_{0h}") and
the estimated logarithm of the normalizing constant. Linear
interpolation is used to estimate
![Z(a\_{0h} \| D\_{0h})](https://latex.codecogs.com/png.latex?Z%28a_%7B0h%7D%20%7C%20D_%7B0h%7D%29 "Z(a_{0h} | D_{0h})")
for values not in the fine grid. Thus, it may be a good idea to conduct
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
#> Chain 2 finished in 0.9 seconds.
#> Chain 3 finished in 0.9 seconds.
#> Chain 1 finished in 1.0 seconds.
#> Chain 4 finished in 0.9 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.9 seconds.
#> Total execution time: 1.1 seconds.
fit.npp[, -c(1, 6)] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 4 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.855  0.851 0.239 0.244  0.463  1.25   1       4680.    5065.
#> 2 z            0.596  0.591 0.27  0.271  0.16   1.05   1.00    5529.    4746.
#> 3 x           -0.797 -0.793 0.164 0.162 -1.07  -0.536  1.00    5059.    4883.
#> 4 a0_hist_1    0.677  0.714 0.226 0.255  0.254  0.976  1       4909.    3632.
```

### Normalized asymptotic power prior (NAPP)

NAPP uses a large sample theory argument to formulate a normal
approximation to the power prior, i.e., the prior is given by

![\beta \| a\_{0h} \sim N(\hat{\beta}\_0, a\_{0h}^{-1} \[I_n(\beta)\]^{-1}),](https://latex.codecogs.com/png.latex?%5Cbeta%20%7C%20a_%7B0h%7D%20%5Csim%20N%28%5Chat%7B%5Cbeta%7D_0%2C%20a_%7B0h%7D%5E%7B-1%7D%20%5BI_n%28%5Cbeta%29%5D%5E%7B-1%7D%29%2C "\beta | a_{0h} \sim N(\hat{\beta}_0, a_{0h}^{-1} [I_n(\beta)]^{-1}),")

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
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: multi_normal_lpdf: Covariance matrix is not symmetric. Covariance matrix[1,2] = -inf, but Covariance matrix[2,1] = -inf (in '/var/folders/_g/wwn3pnf56mx5q1kv0t5574y80000gn/T/RtmpFSR3a9/model-fe8850461191.stan', line 129, column 6 to column 89)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3
#> Chain 1 finished in 0.6 seconds.
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
#> # A tibble: 5 × 10
#>   variable       mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>         <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)   0.839  0.837 0.233 0.233  0.466  1.23   1.00    4035.    4757.
#> 2 z             0.6    0.599 0.268 0.263  0.162  1.04   1.00    4589.    3985.
#> 3 x            -0.787 -0.786 0.162 0.159 -1.06  -0.519  1.00    4323.    4214.
#> 4 a0_hist_1     0.672  0.707 0.227 0.258  0.254  0.972  1       5792.    4259.
#> 5 logit_a0s[1]  1.01   0.883 1.43  1.36  -1.08   3.54   1       5792.    4259.
```

### Latent exchangeability prior (LEAP)

The LEAP assumes that the historical data are generated from a finite
mixture model consisting of
![K \ge 2](https://latex.codecogs.com/png.latex?K%20%5Cge%202 "K \ge 2")
components, with the current data generated from the first component of
this mixture. For single historical data set settings, the posterior
under the LEAP may be expressed hierarchically as

![\begin{align\*}
  y_i \| x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \\\\i = 1, \ldots, n\\
  y\_{0i} \| x\_{0i}, \beta, \beta\_{0k}, \gamma &\sim \gamma_1 \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta) \right) + \sum\_{k=2}^K \gamma_k \text{Bernoulli}\left( \text{logit}^{-1}(x\_{0i}'\beta\_{0k}) \right)
  , \\\\i = 1, \ldots, n_0
  \\
  \beta, \beta\_{0k} &\overset{\text{i.i.d.}}{\sim} N(\mu_0, \Sigma\_{0}), \\k = 2, \ldots, K, \text{ where }\Sigma\_{0} = \text{diag}(\sigma\_{01}^2, \ldots, \sigma\_{0p}^2) \\
  \gamma &\sim \text{Dirichlet}(\alpha\_{0})
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20y_i%20%7C%20x_i%2C%20%5Cbeta%20%26%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%0A%20%20%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n%5C%5C%0A%20%20y_%7B0i%7D%20%7C%20x_%7B0i%7D%2C%20%5Cbeta%2C%20%5Cbeta_%7B0k%7D%2C%20%5Cgamma%20%26%5Csim%20%5Cgamma_1%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta%29%20%5Cright%29%20%2B%20%5Csum_%7Bk%3D2%7D%5EK%20%5Cgamma_k%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_%7B0i%7D%27%5Cbeta_%7B0k%7D%29%20%5Cright%29%0A%20%20%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n_0%0A%20%20%5C%5C%0A%20%20%5Cbeta%2C%20%5Cbeta_%7B0k%7D%20%26%5Coverset%7B%5Ctext%7Bi.i.d.%7D%7D%7B%5Csim%7D%20N%28%5Cmu_0%2C%20%5CSigma_%7B0%7D%29%2C%20%5C%3Ak%20%3D%202%2C%20%5Cldots%2C%20K%2C%20%5Ctext%7B%20where%20%7D%5CSigma_%7B0%7D%20%3D%20%5Ctext%7Bdiag%7D%28%5Csigma_%7B01%7D%5E2%2C%20%5Cldots%2C%20%5Csigma_%7B0p%7D%5E2%29%20%5C%5C%0A%20%20%5Cgamma%20%26%5Csim%20%5Ctext%7BDirichlet%7D%28%5Calpha_%7B0%7D%29%0A%5Cend%7Balign%2A%7D "\begin{align*}
  y_i | x_i, \beta &\sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \ \ i = 1, \ldots, n\\
  y_{0i} | x_{0i}, \beta, \beta_{0k}, \gamma &\sim \gamma_1 \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta) \right) + \sum_{k=2}^K \gamma_k \text{Bernoulli}\left( \text{logit}^{-1}(x_{0i}'\beta_{0k}) \right)
  , \ \ i = 1, \ldots, n_0
  \\
  \beta, \beta_{0k} &\overset{\text{i.i.d.}}{\sim} N(\mu_0, \Sigma_{0}), \:k = 2, \ldots, K, \text{ where }\Sigma_{0} = \text{diag}(\sigma_{01}^2, \ldots, \sigma_{0p}^2) \\
  \gamma &\sim \text{Dirichlet}(\alpha_{0})
\end{align*}")

where
![\gamma = (\gamma_1, \ldots, \gamma_K)'](https://latex.codecogs.com/png.latex?%5Cgamma%20%3D%20%28%5Cgamma_1%2C%20%5Cldots%2C%20%5Cgamma_K%29%27 "\gamma = (\gamma_1, \ldots, \gamma_K)'")
is a vector of mixing probabilities,
![\alpha\_{0} = (\alpha_1, \ldots, \alpha_K)'](https://latex.codecogs.com/png.latex?%5Calpha_%7B0%7D%20%3D%20%28%5Calpha_1%2C%20%5Cldots%2C%20%5Calpha_K%29%27 "\alpha_{0} = (\alpha_1, \ldots, \alpha_K)'")
is a vector of concentration parameters, and
![\mu_0](https://latex.codecogs.com/png.latex?%5Cmu_0 "\mu_0") and
![\Sigma\_{0}](https://latex.codecogs.com/png.latex?%5CSigma_%7B0%7D "\Sigma_{0}")
are respectively the prior mean and covariance matrices for the
![K](https://latex.codecogs.com/png.latex?K "K") regression
coefficients. The defaults in `hdbayes` are

- ![\mu_0 = \textbf{0}\_p](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%20%5Ctextbf%7B0%7D_p "\mu_0 = \textbf{0}_p")
- ![\Sigma\_{0} = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_%7B0%7D%20%3D%20100%20%5Ctimes%20I_p "\Sigma_{0} = 100 \times I_p")
- ![\alpha\_{0} = \textbf{1}\_K](https://latex.codecogs.com/png.latex?%5Calpha_%7B0%7D%20%3D%20%5Ctextbf%7B1%7D_K "\alpha_{0} = \textbf{1}_K"),
  where
  ![\textbf{1}\_K](https://latex.codecogs.com/png.latex?%5Ctextbf%7B1%7D_K "\textbf{1}_K")
  denotes a ![K](https://latex.codecogs.com/png.latex?K "K")-dimensional
  vector of 1s
- ![K = 2](https://latex.codecogs.com/png.latex?K%20%3D%202 "K = 2")

For multiple historical data sets, `hdbayes` assumes that all historical
data sets come from a single finite mixture model.

The LEAP can be fit as follows:

``` r
fit.leap = glm.leap(
  formula = formula, family = family, data.list = data.list,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 2 finished in 4.0 seconds.
#> Chain 3 finished in 4.6 seconds.
#> Chain 4 finished in 4.8 seconds.
#> Chain 1 finished in 5.8 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 4.8 seconds.
#> Total execution time: 5.9 seconds.

suppressWarnings(
 fit.leap[, 2:6] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3)) 
)
#> # A tibble: 5 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.85   0.847 0.245 0.245  0.465  1.26      1    2920.    3210.
#> 2 z            0.602  0.599 0.277 0.277  0.155  1.06      1    3595.    3336.
#> 3 x           -0.789 -0.785 0.17  0.172 -1.07  -0.518     1    2849.    2393.
#> 4 probs[1]     0.863  0.889 0.108 0.101  0.647  0.99      1    1395.    2503.
#> 5 probs[2]     0.137  0.111 0.108 0.101  0.01   0.353     1    1395.    2503.
```

### Normal/half-normal prior

We also include a normal/half-normal prior, where the regression
coefficients are assigned independent normal priors, and, if applicable,
the dispersion parameter is assigned a half-normal prior. For this
example, the normal/half-normal prior takes the form

![\begin{align\*}
  &y_i \| x_i, \beta \sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \\\\i = 1, \ldots, n\\
  &\beta \sim N_p(\mu_0, \Sigma_0), \text{ where }\Sigma_0 = \text{diag}(\sigma\_{01}^2, \ldots, \sigma\_{0p}^2).
\end{align\*}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%2A%7D%0A%20%20%26y_i%20%7C%20x_i%2C%20%5Cbeta%20%5Csim%20%5Ctext%7BBernoulli%7D%5Cleft%28%20%5Ctext%7Blogit%7D%5E%7B-1%7D%28x_i%27%5Cbeta%29%20%5Cright%29%0A%20%20%2C%20%5C%20%5C%20i%20%3D%201%2C%20%5Cldots%2C%20n%5C%5C%0A%20%20%26%5Cbeta%20%5Csim%20N_p%28%5Cmu_0%2C%20%5CSigma_0%29%2C%20%5Ctext%7B%20where%20%7D%5CSigma_0%20%3D%20%5Ctext%7Bdiag%7D%28%5Csigma_%7B01%7D%5E2%2C%20%5Cldots%2C%20%5Csigma_%7B0p%7D%5E2%29.%0A%5Cend%7Balign%2A%7D "\begin{align*}
  &y_i | x_i, \beta \sim \text{Bernoulli}\left( \text{logit}^{-1}(x_i'\beta) \right)
  , \ \ i = 1, \ldots, n\\
  &\beta \sim N_p(\mu_0, \Sigma_0), \text{ where }\Sigma_0 = \text{diag}(\sigma_{01}^2, \ldots, \sigma_{0p}^2).
\end{align*}")

The defaults in `hdbayes` are

- ![\mu_0 = \textbf{0}\_p](https://latex.codecogs.com/png.latex?%5Cmu_0%20%3D%20%5Ctextbf%7B0%7D_p "\mu_0 = \textbf{0}_p")
- ![\Sigma\_{0} = 100 \times I_p](https://latex.codecogs.com/png.latex?%5CSigma_%7B0%7D%20%3D%20100%20%5Ctimes%20I_p "\Sigma_{0} = 100 \times I_p")

The normal/half-normal prior can be fit as follows:

``` r
fit.post = glm.post(
  formula = formula, family = family, data.list = data.list,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
#> Running MCMC with 4 chains, at most 15 in parallel...
#> 
#> Chain 1 finished in 0.5 seconds.
#> Chain 2 finished in 0.4 seconds.
#> Chain 3 finished in 0.4 seconds.
#> Chain 4 finished in 0.5 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.5 seconds.
#> Total execution time: 0.6 seconds.

fit.post[, -1] %>% 
    summarise_draws() %>% 
    mutate(across(where(is.numeric), round, 3))
#> # A tibble: 3 × 10
#>   variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 (Intercept)  0.778  0.772 0.279 0.272  0.328  1.25   1       3631.    3699.
#> 2 z            0.693  0.697 0.304 0.301  0.19   1.20   1.00    4608.    4387.
#> 3 x           -0.772 -0.767 0.186 0.182 -1.08  -0.471  1.00    3524.    3683.
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
                'npp' = fit.npp, 'pp' = fit.pp, 'leap' = fit.leap,
                'normal/half-normal' = fit.post)

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
#> (Intercept)   1.0   0.760    1.036  0.842        0.838     0.849  0.839  0.855
#> z             0.5   0.677    0.313  0.619        0.625     0.617  0.600  0.596
#> x            -1.0  -0.750   -0.856 -0.788       -0.786    -0.791 -0.787 -0.797
#>                 pp   leap normal/half-normal
#> (Intercept)  0.834  0.850              0.778
#> z            0.613  0.602              0.693
#> x           -0.788 -0.789             -0.772

## posterior std dev.
round( post.sd, 3 )
#>             mle.cur mle.hist   bhm commensurate robustmap  napp   npp    pp
#> (Intercept)   0.274    0.377 0.246        0.251     0.253 0.233 0.239 0.251
#> z             0.308    0.442 0.286        0.280     0.280 0.268 0.270 0.280
#> x             0.181    0.266 0.166        0.170     0.169 0.162 0.164 0.165
#>              leap normal/half-normal
#> (Intercept) 0.245              0.279
#> z           0.277              0.304
#> x           0.170              0.186
```
