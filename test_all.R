library(hdbayes)

## obtain number of cores
ncores = max(1, parallel::detectCores() - 1)
warmup  = 1000
total.samples = 10000   ## number of samples post warmup
samples = ceiling(warmup + total.samples / ncores)  ## outputs approx total.samples samples

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

## set formula and family--conduct frequentist analysis of data sets
formula = y ~ z + x
family  = binomial('logit')

fit.mle.cur  = glm(formula, family, data)
fit.mle.hist = glm(formula, family, histdata)


## Bayesian Hierarchical Model
fit.bhm = glm.bhm(
  formula, family, data, histdata,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup
)

## Commensurate prior
fit.commensurate = glm.commensurate(
  formula, family, data, histdata, tau = rep(5, 3),
  cores = ncores, chains = ncores, iter = samples, warmup = warmup
)

## Robust MAP prior
fit.robustmap = glm.robustmap(
  formula, family, data, histdata,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup
)


## normalized asymptotic power prior
fit.napp = glm.napp(
  formula, family, data, histdata,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup
)

##
## NORMALIZED POWER PRIOR
##
  ## parallelize estimation of log normalizing constant
  library(parallel)
  a0     = seq(0, 1, length.out = ncores * 5)

  ## wrapper to obtain log normalizing constant in parallel package
  logncfun = function(a0, ...)
    hdbayes::glm.npp.lognc(formula = formula, family = family, histdata = histdata, a0 = a0, ...)

  cl = makeCluster(ncores)
    clusterSetRNGStream(cl, 123)
    clusterExport(cl, varlist = c('formula', 'family', 'histdata'))
    a0.lognc = parLapply(cl = cl, X = a0, fun = logncfun, iter = 5000, warmup = warmup)
  stopCluster(cl)

  a0.lognc = data.frame( do.call(rbind, a0.lognc) )
  plot(a0.lognc$a0, a0.lognc$lognc)
  min(a0.lognc$min_n_eff)
  max(a0.lognc$max_Rhat)

  ## fit normalized power prior
  fit.npp = glm.npp(
    formula, family, data, histdata, a0 = a0.lognc$a0, lognc = a0.lognc$lognc,
    cores = ncores, chains = ncores, iter = samples, warmup = warmup
  )


##
## POWER PRIOR (fixed a0)
##
fit.pp = glm.pp(
  formula, family, data, histdata, a0 = 0.5,
  cores = ncores, chains = ncores, iter = samples, warmup = warmup
)



##
## COMPARE METHODS
##
fit.list = list('bhm' = fit.bhm, 'commensurate' = fit.commensurate,
                'robustmap' = fit.robustmap, 'napp' = fit.napp,
                'npp' = fit.npp, 'pp' = fit.pp)

post.mean = sapply(fit.list, function(x) summary(x)$summary[names(coef(fit.mle.cur)), 'mean'])
post.sd   = sapply(fit.list, function(x) summary(x)$summary[names(coef(fit.mle.cur)), 'sd'])

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
post.mean

## posterior std dev.
post.sd

