% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aft_bhm.R
\name{aft.bhm}
\alias{aft.bhm}
\title{Posterior of Bayesian hierarchical model (BHM)}
\usage{
aft.bhm(
  formula,
  data.list,
  dist = "weibull",
  meta.mean.mean = NULL,
  meta.mean.sd = NULL,
  meta.sd.mean = NULL,
  meta.sd.sd = NULL,
  scale.mean = NULL,
  scale.sd = NULL,
  get.loglik = FALSE,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  ...
)
}
\arguments{
\item{formula}{a two-sided formula giving the relationship between the response variable and covariates.
The response is a survival object as returned by the \code{survival::Surv(time, event)} function,
where event is a binary indicator for event (0 = no event, 1 = event has occurred). The type of
censoring is assumed to be right-censoring.}

\item{data.list}{a list of \code{data.frame}s. The first element in the list is the current data, and the rest
are the historical data sets. For fitting accelerated failure time (AFT) models, all historical
data sets will be stacked into one historical data set.}

\item{dist}{a character indicating the distribution of survival times. Currently, \code{dist} can be one of the
following values: "weibull", "lognormal", or "loglogistic". Defaults to "weibull".}

\item{meta.mean.mean}{a scalar or a vector whose dimension is equal to the number of regression coefficients giving
the means for the normal hyperpriors on the mean hyperparameters of regression coefficients.
If a scalar is provided, \code{meta.mean.mean} will be a vector of repeated elements of the given
scalar. Defaults to a vector of 0s.}

\item{meta.mean.sd}{a scalar or a vector whose dimension is equal to the number of regression coefficients giving
the sds for the normal hyperpriors on the mean hyperparameters of regression coefficients. If
a scalar is provided, same as for \code{meta.mean.mean}. Defaults to a vector of 10s.}

\item{meta.sd.mean}{a scalar or a vector whose dimension is equal to the number of regression coefficients giving
the means for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
If a scalar is provided, same as for \code{meta.mean.mean}. Defaults to a vector of 0s.}

\item{meta.sd.sd}{a scalar or a vector whose dimension is equal to the number of regression coefficients giving
the sds for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
If a scalar is provided, same as for \code{meta.mean.mean}. Defaults to a vector of 1s.}

\item{scale.mean}{location parameter for the half-normal prior on the scale parameters of current and historical
data models. Defaults to 0.}

\item{scale.sd}{scale parameter for the half-normal prior on the scale parameters of current and historical data
models. Defaults to 10.}

\item{get.loglik}{whether to generate log-likelihood matrix. Defaults to FALSE.}

\item{iter_warmup}{number of warmup iterations to run per chain. Defaults to 1000. See the argument \code{iter_warmup} in
\code{sample()} method in cmdstanr package.}

\item{iter_sampling}{number of post-warmup iterations to run per chain. Defaults to 1000. See the argument \code{iter_sampling}
in \code{sample()} method in cmdstanr package.}

\item{chains}{number of Markov chains to run. Defaults to 4. See the argument \code{chains} in \code{sample()} method in
cmdstanr package.}

\item{...}{arguments passed to \code{sample()} method in cmdstanr package (e.g., \code{seed}, \code{refresh}, \code{init}).}
}
\value{
The function returns an object of class \code{draws_df} giving posterior samples, with an attribute called 'data' which includes
the list of variables specified in the data block of the Stan program.
}
\description{
Sample from the posterior distribution of an accelerated failure time (AFT) model using the Bayesian
hierarchical model (BHM).
}
\details{
The Bayesian hierarchical model (BHM) assumes that the regression coefficients for the historical and
current data are different, but are correlated through a common distribution, whose hyperparameters
(i.e., mean and standard deviation (sd) (the covariance matrix is assumed to have a diagonal structure))
are treated as random. The number of regression coefficients for the current data is assumed to be the
same as that for the historical data.

The hyperpriors on the mean and the sd hyperparameters are independent normal and independent half-normal
distributions, respectively. The scale parameters for both current and historical data models are assumed
to be independent and identically distributed, each assigned a half-normal prior.
}
\examples{
if (instantiate::stan_cmdstan_exists()) {
  if(requireNamespace("survival")){
    library(survival)
    data(E1684)
    data(E1690)
    ## take subset for speed purposes
    E1684 = E1684[1:100, ]
    E1690 = E1690[1:50, ]
    ## replace 0 failure times with 0.50 days
    E1684$failtime[E1684$failtime == 0] = 0.50/365.25
    E1690$failtime[E1690$failtime == 0] = 0.50/365.25
    E1684$cage = as.numeric(scale(E1684$age))
    E1690$cage = as.numeric(scale(E1690$age))
    data_list = list(currdata = E1690, histdata = E1684)
    aft.bhm(
      formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
      data.list = data_list,
      dist = "weibull",
      chains = 1, iter_warmup = 500, iter_sampling = 1000
    )
  }
}
}
