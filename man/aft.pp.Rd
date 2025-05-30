% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aft_pp.R
\name{aft.pp}
\alias{aft.pp}
\title{Posterior of power prior (PP) with fixed \eqn{a_0}}
\usage{
aft.pp(
  formula,
  data.list,
  a0,
  dist = "weibull",
  beta.mean = NULL,
  beta.sd = NULL,
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

\item{a0}{a scalar between 0 and 1 giving the (fixed) power prior parameter for the historical data.}

\item{dist}{a character indicating the distribution of survival times. Currently, \code{dist} can be one of the
following values: "weibull", "lognormal", or "loglogistic". Defaults to "weibull".}

\item{beta.mean}{a scalar or a vector whose dimension is equal to the number of regression coefficients giving
the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
\code{beta.mean} will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.}

\item{beta.sd}{a scalar or a vector whose dimension is equal to the number of regression coefficients giving
the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
same as for \code{beta.mean}. Defaults to a vector of 10s.}

\item{scale.mean}{location parameter for the half-normal prior on the scale parameter of the AFT model. Defaults to 0.}

\item{scale.sd}{scale parameter for the half-normal prior on the scale parameter of the AFT model. Defaults to 10.}

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
Sample from the posterior distribution of an accelerated failure time (AFT) model using the power prior (PP) by
Ibrahim and Chen (2000) \url{doi:10.1214/ss/1009212673}.
}
\details{
The power prior parameters (\eqn{a_0}'s) are treated as fixed. The initial priors on the regression coefficients
are independent normal priors. The current and historical data models are assumed to have a common scale parameter
with a half-normal prior.
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
    aft.pp(
      formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
      data.list = data_list,
      a0 = 0.5,
      dist = "weibull",
      chains = 1, iter_warmup = 500, iter_sampling = 1000
    )
  }
}
}
\references{
Chen, M.-H. and Ibrahim, J. G. (2000). Power prior distributions for Regression Models. Statistical Science, 15(1).
}
