% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curepwe_logml_commensurate.R
\name{curepwe.logml.commensurate}
\alias{curepwe.logml.commensurate}
\title{Log marginal likelihood of a mixture cure rate (CurePWE) under the commensurate prior (CP)}
\usage{
curepwe.logml.commensurate(
  post.samples,
  bridge.args = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  ...
)
}
\arguments{
\item{post.samples}{output from \code{\link[=curepwe.commensurate]{curepwe.commensurate()}} giving posterior samples of a CurePWE model under the
commensurate prior (CP), with an attribute called 'data' which includes the list of variables
specified in the data block of the Stan program.}

\item{bridge.args}{a \code{list} giving arguments (other than \code{samples}, \code{log_posterior}, \code{data}, \code{lb}, and \code{ub}) to
pass onto \code{\link[bridgesampling:bridge_sampler]{bridgesampling::bridge_sampler()}}.}

\item{iter_warmup}{number of warmup iterations to run per chain. Defaults to 1000. See the argument \code{iter_warmup}
in \code{sample()} method in cmdstanr package.}

\item{iter_sampling}{number of post-warmup iterations to run per chain. Defaults to 1000. See the argument \code{iter_sampling}
in \code{sample()} method in cmdstanr package.}

\item{chains}{number of Markov chains to run. Defaults to 4. See the argument \code{chains} in \code{sample()} method
in cmdstanr package.}

\item{...}{arguments passed to \code{sample()} method in cmdstanr package (e.g., \code{seed}, \code{refresh}, \code{init}).}
}
\value{
The function returns a \code{list} with the following objects

\describe{
\item{model}{"Commensurate"}

\item{logml}{the estimated logarithm of the marginal likelihood}

\item{bs}{an object of class \code{bridge} or \code{bridge_list} containing the output from using \code{\link[bridgesampling:bridge_sampler]{bridgesampling::bridge_sampler()}}
to compute the logarithm of the normalizing constant of the commensurate prior (CP) using all data sets}

\item{bs.hist}{an object of class \code{bridge} or \code{bridge_list} containing the output from using
\code{\link[bridgesampling:bridge_sampler]{bridgesampling::bridge_sampler()}} to compute the logarithm of the normalizing constant of the CP using historical
data sets}

\item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}

\item{max_Rhat}{the maximum Rhat}
}
}
\description{
Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
likelihood of a CurePWE model under the commensurate prior (CP).

The arguments related to MCMC sampling are utilized to draw samples from the commensurate prior.
These samples are then used to compute the logarithm of the normalizing constant of the commensurate prior using
historical data sets.
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
    nbreaks = 3
    probs   = 1:nbreaks / nbreaks
    breaks  = as.numeric(
      quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
    )
    breaks  = c(0, breaks)
    breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
    d.cp = curepwe.commensurate(
      formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
      data.list = data_list,
      breaks = breaks,
      p.spike = 0.1,
      chains = 1, iter_warmup = 500, iter_sampling = 1000
    )
    curepwe.logml.commensurate(
      post.samples = d.cp,
      bridge.args = list(silent = TRUE),
      chains = 1, iter_warmup = 500, iter_sampling = 1000
    )
  }
}
}
\references{
Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.

Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
}
