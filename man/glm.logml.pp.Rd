% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_logml_pp.R
\name{glm.logml.pp}
\alias{glm.logml.pp}
\title{Log marginal likelihood of a GLM under power prior (PP)}
\usage{
glm.logml.pp(
  post.samples,
  bridge.args = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  ...
)
}
\arguments{
\item{post.samples}{output from \code{\link[=glm.pp]{glm.pp()}} giving posterior samples of a GLM under the power prior (PP), with
an attribute called 'data' which includes the list of variables specified in the data block
of the Stan program.}

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
If all of the power prior parameters (\eqn{a_0}'s) are equal to zero, or if the posterior samples are obtained from using only
one data set (the current data), then the function will return the same result as the output from \code{\link[=glm.logml.post]{glm.logml.post()}}.

If at least one of the power prior parameters (\eqn{a_0}'s) is non-zero, the function will return a \code{list} with the following objects

\describe{
\item{model}{"PP"}

\item{logml}{the estimated logarithm of the marginal likelihood}

\item{bs}{an object of class \code{bridge} or \code{bridge_list} containing the output from using \code{\link[bridgesampling:bridge_sampler]{bridgesampling::bridge_sampler()}}
to compute the logarithm of the normalizing constant of the power prior (PP) using all data sets}

\item{bs.hist}{an object of class \code{bridge} or \code{bridge_list} containing the output from using
\code{\link[bridgesampling:bridge_sampler]{bridgesampling::bridge_sampler()}} to compute the logarithm of the normalizing constant of the PP using historical
data sets}

\item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}

\item{max_Rhat}{the maximum Rhat}
}
}
\description{
Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
likelihood of a GLM under the power prior (PP).

The arguments related to MCMC sampling are utilized to draw samples from the power prior (PP).
These samples are then used to compute the logarithm of the normalizing constant of the PP using only historical
data sets.
}
\examples{
if (instantiate::stan_cmdstan_exists()) {
  data(actg019)
  data(actg036)
  ## take subset for speed purposes
  actg019 = actg019[1:100, ]
  actg036 = actg036[1:50, ]
  data_list = list(currdata = actg019, histdata = actg036)
  formula = cd4 ~ treatment + age + race
  family = poisson('log')
  a0 = 0.5
  d.pp = glm.pp(
    formula = formula, family = family,
    data.list = data_list,
    a0.vals = a0,
    chains = 1, iter_warmup = 500, iter_sampling = 1000
  )
  glm.logml.pp(
    post.samples = d.pp,
    bridge.args = list(silent = TRUE),
    chains = 1, iter_warmup = 1000, iter_sampling = 2000
  )
}
}
\references{
Chen, M.-H. and Ibrahim, J. G. (2000). Power prior distributions for Regression Models. Statistical Science, 15(1).

Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
}
