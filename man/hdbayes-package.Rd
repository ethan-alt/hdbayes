% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdbayes-package.R
\docType{package}
\name{hdbayes-package}
\alias{hdbayes-package}
\alias{hdbayes}
\title{hdbayes: Bayesian Analysis of Generalized Linear Models with Historical Data}
\description{
User-friendly functions for leveraging (multiple) historical data set(s) for generalized linear models (GLMs) and survival models. The package provides functions for sampling from the posterior distribution under various informative priors, including the prior induced by the Bayesian hierarchical model, power prior by Ibrahim and Chen (2000) \doi{10.1214/ss/1009212673}, normalized power prior by Duan et al. (2006) \doi{10.1002/env.752}, normalized asymptotic power prior by Ibrahim et al. (2015) \doi{10.1002/sim.6728}, commensurate prior by Hobbs et al. (2011) \doi{10.1111/j.1541-0420.2011.01564.x}, robust meta-analytic-predictive prior by Schmidli et al. (2014) \doi{10.1111/biom.12242}, latent exchangeability prior by Alt et al. (2024) \doi{10.1093/biomtc/ujae083}, and a normal (or half-normal) prior. In addition to GLMs, the package supports survival models including: (1) accelerated failure time (AFT) models, (2) piecewise exponential (PWE) models, i.e., proportional hazards models with piecewise constant baseline hazards, and (3) mixture cure rate models that assume a common probability of cure across subjects, paired with a PWE model for the non-cured population. Functions for computing the marginal log-likelihood under each implemented prior are also included. The package compiles all the 'CmdStan' models once during installation using the 'instantiate' package.

Bayesian analysis of generalized linear models using historical data
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ethan-alt/hdbayes}
  \item Report bugs at \url{https://github.com/ethan-alt/hdbayes/issues}
}

}
\author{
\strong{Maintainer}: Ethan M. Alt \email{ethanalt@live.unc.edu} (\href{https://orcid.org/0000-0002-6112-9030}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Xinxin Chen
  \item Luiz M. Carvalho
  \item Joseph G. Ibrahim
}

Other contributors:
\itemize{
  \item Xiuya Chang [contributor]
}

}
\concept{help}
\keyword{internal}
