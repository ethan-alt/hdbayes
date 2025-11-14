# hdbayes 0.2.0

## Enhancements

* Added functions for computing model averaging weights using Bayesian model averaging (BMA), pseudo-BMA, pseudo-BMA+ (pseudo-BMA 
with the Bayesian bootstrap), and stacking. Also added a function for generating samples from the ensemble of posterior distributions 
based on the computed weights.

* Added a vignette demonstrating model averaging methods and ensemble inference.

* Added functions for sampling from the posterior distributions of several survival models, including accelerated failure 
time (AFT) models, the piecewise exponential (PWE) model, and the mixture cure rate model with a PWE component for the non-cured 
popluation (CurePWE). For all survival model implementations, multiple historical data sets are now stacked into a single 
combined data set for model fitting.

* Added functions to compute the log marginal likelihood for the AFT, PWE, and CurePWE models under various prior specifications.

* Added the implementation of stratified power prior.

* Updated the implementation of generalized linear models (GLMs) to support computation of the pointwise log-likelihood matrix,
consistent with the survival model implementations. This enables downstream use in model comparison and averaging procedures, 
such as computing pseudo-BMA and stacking weights, as well as estimating the expected log predictive density (ELPD). 

* Added lower and upper bounds for the probability of being exchangeable in the LEAP implementation for GLMs.

* Modified the output of the `glm.rmap()` function to include the updated mixture weight for the posterior density under the 
meta-analytic predictive (MAP) prior.


# hdbayes 0.1.1

## Enhancements

* Added functions for computing the logarithm of the marginal likelihood of a GLM under all priors implemented in the package.

* Updated the implementation of commensurate prior to be fully Bayesian.

* Updated the implementation of robust meta-analytic predictive prior (RMAP) from using a Gaussian mixture model to approximate the MAP prior to computing the updated mixture weights based on marginal likelihoods. Specifically, we removed `glm.rmap.bhm()` and `glm.rmap.bhm.approx()` functions. The posterior samples from using the RMAP now can be obtained via calling `glm.rmap()` function directly.

* Added function for sampling from the posterior distribution of a GLM under a normal/half-normal prior.

* Added the vignette "AIDS_Progression".

* Updated Stan files for NAPP and NPP by eliciting priors on logit(a0) instead of a0. 

* Added two data sets: E1684 and E1690.


## Bug Fixes

* Fixed bugs in checking input `offset.list` for `glm.leap()`.

* Fixed bugs in computing normalizing constants for normal and Gamma models using `glm.npp.lognc()`.

* Fixed bugs in renaming/reordering output variables in `glm.bhm()`.


# hdbayes 0.0.3

* First CRAN release.
