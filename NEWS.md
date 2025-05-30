# hdbayes 0.2.0 (development version)

## Enhancements

* Added functions for sampling from the posterior distributions of various survival models, including an accelerated failure 
time (AFT) model, a piecewise exponential (PWE) model, and a mixture cure rate (CurePWE) model.

* Added functions to compute the log marginal likelihoods of the AFT model, PWE model, and CurePWE model under various priors.

* Added lower and upper bounds for probability of being exchangeable for LEAP implementation for GLM models.

* Modified the output of the `glm.rmap()` function to include the updated mixture weight for the posterior density under the MAP prior.


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
