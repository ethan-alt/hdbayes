# hdbayes (development version)

## Enhancements

* Added function for sampling from the posterior distribution of a GLM under a normal/half-normal prior.

* Added functions for computing the logarithm of the marginal likelihood of a GLM under all priors implemented in the package except for the robust meta-analytic predictive prior (RMAP).

* Changed the implementation of robust meta-analytic predictive prior (RMAP) from using a Gaussian mixture model to approximate the MAP prior to computing the updated mixture weights based on marginal likelihoods. Specifically, we removed `glm.rmap.bhm()` and `glm.rmap.bhm.approx()` functions. The posterior samples from using the RMAP now can be obtained via calling `glm.rmap()` function directly.

* Updated Stan files for NAPP and NPP by eliciting priors on logit(a0) instead of a0. 

* Added two data sets: E1684 and E1690.


## Bug Fixes

* Fixed bugs in checking input `offset.list` for `glm.leap()`.

* Fixed bugs in computing normalizing constants for normal and Gamma models using `glm.npp.lognc()`.

* Fixed bugs in renaming/reordering output variables in `glm.bhm()`.


# hdbayes 0.0.3

* First CRAN release.
