# DPMM 2.0.2

Performance and correctness pass over model fitting and posterior sampling. No user-facing API changes; verified against the pre-change code across every supported data shape (continuous-only, categorical-only, mixed) and missingness pattern using fixed seeds.

* `runModel()` no longer uses `synthpop::syn()` to build the placeholder row stored in `$dataset`. Only the schema (column names/classes/factor levels) of that row is ever used downstream, not its values, so it's now built directly by sampling each column independently from its own observed values/levels — about 100x faster on a benchmark dataset and removes `synthpop` as a package dependency entirely.
* Removed redundant duplicate density calculations in `posterior_dpmm()`'s posterior-predictive draws (the same per-component log-likelihood was computed twice — once for the normalising constant, once for the conditional weights).
* Replaced the `tidyr`-based reshaping of MCMC samples in `posterior_dpmm()` (`muL`/`tauL`/`phiL`/`w`) with direct index-based array construction, removing a significant per-call overhead. This also fixes a latent bug where `postW`'s component ordering relied on a broken column-name parse that happened to still produce the right order by luck of a stable sort (surfaced as a silent `NAs introduced by coercion` warning).
* Fixed a bug in `posterior_dpmm()`'s mixed-data branch that crashed with `unused argument (phi = ...)` whenever a patient had every continuous variable missing but at least one categorical variable observed.
* Reduced duplicated code in `runModel()` (~1550 → ~1260 lines) by extracting the R-level scaffolding shared across all 8 data-shape/missingness combinations (initial-value draws, the model-build retry loop, and the custom missing-data sampler assignment) into internal helper functions. Each combination's `nimbleCode` model definition is unchanged, since nimble requires a fixed model graph per data shape.

# DPMM 2.0.1

* Added a `seed` argument to `runModel()`, `predict_dpmm_fit()`, and `plot_ggpairs()` so model fitting and posterior sampling can be made reproducible. `posterior_dpmm()` already supported this and is unchanged.
* `predict_dpmm_fit()` now passes its `seed` argument through to `posterior_dpmm()` instead of always calling it with `seed = NULL`.
* Added this `NEWS.md` file to track user-facing changes going forward.

# DPMM 2.0.0

* See git history prior to this file for changes up to this release.
