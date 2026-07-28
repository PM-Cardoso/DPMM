# DPMM 2.0.4

* **Fixed:** `runModel()` seeded every MCMC chain with the *same* value (`setSeed = rep(seed, mcmc_chains)`), rather than giving each chain its own seed. Per nimble's `runMCMC()` documentation, a `setSeed` vector should provide a distinct seed per chain -- reusing one seed for every chain re-initialises each chain's random stream to the same state, undermining the independence between chains that `mcmc_chains > 1` (and downstream chain-comparison diagnostics like `plot_alpha()` and `plot_ggpairs(compare_chains = TRUE)`) rely on. Now uses `seed + seq_len(mcmc_chains) - 1`, giving each chain its own deterministic (and hence still fully reproducible) seed. Single-chain fits (`mcmc_chains = 1`) are unaffected.

# DPMM 2.0.3

Correctness fix and a major speed-up for posterior prediction. Verified with a fixed-seed smoke test across every supported data shape and missingness pattern combination, an old-vs-new speed/shape benchmark, and a full run of both vignettes (all 6 bundled datasets, standardised and non-standardised fits).

* **Fixed:** `predict_dpmm_fit()` now automatically converts continuous predictions back to their original scale when the model was fit with `standardise = TRUE`. Previously it returned predictions in the standardised space, silently, with no way to tell from the output alone -- callers had to know to manually invert `(x - mean) / sd` themselves. `plot_ggpairs()` had a manual (and duplicated) workaround for this internally; that workaround is now dead code and has been removed, since `predict_dpmm_fit()` handles it for every caller.
* **Fixed:** a precision-matrix "symmetrise" step in `posterior_dpmm()`'s continuous value-drawing code referenced the whole list-column of per-iteration matrices instead of the single matrix being processed, due to R's lexical scoping inside a `dplyr::mutate()`/`purrr::map()` closure. This only produced correct output by accident when exactly one posterior sample was requested; anything else silently used garbage-reshaped values. Replaced with a plain, correct `(tau + t(tau)) / 2`.
* **Performance:** rewrote `posterior_dpmm()`'s prediction loop. Patients are now grouped by identical missingness pattern and processed with vectorised linear algebra (matrix-input `dmvnorm`, batched conditional draws sharing one Cholesky factor per iteration/component/group) instead of one `tidyr`/`purrr`-based pass per patient per posterior sample. Component-covariance inversion is computed once per (iteration, component) and shared across every patient, instead of being redone per patient. All output is written into preallocated arrays instead of growing objects in the sampling loop. ~64x faster in a benchmark (150 patients, 200 posterior samples, mixed continuous/categorical data: 18.2s -> 0.28s); scales better the more patients share a missingness pattern. Patients are always returned in the exact order supplied, regardless of internal grouping.
* Removed `condMVNorm` and `tidyr` as package dependencies -- no longer used anywhere after the `posterior_dpmm()` rewrite.

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
