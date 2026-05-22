# Changelog

Notable project changes are listed from newest to oldest.

## 2026-05-21 - Nonlinear Performance Follow-up

- Added a narrow R fast covariance path for internal fixed-dispersion GLM fits with `vce="HC1"` and one-way or no clustering, matching `sandwich::vcovCL()` to numerical precision while reducing `binsglm` logit/probit runtime.
- Trimmed unused R internal GLM metadata and avoided unnecessary confidence-band fit prediction work when only simulation standard errors are needed, preserving existing numerical output.
- Reduced Python quantile-regression fitting overhead by removing an unused per-iteration internal history statistic, preserving fitted coefficients and covariance calculations while lowering memory and arithmetic work in large `binsqreg` fits.
- Added a guarded Python fast path for unweighted binomial-logit generalized linear fits by using `statsmodels`' discrete `Logit` backend when it matches the previous GLM coefficients and robust covariance to machine precision; unsupported cases automatically fall back to the established GLM path.
- Fixed Python `binstest` and `binspwc` generalized-model family construction for non-default links such as `Binomial`/`Logit` and `Binomial`/`Probit`, and reused already-constructed family objects to avoid repeated parsing in GLM-heavy paths.
- Fixed Python `binsglm(..., nolink=True)` confidence-interval and confidence-band paths so linear-predictor fits and standard errors are used directly instead of relying on inverse-link branch temporaries.
- Guarded the R least-squares covariance fast path so generalized linear model objects always use the established GLM covariance backend.

## 2026-05-20 - Release Candidate: Precision Defaults And OLS Performance

- Bumped release metadata for this checkpoint: R `2.1`, Python `3.1.0`, and Stata command/help headers dated `20-MAY-2026` with distribution date `20260520`.
- Kept Stata `precision(double)` as the default so Stata uses double floating precision for internally generated variables whenever appropriate, better matching R and Python's double-precision numerical path.
- Backwards compatibility reference: Stata users who need to reproduce the legacy pre-double-precision numerical path should add `precision(single)` explicitly; this preserves the old Stata behavior that stored internally generated variables as `float`.
- Completed the OLS-focused speed checkpoint across R and Python for `binsreg`, `binsregselect`, `binstest`, and `binspwc`, including fast least-squares/covariance paths, covariance reuse, lower-allocation p-value simulation, and R spline-basis fast paths.
- Added a conservative Stata fast OLS path for unweighted double-precision robust least-squares testing in `binstest` and `binspwc`, with automatic fallback to the established regression path for unsupported cases and for `precision(single)`.
- Benchmarked a pure-R local/sparse basis prototype and left it out because dense BLAS remained faster; future sparse-basis gains should use compiled helpers rather than R-level bookkeeping.
- Revalidated the modernization round with local package checks, replication smoke checks, Stata package checks, and cross-language numerical/speed snapshots.

## 2026-05-19 - Numerical Precision And Cross-Language Alignment

- Added Stata `precision(single)` and `precision(double)` options across the public Stata commands and internal helpers.
- Backwards compatibility reference: use `precision(single)` in Stata to preserve the legacy Stata behavior that stored internally generated variables as `float`; use `precision(double)` to store those variables as `double` and better align Stata output with R and Python. The default is now `precision(double)`, so old replication scripts that need the previous Stata numerical path should add `precision(single)` explicitly.
- Rebuilt the Stata help files and help PDFs to document the new precision option, and regenerated the compiled Mata package artifacts needed by the Stata net-install package.
- Updated Stata numerical smoke fixtures for the new default double-precision path while keeping explicit `precision(single)` available for legacy comparisons.
- Aligned Python quantile calculations with Stata/R type-2 empirical quantiles for quantile-spaced knots, column medians, quantile-regression IQR calculations, and simulation critical values.
- Aligned Python interval/bin assignment at internal knots with Stata/R, assigning observations exactly on a knot to the bin on the left while keeping the right boundary in the final bin.
- Aligned Python `binsreg` plotting output bin IDs with Stata/R's one-based user-facing convention while keeping Python's internal zero-based indexing unchanged.
- Aligned the R spline basis constructor with the Stata/Python local B-spline basis, replacing the previous unsmoothed local power-basis shortcut so the three implementations share the same basis expansion contract.
- Aligned the R selector bias-location helper with Stata/Python knot-side conventions, so bias calculations use the same left-bin assignment at internal knots.
- Aligned R/Python DPI selector internals with Stata's `irecode()` bin assignment at quantile-knot boundaries, including the p=0 candidate bias and variance constants used during degree/smoothness selection.
- Optimized the R/Python spline basis constructors by vectorizing the final expansion from local nonzero basis values into the full design matrix, preserving numerical output while reducing p=2 and p=3 basis-construction overhead.
- Optimized OLS covariance paths without changing numerical output: Python now avoids statsmodels for unweighted and weighted `const`/`HC0`/`HC1`/`HC2`/`HC3` fits plus one-way clustered `HC1` fits, while R now fast-paths unweighted and weighted `const`/`HC1` fits plus one-way clustered `HC1` fits that match `vcov()`/`vcovCL()` exactly. GLM, quantile-regression, multiway-clustered, and non-`HC1` clustered variance paths continue to use the established package backends.
- Added Stata-generated fit/covariance alignment fixtures plus R and Python checkers for OLS coefficients, covariance matrices, fitted values, residuals, ranks, clustered sample counts, and rank-deficient estimable fitted-value variances including clustered designs.
- Optimized Python's internal OLS solver by reusing the rank returned by `numpy.linalg.lstsq()` and using a full-rank covariance inverse path, avoiding separate rank and pseudoinverse decompositions while preserving Stata/statsmodels-aligned numerical output.
- Reduced selector-loop bookkeeping in R/Python by reusing normalized selector inputs and local-mass knots across ROT/DPI calculations; Python also vectorizes the canonical p=0 bin projection with bin counts and sums.
- Reduced selector basis-construction work without changing numerical output: Python now reuses precomputed DPI basis objects within each selector calculation, and R adds a local basis cache for fixed-`nbinsrot` degree/smoothness sweeps with enough candidates to amortize cache overhead. R/Python also now mimic Stata's first-row tie behavior when all candidate distance measures are missing.
- Reduced repeated Python covariance extraction in variance-heavy plotting and testing paths by letting `binsreg_pred()` consume a precomputed covariance matrix, matching the existing R prediction helper pattern. Python `binsreg`, `binsqreg`, `binsglm`, `binstest`, and `binspwc` now reuse fitted-model covariance matrices across CI, CB, simulation, and test-statistic calculations when the underlying fitted model is unchanged.
- Reduced Python simulation memory pressure in confidence-band and testing p-value calculations by chunking `binsreg_pval()` and `binspwc_pval()` draws, matching R's streaming strategy while preserving the previous NumPy random-number draw order and numerical output.
- Reduced repeated Python simulation setup in `binstest` and `binspwc` by centralizing simulated numerator/denominator construction and reusing covariance square-root decompositions when shape, polynomial-model, and external-model tests share the same fitted model.
- Added Stata-generated low-level and command-level alignment fixtures plus R and Python checkers for knots, bin IDs, bin counts, grids, within-bin summaries, basis matrices, saved plotting data, and selected `binsregselect` outputs; regenerated Python and R numerical smoke fixtures for the intentional alignment changes.

## 2026-05-15 - Package Modernization

- Modernized release and repository infrastructure, including the Python 3.0.0 release metadata, PyPI publishing workflow safeguards, and the updated `setuptools` build requirement.
- Added top-level README function descriptions so users can quickly identify the main R, Python, and Stata commands.
- Added the Python package license file and standardized GitHub project maintenance files, including issue templates, security policy, pull request template, and CI workflow configuration.
- Removed retired development-only test scripts, stale fixtures, and obsolete generated artifacts that were no longer part of the active package validation workflow.
