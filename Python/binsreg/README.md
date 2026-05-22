# BINSREG

The `binsreg` package provides tools for statistical analysis using the binscatter methods.

- `binsreg`: implements binscatter least squares regression with robust inference and plots, including curve estimation, pointwise confidence intervals and uniform confidence band.
- `binsqreg`: implements binscatter quantile regression with robust inference and plots, including curve estimation, pointwise confidence intervals and uniform confidence band.
- `binsglm`: implements binscatter generalized linear regression with robust inference and plots, including curve estimation, pointwise confidence intervals and uniform confidence band.
- `binstest`: implements binscatter-based hypothesis testing procedures for parametric specifications of and shape restrictions on the unknown function of interest.
- `binspwc`: implements hypothesis testing procedures for pairwise group comparison of binscatter estimators.
- `binsregselect`: implements data-driven number of bins selectors for binscatter implementation using either quantile-spaced or evenly-spaced binning/partitioning.

All the commands allow for covariate adjustment, smoothness restrictions, and clustering, among other features. See Cattaneo, Crump, Farrell and Feng (2024, 2025, 2026) for references.

Website: [https://nppackages.github.io/](https://nppackages.github.io/).

Source code: [https://github.com/nppackages/binsreg](https://github.com/nppackages/binsreg).

## Authors

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Richard K. Crump (<richard.crump@gmail.com>)

Max H. Farrell (<mhfarrell@gmail.com>)

Yingjie Feng (<fengyingjiepku@gmail.com>)

Ricardo Masini (<ricardo.masini@gmail.com>)


## Installation

To install/update use pip
```
pip install binsreg
```

# Usage
```
from binsreg import binsregselect, binsreg, binsqreg, binsglm, binstest, binspwc
```

- Replication: [binsreg illustration](https://github.com/nppackages/binsreg/blob/main/Python/binsreg_illustration.py), [plot illustration](https://github.com/nppackages/binsreg/blob/main/Python/binsreg_illustration_plot.py), [simulated data](https://github.com/nppackages/binsreg/blob/main/Python/binsreg_sim.csv).


## Dependencies

- numpy
- pandas
- scipy
- statsmodels
- plotnine

## References

For overviews and introductions, see [NP Packages website](https://nppackages.github.io/).

### Software and Implementation

- Cattaneo, Crump, Farrell and Feng (2025): [Binscatter Regressions](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2025_Stata.pdf).<br>
_Stata Journal_ 25(1): 3-50.

### Technical and Methodological

- Cattaneo, Crump, Farrell and Feng (2024): [On Binscatter](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_AER.pdf).<br>
_American Economic Review_ 114(5): 1488-1514.<br>
[Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_AER--Supplemental.pdf)

- Cattaneo, Crump, Farrell and Feng (2026): [Nonlinear Binscatter Methods](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2026_RESTAT.pdf).<br>
_Review of Economics and Statistics_, revise and resubmit.<br>
[Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2026_RESTAT--Supplemental.pdf)
