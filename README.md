# Binscatter Methods

The `binsreg` package implements estimation, inference and graphical procedures using binscatter methods.

- `binsreg`: least squares binscatter regression with robust inference and plots.
- `binsqreg`: quantile binscatter regression with robust inference and plots.
- `binsglm`: generalized linear binscatter regression with robust inference and plots (R and Python).
- `binslogit`: logit binscatter estimation with robust inference and plots (Stata).
- `binsprobit`: probit binscatter estimation with robust inference and plots (Stata).
- `binstest`: binscatter-based tests for parametric specifications and shape restrictions.
- `binspwc`: binscatter-based pairwise group comparison tests.
- `binsregselect`: data-driven selection of the number of bins for binscatter estimation.

## Python Implementation

To install/update in Python type:
```
pip install binsreg
```

- Help: [PyPI repository](https://pypi.org/project/binsreg/).

- Replication: [py-script](Python/binsreg_illustration.py), [plot illustration](Python/binsreg_illustration_plot.py), [data](Python/binsreg_sim.csv).

## R Implementation

To install/update in R type:
```
install.packages('binsreg')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/binsreg/binsreg.pdf), [CRAN repository](https://cran.r-project.org/package=binsreg).

- Replication: [R-script](R/binsreg_illustration.R), [plot illustration](R/binsreg_illustration_plot.R), [data](R/binsreg_sim.csv).

## Stata Implementation

To install/update in Stata type:
```
net install binsreg, from(https://raw.githubusercontent.com/nppackages/binsreg/main/stata) replace
```

- Help: [binsreg](stata/binsreg.pdf), [binslogit](stata/binslogit.pdf), [binsprobit](stata/binsprobit.pdf), [binsqreg](stata/binsqreg.pdf), [binstest](stata/binstest.pdf), [binspwc](stata/binspwc.pdf), [binsregselect](stata/binsregselect.pdf).

- Replication files: [do-file](stata/binsreg_illustration.do), [plot illustration](stata/binsreg_illustration_plot.do), [data](stata/binsreg_simdata.dta), [speed test](stata/binsreg_speedcomparison.do).


## References

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


## Funding

This work was supported in part by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432), and [SES-2241575](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2241575).


<br><br>
