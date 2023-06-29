# BINSREG

The `binsreg` package provides Python, R and Stata implementations of binscatter methods, including partition selection, point estimation, pointwise and uniform inference methods, and graphical procedures.

This work was supported in part by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432), and [SES-2241575](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2241575).

## Website

https://nppackages.github.io/binsreg

## Queries and Requests

Please email: [binsreg@googlegroups.com](mailto:binsreg@googlegroups.com)

## Major Upgrades

This package was first released in Winter 2019, and had one major upgrade in Summer 2021.

- _Summer 2021 new features include_: (i) generalized linear models (logit, Probit, etc.) binscatter; (ii) quantile regression binscatter; (iii) new generic specification and shape restriction hypothesis testing function (now including Lp metrics); (iv) multi-group comparison of binscatter estimators; (v) generic point evaluation of covariate-adjusted binscatter; (vi) speed improvements and optimization. A complete list of upgrades is here: [UPGRADES](https://nppackages.github.io/binsreg/binsreg-0.4_upgrades.txt)


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

- Help: [R Manual](R/binsreg.pdf), [CRAN repository](https://cran.r-project.org/package=binsreg).

- Replication: [R-script](R/binsreg_illustration.R), [plot illustration](R/binsreg_illustration_plot.R), [data](R/binsreg_sim.csv).

## Stata Implementation

To install/update in Stata type:
```
net install binsreg, from(https://raw.githubusercontent.com/nppackages/binsreg/master/stata) replace
```

- Help: [binsreg](stata/binsreg.pdf), [binslogit](stata/binslogit.pdf), [binsprobit](stata/binsprobit.pdf), [binsqreg](stata/binsqreg.pdf), [binstest](stata/binstest.pdf), [binspwc](stata/binspwc.pdf), [binsregselect](stata/binsregselect.pdf).

- Replication files: [do-file](stata/binsreg_illustration.do), [plot illustration](stata/binsreg_illustration_plot.do), [data](stata/binsreg_simdata.dta), [speed test](stata/binsreg_speedcomparison.do).


## References

### Software and Implementation

- Cattaneo, Crump, Farrell and Feng (2023): [Binscatter Regressions](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_Stata.pdf).<br>
Working paper, prepared for _Stata Journal_.

### Technical and Methodological

- Cattaneo, Crump, Farrell and Feng (2023): [On Binscatter](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_AER.pdf).<br>
Working paper.<br>
[Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_AER--Supplemental.pdf)

- Cattaneo, Crump, Farrell and Feng (2023): [Nonlinear Binscatter Methods](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_NonlinearBinscatter.pdf).<br>
Working paper.<br>
[Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_NonlinearBinscatter--Supplemental.pdf)

<br><br>
