# BINSREG

The `binsreg` package provides Stata and R implementations of binscatter methods, including partition selection, point estimation, pointwise and uniform inference methods, and graphical procedures.

This work was supported in part by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805) and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).

## Queries and Requests

Please email: binsreg.package@gmail.com

## Major Upgrades coming in Fall 2020:

- L2 and other metrics for hypothesis testing.

- New command/function binsreglincom for testing of linear combinations across subgroups (e.g., H0: mu1(x)=mu2(x) for all x). For now, see option by() for joint plotting of marginal confidence bands.

- New command/function binsxtreg for panel data estimation, inference and binned scatter plots. For now, in Stata use the command i. or ib(). for incorporating fixed effects (as with the regress command).

- Handling of formulas in R package.

- Recentering of binscatter estimate of \mu(x) to account for additional covariates. For now, the package sets additional covariates at zero.

- Backwards compatibility with Stata 13. For now, Stata 14 or better is needed.

## Stata Implementation

To install/update in Stata type:
```
net install binsreg, from(https://raw.githubusercontent.com/nppackages/binsreg/master/stata) replace
```

- Help: [binsreg](stata/binsreg.pdf), [binsregtest](stata/binsregtest.pdf), [binsregselect](stata/binsregselect.pdf).

- Replication files: [do-file](stata/binsreg_illustration.do), [data](stata/binsreg_simdata.dta), [speed test](stata/binsreg_speedtest.do).

## R Implementation

To install/update in R type:
```
install.packages('binsreg')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/binsreg/binsreg.pdf), [CRAN repository](https://cran.r-project.org/package=binsreg).

- Replication: [R-script](R/binsreg_R_illustration.R), [data](R/binsreg_sim.csv).


## References

### Software and Implementation

- Cattaneo, Crump, Farrell and Feng (2019): [Binscatter Regressions](references/Cattaneo-Crump-Farrell-Feng_2019_Stata.pdf), working paper.

### Technical and Methodological

- Cattaneo, Crump, Farrell and Feng (2019): [On Binscatter](references/Cattaneo-Crump-Farrell-Feng_2019_wp.pdf), working paper. [Supplemental Appendix](references/Cattaneo-Crump-Farrell-Feng_2019_wp--Supplemental.pdf).

<br><br>
