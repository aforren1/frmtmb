# Generics frmtmb re-exports from the package that owns them.

nlme is Recommended and ships with R; generics declares only base
packages, [`library(generics)`](https://generics.r-lib.org) costing
0.000 s at the minimum of 12 replicates against 0.860 s for lme4
(`dev/generics-loadcost.R`). So these four are free. Three of them need
nothing further; `refit` also joins the run-time table, because lme4
defines a rival `refit` that this import cannot reach.
