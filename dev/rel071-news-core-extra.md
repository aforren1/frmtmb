## Documentation

* `?brmsfamily` and `?frmtmb-families` ("Conway-Maxwell-Poisson
  counts") say that the `mu` of `compois()`, which
  `brmsfamily("com_poisson")` gives, is the mean of the count. brms
  2.23.0 fits `mu` as the rate of the density, so the coefficients of
  the two fits are on different scales unless `nu = 1`
  (`dev/remeasure-findings.md`, R6).
