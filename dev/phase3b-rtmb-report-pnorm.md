# DRAFT, NOT FILED: derivatives of pnorm(x, log.p = TRUE) lose accuracy as x decreases

Status: draft for the user to review. Do not file without their say so.
Address: RTMB (github.com/kaskr/RTMB).

## Summary

In RTMB 2.0 (TMB 1.9.25, R 4.6.1, Windows), the value of
`pnorm(x, log.p = TRUE)` on a tape is correct far into the lower tail,
but its derivatives are not. The relative error of the first
derivative grows with x^2, at about x^2 eps / 2: 3.2e-7 at x = -1e5,
9.1e-4 at -1e7, and 1.43 at -2e8 (485165195 against the true 2e8). The
second derivative is wrong sooner: 3170 at x = -1e5, where the true
value is -1.

Where the defect is. For an `advector` argument with `log.p = TRUE`,
RTMB's `pnorm` method calls RTMB's `distr_log_pnorm(q, lower.tail)`
(read from the method body). That routine does no derivative work: it
loops over the elements and calls `log_pnorm_both(tmp)[tail]` for each
(read from its source). `log_pnorm_both` is in neither TMB 1.9.25's
headers nor RTMB's installed headers, so it is RTMB's own atomic. Its
VALUE is right; its DERIVATIVE loses accuracy at about x^2 eps / 2.

INFERRED, not read: the error growth fits a derivative formed as the
ratio dnorm(x) / pnorm(x) by subtracting two logs,
`exp(log dnorm(x) - log pnorm(x))`. Both logs are about -x^2 / 2, so
their difference keeps an absolute error of about x^2 eps / 2, which
becomes the relative error of the derivative.

## Reproduction

```r
library(RTMB)
tp <- MakeTape(function(x) pnorm(x, log.p = TRUE), 0)
hp <- tp$jacfun()
# truth: the Mills-ratio expansion of phi(x) / Phi(x) for x << 0
d1 <- function(x) -x - 1 / x + 2 / x^3
d2 <- function(x) -1 + 1 / x^2 - 6 / x^4
for (x in c(-1e3, -1e5, -1e7, -2e8)) {
  cat(x, tp$jacobian(x), d1(x), hp$jacobian(x), d2(x), "\n")
}
```

Output on RTMB 2.0 (`dev/phase3b-log/rtmb-repro.txt`, from
`dev/phase3b-rtmb-repro.R`):

| x | d, tape | d, truth | relative error | d2, tape | d2, truth |
|---|---|---|---|---|---|
| -1e3 | 1000.001 | 1000.001 | 4.8e-11 | -0.999951 | -0.999999 |
| -1e5 | 99999.9683 | 100000 | 3.2e-7 | 3169.6 | -1 |
| -1e6 | 999980.165 | 1000000 | 2.0e-5 | 1.98e7 | -1 |
| -1e7 | 9990922.61 | 10000000 | 9.1e-4 | 9.07e10 | -1 |
| -2e8 | 485165195 | 200000000 | 1.43 | -1.38e17 | -1 |

A truth computed in plain R as
`exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))` agrees with the
tape, because it has the same cancellation. Use the expansion above,
or Rmpfr.

## Why it matters

A Laplace approximation uses the second derivative. In frmtmb.eam the
Wiener distribution function evaluates `pnorm(log.p = TRUE)` at large
negative arguments near the decision time, and the inner Hessian was
wrong there. frmtmb.eam now holds the time away from zero
(`ddm_rt_u_floor = 1e-10` in `R/wiener-rtcdf.R`), which keeps the
argument above about -1e5.

## Suggested fix

In `log_pnorm_both`'s derivative, compute the lower-tail derivative
phi(x) / Phi(x) by the Mills ratio for x below about -5: a continued
fraction, or the expansion -x - 1/x + 2/x^3 - ..., and not the
difference of two logs. The same applies to the upper tail by
symmetry.

## Withdrawn

An earlier note said the tape derivative is NaN at x = -8.3e9. That did
not reproduce on RTMB 2.0 and is not part of this report.
