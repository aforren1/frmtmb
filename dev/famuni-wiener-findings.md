# Lane famuni: frmtmb.eam `wiener()` against brms 2.23.0 `wiener()`

Helper audit, 2026-10-13. Base: rellib-r8 (frmtmb 0.70.0, eam 0.11.3).
Every eam run below used `pkgload::load_all()` of the worktree's eam
over rellib-r8's core. No package was installed. No git operation.
Logs are in `dev/famuni-out/` (gitignored); the numbers here are
copied from them.

## 1. What brms 2.23.0 does (run, not recalled)

Script `dev/famuni-wiener-brms.R`, log `wiener-brms.txt`, seed
20261013, n = 80.

- `formals(brms::wiener)`: `link = "identity"`, `link_bs = "log"`,
  `link_ndt = "log"`, `link_bias = "logit"`.
- `.family_wiener()`: links `identity, log, softplus, squareplus`;
  dpars `mu, bs, ndt, bias`; `ad = weights, subset, dec, index`;
  include `fun_wiener_diffusion.stan`.
- Links accepted, by calling `wiener()` with each: `link`, `link_bs`
  and `link_ndt` take `identity, log, softplus, squareplus`;
  `link_bias` takes `identity, logit`. Every other name tried (`logit,
  inverse, probit, cloglog, sqrt, tan_half, logm1`) is refused.
- Density: rstan 2.32.7 compiles against Stan 2.39.0, but brms emits
  only the 4-parameter `wiener_lpdf(y | alpha, tau, beta, delta)`
  through `wiener_diffusion_lpdf`, with dec == 0 reflected as
  `(1 - beta, -delta)`. No `sv`, `sw` or `st0`; the 5- and 7-parameter
  Stan densities are not used.
- `cens()`, `trunc()` and `mi()` are refused for wiener; `weights()`
  works; `mixture(wiener, wiener)` works. A missing `dec()` is refused.
  A factor or `"lower"/"upper"` character `dec()` maps to 0/1 by level,
  as eam's does.
- Scalar spelling: parameter block `real<lower=0> bs`,
  `real<lower=0,upper=min_Y> ndt`, `real<lower=0,upper=1> bias`.
  Formula spelling: `bs = exp(bs)`, `ndt = exp(ndt)`,
  `bias = inv_logit(bias)`, no bound on ndt.
- `default_prior()`:

  | slot | scalar spelling | formula spelling (`~ 1`) |
  |---|---|---|
  | mu Intercept | `student_t(3, 0.5, 2.5)` | same |
  | bs | `gamma(1, 1)` | Intercept `normal(-0.6, 1.3)` |
  | ndt | `uniform(0, min_Y)` | Intercept flat |
  | bias | `beta(1, 1)` | Intercept `logistic(0, 1)` |

  The mu location is the rounded median of the response
  (`student_t(3, 0, 2.5)` under softplus).
- `log_lik_wiener()` calls `brms::dwiener()` (RWiener backend) at the
  observed `dec`. `posterior_predict_wiener()` draws `(q, resp)` with
  `RWiener::rwiener()` and returns `q` only, or `q` signed by the DRAWN
  `resp` when `negative_rt = TRUE`. `posterior_epred_wiener()` is
  `ndt - bias / mu + bs / mu * (exp(-2 mu bias) - 1) /
  (exp(-2 mu bs) - 1)`.

## 2. Default priors on the sampling route

Scripts `dev/famuni-wiener-eamprobe.R` (log `wiener-eamprobe.txt`) and
`dev/famuni-wiener-sampleprior.R` (log `wiener-sampleprior.txt`).

`frmtmb::default_prior(..., route = "sample")` and frmtmb.sample's
`default_priors_for()` give only the mu intercept,
`student_t(3, 0.5, 2.5)`, which matches brms. `bs`, `ndt` and `bias`
are flat in both spellings, and `default_prior_notes()` returns NONE
for either model, so `frm_sample()`'s disclosure passes the three in
silence. That contradicts the notes' own purpose ("a family that gets
no defaults never passes in silence").

What would close it, NOT done (frmtmb.sample is out of my scope and it
is not one line):

- bias, both spellings: brms's `beta(1, 1)` on the natural scale and
  `logistic(0, 1)` on the logit intercept are the same density (an
  identity: the logit of a uniform is standard logistic).
- ndt, scalar spelling: `uniform(0, min_Y)` is exactly
  `logistic(0, 1)` on eam's default scaled-logit intercept, because
  that link's bound is `min(y)` (same identity). Under `ndt_group()`
  or `link_ndt` there is no such equivalence; brms has flat there.
- bs: formula spelling `normal(-0.6, 1.3)` on the log intercept is
  expressible; scalar `gamma(1, 1)` on the natural scale is in the class
  of gamma defaults the notes already say `set_prior()` cannot carry.
- At minimum, `default_prior_notes()` should name `bs`, `ndt` and
  `bias` for a wiener model. These rules are generic in brms
  (`def_dprior()` keys them on the dpar name: `bs`, `bias`, `zi`, `hu`,
  `zoi`, `coi`, `quantile`, `xi`, `alpha`, ...), so frmtmb.sample likely
  has the same silent gap for every such dpar, not only wiener's. I did
  not measure the other families.

## 3. Density agreement

Script `dev/famuni-wiener-density.R`, logs `wiener-density.txt` and
`wiener-density2.txt` (`dev/famuni-wiener-density2.R`). Grid: decision
time {1e-4, 1e-3, 0.01, 0.05, 0.1, 0.3, 1, 3, 10} x bs {0.5, 1, 2, 4}
x bias {0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999} x drift {-20, -10,
-5, -1, 0, 1, 5, 10, 20} x both boundaries = 4536 points, ndt 0.1.
Reference: 1200-bit Rmpfr, small-time series for u < 1 and large-time
series otherwise, 60 terms each; finite at 4536 of 4536, range
-79905.7 to 7.97. brms's Stan density: brms's own generated functions
block, compiled with `rstan::expose_stan_functions()` (32.4 s).

Max abs log-density error against the reference:

| implementation | all 4536 | where ref > -700 (3964) | non-finite |
|---|---|---|---|
| frmtmb.eam `ddm_lpdf_both()` | 1.455e-11 | 1.137e-13 | 0 |
| brms Stan `wiener_diffusion_lpdf` | 1.829e-05 (finite rows) | 1.829e-05 | 450 |
| brms `dwiener()` (RWiener) | 4.547e-13 (finite rows) | 1.421e-13 | 450 |
| RWiener `dwiener()` | 4.547e-13 (finite rows) | 1.421e-13 | 450 |
| WienR `dWDM(precision = 12)` | 2.793e-01 | 2.793e-01 | 0 |
| rtdists `ddiffusion()` | 1.491e+02 | 1.491e+02 | 562 |

eam's worst relative error is 1.047e-14 (scaled by max(1, |ref|)).
eam against Stan directly, where both finite: 1.829e-05 (n = 4086).

- Stan's density is `-Inf` on 450 points, 448 of them with the true
  log density below -700 (it forms the density, then takes the log,
  so it underflows) and 2 at -699.51 (t = 0.01, bs = 4, bias
  0.99 or 0.01, |v| = 20). Where finite its error is above 1e-8 on 54
  of 4086 points; the ten worst are all at bias 0.999/0.001 and
  normalized time 0.625, each off by 1.829e-05. This is brms's (Stan's)
  inaccuracy, not eam's: eam, RWiener and the Rmpfr series agree
  there.
- WienR's and rtdists's large errors are at bias 0.999/0.001 with
  normalized time 0.1875 (the six worst WienR rows); they
  are references with their own defects at extreme bias, not evidence
  against eam.
- At a decision time of 0 or -0.001, eam returns `-Inf` and Stan throws
  "Random variable ... must be greater than nondecision time".

## 4. Fit agreement against brms

Script `dev/famuni-wiener-fit2.R`, log `wiener-fit2.txt`. Data: n = 400,
`RWiener::rwiener()` per row, seed 20261013, mu = 0.6 + 0.8 x,
bs = exp(0.25 + 0.2 x), ndt 0.3, bias 0.5; min rt 0.3316, 75.25
percent upper. eam ML `bf(rt | dec(d) ~ x, bs ~ x)`, default links,
0.24 s, convergence 0. brms: same formula, default priors (note the
scalar spelling's `ndt` bound `(0, min_Y)` and `uniform(0, min_Y)`
prior, `bs` gamma(1, 1), `bias` beta(1, 1), mu intercept
`student_t(3, 0.6, 2.5)`), rstan, 4 chains x 1000 post-warmup, seed
20261013, 88.9 s; max Rhat 1.0014, min bulk ESS 1733, 0 divergences.

| par | truth | ML | SE | post. median | post. SD | z shift | SE/SD |
|---|---|---|---|---|---|---|---|
| mu Intercept | 0.60 | 0.8712 | 0.10699 | 0.8777 | 0.10660 | -0.0604 | 1.0036 |
| mu x | 0.80 | 0.6814 | 0.09003 | 0.6800 | 0.08932 | 0.0152 | 1.0080 |
| bs Intercept | 0.25 | 0.2683 | 0.02709 | 0.2696 | 0.02691 | -0.0461 | 1.0065 |
| bs x | 0.20 | 0.2114 | 0.02253 | 0.2111 | 0.02276 | 0.0140 | 0.9896 |
| ndt | 0.30 | 0.3001 | 0.005178 | 0.2995 | 0.005309 | 0.1174 | 0.9754 |
| bias | 0.50 | 0.4922 | 0.01806 | 0.4911 | 0.01793 | 0.0570 | 1.0073 |

max |z shift| 0.117. ndt and bias are on the natural scale (eam's
`summary()$spec_pars`, delta-method SE).

Second brms fit, `dev/famuni-wiener-fit3.R`, log `wiener-fit3.txt`:
the same data with `ndt ~ 1` (brms's log link, flat intercept prior)
against eam `wiener(link_ndt = "log")`; brms 154.8 s, max Rhat 1.0011,
min bulk ESS 2075, 0 divergences. eam's log-likelihood equals the
default-link fit's to -3.7e-11.

| par | ML | SE | post. median | post. SD | z shift | SE/SD |
|---|---|---|---|---|---|---|
| Intercept | 0.8712 | 0.10699 | 0.8734 | 0.10616 | -0.0202 | 1.0078 |
| x | 0.6814 | 0.09003 | 0.6821 | 0.09067 | -0.0083 | 0.9929 |
| bs Intercept | 0.2683 | 0.02709 | 0.2715 | 0.02794 | -0.1116 | 0.9694 |
| bs x | 0.2114 | 0.02253 | 0.2117 | 0.02231 | -0.0130 | 1.0096 |
| ndt Intercept (log) | -1.2036 | 0.01725 | -1.2062 | 0.01802 | 0.1484 | 0.9572 |
| bias | 0.4922 | 0.01806 | 0.4921 | 0.01756 | 0.0033 | 1.0284 |

max |z shift| 0.148.

## 5. log_lik(), posterior_predict(), posterior_epred()

Same script and log as section 4.

- log_lik: eam's density at 50 brms draws against
  `brms::log_lik(draw_ids = )`: max abs 1.088e-14. (brms returns the
  draws of `draw_ids` in ascending order; an unsorted id vector gave a
  spurious 2.9 difference, `dev/famuni-wiener-lldbg.R`.) At the ML,
  eam's density summed is -118.6265715114, equal to `logLik()`, and
  against `brms::dwiener()` at the ML parameters: 3.775e-15.
  `frmtmb::log_lik()` on an ML fit is refused by design (per draw).
  Same semantics: the log density at the observed boundary.
- posterior_predict: brms's default returns no negative values;
  with `negative_rt = TRUE` 25.12 percent of draws are lower-boundary
  against 24.75 percent observed, and rows observed at the UPPER
  boundary get 21.06 percent lower-boundary draws (rows observed lower:
  37.46 percent). So brms draws the boundary marginally and drops it.
  eam's `simulate()` conditions on each row's observed boundary (no
  negative values). Mean simulated rt on upper rows: eam 0.69139, brms
  0.69056; on lower rows: eam 0.61231, brms 0.61739. Core's
  `predict(negative_rt = TRUE)` and frmtmb.sample's
  `posterior_predict(negative_rt = TRUE)` refuse. Deliberate divergence,
  now documented in `?wiener`.
- posterior_epred: brms's mean is marginal over both boundaries (upper
  rows 0.72736, lower 0.62517); eam's `fitted()` is conditional (0.69127,
  0.61360). DEFECT IN BRMS 2.23.0, measured by calling
  `brms:::posterior_epred_wiener()` on a 48-cell grid against the exact
  unconditional mean (eam's conditional means weighted by the boundary
  probabilities): relative error at bs = 1 at most 4.07e-16; at
  bs != 1 up to 0.637 (bs 2.5, mu 1.5, bias 0.7); with the absolute
  start `z = bias * bs` substituted the formula agrees to 1.04e-15. A
  20000-draw RWiener Monte Carlo agrees with the exact value on four
  cells (0.44867 vs 0.44902; 0.79461 vs 0.79217; 1.68820 vs 1.68650;
  0.79248 vs 0.79217). brms's formula is also NaN at mu = 0. Not
  reported upstream by me.

## 6. Changes made (frmtmb.eam only)

Files: `extensions/frmtmb.eam/R/wiener-family.R`,
`extensions/frmtmb.eam/R/ddm-shared.R`,
`extensions/frmtmb.eam/tests/testthat/test-wiener-links.R` (new),
`extensions/frmtmb.eam/man/wiener.Rd` (roxygenised in the worktree),
`extensions/frmtmb.eam/NEWS.md`.

1. `wiener()` arguments now start `link, link_bs, link_ndt, link_bias`
   (brms's order), then `max_ndt, variability, nodes,
   allow_unreachable, contaminant, contaminant_range`. BREAKING for a
   positional `max_ndt`; no positional call exists in the repository
   (grep of `wiener([0-9.`).
2. `link` is refused outside brms's set. On rellib-r8,
   `wiener(link = "logit")`, `"probit"` and `"sqrt"` fitted
   (logLik -46.32, -46.32, -54.56 on the probe data; `"inverse"`
   failed with a NaN gradient), `wiener-eamprobe.txt`. An unquoted
   name in the set reads as its name (`wiener(log)`), as brms does; a
   custom link list still passes.
3. `link_bs` and `link_bias` added with brms's sets. An identity link
   gets a support guard, `0 / (x + |x|)`, added to the density and
   distribution-function slots, so `bs <= 0`, `bias` outside (0, 1) and
   `ndt <= 0` give NaN. Without it the density at bias 1.2 is finite
   (-2.64, -3.22, -3.14 on three rows) because `ddm_pos()` folds a
   negative series sum (`wiener-links.txt`).
4. `link_ndt` added: NULL (default) is the bounded link as before;
   `"log"`, `"identity"`, `"softplus"`, `"squareplus"` are brms's. The
   family sets `ndt_link_fixed`, and `ddm_ndt_install()` then leaves the
   link alone; the bound record stays (`ndt_time()` works). Refused with
   `max_ndt`, `ndt_group()`, `variability`, `contaminant = TRUE` and
   `allow_unreachable = TRUE`. Why the last three: they score a row past
   ndt at a floored decision time, and with no bound every row can sit
   there at once. Measured before the refusal existed:
   `mixture(wiener(link_ndt = "log", allow_unreachable = TRUE),
   lognormal())` on data with no fast guesses ran to ndt 1.8e+308,
   bs 7.05e-05, mu 1169 and logLik +7713.9 on 400 rows (the Wiener
   component's log density 19.2 on every row), `wiener-mix2.txt`; the
   documented `max_ndt = 0.4` recipe on the same data stops with a NaN
   gradient, and on data with 20 fast guesses the two recipes reach the
   same logLik, -64.022 (`wiener-mix.txt`).
5. `?wiener`: new sections "Links, and the one difference from brms"
   and "Post-fit quantities against brms", parameter docs, Grasman et
   al. (2009) reference. Rendered with `tools::Rd2txt` and read;
   `tools::checkRd()` clean.

Reparameterization invariance on a scalar dpar (`wiener-links.txt`,
seed 20261013, n 400): every arm converged; |dlogLik| at most 1.334e-09
(bs squareplus), max |dEst|/SE at most 3.41e-05, SE ratio
0.999359 to 1.000003, for bs softplus/squareplus/identity, bias
identity, ndt log/softplus/squareplus/identity and all-identity. Right
censoring at 1.2 s (33 rows): default and `link_ndt = "log"` differ by
-5.7e-13 in logLik.

### Optimizer under link_ndt = "log" (`dev/famuni-wiener-ndtrep.R`,
summary `dev/famuni-wiener-ndtrep-sum.R`, log `wiener-ndtrep-sum.txt`)

- Scalar arm, seeds 20261101 to 20261200, n 30 and 400, mu = 0.8 +
  0.5 x, bs 1.3, ndt 0.3: 400 fits, convergence 0 on 100 of 100 per
  cell for both links, 0 warnings. log minus default logLik: -1e-10 to
  1e-10 at n 30, -6e-10 to 4e-10 at n 400; max |dndt|/SE 5.7e-06 and
  1.56e-05. Median time 0.03 s vs 0.03 s and 0.21 s vs 0.23 s.
- Hierarchical arm, seeds 20261101 to 20261120, 20 subjects x 50
  trials, log ndt_s ~ N(log 0.25, 0.12), `ndt ~ 1 + (1 | s)`:

  | arm | conv 0 | fits warning | median RMSE ms (range) | median sd |
  |---|---|---|---|---|
  | global bound (default) | 12/20 | 8 | 24.43 (16.25 to 35.98) | 1.75 (logit scale) |
  | `ndt_group(s)` | 20/20 | 0 | 16.38 (10.70 to 20.42) | 3.0e-05 |
  | `link_ndt = "log"` | 20/20 | 0 | 13.81 (8.74 to 17.51) | 0.115 (truth 0.12) |

  Paired: the log link's RMSE was below ndt_group's on 20 of 20 and
  below the global bound's on 20 of 20. Its sd ranged 0.0887 to 0.162.
  Each of the eight global-bound warnings begins "Optimizer did not
  report convergence: fa" (the one printed in full is "false
  convergence (8)").

### Tests run (one file per process, `dev/famuni-wiener-runtest.R`)

- New `test-wiener-links.R` against rellib-r8's installed eam 0.11.3
  (`test-links-rellib.txt`): `tests=8 expectations=7 failed=7 error=7`.
  Behavioral failures seen: formals order (2 failures); `wiener(link =
  lk)` did not throw for logit, probit, inverse, sqrt, cloglog (5);
  `wiener(log)` errored as a bad `max_ndt`. The rest are missing-symbol
  or unused-argument errors (the weak form) for the new arguments.
- Against the worktree (`test-links-worktree.txt`): `tests=8
  expectations=55 failed=0 error=0 skipped=0 warning=0`.
- Whole eam suite against the worktree, 30 files
  (`dev/famuni-wiener-suite.sh`, `suite-worktree.txt`): 0 failed, 0
  error, 0 warning in all 30; the only skips are `test-scale.R`'s 3
  gated tests (`skip_unless_scale()`, pre-existing). Lib line: frmtmb
  0.70.0 from rellib-r8 (the lane owner's uncommitted core edits were
  NOT loaded by these runs).

## 7. Deliberate divergences kept

- `ndt` default link: bounded (scaled logit onto `(0, min(y))`, or a
  logit on a fraction under `ndt_group()`), against brms's log link.
  brms's scalar spelling bounds ndt in its parameter block; frmtmb has
  no parameter block, so the link carries the bound.
- `fitted()` / `simulate()` condition on the observed boundary; brms's
  `posterior_epred()` / `posterior_predict()` are marginal over it.
- eam supports `cens()`, `trunc()`, variability and a contaminant,
  which brms refuses or lacks.

## 8. Recommendation on link_ndt = "log"

Offer it, as implemented, and keep the bounded link as the default.

Reasons for offering it: it is brms's parameterization, so `ndt ~ ...`
coefficients and priors port one to one (section 4, second fit, max
|z| 0.148); on a hierarchical ndt it is the best of the three
constraints measured (20 of 20 converged, lowest per-subject error on
20 of 20 data sets, sd recovered at 0.115 against 0.12), because each
row's own density bounds that row, which is per subject without
`ndt_group()`'s floor artifact (sd driven to zero); and it costs the
optimizer nothing on the plain model (400 of 400 converged, same
maxima).

What the optimizer needs, and has: a feasible start (half the
fastest response, already the default `init`); a density that is not
finite past the wall (eam's is exactly `-Inf` at a decision time at or
below 0, so the objective is `+Inf` and `nlminb` rejects the step, as
core's `shifted_lognormal` relies on NaN); and no floored density
anywhere, which is why variability, the contaminant and
`allow_unreachable` are refused with it.

Reasons to keep the bounded default: the bound is structural, needs no
rejected steps, works inside `mixture()` and with the contaminant and
variability, where the unbounded link is unsafe (section 6, item 4).

## 9. Defects found, not fixed

1. brms 2.23.0 `posterior_epred()` for wiener uses `bias` where the
   formula needs `bias * bs` (section 5). Upstream bug.
2. Stan 2.39's 4-parameter `wiener_lpdf` (as brms calls it) returns
   `-Inf` below about -700 and is off by 1.83e-05 at bias 0.999/0.001
   (section 3). Upstream.
3. frmtmb.sample: no defaults for wiener's `bs`, `ndt`, `bias`, and
   `default_prior_notes()` does not name them (section 2). Not a
   one-line fix; for the lane owner.
4. Under an unbounded `ndt` link, any density that floors the decision
   time admits an unbounded likelihood (section 6, item 4). Guarded by
   refusal here; the same mechanism under `max_ndt` above min(y) with
   `allow_unreachable` (the documented mixture recipe) is bounded only
   because ndt is, and on data with no contaminant it ends in a NaN
   gradient rather than a degenerate fit (`wiener-mix.txt`).

## 10. Not done

- No `R CMD check --as-cran` of eam (the brief forbids installs; the
  lane owner runs it). No OpenBLAS pass: the change adds no rounding-
  decided branch, and the invariance tests compare against ratios to
  the run's own standard errors and log-likelihood.
- No brms-suite ledger rows checked; `test-brms-parity.R` still passes.
- Bump: minor for frmtmb.eam (new arguments, BREAKING argument order).
  Floor: frmtmb 0.70.0 (no new core API used).
