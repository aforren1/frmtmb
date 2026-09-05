# brms likelihood-identity tests: plan

Status: partly implemented on branch `wt-brmslp`, 2026-09-05. The
harness, the translator with its round-trip unit tests, the CI
workflow and 19 tests covering plan rows 1, 2, 3, 5, 7, 12, 13, 14,
15, 16, 17, 20 and 21 plus three check C rows are in place and green.
Rows 4, 6, 8, 9, 10, 11, 18 and 19 are not written. One real
divergence between the packages was found and is recorded below, and
the exemption list it creates has exactly that one entry. Results,
deviations, findings and the remaining work are under
"Implementation log" at the end of this document.

## Claim under test

For every model shape brms and frmtmb both spell, frmtmb's objective is
the same function of the parameters as the Stan program brms generates
with flat priors. Estimates are never compared to posterior summaries;
the estimand differs and the tolerance would have to hide that. The
density is compared at a point, exactly.

The existing tier in `test-brms-agreement.R` already does this twice
(the flat-prior `optimizing` test and the `mo()` `log_prob` test). This
plan generalizes those two and retires the posterior-mean comparison.

## The three checks

Let `fit` be the frmtmb fit and `sf` an uncompiled-data Stan fit object
(`rstan::sampling(mod, data = sdat, chains = 0)`), built from
`brms::make_stancode()` / `make_standata()` with every prior set to
`prior("")`.

**A. Value identity.** `log_prob(sf, upars, adjust_transform = FALSE)`
at frmtmb's estimates equals `logLik(fit)` plus a known constant.
Tolerance `1e-6 * max(1, |logLik|)`. The constant is zero unless brms
keeps an implicit prior the empty string cannot remove; the two known
cases are the Dirichlet on a `mo()` simplex (`lgamma(D)`) and the
Dirichlet on mixture weights. Any other nonzero constant is a finding.

**B. Stationarity.** `grad_log_prob(sf, upars, adjust_transform =
FALSE)` at the same point has max absolute entry below `1e-3`. This
proves frmtmb's estimate is brms's optimum without running Stan's
optimizer, and it is the check that catches a parameter-mapping error,
because a wrong map lands off the optimum.

**C. Joint identity, random-effect models.** brms has no marginal
likelihood, so A and B run on the joint density instead. The point is
frmtmb's outer estimates plus its conditional modes `b_hat`, mapped to
brms's non-centered `z` (below). Check A compares `log_prob` to minus
RTMB's inner objective, `obj$env$f(c(theta_hat, b_hat))`, plus the
map's log-Jacobian. Check B requires zero gradient with respect to the
`z` block only; the outer gradient is not zero at the joint point and
is not asserted. The Laplace step itself is validated in-house by the
tmbstan tier and is out of scope here.

## Parameter map

One translator, `stan_pars_from_fit(fit, sdat)`, returning the
constrained parameter list `rstan::unconstrain_pars()` accepts. Rules:

- `Intercept`: brms centers `X`, so `Intercept = b0 + colMeans(X) %*%
  b`, per dpar, unless `center = FALSE` was set. `Xmo`, `Xs`, `Xcs`
  are not centered.
- `b_<dpar>`, `bsp`, `bs`, `bcs`, `simo_*`: straight copies; simplexes
  from frmtmb's `zeta` softmax.
- `sd_<g>`, `L_<g>`: from `VarCorr(fit)`; `L` is the Cholesky factor of
  the correlation matrix, lower triangular.
- `z_<g>`: brms builds `r = (diag(sd) L z)^T` per group, so `z =
  solve(diag(sd) %*% L, t(r))`, with `r` as levels-by-coefficients in
  brms's level order (`J_<g>` in standata gives it). frmtmb's `b`
  order follows its `Zt` column order, which is NOT brms's; take it
  from `ranef(fit)`, not from the raw vector. The log-Jacobian added
  to check C is `n_levels * (sum(log(sd)) + sum(log(diag(L))))` per
  group.
- `sigma`, `shape`, `nu`, `phi`, `zi`, `hu`: natural scale, from
  frmtmb's link-scale estimate through the inverse link.
- `rescor` `Lrescor`, autocor `ar`/`ma`/`cosy`, `cortime`, `sdgp`,
  `lscale`, `zgp`: same pattern; each one added when its model shape
  joins the matrix, with a unit test on the map alone (round trip
  through `rstan::constrain_pars`).

## Model matrix

Reference-free shapes first, because they are the ones with no other
reference. One row each, small data, fixed seed.

1. `y ~ x + z, sigma ~ x` (exists; A, B)
2. `y ~ mo(inc) + z` (exists; A, B)
3. `y ~ mo(inc) * z` and `mo(inc):z`
4. `y | mi() ~ x`, `x | mi(sdx)` measurement error
5. nonlinear `bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE)`
6. nonlinear with random effects `a ~ 1 + (1 | g)` (C)
7. `(1 | ID | g)` across `mu` and `sigma` (C)
8. `mvbf(y1 ~ x, y2 ~ x) + set_rescor(TRUE)`
9. multivariate with `(1 | ID | g)` in both responses (C)
10. `gp(x)` exact and `gp(x, k = 10)` HSGP
11. `s(x)` and `t2(x, z)` (C: the wiggly part is a random effect)
12. ordinal: cumulative, sratio, cratio, acat, with `cs(x)`
13. categorical and multinomial
14. `y | cens(c) ~ x`, `y | trunc(lb = 0) ~ x`, `y | se(s) ~ 1`
15. `y | trials(n) ~ x` binomial, beta-binomial
16. zero-inflated and hurdle families with `zi ~ x`
17. `mixture(gaussian, gaussian)` with `theta1 ~ x`
18. `ar(time, gr = g, p = 1)`, `cosy`, `unstr` (gaussian)
19. `car(M, gr = g)` each type
20. `weights` and `subset` addition terms
21. every remaining family in the shared roster, as `y ~ x` (A, B only)

Extend with one `(1 | g)` variant of rows 3, 4, 8, 10, 12, 16 for C.

## Harness

- New file `tests/testthat/test-brms-likelihood.R`, gated by the
  existing `skip_unless_brms_fit()`.
- One helper `brms_lp_check(bform, family, data, fit, joint = FALSE)`
  that builds the flat prior set from `get_prior()` (every row set to
  `""`), compiles, translates, and runs A and B (and C when `joint`).
  A test is one call per row of the matrix.
- Compiled models cached on disk keyed by the hash of the Stan code
  plus the rstan version, directory from `FRMTMB_STAN_CACHE`. The code
  depends on the formula and family, not the data, so one compile
  serves every data variant of a row.
- CI: a scheduled weekly workflow with rstan installed, the cache
  restored by `actions/cache` on the same key, and a manual trigger.
  Not on push.
- Failures are failures. No known-divergence list for numeric
  mismatches. The only admissible exemption is a documented design
  choice, recorded in the test with the reason, and today that list is
  empty for flat priors.

## Follow-on, priors

**Done, 2026-09-05, on branch `wt-brms-priors`.** Checks A and B were
repeated with brms's own `get_prior()` defaults on seven shapes
(random intercept, distributional, correlated slope, ordinal,
nonlinear, mixture, and a fixed-effect `sleepstudy` added to isolate
one question), each with three Stan programs and three frmtmb fits.
The measurement is in `dev/brms-priors-findings.md` and is pinned by
`tests/testthat/test-brms-priors.R` with
`tests/testthat/helper-brms-priors.R`. Nothing under `R/` changed; the
maintainer decides what does.

**The prediction this section made was wrong, and the reason is worth
keeping.** It said the link-scale versus natural-scale placement of
`sd` priors would show as a value mismatch. It does not: `R/priors.R`
`prior_logdens()` already evaluates an `"sd"`-scaled entry at
`exp(theta)` and adds `theta` as the log-Jacobian, which is brms's
placement under `adjust_transform = TRUE`, and `man/set_prior.Rd`
already documents it. Measured on `sleepstudy`, frmtmb's `sd` entry
equals brms's half-t statement plus `log(sd)` minus `log(2)` to eight
figures. The paragraph was out of date, and only running it found that
out.

What the measurement turned up instead, in descending order of size:

1. **brms's `get_prior()` defaults reach `frm(prior = )` as nothing at
   all.** Every row carries `source == "default"` and `as_priorlist()`
   drops exactly those, so the fit is unpenalized and `fit$prior` is
   `NULL`. True on all seven shapes.
2. **`class = "Intercept"` is on the intercept at zero; brms's is on
   the intercept at the mean of the predictors.** Same density, same
   scale, different argument. On `Reaction ~ Days` this biases the
   SLOPE by 0.069 standard errors, where brms's identical prior string
   moves it by 3.5e-05. It is the only one of these that touches a
   regression coefficient.
3. **A distributional parameter without a linear predictor is the real
   placement difference**, and frmtmb refuses the row rather than
   mistranslating it. The nearest spelling its own error message
   suggests sits on the link scale and captures between -0.5% and 35%
   of what brms's prior intends, because brms's default scales are
   derived from the spread of the response and are natural-scale
   quantities. The natural placement exists internally
   (`natural = TRUE`, set only by `frm_sample()`), is unreachable from
   `frm(prior = )`, and reproduces brms's mode to nine figures.
4. **Three prior kinds have no spelling at all**: the ordinal
   threshold prior (accepted class, no target), brms's mixture
   `theta`, and any density outside
   normal/student_t/cauchy/exponential/lkj (`logistic`, `dirichlet`,
   `gamma`).

Check B answered its question cleanly on the nonlinear row, which has
no confounder: with the prior on the natural scale,
`log_prob(..., adjust_transform = TRUE)` minus frmtmb's penalized
objective is `log(2)` to eight figures and that gradient is 3.5e-05,
so frmtmb maximizes exactly the density brms samples. With the link
spelling neither gradient vanishes.

One methodological note for anyone extending this: the mode-distance
question is not answerable on a random-effect shape. Stan's optimum
over the joint density runs the standard deviations away (on
`(Days | Subject)` it gives `sd_1 = (110.2, 66.1)` against frmtmb's
marginal `(24.91, 5.99)`), which is the same reason this plan compares
densities at a point rather than optimizers.

## Done when

Every row of the matrix passes A and B; every row marked C passes C;
the sleepstudy posterior-mean test is deleted; the README's validation
paragraph says "the log-likelihood equals the Stan program's log
density at the estimate" instead of "verifies that our estimates equal
the mode".

## Implementation log

### Toolchain finding: rstan cannot compile on this machine

Before any row could run, every `rstan::stan_model()` call failed with
`Compilation ERROR, function(s)/method(s) not created!`, followed by
`Error in sink(type = "output"): invalid connection` when rstan tried to
unwind its output sink. The second error is a symptom of the first and
hides it, which is why the cause took a manual compile to find.

The installed pair is rstan 2.32.7 with StanHeaders 2.39.1 on R 4.6.1.
rstan asks for the C++17 toolchain, which is where the leading
`-std=gnu++17` comes from: it is rstan's own `CXX17STD`, not an R
default, since R 4.6.1 defaults to `gnu++20`. rstan then emits its
generated model with `// [[Rcpp::plugins(cpp14)]]`, and the inline
plugin turns that into `-std=c++1y` inside `PKG_CXXFLAGS`, which
Makeconf places AFTER `CXX17STD`. The last `-std` wins, so the unit
compiles as C++14, and StanHeaders 2.39.1 requires C++17: the first
error is `'void_t' is not a member of 'std'`. Extracting the generated
translation unit and compiling it by hand with only `-std=gnu++17`
succeeds with zero errors, which is the measurement that identifies the
flag and not the headers as the fault.

This is an environment fault, not a frmtmb fault, and it is outside this
lane's surface, so nothing in the package works around it. Local runs
set `R_MAKEVARS_USER` to a file whose `CXX17FLAGS` ends in
`-std=gnu++17`; that flag lands last and restores C++17. The alternative
repair is to install the StanHeaders 2.32.x that rstan 2.32.7 was built
against.

**CI is affected too**, which an earlier draft of this section wrongly
left open. RSPM serves rstan 2.32.7 with StanHeaders 2.39.1 on both
noble and jammy, the same two version numbers, and none of the three
ingredients is platform-guarded: `rstan:::rstanplugin()` returns the
`cpp14` plugin from BOTH arms of its Windows branch,
`rstan:::cxxfunctionplus()` sets `USE_CXX17` unconditionally, and the
`PKG_CXXFLAGS` after `CXX17STD` ordering is generic Makeconf. The
workflow therefore writes `~/.R/Makevars` with
`CXX17FLAGS = -O2 -Wall -std=gnu++17` before running the tests.
`CXX17STD` is not a usable lever: it lands before `PKG_CXXFLAGS` and so
loses to the plugin. The step is a no-op once the pair is fixed
upstream. This remains a prediction rather than an executed result,
because no Linux runner was available here.

Measured cold compile of a one-parameter Stan model after the repair:
97 s.

### Deviations from the plan

1. `stan_pars_from_fit(fit, sdat)` gains two arguments and becomes
   `stan_pars_from_fit(fit, sdat, code, rtab)`. The set of parameters to
   translate is only in the Stan code, and standata carries neither the
   level labels that `J_<id>` indexes nor the coefficient order inside a
   group. Driving the translator off the declared parameter block, and
   not off frmtmb's parameter vector, is what makes an untranslated
   shape an error that names the parameter instead of a silently
   omitted term.

2. `rtab` is `brms::brm(..., empty = TRUE)$ranef`, brms's own
   group-level table. It costs 0.34 s, compiles nothing, and is the only
   public source of the two facts standata omits.

3. The workflow runs on push and pull_request as well as on the weekly
   schedule and on dispatch, per the maintainer. The plan said "Not on
   push". A path filter on `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE`
   and the workflow file keeps it off pushes that cannot affect it.

### The parameter map as implemented

Each rule is exercised by a round trip through `rstan::constrain_pars`
in the first two tests of the file, and then by checks A and B on every
row that reaches it.

- `Intercept`, `Intercept_<dpar>`: `b0 + sum(colMeans(X)[-1] * b)`.
  brms passes the uncentered `X` with its intercept column and centers
  inside the Stan program, so `colMeans(X)[-1]` IS the `means_X` the
  generated quantity `b_Intercept = Intercept - dot_product(means_X, b)`
  subtracts. Verified directly against that identity.
- `b`, `b_<dpar>`: frmtmb's `fixef()` indexed BY brms's column names,
  so the order is brms's. `Intercept` maps to frmtmb's `(Intercept)`.
- `bsp`, `bsp_<dpar>`: the coefficients of that predictor that are not
  columns of `X`, in frmtmb's order.
- `simo_<j>`: `softmax(c(0, zeta_j))` from frmtmb's `zeta<j>`.
- `sd_<i>`, `L_<i>`: from `VarCorr(fit)`, subset and ordered to brms's
  coefficient order. `L` is the lower Cholesky factor of the
  correlation matrix.
- `z_<i>`: `solve(diag(sd) %*% L, t(r))` with `r` from `ranef(fit)`
  reindexed to brms's level labels. Taking `r` from the raw `b` vector
  would follow frmtmb's `Zt` order, which is not brms's.
- bare dpar names (`sigma`, `shape`, `phi`, ...): brms declares these on
  the natural scale, so frmtmb's link-scale intercept goes through
  `family(fit)$links[[dpar]]$linkinv`.
- log-Jacobian for check C: `n_levels * (sum(log(sd)) +
  sum(log(diag(L))))` per group, summed.

### Finding: the coefficient order inside a group is not in get_prior

The first implementation read brms's within-group coefficient order off
`get_prior()`'s `sd` rows. Those rows are sorted alphabetically, so
`(Days | Subject)` came back as Days then Intercept while Stan's
`Z_1_1` is the intercept and `Z_1_2` is Days. Check A was wrong by 1149
nats and check B put 295 on the z block where zero belongs. This is the
mapping error check B exists to catch, and it caught it on the first
random-effect row. The order now comes from brms's ranef table, and the
test asserts it explicitly.

### Divergence: mo() simplexes are shared per variable, not per term

Row 3, `y ~ mo(inc) * z`, is not the same model in the two packages,
and no parameter map can make it one.

brms builds one simplex per SPECIAL TERM. For `mo(inc) * z` its
standata carries `Xmo_1` and `Xmo_2`, `Imo = 2`, `Jmo = c(3, 3)` and
`con_simo_1`, `con_simo_2`, and its parameter block declares `simo_1`
and `simo_2` as independent simplexes. The monotonic shape of the main
effect and the monotonic shape inside the interaction are free to
differ, which is the documented point of the feature.

frmtmb builds one simplex per mo() VARIABLE. `R/frame.R:1865` says so
in as many words, and the key is `deparse1(mexpr)`, the variable
expression, so `mo(inc)` and `mo(inc):z` share `zeta1`. frmtmb's model
therefore has two fewer free parameters than brms's and constrains the
interaction to the shape of the main effect.

Which side is right: brms's. Sharing is the stricter model, and a user
who writes `mo(inc) * z` after reading brms's monotonic vignette is
asking for a shape that may vary with `z`. frmtmb silently fits
something else. This is a feature gap, not an arithmetic error, and it
is only reachable when one mo() variable appears in more than one term;
`y ~ mo(inc) + z` and `y ~ mo(inc):z` each use a single simplex in both
packages and agree exactly.

This is the first and only entry on the exemption list the plan keeps.
It is admitted because it is a design choice frmtmb states in its own
source, not because a tolerance was widened. The test records it by
asserting the structural difference, so that it fails loudly if either
package changes, rather than by skipping the row.

**This is the second independent report, not the first.**
`dev/brms-vignette-audit.md:513-517` already recorded it: "`mo(income)
* age` makes the interaction SHARE the main effect's simplex. brms fits
two. The shape of the monotonic effect therefore cannot vary with age",
and that audit lists the shared simplex at its line 680 among the items
new to it. The evidence here is stronger only in that it has the Stan
parameter block; the conclusion is the same one, reached twice.

**Why it survived both.** `NEWS.md:2364`, in the v0.18 entry, gives the
rationale as "`mo()` interactions share their variable's simplex (brms
convention)". That parenthetical is false: one simplex per term IS the
brms convention. The design was chosen to match brms and does not match
brms, so a maintainer reading only NEWS would conclude the behavior is
correct and compatible. The mechanism is deliberate; the premise is
mistaken. `dev/feature-gaps.md:421` and `R/parse.R:878` record the
choice without repeating the false justification.

Decision for the maintainer, which this lane does not make: either
frmtmb keys `mo_zetas` on the term rather than on `deparse1(mexpr)` at
`R/frame.R`, or the sharing is announced where users of `mo()` will
meet it. The one thing that should change either way is the NEWS
rationale, since it is what made the divergence invisible.

### Cache design

`brms_stan_model(code)` is content addressed on the Stan program.

- Key: md5 of the generated Stan code with the rstan version appended
  as a trailing comment line. The rstan version belongs in the key
  because the cached object carries a compiled DSO, and a DSO built by
  one rstan is not loadable by another.
- Store: one `<key>.rds` per program under `FRMTMB_STAN_CACHE`,
  defaulting to `dev/stan-cache`. `rstan::stan_model(save_dso = TRUE)`
  is what makes the saved object usable from a later session; without
  it the RDS reloads into a recompile.
- The code depends on the SHAPE of the model, not on the values in the
  data, so one compile serves every data variant of a row and the two
  translator tests share their programs with the rows that follow them.
  Distinct formulas can still land on the same key when they generate
  byte-identical programs (`y ~ x` and `y ~ x + I(x^2)` differ only in
  standata), which is a saving rather than a collision: the key is the
  program, and the program is what is being cached.
- A session cache, `.brms_stan_models`, sits in front of the disk
  cache, and it is a correctness fix rather than a speed one. A model
  object compiled in THIS session and then re-read from its own RDS
  comes back with a DSO that will not initialize ("NULL value passed
  for DllInfo"), so the first request for a program worked and every
  later one in the same process failed. An RDS written by an EARLIER
  session reads back fine, which is why a warm run never saw it. The
  file uses two programs twice, so from an empty cache the old code
  errored on rows 1 and 2: 93 assertions with 2 errors instead of 97
  with none. Found by review, not by this lane, which had only ever run
  warm. Every program is now read from disk at most once per process.
- A corrupt or unreadable RDS falls through to a recompile rather than
  failing the run, because a cache is not allowed to be the thing that
  breaks a test.
- `dev/stan-cache` is in `.gitignore`. It needs no `.Rbuildignore`
  entry: `^dev$` already excludes the whole directory from the tarball.

In CI `actions/cache` stores the same directory. The key is the rstan
version plus a hash of the two test files, because those two files
determine the set of programs a run needs; `restore-keys` falls back to
the version prefix, so a changed test file still restores every program
that did not change and only the new ones compile.

### The harness as built

`tests/testthat/helper-brms.R` gains, in order:

- `brms_stan_cache_dir()`, `brms_stan_cache_key()`, `brms_stan_model()`
  and `.brms_stan_models`: the content-addressed compile cache and the
  session cache in front of it.
- `brms_flat_prior()`: `get_prior()` with every `prior` cell set to `""`.
- `brms_stan_par_names()`: the names declared in the Stan `parameters`
  block, in declaration order. The translator walks these, so a shape
  with no rule is an error that names the parameter.
- `brms_Xc_cols()`, `brms_fe_of()`, `brms_X_of()`, `brms_coef_to_frm()`:
  the population-level plumbing.
- `brms_ranef_table()`, `brms_group_info()`, `brms_ranef_block()`,
  `brms_group_pars()`, `brms_blockdiag()`: the group-level plumbing.
- `stan_pars_from_fit()`: the translator.
- `brms_z_index()`: the z slots of the unconstrained vector, found by
  perturbing z and diffing, so it needs no knowledge of declaration
  order or of any other parameter's unconstrained size.
- `expect_par_roundtrip()`: the map's own unit test, through
  `rstan::unconstrain_pars()` and back through `constrain_pars()`.
- `brms_lp_check()`: checks A and B, and C when `joint = TRUE`.

`brms_lp_check()` takes `const` (the admitted constant, default zero)
and returns the measured constant and gradient invisibly, so a row that
drifts reports a number rather than only a pass or a fail.

The ordinal threshold rule calls `ord_tau_from_raw()`, which is
internal to frmtmb. testthat runs in the package namespace, so this
needs no export and no change under `R/`. Copying the two-line
transform into the test instead would have made the test stale rather
than failing if frmtmb ever changed it.

### Results

Measured on Windows 11, R 4.6.1, rstan 2.32.7, brms 2.23.0, with the
C++17 repair described above. "constant" is `log_prob` minus frmtmb's
log density, so zero is the claim and anything else is either an
admitted Dirichlet or a finding. "max abs gradient" is check B, over
every parameter for A and B rows and over the z block for C rows; the
threshold is 1e-3, and what is left is frmtmb's own convergence
tolerance, not a mapping residual.

Times are wall clock for the whole row, including the Stan compile on
its first run. A cached row costs about 4 s, which is what the C0 line
shows.

The file as a whole, one file per process, with
`FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`:

| run | tests | assertions | failed | skipped | errors | wall | programs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cold, empty cache | 19 | 97 | 0 | 0 | 0 | 1481.7 s | 24 compiled |
| warm | 19 | 97 | 0 | 0 | 0 | 15.3 s | 24 restored |

The cold run is the one that matters, because it is what the first run
of the CI workflow does. Before the session cache went in it was 93
assertions with 2 errors, and there was no configuration in which the
workflow's debut could be green.

Three rows read FAIL and none of them is a numeric mismatch.

`3a` and `C3` are the one documented exemption. Both stop inside rstan
because brms declares a `simo_2` that frmtmb has no parameter for; see
the divergence section above. The spellings the two packages do agree
on, `3b` and `C3b`, are in the table and are green.

`8` is simply not reached. The translator has no rule yet for the
multivariate response suffixes or for `Lrescor`, and an unwritten rule
is an error naming what is missing. It is kept in the table so the gap
stays visible, and it is not in the test file.

| row | result | constant | max abs gradient | s |
| --- | --- | --- | --- | --- |
| 1 gaussian, sigma ~ x | pass | 0 (to 1e-9) | 0.000264 | 121 |
| 12a cumulative(y ~ x) | pass | 0 (to 1e-9) | 0.000163 | 2 |
| 12b sratio(y ~ x) | pass | 0 (to 1e-9) | 0.000264 | 0 |
| 12c cratio(y ~ x) | pass | 0 (to 1e-9) | 3.14e-05 | 0 |
| 12d acat(y ~ x) | pass | 0 (to 1e-9) | 4.36e-06 | 1 |
| 12e sratio with cs(z) | pass | 0 (to 1e-9) | 0.00035 | 6 |
| 13 categorical(y ~ x) | pass | 0 (to 1e-9) | 3.42e-05 | 101 |
| 14a right censoring | pass | 0 (to 1e-9) | 0.000145 | 65 |
| 14c known standard errors | pass | 0 (to 1e-9) | 0.000169 | 95 |
| 15 binomial with trials(n) | pass | 0 (to 1e-9) | 2.85e-09 | 92 |
| 16 zero-inflated poisson, zi ~ x | pass | 0 (to 1e-9) | 2e-05 | 61 |
| 17 mixture(gaussian, gaussian), theta1 ~ x | pass | 0 (to 1e-9) | 0.000318 | 4 |
| 2 mo(inc) + z | pass | 0.693147181 | 2e-05 | 111 |
| 20 weights(w) | pass | 0 (to 1e-9) | 0.000146 | 116 |
| 21a poisson | pass | 0 (to 1e-9) | 7.96e-05 | 78 |
| 21b Gamma(log) | pass | 0 (to 1e-9) | 0.000408 | 66 |
| 21c negbinomial | pass | 0 (to 1e-9) | 3.02e-07 | 141 |
| 21d bernoulli | pass | 0 (to 1e-9) | 2.77e-05 | 74 |
| 3a mo(inc) * z | FAIL | Exception: mismatch in number dimensions declared and found in context | | 106 |
| 3b mo(inc):z | pass | 0.693147181 | 5.82e-07 | 128 |
| 5 nonlinear a * exp(-b * x) | pass | 0 (to 1e-9) | 1.55e-06 | 93 |
| 7 (1 \| q \| g) across mu and sigma, joint | pass | 0 (to 1e-9) | 1.5e-10 | 4 |
| 8 mvbf with rescor | FAIL | frmtmb fit has no linear predictor named y1 | | 97 |
| C0 sleepstudy (Days \| Subject) joint | pass | 0 (to 1e-9) | 1.04e-14 | 4 |
| C16 zero-inflated poisson + (1 \| g) joint | pass | 0 (to 1e-9) | 1.63e-19 | 108 |
| C3 mo(inc) * z + (1 \| g) joint | FAIL | Exception: mismatch in number dimensions declared and found in context | | 155 |
| C3b mo(inc):z + (1 \| g) joint | pass | 0.693147181 | 8.88e-16 | 161 |

### Coverage against the plan's matrix

Green and carried into `tests/testthat/test-brms-likelihood.R`:

- 1, distributional gaussian.
- 2, `mo(inc) + z`, constant `lgamma(3)`.
- 3, `mo(inc):z`, constant `lgamma(3)`.
- 5, nonlinear.
- 7, `(1 | q | g)` merged across mu and sigma (C).
- 12, all four ordinal families and `cs()`.
- 13, categorical.
- 14, `cens()` and `se()`.
- 15, binomial with `trials()`.
- 16, zero-inflated poisson with `zi ~ x`.
- 17, `mixture(gaussian, gaussian)` with `theta1 ~ x`.
- 20, `weights()`.
- 21, poisson, Gamma(log), negbinomial, bernoulli.
- C anchors: sleepstudy `(Days | Subject)` with `sigma ~ Days`, row
  3's `(1 | g)` variant, row 16's `(1 | g)` variant.

Row 3's `mo(inc) * z` spelling is the documented exemption above, and
it is asserted structurally rather than skipped.

Rows not reached, and what each needs beyond what is written. Every
one of these is a branch in `stan_pars_from_fit()` plus a row; an
unwritten branch is an error that names the Stan parameter, so the
remaining work is enumerable rather than exploratory.

- 4, `mi()` measurement error: a rule for the latent response block.
- 6, nonlinear with random effects: nothing new beyond rows 5 and 7.
- 8 and 9, multivariate with `set_rescor(TRUE)`: `Lrescor` and the
  per-response naming, measured in the section below.
- 10, `gp()` exact and HSGP: `sdgp`, `lscale`, `zgp`.
- 11, `s()` and `t2()`: `bs`, `sds_*`, `zs_*`, and C on the wiggly
  part.
- 12, the multinomial half of row 13, and `trunc()` from row 14.
- 15's beta-binomial, and the hurdle families from row 16.
- 18, `ar`, `cosy`, `unstr`: `ar`, `ma`, `cosy`, `Lcortime`, `sderr`.
- 19, `car`: `rcar`, `sdcar`.
- 21, the rest of the family roster.
- the C variants of rows 4, 8, 10, 12.

### Known limits of the translator as written

These are not failures. They are places where a rule is correct for
every row now in the file and would need widening before a row that
reaches them is added.

- `simo_<j>` maps to frmtmb's `zeta<j>` by the trailing index alone. A
  model with monotonic terms under more than one dpar would need the
  dpar to take part in the match.
- `bsp` takes the special-term coefficients in frmtmb's own order,
  since brms's standata does not carry their names. With one special
  term per predictor this is unambiguous; with several, check B is what
  would catch a wrong order, as it did for the group-level coefficients.
- `brms_ranef_block()` pools every frmtmb ranef block on a grouping
  factor and then selects by coefficient name. A merged `(1 | ID | g)`
  block spans two dpars and repeats the coefficient name across them,
  which row 7 exposed: the selection now goes through
  `brms_frm_coef()`, which falls back to frmtmb's prefixed spelling
  (`y.sigma:(Intercept)`) when the bare name is not unique. Two blocks
  on the same factor under the same dpar would still be ambiguous.
- `bcs` is filled as `matrix(v, nrow = Kcs, ncol = nthres,
  byrow = TRUE)` from frmtmb's `bcs<j>` vectors. Row 12 exercises it
  only at `Kcs = 1`, where `byrow` cannot be distinguished from its
  opposite, so a second category-specific covariate is the case that
  would first test the assumption. The length check fires either way.
- frmtmb labels a merged block `1 | g + sigma: 1 | g [ID]`, so
  `brms_block_group()` strips the trailing `[ID]` marker before reading
  the grouping factor off the last bar-separated segment. The plain
  spelling, `Days | Subject`, needs no stripping.

### Finding: frmtmb stores ordinal thresholds two different ways

Rows 12c and 12d failed check A, in value and in gradient. The cause
was in the translator, not in frmtmb, and the check found it the way it
was designed to.

Measured by reverting `brms_ord_thresholds()` to pass `ordered = TRUE`
unconditionally, on exactly the data the test file's row 12 uses:

| family | reverted miss | reverted grad | current miss | current grad |
| --- | --- | --- | --- | --- |
| cumulative | 0.0000 | 0.0002 | 0 | 1.63e-04 |
| sratio | 0.0000 | 0.0001 | 0 | 5.42e-05 |
| cratio | -0.0234 | 1.5295 | 0 | 1.75e-04 |
| acat | -17.2373 | 44.8362 | 0 | 2.78e-05 |

Exactly the two families that store raw thresholds break, and the two
that store log increments are untouched. That is the conditional the
fix introduces, and both break in value AND in gradient, so check B
catches the map error independently of check A.

An earlier draft recorded 0.83 and 9.6 nats here. Those came from the
scratch driver, whose data frame draws an extra `rnorm()` column before
`y`, so its `y` is not the test file's. The table above is the
reproducible measurement and agrees with an independent review that
obtained the same two figures by two routes.

`R/families.R` says it plainly: cumulative and sratio hold
`(tau_1, log increments)` and pass `ordered = TRUE`; cratio and acat
hold the thresholds themselves and pass `ordered = FALSE`. The first
translator applied the ordered transform unconditionally, so for cratio
and acat it exponentiated numbers that were already thresholds. The
family object does not carry the flag, so `brms_ord_thresholds()` reads
the family name and duplicates that one fact on purpose. If frmtmb ever
changes a convention, these four rows fail.

Worth recording alongside it: frmtmb's sratio and cratio fits reach the
same maximized log-likelihood on the same data, to eight decimals, and
both now match brms exactly. Under the logit link brms's own
`sratio_logit_lpmf` and `cratio_logit_lpmf` are the same function
written two ways, so this is agreement, not a collision.

### Notes for row 8, multivariate with rescor

The row stops at `frmtmb fit has no linear predictor named y1`. The
facts needed to finish it, measured rather than guessed:

- brms declares `b_y1`, `Intercept_y1`, `sigma_y1`, `b_y2`,
  `Intercept_y2`, `sigma_y2`, `Lrescor`.
- frmtmb's `fixef()` keys the same quantities as `y1_mu`, `y1_sigma`,
  `y2_mu`, `y2_sigma`, so brms's suffix and frmtmb's prefix are the
  same name in the opposite order, and brms's `sigma_y1` is frmtmb's
  `y1_sigma`.
- standata carries `X_y1`, `Kc_y1`, `N_y1` per response.
- `family(fit)$links` is EMPTY for a multivariate fit, so the
  natural-scale rule cannot find its link there and needs the
  per-response family instead.
- the residual correlation lives in `fit$estimates$thetar`, and
  `Lrescor` wants its Cholesky factor.
