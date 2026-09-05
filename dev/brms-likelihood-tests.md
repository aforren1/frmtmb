# brms likelihood-identity tests: plan

Status: implemented in two rounds, on branches `wt-brmslp` and
`wt-brms-rows`, 2026-09-05. The harness, the translator with its
round-trip unit tests and the CI workflow came from the first round,
with rows 1, 2, 3, 5, 7, 12, 13, 14, 15, 16, 17, 20 and 21. The second
round added rows 4, 6, 8, 9, 10, 11, 18 and 19, which is every
remaining row of the matrix. The first round's one exemption, row 3's
`mo(inc) * z`, was closed on branch `wt-mo-terms`: frmtmb now builds
one simplex per mo() TERM and row 3b asserts the identity with the flat
Dirichlet admitted once per simplex. Four real divergences between the
packages remain, over three rows, and each is asserted rather than
skipped. Results, deviations, findings and the remaining work are under
"Implementation log" at the end of this document, and the second round
has its own section there.

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
  existing `skip_unless_brms_fit()`. That gate is TWO variables, not
  one: it calls `skip_unless_brms()`, which calls `skip_on_cran()`, so
  outside `R CMD check` a run needs `NOT_CRAN = "true"` beside
  `FRMTMB_BRMS_FIT_TESTS = "true"` or the whole file skips silently.
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

Once the flat-prior matrix is green, repeat A and B with the priors
brms's own `get_prior()` defaults supply, against `frm(prior =)`. This
is where the link-scale versus natural-scale placement of `sd` priors
will show as a value mismatch. That result decides whether the
placement changes or the documentation does; either way it stops being
a paragraph and becomes a measurement.

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

RESOLVED on branch `wt-mo-terms`: frmtmb builds one simplex per mo()
term, in brms's `stats::terms()` order, and row 3b now asserts the
identity at a constant of `2 * lgamma(3)`. What follows describes the
state before that change and is kept because it is the argument for
which side was right.

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
  order or of any other parameter's unconstrained size. The second
  round renamed it `brms_inner_index()` and widened its pattern, since
  a smooth, a GP, an `mi()` response and a CAR field all put an inner
  parameter where check C wants a zero gradient.
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

Row 3's `mo(inc) * z` spelling was the first round's one exemption and
is now an identity, closed on branch `wt-mo-terms`; see the RESOLVED
note above.

Rows not reached BY THE FIRST ROUND, and what each needed beyond what
was written then. Every one of these was a branch in
`stan_pars_from_fit()` plus a row; an unwritten branch is an error that
names the Stan parameter, so the remaining work was enumerable rather
than exploratory. All eight were written in the second round, whose
section at the end of this document supersedes this list.

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

Of these, the estimate of what row 6 needed was right and the estimate
of what rows 10, 11, 18 and 19 needed was not: each turned out to need
a scale or a basis correction as well as a name, and each of those
corrections is a finding in the second round's section.

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

## Second round: rows 4, 6, 8, 9, 10, 11, 18 and 19

Written on branch `wt-brms-rows` in the plan's order, which puts the
two multivariate rows first because their translator work unblocks the
rest. Reviewed independently afterwards, which reproduced every number
here and corrected three claims: the length-scale section said something
false about frmtmb (below), the `ar()` sigma section called a
documented divergence a discovery, and the esicar entry described a
softly constrained field as an unconstrained one. All three are fixed in
place, and the review added the measurement that at brms's nugget the
exact GP does not converge at all. Every row carries its translator
rules and a round trip through
`rstan::constrain_pars()` before the identity checks, and every rule is
driven off brms's declared parameter block, so a shape with no rule is
still an error that names the parameter.

### The parameter map, second round

- **multivariate responses.** brms suffixes every parameter of a
  response with the response name and puts the dpar FIRST (`b_sigma_y1`,
  `sigma_y1`); frmtmb keys the same linear predictor the other way round
  (`y1_sigma`). `brms_lp_parts()` rewrites one into the other and is the
  identity for a univariate fit, so every population-level rule went
  through it rather than being duplicated.
- `family(fit)$links` is EMPTY for a multivariate fit, because each
  response carries its own family object. `brms_link_of()` picks the
  family by response, and the natural-scale dpar rule reads the link
  from there.
- `Lrescor`: the Cholesky factor of `rescor_matrix(fit)`, reordered into
  brms's response order, which is the order the responses appear in
  standata and not an alphabetical one.
- a group-level block that spans two RESPONSES repeats the coefficient
  name once per response, so `brms_frm_coef()` now takes the response as
  well as the dpar; row 7's dpar fallback alone cannot separate
  `y1.mu:(Intercept)` from `y2.mu:(Intercept)`.
- `Ymi_<r>`, `Yl_<r>`: frmtmb keeps a plain `mi()` response's
  imputations and a measurement-error response's latent values in ONE
  `miss` vector, per response in formula order, and
  `fit$frame$mi_map[[r]]` gives the rows and the slots. The map to brms
  is the identity, so these add nothing to the log-Jacobian, and the
  rule asserts that `Jmi_<r>` names the same rows frmtmb imputes.
- `bs`: the unpenalized columns of a smooth, which brms puts in `Xs` and
  does NOT center; frmtmb names them `<term>.fx<j>` among the
  population-level effects. The two spellings span the same columns on a
  different scale, so the coefficients go through a change of basis read
  off the two design matrices. `bs` carries no prior, so any invertible
  map is admissible there.
- `sds_<i>`, `zs_<i>_<j>`: brms gives a smooth TERM one SD vector of
  length `nb_<i>` and one standardized vector per basis. frmtmb makes
  each basis its own homogeneous-diagonal block and labels all of them
  with the term, so the term is the run of blocks that share a label and
  the basis is its position in that run. The coefficients go through the
  same change of basis, which for a PENALIZED block must be orthogonal
  (below).
- `sdgp_<i>`, `lscale_<i>`, `zgp_<i>`: `exp(theta)` for the two scales
  and `solve(chol_lower(Sigma), f)` for the latent variables, with
  `Sigma` from the same covstruct registry entry the objective uses.
  `VarCorr()` cannot serve here: it reports only the marginal SD for a
  smooth or a GP block.
- `ar`, `cosy`, `Lcortime`: through frmtmb's own
  `autocor_natural(theta, ac)`, which already returns these on brms's
  scales under brms's names, and `us_chol_L()` for the unstructured
  factor, whose row-normalized lower-triangular form IS a Cholesky
  correlation factor.
- `sigma` under `ar(cov = TRUE)`: times `sqrt(1 - ar^2)` (below).
- `sdcar`, `car`, `rcar`, `zcar`: `exp(theta[1])` and
  `car_rho(theta[2])`; the proper form (escar) declares the field
  itself, so its map is the identity, and the intrinsic one (icar)
  declares `zcar = rcar / sdcar` and so carries `Nloc * log(sdcar)`.

### Divergence: the exact gp() nugget

`gp(x)` without `k` is not the same joint density in the two packages,
and no parameter map can make it one.

Both stabilize the kernel with a nugget and the two nuggets differ by
six orders of magnitude. brms's generated `gp_exp_quad()` writes
`cov[n, n] += 1e-12`, an absolute floor. frmtmb's `gp_corr()`
(`R/covstruct.R:1421`) returns `exp(-Q) + diag(1e-6, n)`, a nugget on
the CORRELATION, which is `1e-6 * sdgp^2` on the covariance.

Everything else agrees. On the row's design, 80 points of `sin(x)` on
[0, 10], the two kernels are equal off the diagonal to 3.9e-16 once the
length-scale is put on brms's scale, and the diagonals differ by exactly
`1e-6 * sdgp^2 - 1e-12`. That difference is not small in the joint
density: 30 of the 80 eigenvalues of the correlation lie below 1e-6 and
none lies below 1e-12, so the log density of the SAME field differs by
216.6 nats. The two halves are a log-determinant contribution of
+429.72, half the raw log-determinant gap of 859.4 (frmtmb -967.1
against brms -1826.6), and a quadratic-form contribution of -213.11
(7.1 against 433.3).

The marginal likelihood barely moves, which is why the two fits agree
and the two joint densities do not: the nugget enters the marginal
covariance as `Z (K + 1e-6 sdgp^2 I) Z' + sigma^2 I`, where it is 1e-5
of `sigma^2` on this design. Measured by refitting with `gp_corr()`
replaced in the installed namespace, on the row's design and on a
long-length-scale design built to be hostile (`2 sin(x/5) + 0.3x`,
fitted length-scale 12.2, 75 of 80 eigenvalues below 1e-6), moving the
nugget from 1e-6 to 1e-12 shifts `sdgp` by 0.22 and 3.60 percent, the
length-scale by 0.25 and 1.86 percent, `sigma` in the fourth decimal,
and the log-likelihood by less than 6e-4 nats against the 1.92 nats of a
one-parameter 95 percent confidence boundary. The nugget does not
measurably bias the fit.

Which side is right: brms's, for the density both packages advertise. A
squared exponential kernel is what both document, and a relative nugget
of 1e-6 is a white-noise term with SD `1e-3 * sdgp` that frmtmb does not
mention. frmtmb's own comment says why it is there ("keeps the
factorization stable at long lengthscales"), and that is not a hedge but
a measured fact: at 1e-12 the fit does not converge on EITHER design.
Both raise "false convergence (8)" with a large maximum absolute
gradient at the optimum, 0.117 on the row's design and 0.157 on the long
one, where 1e-6 raises neither. Part of the small shift tabulated above
is therefore optimizer failure rather than bias, and the true bias is
smaller still.

The decision for the maintainer, which this lane does not make, is
whether to document the nugget as part of the model, to keep a smaller
relative one, or to match brms's constant. That last option is listed
last on purpose: as the code stands it is the one that demonstrably does
not work, because frmtmb cannot fit at 1e-12 on either design tried.

The row is not skipped. The test asserts brms's constant by matching it
in the generated program, asserts frmtmb's by comparing the two kernels
entry by entry, and measures the gap, so it fails loudly if either
package changes its mind. The Hilbert-space spelling, `gp(x, k = 10)`,
has no nugget on either side and passes checks A, B and C exactly.

### Finding: the two GP forms carry different INTERNAL length-scales

Not a divergence, because it is a reparameterization the map absorbs,
but it is the fact the translator's `dmax` division rests on.

brms scales a `gp()` term's coordinates to unit maximum pairwise
distance before passing them to Stan (`Xgp_1` has maximum squared
distance 1 and `dmax_1` carries the divisor), so its `lscale` is on that
scale. frmtmb's exact form keeps the raw coordinates in `aux_D2`
(`R/covstruct.R:1415-1422`), so its internal `theta` is a data-scale
length-scale, while the Hilbert-space form divides by `dmax` at frame
construction (`R/frame.R:1763-1772`). On the row's design,
`dmax = 9.768956`, the same kernel is `exp(theta[2]) = 1.731704` under
`gp(x)` and `0.177266` under `gp(x, k = )`. The translator divides the
exact form by `dmax_<i>`, and with that one correction the kernels agree
to 4e-16; the Hilbert-space form needs no correction, which the row's
spectral-density check confirms to 2.9e-9.

The inconsistency is internal only. `confint_varcorr()` shifts the
Hilbert-space length-scale back to data units by `log(gp_dmax)`
(`R/confint.R:498-510`), so both forms REPORT `range(gp)` in data units,
1.731704 and 1.731698 on the row's design, and the vignette's promise
that "reported lengthscales stay in data units"
(`vignettes/frmtmb.Rmd:133`) holds. What differs is the internal theta,
which the translator absorbs. brms rescales for both of its forms too,
so neither of frmtmb's reported ranges is directly comparable with a
brms `lscale` without multiplying by `dmax`.

An earlier draft of this section said the two forms report on two
different scales and that only the Hilbert-space one is comparable with
brms. Both clauses were wrong, and neither was measured: the reported
ranges agree to 3e-6 relative. Recorded here rather than deleted,
because the wrong version named frmtmb as the inconsistent party.

### Finding: s() and t2() do not use the same basis convention

Already recorded structurally in `test-brms-agreement.R` ("brms calls
`mgcv::smoothCon()` with `diagonal.penalty = TRUE`, we do not ... the
individual basis columns are a reparameterization of each other"), but
what that reparameterization IS was not measured, and the identity
needs it.

Measured on row 11's design, the change of basis between the two `s(x)`
penalized bases is a PERMUTATION: the first column is shared and the
remaining seven are reversed, `A'A = I` to 2.3e-14 and `det A = -1`.
Because it is orthogonal, the i.i.d. prior on the coefficients survives
it and the two packages fit the same model in two orders. The
unpenalized column is the same direction on a different scale, frmtmb's
being 5.908294 times brms's on this design; that column carries no
prior, so the scale is free.

`t2(x, z)` uses identical bases element for element and its map is the
identity, which is why the same rule serves both.

The translator therefore reads the change of basis off the two design
matrices rather than assuming one, and REFUSES a penalized map that is
not orthogonal, because a non-orthogonal one would mean the two i.i.d.
priors are not the same prior and the models genuinely differ. Without
this, check A on `s(x)` missed by 77.58 nats and the z-block gradient
was 113.3, which is what a wrong basis looks like.

### Confirmed: brms's sigma under ar(cov = TRUE) is the innovation SD

Not a discovery. frmtmb had this divergence written down before the lane
began, in a source comment (`R/autocor.R:16-30`) and in an exported
`@section Divergence from brms` that carries it to the user
(`R/autocor.R:161-169`), down to the conversion
`sigma_marginal = sigma_innovation / sqrt(1 - phi^2)`. What the lane
adds is the first NUMERIC verification of that documented claim, which
is worth more than a discovery: the factor is not a reading of brms's
source but a term without which the two log densities differ by 5.2834
nats with a gradient of 58.17.

brms's generated `cholesky_cor_ar1()` ends with

    return cholesky_decompose(mat ./ (1 - ar^2));

so what it returns is not a correlation factor: it is the Cholesky
factor of the stationary AR(1) covariance with unit innovation
variance. `normal_time_hom_lpdf()` then multiplies it by `sigma`, which
makes brms's `sigma` the INNOVATION standard deviation. brms's own
comment above that function calls it "the cholesky factor of an AR1
correlation matrix", which it is not: a correlation matrix has a unit
diagonal and this one has `1 / (1 - ar^2)`. frmtmb builds a proper
correlation matrix and its `sigma` is the marginal residual SD, as it is
in nlme and everywhere else in the package.

Same model, two scales, related exactly by `sqrt(1 - ar^2)`. The
translator applies that factor and keys the rule on the function brms
actually emitted rather than on the family, so a program that stops
dividing stops being corrected. Without it check A missed by 5.2834
nats on this design and the gradient was 58.17.

The user-facing consequence is worth stating plainly: for
`ar(cov = TRUE)` frmtmb's reported `sigma` is not comparable with brms's
posterior for `sigma` on the same model. `cosy()` and `unstr()` carry
proper correlation matrices and need no correction.

### Divergence: ar() defaults to a different likelihood

brms's `ar()` defaults to `cov = FALSE`, the residual REGRESSION form,
which conditions on the first observations of each group instead of
giving them their stationary distribution. frmtmb implements only the
marginal residual-covariance form and refuses the other spelling by
name:

    ar(): only the residual-covariance formulation is implemented, so
    the call needs cov = TRUE: ar(time, g, cov = TRUE)

The two are different likelihoods on the same data, not two
parameterizations of one, so there is no map between them: brms's
default program declares no correlation factor at all and carries
`J_lag` instead. The row therefore runs on `cov = TRUE`, which both
packages implement and which agrees exactly, and the refusal is asserted
so that it stays a refusal rather than becoming a silent substitution.
This is a feature gap with a loud error, which is the right shape for
one, and it is the only one of the four exemptions that a user cannot
walk into by accident.

### Divergence: esicar and bym2 keep different latent variables

Two shapes where the joint densities are functions of different
arguments.

**esicar.** brms imposes the sum-to-zero constraint HARD: it declares
`Nloc - 1` free values, sets `rcar[Nloc] = -sum(zcar)`, and normalizes
by `(Nloc - 1) log(tau)`. frmtmb imposes it softly, folding
`1 / (con_sd n_j)^2` into the precision as a rank-`c` update
(`R/covstruct.R:1231-1245`), keeps all `Nloc` values, and normalizes by
`Nloc log(tau)` plus the log determinant of that precision. The gap
between the two log densities is
`log(sdcar) - 0.5 ldet_K + 0.5 Nloc log(2 pi)`, which MOVES WITH
`sdcar`: it is not a constant and no admitted constant can absorb it.

The soft constraint is a documented choice, not an omission.
`R/covstruct.R:916-947` argues it, names brms's own non-centered `zcar`
parameterization as the precedent, and tabulates the bias: at the 1e-3
default the log-likelihood is 4.7e-4 off the hard-constrained reference
and `sdcar` 3.6e-5 relative, while tightening `con_sd` costs optimizer
robustness (nlminb false convergence once in 25 refits at 1e-4 and six
times at 1e-5). The vignette says the same to users
(`vignettes/frmtmb.Rmd:178-183`).

**bym2.** brms keeps the spatial and the non-spatial parts as SEPARATE
latent vectors, `2 Nloc` values plus `rhocar`, and frmtmb integrates the
mixture into one dense marginal covariance over `Nloc` values. A point
identity between densities over different numbers of latent variables
does not exist.

Both are asserted structurally, from brms's generated parameter block
against frmtmb's block dimension, and the translator refuses them by
name rather than producing a map that cannot exist.

**Alongside it: frmtmb's esicar IS its icar, bit for bit.** `car_aux()`
(`R/covstruct.R:1204`) has no `esicar` branch at all: `escar` returns
early (`R/covstruct.R:1209-1230`) and everything else falls into the
shared intrinsic path (`R/covstruct.R:1231-1245`). On the row's design
the two fits agree to the last digit, log-likelihood -90.316306819585
and `sdcar` 1.413415 for both, and every field of `aux_car` is equal
except the `type` string. The vignette calls `esicar` an alias for
`icar` (`vignettes/frmtmb.Rmd:150-151`), so this is a documented choice
rather than an accident, and frmtmb's field IS constrained, softly.

What is left is a labeling contradiction inside one vignette section.
The same section opens by saying `car()` fits brms's spelling "with all
four of its types" while the package has three distinct densities under
four spellings, and in brms `esicar` is NOT an alias for `icar`: it is
the hard-constrained one, which this row's test asserts directly from
brms's generated program. A user porting a brms `esicar` model gets a
different, softly constrained likelihood under the same call and is told
nothing at the call site. That is one sentence of vignette text, not an
arithmetic defect. The test asserts the equality of the two frmtmb fits,
so it fails if the aliasing ever changes.

### The CAR normalizers brms omits, and why they are admitted constants

brms writes both of its CAR densities unnormalized in the field and
frmtmb keeps a proper density, so the two differ by a closed form in the
DATA alone. These are admitted the way the `mo()` Dirichlet is: not a
tolerance, an exact expression, and the test builds it from brms's OWN
standata so that it is a check on frmtmb's precision matrix rather than
a restatement of it.

- escar: `0.5 (Nloc log(2 pi) - sum(log(Nneigh)))`. brms's
  `sparse_car_lpdf` drops the `-Nloc/2 log(2 pi)` and the
  `0.5 log det D` of a proper CAR and keeps the rest term for term.
  Measured 6.149684293 on the row's 4 x 4 rook lattice, against a
  closed form of the same value.
- icar: `0.5 (Nloc - 1) log(2 pi) - log(0.001 Nloc) - 0.5 log det K`,
  with `K` the graph Laplacian built from brms's own edge list plus
  `1 / (0.001 Nloc)^2` in every cell. brms drops the ICAR normalizer
  entirely and adds a normalized soft sum-to-zero term,
  `normal(sum(zcar) | 0, 0.001 Nloc)`, whose precision is exactly the
  rank-one term frmtmb folds into its precision: frmtmb's `con_sd`
  defaults to 0.001 and its term is `1 / (con_sd n)^2`, the same number.
  Measured 5.253269634.

That the two quadratic forms agree at all is the real result: frmtmb's
`K` is brms's `D - W` plus brms's own constraint precision, checked
against brms's edge list rather than against frmtmb's adjacency matrix.

### Results, second round

Same conventions as the first round's table. "constant" is `log_prob`
minus frmtmb's log density; zero is the claim, and the only nonzero
entries are the two CAR normalizers above, each equal to its closed form
to the 1e-6 relative tolerance the check uses. "max abs gradient" is
check B, over the inner block for a C row and over every parameter
otherwise. Times are warm, one test per process.

| row | result | constant | max abs gradient | s |
| --- | --- | --- | --- | --- |
| 8 mvbf(y1, y2) + rescor | pass | -1.14e-13 | 0.000318 | 20.1 |
| 9 mv with (1 \| p \| g), joint | pass | 5.68e-14 | 1.78e-15 | 26.6 |
| 4 y \| mi() with x \| mi(sdx), joint | pass | 5.68e-14 | 7.55e-15 | 8.8 |
| 6 nonlinear, a ~ 1 + (1 \| g), joint | pass | -2.84e-14 | 8.66e-15 | 4.5 |
| 11a s(x), joint | pass | 0 | 6.51e-14 | 3.2 |
| 11b t2(x, z), joint | pass | 2.84e-14 | 2.44e-14 | 4.7 |
| 10b gp(x, k = 10), joint | pass | -7.11e-15 | 2.04e-14 | 2.7 |
| 10a gp(x) exact | EXEMPT | 216.6 nats, structural | | 1.1 |
| 18a ar(time, gr = g, p = 1, cov = TRUE) | pass | -2.56e-13 | 0.000459 | 7.7 |
| 18b cosy(time, gr = g) | pass | 2.56e-13 | 3.26e-05 | (same) |
| 18c unstr(time, gr = g) | pass | -1.99e-13 | 0.000232 | (same) |
| 18d ar(cov = FALSE) | EXEMPT | refused by frmtmb | | 0.7 |
| 19a car(escar), joint | pass | 6.149684293 | 3.77e-15 | 6.1 |
| 19b car(icar), joint | pass | 5.253269634 | 2.50e-12 | (same) |
| 19c car(esicar), car(bym2) | EXEMPT | different latent set | | 3.1 |

Three EXEMPT rows carrying four divergences: row 19c holds two, esicar
and bym2, which is why the count of divergences and the count of rows
that stand aside are not the same number.

The gradients that are not at machine zero are the A-and-B rows, where
check B is over every parameter and what is left is frmtmb's own
convergence tolerance. The C rows are asserted on the inner block only,
where frmtmb's inner solve puts the gradient at 1e-12 or below.

The file as a whole, one file per process, with
`FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`:

| run | tests | assertions | failed | skipped | errors | wall | programs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cold, empty cache | 31 | 347 | 0 | 0 | 0 | 4492.0 s | 36 compiled |
| warm | 31 | 347 | 0 | 0 | 0 | 51.0 s | 36 restored |

Twelve programs joined the 24 the first round needed. The cold figure
is an UPPER bound and is not comparable with the first round's
1481.7 s for 24: that one was measured on a quiet machine at about 62 s
a program, and this one ran with five sibling worktrees compiling and
fitting at the same time, which works out at about 125 s a program. The
CI workflow's 90-minute budget still covers it, but only just, and the
honest reading is that a cold run costs between 25 and 75 minutes
depending on what else has the cores.

One earlier attempt at the same cold run was killed at 25 of the 36
programs with no R output and exit 255, under the same contention. The
run was resumed against the 25 cached programs and finished green, and
the clean cold run above was then done from an empty directory, so the
figure is a measurement rather than a reconstruction.

### Missing accessors

Two values had no public accessor and were reached through the fit
object instead, per this round's rule that `R/` belongs to other lanes.
Both are recorded here rather than patched.

1. **A block's covariance at the fitted `theta`.** `VarCorr()` collapses
   a smooth or a GP block to its marginal SD on purpose ("the k x k
   identity blowup is noise"), so it cannot give the matrix the z rule
   has to invert. The test calls
   `covstruct_registry[[bk$covstruct]]$vcov(theta[bk$theta_idx], bk)`,
   which is exactly what `VarCorr.frmtmb_fit()` calls for every other
   structure. A public accessor for one block's fitted covariance would
   serve the tier and anyone writing a predictor by hand.
2. **The mi() row map.** `fit$frame$mi_map[[resp]]` carries the rows a
   response imputes and their slots in `estimates$miss`; nothing exports
   it, and `estimates$miss` alone cannot say which rows it belongs to
   when a model has more than one `mi()` response.

Reached without an accessor but worth noting as internals the tier now
depends on, alongside `ord_tau_from_raw()` from the first round:
`autocor_natural()`, `us_chol_L()`, `car_rho()`, `covstruct_registry`,
and the frame fields `re_blocks` (for `b_idx`, `c_idx`, `theta_idx`,
`dim`, `term_label`, `aux_car`), `linpreds` (for `X` and `Z`) and
`autocor`. testthat runs in the package namespace, so none of these
needs an export, but each is a place where a rename in `R/` breaks the
tier rather than the package.

### Coverage after the second round

Every row of the plan's matrix is now written. What remains inside rows
that are otherwise green, unchanged from the first round's list:

- 12, the multinomial half of row 13, and `trunc()` from row 14.
- 15's beta-binomial, and the hurdle families from row 16.
- 21, the rest of the family roster.
- the C variants of rows 4, 8, 10 and 12 that the plan asks for as
  extensions; rows 4, 9, 10b and 11 are already C rows as written, and
  row 8's C variant is row 9.
- `ma()` and `arma()` from row 18, which need the same rule as `ar` and
  a check on whether brms's `cholesky_cor_ma` divides the way
  `cholesky_cor_ar1` does.
- `gp(..., by =)`, where `Kgp_<i>` exceeds one; the translator refuses
  it by name.

### Known limits of the translator, second round

- `brms_autocor_of()` and `brms_car_block()` each require exactly ONE
  block of their kind, because brms declares `ar`, `cosy`, `Lcortime`,
  `rcar` and `sdcar` without an index. A model with two `car()` terms
  would need brms's own naming to change first.
- the smooth-term index maps to frmtmb's blocks by position within a
  linear predictor. Two smooth terms under one dpar are ordered the same
  way in both packages, and the basis dimensions are asserted against
  `knots_<i>`, so a reordering fails rather than passing wrongly.
- `brms_sigma_scale()` covers `cholesky_cor_ar1` only. An `ma()` or
  `arma()` row would have to check its own generated function.
- the mi() rule assumes brms's `Jmi_<r>` is sorted, which is what
  `which(is.na(y))` gives on both sides; it asserts the equality rather
  than trusting it.
- `brms_car_const()` hardcodes `s <- 0.001 * Nloc` rather than reading
  `con_sd` off the fit, because brms hardcodes the same number in its
  generated program and the constant is derived from BRMS's side. The
  row asserts frmtmb's `con_sd` is 0.001 next to it, so a changed
  default fails loudly rather than drifting; a `con_sd = ` argument on
  the row would need both halves to move together.
