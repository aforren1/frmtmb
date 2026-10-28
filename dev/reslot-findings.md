# Lane reslot: one prior slot per parameter, brms's (flat)

Base 435e7ce6 (frmtmb 0.76.0), worktree `frmtmb-wt-reslot`, branch
`wt-reslot`, library `C:/Users/adf44/source/r/wt-reslot-lib`. "Before"
is rellib-r14 throughout. Nothing is committed. Logs are in
`dev/reslot-out/` (gitignored); every number below is pasted from them.

The items are the four of the lane brief; the user decided all four in
favor of brms. Status: DONE, all four items (this file is the record if the
report does not arrive).

## What brms 2.23.0 does (run, not recalled)

`dev/reslot-brms.R` (log `brms.txt`), `dev/reslot-brms2.R`
(`brms2.txt`), `dev/reslot-brms-corlb.R` (`brms-corlb.txt`);
`default_prior()`, `validate_prior()`, `stancode()`, `standata()` and
`rstan::stanc()` only. Seed 1, 8 groups of 10.

- brms has no class `theta` outside a mixture. On `(1 + x + z | g)`,
  `validate_prior()` refuses `set_prior("normal(0, 1)", class = "theta",
  coef = "theta_4")` ("theta_theta_4 ~ normal(0, 1)" does not
  correspond to any model parameter) and the class-wide `theta` too.
- brms has no prior on one correlation: `class = "cor", coef = "x"` is
  refused ("L_x ~ lkj(2)"). `lb`/`ub` on class `cor` are accepted by
  `validate_prior()`, and `stancode()` then writes
  `cholesky_factor_corr<lower=0>[M_1] L_1` and
  `lkj_corr_cholesky_lccdf(0 | 2)`, which `rstan::stanc()` rejects:
  "Syntax error in 'string', line 46, column 22" for `lb = 0` and for
  `lb = -0.5, ub = 0.5`; the same model without the bound parses. So
  brms has no working bound on a correlation either.
- `default_prior()` on `(1 + x + z | g)`: `b`, `b x`, `b z`, `cor`, `cor
  g`, `Intercept`, `sd`, `sd g`, `sd Intercept g`, `sd x g`, `sd z g`,
  `sigma`: no theta rows.
- The ported rows: `priors:5`'s classes are `b` x4, `cor` x2,
  `Intercept`, `sd` x6 (13 rows); `priors:123`'s validated table has 9
  rows.
- `""` in a slot (brms's `(flat)`) leaves no `lprior` term: with `""`
  on Intercept, sd, sigma and cor the program declares
  `vector<lower=0>[M_1] sd_1`, `real<lower=0> sigma`,
  `cholesky_factor_corr[M_1] L_1` and has no lkj line, so the density
  is Stan's implicit uniform on L's free entries. `""` on rescor: no
  `Lrescor` term. `ar(cov = TRUE)` with `ar` `""`:
  `vector<lower=-1,upper=1>[Kar] ar` and no term. Under `cov = TRUE`
  brms's `sigma` is the innovation sd (`normal_time_hom_lpdf(Y | mu,
  sigma, Lcortime, ...)` with `cholesky_cor_ar1()` dividing by
  `1 - ar^2`).
- `""` on a simplex keeps the Dirichlet with concentration 1:
  `con_simo_1` is `1 1 1` under `""` and under the default; a mixture's
  `con_theta` is `1 1` both ways.
- `delta` of an equidistant `cumulative()` is `(flat)` with `lb = 0` by
  default, declared `real<lower=0> delta`, no `lprior` term.
- `default_prior()` prints `(flat)` with `lb = -1, ub = 1` for
  `ar(cov = TRUE)`, `ma(cov = TRUE)` and both classes of
  `arma(cov = TRUE)`; `(flat)` with no bound for `ar()`; `(flat)`
  `lb = 0, ub = 1` for `cosy()`.
- Other spellings brms refuses: `class = "b", coef = "Intercept"` on a
  centered model ("b_Intercept"), `class = "Intercept", nlpar = "a"`
  ("Intercept_a"), `class = "sigma"` where sigma has a formula, and
  `class = "Intercept", dpar = "sigma"` where it has none.

## 1. One prior slot per random-effect and rescor parameter

### The defect, measured (`dev/reslot-double.R before`, `double-before.txt`)

Seed 1, `y ~ x + z + (1 + x + z | g)`, 8 groups of 10; theta_1 to 3 are
the log sds, theta_4 to 6 the correlation parameters.

- `frm()` resolution: `cor lkj(2)` alone is one entry `theta 4,5,6 lkj`;
  `cor lkj(2)` then `theta theta_4 normal(0, 1)` is one entry
  `theta 4 normal`, with no density on theta_5 and theta_6; in the
  other order, the LKJ alone. -log prior at the ML estimate:
  1.359461736187 (lkj alone) against 1.051441018731 (lkj + theta_4).
- MAP fits: `lkj + theta_4` stops at objective 86.8517134490 with theta
  = -0.331400 -0.576807 -1.654148 -0.019564 558.850086 1740.958672 (the
  two correlation parameters with no density run off); `theta_4 + lkj`
  at 88.6323489294 with theta_5, theta_6 = 0.196904, 0.354063.
- `frm_sample()` stack: the defaults alone end in `theta 4,5,6 lkj`;
  one `theta theta_4 normal(0, 1)` replaces it by `theta 4 normal`, so
  the default lkj(1) on theta_5 and theta_6 is gone. `sd` with group
  g, coef x plus `theta theta_2` keeps `theta 2 normal` and no sd entry.
- rescor (`bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE)`): `rescor lkj(2)`
  then `theta thetar_1 normal(0, 1)` keeps the normal; the other order
  the LKJ; in frm_sample one thetar_1 prior retires the rescor lkj(1).
- default_prior(): 7 theta rows (class-wide plus theta_1..6) on
  `(1 + x + z | g)`, and a thetar_1 row on the rescor model.

After (`double-after.txt`): every one of these is refused with the
class to write, in either order, on both routes; the class spellings
alone resolve as before; the tables list no theta row.

### What changed (R/priors.R)

- `cov_class_coefs(spec, frame)` lists every `theta` and `thetar`
  position a brms class reaches, with its slot (class, group, coef,
  resp, dpar, nlpar): `sd` over `block_sd_idx()` of every non-gp block
  (registered structures included, so car()'s sd), `sdgp` and `lscale`
  over a gp() block's positions, `cor` over `block_cor_spec()$idx` of
  every block whose `block_cor_prior()` is "lkj" (us, gr(cov), cs,
  homcs, ar1, hetar1), and `rescor` over every thetar position.
- `prior_table()` lists no class theta row for those; the class-wide
  theta row appears only when some theta position has no class.
- The theta branch of `resolve_priorlist()` refuses a coef (or the bare
  position "2") on such a position with `cov_class_refusal()`, which
  names the class and the set_prior() fields (for example
  `set_prior(..., class = "sd", group = "g", coef = "x")`, with
  `resp = "y1"` in a multivariate model) and states the class's scale.
  A class-wide theta covers only the positions no class reaches, and is
  refused where there are none.
- The lb/ub refusal on `cor`/`rescor`/`lncor` no longer points at class
  theta.

Guard-absent cases, all accepted (`test-prior-one-slot.R`,
`dev/reslot-car.R`, `car-before.txt`, `car-after.txt`): a toep()
correlation (theta_11, and a box on theta_12; the class-wide theta
covers exactly theta_11..19 of `toep(0 + tf | g)` with 10 time points);
a registered structure's own parameter (`xrho`, registered in the test
under `local_frmtmb_registries()`); frmtmb.spatial's car(): before, the
table listed theta, theta_1 and theta_2, now theta and theta_2;
`theta theta_2` and the class-wide theta resolve to `theta 2`, and
theta_1 is refused naming `class = "sd", group = "loc"`; the thetaac
route of a registered residual term without an interval is unchanged
(test-autocor-seam.R 80 of 80).

brms has no box on one correlation: `coef` on `cor` is refused, and
`lb`/`ub` on `cor` give a program that does not parse (above). So the
box on a single random-effect correlation is removed.

## 2. The sampling table prints (flat) with lb/ub

- Core `registered_prior_defaults()` prints `(flat)` for a default
  marked `flat_interval` (not `flat_region`). frmtmb.sample's
  `default_priors_for()` writes those defaults with prior `""` and the
  interval in lb/ub, so `prior_summary()` prints brms's row too; the
  announcement prints `(flat), lb = -1, ub = 1  [flat on its interval,
  brms's (flat)]`. The density is unchanged (the uniform on the
  interval).
- A higher-order `cov = TRUE` ar or ma (no brms term, `flat_region`)
  keeps its `uniform(-2, 2)` text: it is flat on a box, not an interval.
- `resolve_ac_class()` no longer writes a bound at an end of the
  interval; it bounded nothing already (`-Inf`/`Inf`).

Round trip (`dev/reslot-roundtrip.R`, data seed 4, points seed 7;
`roundtrip-before.txt`, `roundtrip-after.txt`): for ar(cov = TRUE),
ma(cov = TRUE), arma(cov = TRUE), ar(p = 2, cov = TRUE), cosy and sar,
the table passed back as `prior =` to frm_sample's resolver gives
identical bounds, the same charts, and over 300 points a largest
relative fn difference of 1.55e-16, 2.05e-16, 2.21e-16, 1.68e-16,
1.56e-16 and 2.09e-16, gradient difference 0 (the same on rellib-r14
with the old text). Table rows after: `(flat) ar -1 1`, `(flat) ma -1
1`, `(flat) ar -1 1` with `(flat) ma -1 1`, `uniform(-2, 2) ar NA NA`,
`(flat) cosy 0 1`, `(flat) lagsar -0.5211086 0.5211086`.

A consequence, and a decision: `frm(prior = default_prior(..., route =
"sample"))` on an ar(cov = TRUE) model is no longer the mode of the
sampled density. On frm()'s route the `(flat)` row with lb/ub is a
bounds-only specification, a box with no density, as every bounds-only
specification is there; the old `uniform(-1, 1)` row added the uniform
with the map's Jacobian. The test "frm()'s MAP under the sampling table
is the sampler's mode" (test-default-priors-autocor.R) now passes the
run's prior object (`prior_summary()`'s, which carries the densities),
and its comment says why. A brms (flat) row adds nothing to an
optimizer's objective either.

## 3. prior = "flat" means brms's (flat)

Decided by running brms (above): brms has no opt-out; its nearest model
is `""` in every slot, which is no lprior term, so Stan's implicit
uniform on each declared parameter within its declared bounds. That
has a brms meaning worth keeping, so the opt-out stays and means that.

### What changed

- Core: a density kind `"flat"` (`prior_base_logdens()` returns
  `0 * x`), which the resolver carries through the class's map, so the
  sampled density on the internal parameter is the map's Jacobian
  (`scale = "sd"` adds theta, `"natural"` the link's log|mu_eta|,
  `"ordthres"` the ordered map's, `"trans"` the residual map's, and the
  innovation-sd factor under cov = TRUE through `lc`).
  `lkj_dist(flat = TRUE)` is flat on the Cholesky factor L brms
  declares: row i of frmtmb's chart maps its i - 1 parameters t to L's
  free entries with Jacobian determinant (1 + |t|^2)^(-(i + 1)/2), so
  the row powers are (i + 1)/2, where lkj(1) has (d + 1)/2 on every row
  (the two agree at d = 2 only). On the ar1 and cs maps flat on the one
  correlation is lkj(1). `format_prior_dist()` and
  `prior_dist_string()` print it `(flat)` and `""`.
  `prior_entry_verdict()` (the weights check) reads it as order 0 at 0
  and as free on a link scale.
- frmtmb.sample: `default_priors_for(fit, flat = TRUE)` replaces every
  default density by `flat_spec()`: kind flat with the bounds kept, an
  LKJ flat on L; the flat-interval residual defaults and the
  dirichlet(1) simplexes stay (already flat on brms's parameter); a
  coefficient, a location intercept and a dpar intercept with a formula
  get nothing (the internal parameter is brms's).
  `sample_resolve_priors(fit, "flat")` resolves that stack and keeps
  the mixture weights' dirichlet(1) (brms's con_theta is 1 under "").
  check_laplace() is unchanged (`.diagnostic`, no prior at all).
- Also: a bounds-only user specification that replaces a default is
  brms's (flat) on that slot (`keep_default_bounds()`), and `delta` of
  an ordered equidistant family gets brms's default, `(flat)` with lb
  0, flat on delta (frmtmb holds log delta), by default and under flat.
- Messages and docs: the propriety warning now says the sd is flat on
  its own scale; the non-centering reasons say "carries no density"; a
  flat run non-centers (flat on the sd vanishes at sd = 0 in theta's
  coordinates, and flat on L is proper).

### Measured

`dev/reslot-flat-probe.R` (`flat-probe-before.txt`, `-after.txt`): on
rellib-r14 every model resolves "flat" to no entries. After, for
example, `(1 + x + z | g)`: `betad 1 flat sd`, `theta 1..3 flat sd`,
`theta 4,5,6 lkj internal`; cumulative equidistant: `tau_raw 1 flat
internal` and `tau_raw 2 flat sd`.

Against brms, `dev/reslot-flatpost.R` (data seeds 11, 12, 13, sampler
seed 1, 4 chains of 4000 with 1000 warmup; brms with `""` in every
slot; logs `flatpost-brms.txt`, `flatpost-before.txt`,
`flatpost-after.txt`):

| model | quantity | brms | after | before (rellib-r14) |
|---|---|---|---|---|
| `y ~ 1 + (1 \| g)`, 6 x 5 | sd_g__Intercept mean (sd) | 1.4581 (0.7107) | 1.4554 (0.7505) | 1.1969 (0.5902) |
| | sigma mean | 0.9728 | 0.9759 | 0.9570 |
| | divergent | 6 | 6 | 1 |
| `y ~ x + (1 + x \| g)`, 8 x 6 | cor mean (sd) | -0.2831 (0.5198) | -0.2646 (0.5265) | -1.0000 (0.0000) |
| | sd_g__Intercept mean | 0.4872 | 0.4858 | 0.3421 |
| | divergent | 8 | 1 | 529 |
| `y ~ 1 + ar(cov = TRUE)`, 6 x 5 | ar[1] mean (sd) | 0.4812 (0.1637) | 0.4797 (0.1621) | 0.5414 (0.1698) |
| | innovation sigma mean | 0.9654 | 0.9669 | 0.9410 |

frmtmb's `sigma` column is the marginal sd at 0.76.0, so the script
computes the innovation sd per draw as `sigma * sqrt(1 - ar^2)`.

Independent check of flat on L (test-flat-prior.R): the log
determinant of t to L's free entries, by central differences through
the exported `us_sigma()`, against the resolved density at 6 points:
spread below 1e-6 of the determinant's size, and lkj(1)'s spread
against it is more than 1000 times larger.

## 4. Other double slots

Every class pair that reaches one internal parameter was probed
(`dev/reslot-other-slots.R`, `other-before.txt`, `other-after.txt`):

| parameter | spellings | before | brms | now |
|---|---|---|---|---|
| centered intercept | `Intercept`; `b` coef `Intercept` | both reach beta_1, the later one applies, the b spelling without the centering | refuses `b_Intercept` | b spelling refused, naming class Intercept (with dpar) |
| nonlinear intercept | `b` nlpar a coef Intercept; `Intercept` nlpar a | two slots, the later applies (location 9 or 7 by order); frm() fitted with both | refuses `Intercept_a` | `Intercept` + nlpar is stored as `b` coef Intercept; both together are refused as a duplicate. The spelling stays accepted (user decision of 2026-09-24, test-prior-location-dpar.R) |
| gp() sd, length scale | `sdgp`, `lscale`; theta_1, theta_2 | later wins | no theta | theta refused (item 1) |
| smooth sd | `sd` group "s(x)"; theta_1 | later wins | sds | theta refused (item 1) |
| ar1(), cs() correlation | `cor`; theta_2 | later wins | no such structure | theta refused (item 1) |
| sigma without a formula | `sigma`; `Intercept` dpar sigma | already one slot: each is refused where the other applies | the same | unchanged |

Not changed, reported:

- The named-list prior (`prior = list(theta_1 = prior_normal())`)
  addresses internal parameters. In frm() it cannot be combined with a
  set_prior() specification; in frm_sample() it is merged after the
  class stack and takes over the parameters it names, whatever the
  order (a defined precedence, not a later-wins pair). brms has no
  such spelling. Whether to keep it is a user decision.
- frmtmb.spatial's car() dependence parameter has no brms class here
  (`car` is refused by name with a pointer to class theta) and is flat
  on its logit under frm_sample's defaults and under `prior = "flat"`,
  where brms declares it on (0, 1) and leaves it (flat). A fix needs a
  covariance-structure seam slot for a natural map and its Jacobian.
  Filed, not done.

## Models the change must not touch

`dev/reslot-identical.R before|after|compare` (seed 1;
`identical-compare.txt`): gaussian, poisson, ar(cov = TRUE), cosy,
flexible cumulative, student, Gamma, `(1 | g)`, `(1 + x | g)`, rescor
and equidistant cumulative models. For each: ML estimates, standard
errors, logLik, a MAP fit under class priors (estimates and objective),
both default_prior() tables, frm_sample's default -log prior at the
estimate, and its prior summary. Largest relative difference 0 on
every number of every model. The differences are the intended ones:
the ar and cosy sample tables and summaries print `(flat)`; the
`(1 | g)`, `(1 + x | g)` and rescor tables lose their theta rows; the
equidistant cumulative model gains the delta default: -log prior at
the estimate changes by -0.752901, which is -log(delta) there
(tau_raw_2 = log(delta) = 0.752901; the compare script prints it as
0.336, a difference relative to the larger value).

## Tests

### Seen to fail on rellib-r14 first

- `tests/testthat/test-prior-one-slot.R` (new, core): before 11 pass,
  31 fail (`dev/reslot-out/t-oneslot-before.txt`); the passes are the
  guard-absent cases (toep, the registered `xrho`, the ar1 cor reach,
  the accepted intercept spellings), which must hold on both. After:
  41 pass, 0 fail.
- `extensions/frmtmb.sample/tests/testthat/test-flat-prior.R` (new):
  before 53 pass, 17 fail, 2 error (`t-flat-before.txt`), the failures
  behavioral (for example -log prior + theta_1 + betad_1 = 3.803 where
  0 is expected; the flat-on-L spread 2.106 against a bound of 3e-6).
  After: 75 pass.
- `test-brms-suite-priors.R` (gated) with priors:5 and priors:123 marked
  pass: on rellib-r14, fail 2 (`reslot-g0bp`: priors:5 has 5 extra
  theta classes, `actual[11:18]` "sd" "sd" "sd" "theta" x5; priors:123
  has 11 rows against brms's 9). After: 36 pass, 0 fail.

### Existing tests rewritten (they pinned the old behavior)

Core: test-brms-formula-priors.R (the class set; "class theta rows name
real parameters" becomes "a standard deviation has one slot, class sd"
with `lognormal(-3, 0.05)` on class sd, which is `normal(-3, 0.05)` on
its log), test-lkj.R ("later wins between class cor and the class theta
hatch" becomes a refusal in either order), test-prior-mv-resp.R,
test-priors-autocor-classes.R (the refusal of thetar_1 and theta_1; an
end of an interval writes no bound; the lkj classes' lb/ub message),
test-priors-bounds-grcov.R (class sd bounds instead of class theta),
test-reml-nongaussian.R, test-review-fixes.R and test-tmb-examples.R
(pinning an sd through class sd, on the sd scale: `lb = ub = exp(th)`,
`lb = exp(1)`, `lb = ub = 1`), test-v19.R, test-brms-suite-priors.R
(priors:5 and :123 marked pass).

frmtmb.sample: test-default-priors-autocor.R ("(flat)" in the table and
the message; "frm()'s MAP under the sampling run's prior is its mode"
now passes the run's prior object, see item 2), test-autocor-seam-draws.R
(the announced interval), test-default-priors-dpar.R (flat keeps the
mixture weights' dirichlet(1)), test-laplace-draws.R (the REML refusal
is reached through `.diagnostic`, as check_laplace() reaches it),
test-reparam.R, test-sample-direct.R and test-sampling-ported.R (a flat
run non-centers and its prior summary lists (flat) rows; the boundary
fit's flat run has a median sd of 0.106 on seed 2, where the old flat
run sat below 0.01).

### Every suite, one file per process (`dev/reslot-tier.sh`, N = 6)

Ungated, tier reslot-u4 (`dev/release/reslot-u4.log`), 457 of 457 files,
every `lib:` line showing frmtmb from the lane library:

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 282 | 20109 | 0 | 0 | 204 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 32 | 1838 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 681 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 430 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 71 | 3344 | 0 | 0 | 14 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

Gated (`FRMTMB_BRMS_FIT_TESTS=true`), tier reslot-g1, 457 of 457:

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 282 | 22545 | 0 | 0 | 15 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 32 | 1844 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 12 | 681 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 501 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 71 | 3680 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 8 | 561 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

Every gated skip is another opt-in tier: 15 "set FRMTMB_SCALE_TESTS",
13 "set FRMTMB_DRMTMB_FIT_TESTS", 1 "set FRMTMB_FUZZ".

### OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`)

The new and rewritten files with a numerical assertion, one per process
under the OpenBLAS copy of R (`La_version()` 3.12.0; logs
`dev/reslot-out/ob-*.txt`): test-prior-one-slot.R 41, test-lkj.R 213,
test-priors-autocor-classes.R 169, test-priors-bounds-grcov.R 50,
test-reml-nongaussian.R 9, test-tmb-examples.R 27, test-flat-prior.R
75, test-default-priors-autocor.R 89; all fail 0, error 0. The new
assertions compare an exact identity (to 1e-12 of the scale, where the
two sides are computed by the same arithmetic) or a finite-difference
determinant (to 1e-6), so no rounding decides one; 0.3.32 was not run.
`dev/ciharden-out/` is deleted.

### Every suite in one session, three orders (`dev/reslot-sessions.sh`)

Tier reslot-s1 (`dev/release/reslot-s1.log`), N = 6, all nine
packages, ungated, frmtmb from the lane library:

    SESSION frmtmb forward files=272/282 pass=20109 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=272/282 pass=20352 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb reverse files=272/282 pass=20259 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb.coupling forward/preload/reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0 (each)
    SESSION frmtmb.eam forward/preload/reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0 (each)
    SESSION frmtmb.latent forward/preload/reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0 (each)
    SESSION frmtmb.learn forward/preload/reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0 (each)
    SESSION frmtmb.ode forward/preload/reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0 (each)
    SESSION frmtmb.sample forward/preload/reverse files=66/71 pass=3344 fail=0 err=0 skip=14 warn=0 (each)
    SESSION frmtmb.spatial forward/preload/reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0 (each)
    SESSION frmtmb.survival forward/preload/reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0 (each)

46 LEAK lines, every one tagged `[with load of ...]`; no untagged LEAK
line. (The files=272/282 and 66/71 counts have the same 10 and 5
uncounted files as earlier rounds' session lines.)

## Ported brms-suite rows (do not regenerate here)

`priors:5` and `priors:123` flip from divergence to pass. In
tests/testthat/test-brms-suite-priors.R both are now `brms_port(...,
"pass", "")`; on rellib-r14 that file fails exactly these two (above),
and on the lane library it passes 36 of 36 (gated). At consolidation,
`dev/brmsport-verdicts.tsv` and `dev/brmsport-verdicts-manual.tsv` need
the same two rows set to pass (their reasons were already stale:
per-coefficient sd rows have existed since 0.72.0, and priors:123 had
11 rows on rellib-r14, not 10). No other verdict row reads a class
theta row or `prior = "flat"`; `priors:85` and bin-2's `priors:130`
are about the fit route's flat table and keep their outcome.

## Not done, and why

- The named-list prior spelling and car()'s dependence parameter
  (item 4, above): reported, not changed.
- The sampling table's `lb` column: brms prints `lb = 0` on every `sd`
  and `sigma` row; frmtmb prints NA there (pre-existing, on both
  routes; not a slot question). Filed.
- `vignettes/brms-migration.Rmd` (the paragraph on prior classes kept
  elsewhere) still says sdgp and lscale are refused by name and that a
  `coef` on a `sd` row is refused; both have been wrong since 0.72.0 or
  earlier. Not this lane's text; filed.
- A user-facing `"flat"` density string (brms has none: its `""` is
  (flat)) was not added; frmtmb.sample writes the density object.
- brms's invalid program for `lb` on class `cor` is a brms defect; not
  posted upstream.

## Overlaps with the other lanes

- sigmainn: `R/priors.R` (my edits: `set_prior()` docs, the nlpar
  normalization in `set_prior_one()`, the lkj lb/ub message,
  `target_coefs()`'s centered-intercept refusal, `prior_table()`'s
  theta rows, new `cov_class_coefs()`/`cov_class_refusal()` after
  `theta_coef_target()`, the theta branch, the finite-bound guard in
  `resolve_ac_class()`, the lkj flat flag at four `lkj_dist()` calls,
  `lkj_dist()`/`lkj_logdens()`, `prior_base_logdens()`'s `flat`).
  I did not touch `innov_logc()` or `prior_logdens(lc =)`; the flat
  density on class sigma under cov = TRUE goes through them, so if
  sigmainn removes them, a flat sigma there must still be read at
  brms's sigma (it is the innovation sd in brms).
- priordraws: frmtmb.sample `R/sample.R` (`default_priors_for()`,
  `keep_default_bounds()`, `sample_resolve_priors()`'s flat branch,
  `announce_default_priors()`, docs). `prior = "flat"` now leaves
  entries of kind `"flat"` (and an lkj with `flat = TRUE`) in the
  effective prior, which `prior_draws()` and `sample_prior = "only"`
  cannot draw from (improper); they must refuse or skip them.
- mointerop: none known.
- NEWS.md of core and frmtmb.sample (new top sections).

## Bump and floor

- frmtmb: minor (BREAKING prior refusals and table rows; a new density
  kind read by the resolver). No new export.
- frmtmb.sample: minor (BREAKING `prior = "flat"`, the delta default,
  bounds-only specifications). Floor: the new frmtmb (it writes the
  `"flat"` density kind and the lkj `flat` flag, which frmtmb 0.76.0's
  resolver does not read, and its tests need the one-slot refusals).
- Other extensions: unchanged; all suites pass on the new core with
  their rellib-r14 builds.

## R CMD check --as-cran

Built and checked in `dev/reslot-check/<pkg>` (deleted afterwards; the
logs are `dev/reslot-out/check-frmtmb.log` and
`check-frmtmb.sample.log`), `dev/reslot-check.sh`, with
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH, and
`R_LIBS` = lane library, rellib-r14, user library:

    frmtmb:        Status: 1 NOTE   (tests 22 min, OK)
    frmtmb.sample: Status: 1 NOTE

Both NOTEs are the environmental "Skipping checking math rendering:
package 'V8' unavailable" on the HTML manual.

## Files touched

Code: `R/priors.R`, `R/sampling-api.R`, `R/weights-unbounded.R`,
`extensions/frmtmb.sample/R/sample.R`. Docs: `R/fit.R` and
`R/par-template.R` (roxygen only), `man/set_prior.Rd`,
`man/default_prior.Rd`, `man/frm.Rd`, `man/par_template.Rd`,
`extensions/frmtmb.sample/man/frm_sample.Rd`,
`extensions/frmtmb.sample/man/check_laplace.Rd`, the frmtmb.sample
vignettes `sampling.Rmd`, `brms-posterior.Rmd`,
`posterior-diagnostics.Rmd`, `NEWS.md`,
`extensions/frmtmb.sample/NEWS.md`. Tests: new
`tests/testthat/test-prior-one-slot.R` and
`extensions/frmtmb.sample/tests/testthat/test-flat-prior.R`; rewritten
as listed under Tests. Lane scripts: `dev/reslot-*.R`,
`dev/reslot-*.sh`, this file.

## Fix round 1 (review `dev/reviews/2026-10-27-reslot.md`)

Decisions given with the review: (a) keep the box reading of a `(flat)`
row on frm()'s route, document it and pin it in a test; (b) remove the
named-list prior spelling.

### 1. Finding 1 (major): the named list was a second slot; removed

- Core: `resolve_prior_input()` resolves a priorlist only;
  `resolve_priors()` is deleted. `as_priorlist()` refuses anything that
  is not NULL, a priorlist, a `brmsprior` or a prior table, with
  `prior_input_refusal()`: a named list of prior objects ("prior =
  list(theta_4 = ) put densities on internal parameters by name, ...
  Write each density with set_prior() on its class, e.g.
  set_prior(prior_normal(0, 5), class = "b")"), a plain data frame
  (brms refuses one too, "Argument 'prior' must be a 'brmsprior'
  object", run in `dev/reslot-brms-df.R`; the message points at
  `as.brmsprior()`, which reads brms's columns), anything else. frm()'s
  own `is.list()` check is gone (the refusal covers it); update()'s
  named-list merge and `is_named_prior_list()` are gone; the priorlist
  printer no longer prints `overrides`. The three messages that pointed
  at `prior = list(tau_raw = )` (one threshold, the whole vector, the
  sum-to-zero thresholds) and the compat row no longer offer it; brms
  has no per-threshold bound or a prior on sum-to-zero thresholds'
  free directions, so those go. `?frmtmb-priors` now describes the
  constructors as densities set_prior() takes; they stay.
- frmtmb.sample: `prior_stack()` has no legacy split, the merge in
  `sample_resolve_priors()` and `merge_prior_inputs()` are deleted, and
  `frm_sample()`'s docs and refusal no longer name the list. The
  `overrides` checks in `er_entries()` (evidence-ratio.R) and
  `prior_summary.frmtmb_draws()` are left as they are: they only read
  an attribute old saved draws may carry, and priordraws edits
  evidence-ratio.R.
- Test (test-prior-one-slot.R, "no spelling left reaches a parameter a
  class reaches"): on `(1 + x + z | g)` and a two-response rescor
  model, every position `cov_class_coefs()` lists is refused under
  class theta with its name as coef and as a named list, and the
  class-wide theta reaches none of them.
- Tests rewritten off the named list (17 uses in 7 files, plus
  test-conditions.R's message): test-brms-priors.R, test-lkj.R,
  test-priors-bounds-grcov.R (now: a named list and a plain data frame
  are refused, as.brmsprior() reads the frame and fits as set_prior()
  does), test-update-prior-drop.R (update() refuses the list),
  test-reparam.R (`lognormal(0, s)` on class sd, which is `normal(0, s)`
  on theta_1), test-sample-direct.R, test-sampling-ported.R.

### 2. Finding 2 (major): a flat run's prior_summary() round-trips

- `default_priors_for(flat = TRUE)` lists the location intercepts and
  each intercept of a dpar with a formula as `set_prior("", class =
  "Intercept", ...)`, flattened: kind `"flat"` on the internal
  parameter, which is brms's, so no density, but the slot is in the
  summary.
- `set_prior("", class = ...)` with no bound is accepted (brms's
  `(flat)` on the slot), on every class including the matrix and
  simplex classes. On frm()'s route it applies nothing:
  `drop_flat_slots()` removes such specifications before the fit, so a
  prior made only of them leaves `fit$prior` NULL and the ML fit
  bitwise unchanged (test). The resolver skips a density-free spec on
  the LKJ and simplex classes.
- A prior-table row (frmtmb's or brms's) with an empty prior and no
  bound is read as that `""` specification, where it was dropped.
- frmtmb.sample's `keep_default_bounds()`: a `""` specification that
  replaces a default becomes brms's (flat) on the default's own
  parameter: the uniform on an interval as before, flat on the Cholesky
  factor where the default is an LKJ, dirichlet(1) where it is a
  simplex, kind `"flat"` otherwise.
- Tests (test-flat-prior.R): the flat run's `prior_summary()` object,
  and its table read back with `as.brmsprior()`, give the run's
  objective: log density and gradient equal to 1e-12 relative at 40
  points (seed 4), through `prior_augmented_obj()`, not draws
  (rtmb-pitfalls item 23). `set_prior("", class = "cor")` gives the
  same flat-on-L entry `prior = "flat"` gives and leaves the sd default;
  `""` on the intercept replaces its student_t by kind flat.

### 3. Finding 3 (minor): the box reading, documented and pinned

`dev/reslot-boxfit.R` (the review's data, seed 13, 6 series of 8;
`boxfit-before.txt`, `boxfit-after.txt`), estimates (beta, betad,
thetaac):

| frm() prior | rellib-r14 | lane |
|---|---|---|
| none (ML) | 0.743156338 0.055843918 0.808909090 | the same |
| `default_prior(route = "sample")` | 0.755039925 0.030632434 0.724898994 | 0.745970184 0.063311362 0.805267000 |
| the run's prior object | 0.755039925 0.030632434 0.724898994 | the same |
| `uniform(-1, 1)` on ar | 0.755039925 0.030632434 0.724898994 | the same |

Core NEWS has a bullet with these numbers; `?default_prior` says the
sample table is the sampled density on frm_sample()'s route and a box
on frm()'s. Test (test-default-priors-autocor.R, "frm() reads the
sampling table's (flat) ar row as a box"): `frm(prior = tab)` equals
`frm(prior = tab without its ar row)` bitwise, and `uniform(-1, 1)` on
ar moves the ar estimate by more than 1 percent of it.

### 4. Finding 4 (minor): messages

- frmtmb.sample's lower/upper refusal names class `"sd"` and
  default_prior(), not theta_1.
- brms's `car` vocabulary refusal names the coef default_prior() lists
  under class theta and says the sd is class sd.
- The banner reads "prior = "flat" puts brms's (flat) in every slot".
- The fourth item (bound one threshold with class theta) is the
  ordinal message rewritten under finding 1.

### 5. Finding 5 (minor): car() is the exception

`?frm_sample` (prior = "flat") and frmtmb.sample NEWS say car()'s
dependence parameter is the one brms class that is not yet brms's
(flat): flat on its logit, not integrable at 0, improper under
`prior = "flat"`; give it a class theta prior. The seam slot stays
filed.

### 6. Nits

- The cor refusal with `resp` says resp on class cor is frmtmb's own
  spelling and that brms refuses it (test).
- `Intercept` with `nlpar` keeps the spelling written (`spec$written`)
  for the "Prior target not found" refusal: `(class=Intercept,
  nlpar=a)` (test on `a ~ 0 + x`).
- No test comment cites a lane script or log; numbers are inline.
- Thresholds as ratios: test-sample-direct.R compares the flat and
  default runs' median sd with 5 times the ML mode's sd, and the flat
  ESS with a fifth of the default's; test-flat-prior.R scales by
  `abs(f0)` and `max(abs(g0))`.
- A class-wide theta bound on toep() positions is tested again
  (test-priors-bounds-grcov.R: `lb = -2, ub = 2` lands on theta_4 and
  theta_5 of `toep(0 + tf | g)` with three levels, and on no sd).
- Filed, pre-existing on both arms (the review's
  `dev/review-reslot-p3b.R`): `default_prior()` on
  `y ~ x + rr(0 + tf | g, d = 1)` stops with "arguments imply
  differing number of rows: 1, 0".

### Seen to fail on rellib-r14 (tier reslot-fb1)

test-prior-one-slot.R 11 pass, 33 fail, 2 error; test-flat-prior.R 53
pass, 17 fail, 4 error; test-default-priors-autocor.R 87 pass, 5 fail
(the box test's 3, and the table text's 2). Every new test fails there
behaviorally except the guard-absent expectations, which pass on both.

### For consolidation

- sigmainn removes `innov_logc()`, `autocor_innov_sd()`,
  `prior_logdens(lc =)` and the entries' `innov` field. Then a flat
  sigma under cov = TRUE is plain flat on the exp map, which is right,
  but test-flat-prior.R's "prior = "flat" keeps the flat defaults on an
  interval" asserts `expect_false(is.null(sg[[1L]]$innov))` and must
  assert the plain flat entry instead. `flat_spec()`'s
  `s$innov_sd <- NULL` is harmless.
- priordraws: under `prior = "flat"` the effective prior holds kind
  `"flat"`, lkj with `flat = TRUE`, and now flat Intercept rows;
  `prior_draws()` / `sample_prior = "only"` must refuse or skip them,
  and draw the kept mixture dirichlet(1).
- `dev/brmsport-verdicts.tsv` and `-manual.tsv`: priors:5 and
  priors:123 to pass.

### Runs of fix round 1

Ungated, all nine suites one file per process (tier reslot-u5, 457 of
457, frmtmb from the lane library on every `lib:` line):

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 282 | 20136 | 0 | 1 | 204 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 32 | 1838 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 681 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 430 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 71 | 3354 | 0 | 0 | 14 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

The one error was test-input-validation.R's `frm(prior = 5)`, which
expects the argument named as `` `prior` ``; the generic refusal now
starts "`prior` must be a set_prior() specification", and the file
passes 43 of 43 in the gated rerun below.

Gated (`FRMTMB_BRMS_FIT_TESTS=true`), every touched file and the
review's list, one file per process (tier reslot-g2, 41 of 41, fail 0,
error 0, skip 0, warn 0 on every file): core test-brms-formula-priors
211, test-brms-priors 102, test-brms-suite-priors 36, test-conditions
150, test-get-prior-route 78, test-input-validation 43, test-lkj 213,
test-message-uniqueness 6, test-par-template 51, test-prior-compat 199,
test-prior-location-dpar 55, test-prior-mv-resp 36,
test-prior-one-slot 70, test-priors-autocor-classes 169,
test-priors-bounds-grcov 53, test-reml-nongaussian 9,
test-review-fixes 18, test-setprior 28, test-simulate-ergonomics 52,
test-thres-refit 118, test-tmb-examples 27, test-update-prior-drop 28,
test-v19 23, test-weights-unbounded 75; frmtmb.sample
test-autocor-seam-draws 19, test-default-priors-autocor 92,
test-default-priors-bounds 21, test-default-priors-brms 40,
test-default-priors-dpar 57, test-evidence-ratio 36,
test-flat-prior 85, test-laplace-draws 93, test-prior-route 9,
test-prior-update 7, test-reparam 373, test-sample-direct 135,
test-sampling-ported 214; frmtmb.spatial test-car 132,
test-car-postfit 14, test-car-v28 25, test-sar 126.

One session, three orders (tier reslot-s2, core and frmtmb.sample):

    SESSION frmtmb forward files=272/282 pass=20142 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=272/282 pass=20292 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=272/282 pass=20385 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.sample forward/reverse/preload files=66/71 pass=3354 fail=0 err=0 skip=14 warn=0 (each)

34 LEAK lines, all tagged `[with load of ...]`.

OpenBLAS 0.3.26 with LAPACK (`La_version()` 3.12.0), the files with
new numerical assertions (`dev/reslot-out/ob2-*.txt`):
test-prior-one-slot 70, test-priors-bounds-grcov 53, test-flat-prior
85, test-default-priors-autocor 92; fail 0, error 0. The copy is
deleted.

R CMD check --as-cran (`dev/reslot-check.sh`, one at a time; logs
`dev/reslot-out/check2-frmtmb.log`, `check2-frmtmb.sample.log`):

    frmtmb:        Status: 1 NOTE   (tests 19 min, OK)
    frmtmb.sample: Status: 1 NOTE

Both NOTEs are "Skipping checking math rendering: package 'V8'
unavailable".

### Bump and floor after fix round 1

Unchanged: frmtmb minor; frmtmb.sample minor, floor the new frmtmb
(it also relies on `as_priorlist()` refusing the named list and on
`set_prior("")`). No new export; `resolve_priors()` was internal.

## Fix round 2 (review "Short check of fix round 1")

### S1 (major): frm()'s table is refused by frm_sample()

`default_prior()` defaults to `route = "fit"`, whose `(flat)` rows are
frm()'s defaults. Since fix round 1 an empty row is brms's `(flat)` on
its slot, so a brms user's workflow (take `default_prior(f, data)`,
edit one row, call `frm_sample(prior = )`) sampled flat on every
unedited Intercept, sigma, sd and cor slot, silently; brms keeps its
defaults there.

- frmtmb.sample's `sample_resolve_priors()` refuses a `prior` whose
  `route` attribute is "fit": "prior = is a default_prior() or
  validate_prior() table of route = "fit", whose (flat) rows are
  frm()'s defaults ...", pointing to `default_prior(..., route =
  "sample")` and to `prior = "flat"`. It reaches frm_sample() before
  anything is sampled (frm_sample() keeps the mark through its own
  `as_priorlist()`).
- Core: `with_prior_route()` marks a priorlist read from a fit-route
  table (`as_priorlist()`, `as.brmsprior()`, `$<-`), and `+` and `c()`
  carry the mark when any operand has it (`combined_prior_route()`).
  `frm()` drops the mark from the prior it keeps (`drop_flat_slots()`),
  and frm_sample() drops it once checked, so a fit's own prior and a
  run's `prior_summary()` pass back.
- Survival of the mark, `dev/reslot-route-attr.R`
  (`route-attr-after2.txt`): kept by a column edit (`$<-`), a cell edit
  (`[<-`), a row subset, `rbind()`, `validate_prior()`'s table,
  `as.brmsprior()`, and `set_prior() + as.brmsprior(p)`. Not kept by
  `as.data.frame(p)` (a plain data frame, which `prior =` refuses with
  the pointer to `as.brmsprior()`) or `c(p, set_prior())` (a plain list,
  refused generically); `p + set_prior()` is an R error ("Incompatible
  methods"), before and after. A table with no mark, built by hand
  (`as.brmsprior(as.data.frame(p))`) or brms's own, is read as written,
  its `(flat)` rows as brms's `(flat)`, which is brms's meaning.
- frm() given a sample-route table: no refusal. Its rows with a density
  are brms's defaults, which frm() applies as a MAP penalty, the fit a
  user asks for by passing them; its `(flat)` rows are boxes or nothing
  on frm()'s route (decision (a)), which is also what `(flat)` means in
  brms. Nothing is read against the user's intent, so the asymmetry is
  harmless.
- Tests: test-flat-prior.R "frm()'s prior table is refused, however it
  was edited" (the edited table, a row subset, `as.brmsprior()`, `+`,
  `c()`, `validate_prior()`'s table, and `frm_sample()` itself; the
  edited sample-route table resolves to the defaults plus the edited
  row; a hand-built table is read as written); test-prior-one-slot.R
  "frm()'s table carries its route into a priorlist" (the mark through
  each operation, none on the fit's prior, and the fit equals the
  edited row's alone). On rellib-r14 (tier reslot-fb2) the first fails
  8 expectations and the second 5.

### S2 (minor): c() and a list refuse what is not a specification

`c.frmtmb_priorlist()` refuses an element that is not a specification
with `prior_input_refusal()`: a named list of densities gets the
named-list message ("prior = list(theta_4 = ) put densities on
internal parameters by name ..."), a bare density object "a density
object alone names no parameter; put it on a class with set_prior()".
`prior_input_refusal()` now finds named densities anywhere in a list,
so `prior = list(set_prior(...), theta_4 = prior_normal())` gets the
named-list message too. Test in test-prior-one-slot.R ("c() and a list
refuse what is not a specification"); on rellib-r14 it fails 2
expectations and errors once (vapply's "values must be length 1").

### Runs of fix round 2

Gated, one file per process (tier reslot-g3, 18 of 18, fail 0, error
0, skip 0, warn 0): core test-prior-one-slot 80, test-input-validation
43, test-message-uniqueness 6, test-setprior 28, test-prior-compat 199,
test-brms-formula-priors 211, test-update-prior-drop 28,
test-priors-bounds-grcov 53, test-lkj 213, test-brms-suite-priors 36;
frmtmb.sample test-flat-prior 95, test-message-uniqueness 6,
test-default-priors-autocor 92, test-sample-direct 135, test-reparam
373, test-sampling-ported 214, test-prior-route 9, test-prior-update 7.

### For consolidation (added in fix round 2)

- Do not remove the `overrides` checks in `er_entries()`
  (evidence-ratio.R) and `prior_summary.frmtmb_draws()`
  (methods-draws.R) in this lane; nothing sets the attribute any more,
  so they are dead code, to be removed at consolidation after
  priordraws' edits to evidence-ratio.R are merged.
- Under `prior = "flat"` a run's prior is no longer NULL (it holds
  kind "flat" entries, lkj with `flat = TRUE`, and the flat Intercept
  rows), so `er_entries()`'s "sampled flat" branch is never reached by
  such a run; the evidence ratio must refuse flat entries as improper
  (priordraws' code).

One session, three orders, core and frmtmb.sample (tier reslot-s3):

    SESSION frmtmb forward files=272/282 pass=20152 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=272/282 pass=20302 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=272/282 pass=20395 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.sample forward/reverse/preload files=66/71 pass=3364 fail=0 err=0 skip=14 warn=0 (each)

34 LEAK lines, all tagged `[with load of ...]`.

R CMD check --as-cran, one at a time (`dev/reslot-out/check3-*.log`):

    frmtmb:        Status: 1 NOTE   (tests 19 min, OK)
    frmtmb.sample: Status: 1 NOTE   (tests 350 s, OK)

Both NOTEs: "Skipping checking math rendering: package 'V8'
unavailable". No OpenBLAS run: this round adds no numerical code and
no numerical assertion.

Bump and floor: unchanged (frmtmb minor; frmtmb.sample minor, floor the
new frmtmb, which now also carries the route mark it reads).
