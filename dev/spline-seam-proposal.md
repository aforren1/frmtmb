# A core API for splines the objective consumes nonlinearly

Written from `extensions/frmtmb.spline` on branch `wt-spline`, base
commit 5dfdd84 (core 0.50.0). **Nothing proposed here was implemented.**
This document says what core would have to grow, why, and what each
addition would oblige the rest of core to say.

There are two proposals, and they are independent. The first is small
and pays for itself immediately. The second is the real one.

---

## Part 1. The covariance of a grid prediction (small, immediate)

### What is missing

Three functions in `frmtmb.spline` need one object: the covariance of a
whole grid prediction, `Sigma = C V C'`, where `C` maps the coefficient
vector to the curve and `V` is the joint covariance of the fixed AND
random coefficients. A penalized smooth's wiggly part is a
random-effect block even when the smooth is a population term, so `V`
must span `b`.

Nothing exported returns either piece.

* `vcov(fit, full = TRUE)` (`R/methods-fit.R:377`) returns the OUTER
  parameter vector's covariance. Under ML that is `cov.fixed`, whose
  rows are `beta`, `betad`, `theta`, `thetaac`, `thetar`; `b` is
  integrated out and never was an outer parameter. Under REML or
  `control$profile` the blocks come from `jointPrecision` and the code
  then drops `b` by name. Neither branch can be changed: the documented
  invariant, asserted in `tests/testthat/test-methods-audit.R`, is that
  `full = TRUE`'s row names are exactly `confint()`'s.
* `predict(se.fit = TRUE)` computes exactly the right thing and throws
  away all but the diagonal. `predict.frmtmb_fit` starts at
  `R/predict.R:1019`; the line that discards `Sigma` is
  **`R/predict.R:1235`**:

      jc <- get_joint_cov(object)                # R/predict.R:12, @noRd
      da <- lp_delta_A(object, lp, ed, newdata, ...)   # @noRd
      V  <- jc$V[da$coef_pos, da$coef_pos, drop = FALSE]
      var_eta <- pmax(rowSums((A %*% V) * A), 0)   # R/predict.R:1235

  `rowSums()` is where `Sigma` dies.
* `ranef(condVar = TRUE)` gives `sqrt(diag.cov.random)`: the diagonal of
  the `b` block, with no cross-block against `beta`.
* `emmeans::emm_basis` returns `V = vcov(object)[idx, idx]`, the fixed
  block, so a smooth's wiggly columns are absent from it too.
* `hypothesis()` delta-methods over `beta, betad, theta, thetaac,
  thetar` (`hyp_par_cov()`, `R/confint.R:1660`). `b` is excluded there
  as well.
* `marginaleffects::set_coef` sets `beta` and `betad` only, so the
  design cannot be recovered by unit perturbation through that route.
* `frm_bootstrap(FUN = )` is the ONE exported route that reaches a full
  curve covariance, by refitting `nsim` times. It works and it costs
  `nsim` fits.

### Core already caches both halves, and the extension now reads them

This was missing from the first draft of this document and it changes
the cost argument. Core memoizes the sdreport and the inverted joint
precision on the fit:

    sdr_of(fit)        # R/fit.R:1209-1217 - caches fit$cache$sdr
    get_joint_cov(fit) # R/predict.R:12    - caches fit$cache$Vjoint

`fit$cache` is a `new.env(parent = emptyenv())` (`R/fit.R:938`), so the
memo survives the copy an extension makes. The extension's own
`predict(se.fit = TRUE)` check call populates `Vjoint` on the way past,
so by the time a curve needs the covariance it already exists.

Reading it rather than recomputing it fixed three things at once, all of
which the review measured:

* a second `RTMB::sdreport()` per call, gone;
* a Schur complement over the coefficients the curve does NOT touch,
  taken on a matrix densified out of the sparse one TMB produced, gone.
  That was the scaling wall: 114 s and 2.1 GB at 8000 random
  coefficients. Reading the cache instead is 8.2 s and 1.2 GB, and the
  covariance subset itself is 0.00 s;
* `autoscale_sdreport()` (`R/autoscale.R:134`) was being bypassed, so a
  fit with `par_units` set got a covariance built on the unscaled
  Hessian. Core's cached object is the autoscaled one.

So `frm_lp_basis()` is worth exporting for REACH, not for cost. The
cost is core's own joint-precision solve, and an extension can already
avoid paying for it twice.

### What the extension does instead, and what it costs

`eta` is linear in `(beta, b)`, so the difference between a prediction
and the same prediction with one coefficient raised by one IS that
coefficient's design column, exactly. The extension therefore rebuilds
`C` with one `predict()` call per contributing coefficient, and takes
`V` from the fit's own cache (above), falling back to
`RTMB::sdreport(fit$obj, getJointPrecision = TRUE)` only when core never
formed one.

Neither `fit$obj` nor `fit$estimates` is a `:::` reach:
`frmtmb.sample::frm_sample()` reads `fit$obj` (`sample.R:1271`) and the
documented `frmtmb-extension-api` example reads `fit$frame` and
`fit$estimates`. But both are undocumented slots, and the reconstruction
is a re-implementation of `lp_eta_design()` and `lp_delta_A()` by
numerical means.

It is checked rather than trusted: every call recomputes
`sqrt(diag(Sigma))` and compares it against `predict(se.fit = TRUE)`,
and refuses above `tol`. Measured agreement is 3e-15 to 9e-15 relative,
which is machine precision.

**The cost is not the call count, and an earlier draft of this document
said it was.** Measured at `re.form = NA` on a 20-point grid, one
process each, split by stage:

| model | random coefs | design rebuild | `predict(se.fit = TRUE)` | covariance subset |
|---|---:|---:|---:|---:|
| `s(x, k = 10)` | 8 | 0.01 s (11 calls) | 0.29 s | 0.00 s |
| `s(t, k = 8) + (1 + t \| subject)`, 1000 subjects | 2006 | 0.07 s (101 calls) | 0.98 s | 0.00 s |
| the same, 4000 subjects | 8006 | 0.28 s (351 calls) | 6.87 s | 0.00 s |

The design rebuild, which is the term the call count measures and the
term `frm_lp_basis()` would remove, is a TENTH of the cost at every
size. What grows is core's own joint-precision solve inside
`predict(se.fit = TRUE)`. The covariance subset is free because it is
read from the cache.

So the seam is worth exporting for correctness and reach, and the cost
argument for it is weak. The strong cost argument is against the version
of this extension that recomputed the covariance itself: that was 114 s
and 2.1 GB at 8000 coefficients, and it is now 8.2 s and 1.2 GB.

The call count still matters for a different reason, which is that it
does not grow with the number of LEVELS of a grouping factor: a block
that contributes nothing under `re.form = NA` is skipped in chunks of
24 rather than one coefficient at a time. On the model
`vignette("curve-inference")` fits, `v ~ s(t, k = 12) +
s(t, subject, bs = "fs", k = 5)` over 20 subjects, that is 32 calls
against 110 random coefficients. A factor-smooth model with NO
population smooth needs far fewer, because at `re.form = NA` the `fs`
term contributes nothing and the population curve is a constant; both
are pinned in `tests/testthat/test-curve.R`.

### The proposal: TWO asks, and the small one is the urgent one

They are different asks and they should not be traded against each
other. Part 1a legitimizes a dependency the extension already has and is
five lines. Part 1b adds reach the extension cannot get any other way
and is a larger job. The extension ships today; 1a is what would stop it
reading an internal.

#### Part 1a. An accessor for the cached joint covariance (five lines)

```r
frm_joint_cov(object)     # or simply export get_joint_cov()
#> list(V = <p x p>, names = character(p))
```

`get_joint_cov()` (`R/predict.R:12`) already returns exactly this and
already memoizes it on the fit. The ask is to export it, or to export a
thin wrapper, and to document the `names` element as the component
labels that index `V`.

**Why this one first.** `frmtmb.spline` reads `fit$cache$Vjoint`
directly. That is an internal: `get_joint_cov()` is `@noRd`, and neither
`fit$cache` nor the `list(V =, names =)` shape is documented anywhere in
frmtmb, so unlike `fit$obj` and `fit$estimates` there is no precedent to
point at. It is **the only unsanctioned dependency in the package**, and
exporting the accessor removes it outright.

The reach is defensible in the meantime, and the extension says so in
its own documentation (`?frm_curve`, "The one internal this reaches
into") rather than only here. Its failure modes are bounded: a renamed
or reshaped slot reads as `NULL` and falls back to a slower sparse
route; a changed MEANING is caught by the covariance check every call
makes; absence is the ordinary case on a fresh fit and is handled by
warming the cache first. None of them yields a wrong number. But
"bounded failure modes on an undocumented slot" is not a contract, and a
five-line export would make it one.

It also fixes something the fallback cannot. The fallback is a fresh
`RTMB::sdreport()`, which goes round `autoscale_sdreport()`
(`R/autoscale.R:134`) and would return an unscaled covariance on an
autoscaled fit; there is no exported route to the AUTOSCALED joint
precision at all, so the extension refuses that combination rather than
guessing. An accessor closes that hole too.

#### Part 1b. `frm_lp_basis()` (for reach, not for cost)

One exported function, and it already exists as two private ones.

```r
frm_lp_basis(object, newdata = NULL, dpar = NULL, resp = NULL,
             re.form = NULL, allow_new_levels = FALSE)
```

returning

```
list(eta        = numeric(n),        # what predict(type = "link") gives
     A          = <n x p matrix>,    # the design over the coefficients below
     coef_pos   = integer(p),        # rows of V, in V's own order
     V          = <p x p matrix>,    # the joint covariance at those rows
     coef_names = character(p),      # beta.<name> / b.<label>
     extra_var  = numeric(n))        # new-level and exact-gp variance
```

The name follows `emmeans`'s `emm_basis()`, which is the same idea for
the fixed block, and the shape is what `predict(se.fit = TRUE)` already
assembles: `A`, `coef_pos` and `V` come straight out of `lp_delta_A()`
and `get_joint_cov()`, and `extra_var` out of `lp_extra_var()`. The body
is a return statement.

`predict(se.fit = TRUE)` then becomes a two-line consumer of it, which
is the test that the seam is the right one: if `predict()` cannot be
written in terms of `frm_lp_basis()`, the shape is wrong.

**The cost argument for this one is gone, and it should not be made.**
An earlier draft of this document argued for `frm_lp_basis()` on speed.
Measured on the shipped package, the design rebuild it would replace is
**0.28 s of a 5.91 s call** at 8006 random coefficients: 351
`predict()` calls against core's own joint-precision solve, which is the
rest. Removing the design rebuild entirely would save under five percent
of that call. Part 1a is where the value is on cost, and it is already
realized.

**What this buys is REACH, and nothing else buys it.** Any delta-method
quantity over a fitted curve: contrasts between two grids, an average
marginal effect with a correct standard error, a simultaneous band, a
derivative, an implicit-function feature. All of them are `g(A c)` for
some `g`. The extension can build `A` today for a predictor that is
LINEAR in `(beta, b)`, by unit perturbation, and it cannot build it at
all for the two cases that matter most to Part 2: a nonlinear (`nl`)
body, where `A` is a Jacobian rather than a design, and a reduced-rank
block at `re.form = NULL`, where the loadings live in `theta` and the
perturbation cannot see them. Those two are why `frm_lp_basis()` is
worth exporting.

**What it must say.**

* `emmeans`: nothing changes, but `emm_basis()` should eventually be
  written in terms of it so that a smooth term stops being invisible to
  `emmeans` on a frmtmb fit.
* `simulate()`, `importance`: nothing. This is a post-fit accessor.
* The compat table gains one `method` feature, `frm_lp_basis`, whose
  rows against `s()`, `t2()`, `gp()`, `rr` and `nl` are exactly the rows
  `frmtmb.spline` registers for `frm_curve` today, because they are the
  same rows. Two of them are worth naming. `rr` is `conditional`, not
  `refused`: a reduced-rank block's loadings live in `theta`, so `A`
  over `(beta, b)` is well defined and INCOMPLETE, and core's own delta
  method carries the missing piece separately as `rr_jacobians()`; an
  exported `frm_lp_basis()` should return that piece rather than leave a
  caller to discover it is absent. `gp` is the other: an exact `gp()`
  contributes a kriging variance that is not coefficient uncertainty at
  all, which is why `extra_var` is a separate element rather than folded
  into `A V A'`.

#### Part 1c. An `lccdf` slot, found on the way

A family can declare `lpdf` and `lcdf`, and core forms the
right-censored contribution as `log(Fub - F(y))`, which without
truncation is `log(1 - F)`. So a family can never contribute an accurate
`log S` once `F` rounds to 1, and `log(0)` differentiates to `NaN`.

Measured: on flexsurv's `bc` data, a Royston-Parmar slope of 30 in log
time puts 572 of 686 rows at `F == 1` exactly. `flexsurv` computes
`log(S)` directly and has no such region. `frmtmb.spline` works around it
by squeezing `F` into `(0, 1)` with a smooth positive part, which is
exact in double precision wherever the survival probability exceeds
6.7e-08 and returns 5e-16 rather than 0 below it. That is a workaround,
not a fix: the family is being asked for a number whose complement
cannot be represented.

**Proposal:** an optional `lccdf` slot, the LOG survival function. When
a family declares it, `R/objective.R:100` uses `lccdf(y)` for a
right-censored row instead of `log(Fub - Fv)`, and falls back to the
present form when it does not.

**It is five sites, not one.** The dispatch is one line; the rest is:

1. a `custom_family()` / `frmtmb_family()` argument and its validation;
2. `frmtmb_ad_overload()` wrapping alongside `lcdf`'s at
   `R/families.R:259`;
3. an arity shim beside `fam_lcdf()` (`R/families.R:3898-3903`), which
   already carries one for the four-argument categorical form;
4. the use site at `R/objective.R:100`;
5. the frame guard at `R/frame.R:1249`, which refuses `cens()` when
   `lcdf` is absent and must learn that `lccdf` alone is also enough.

**And it does not close the class.** It fixes right censoring. Left
censoring keeps `log(Fv - Flb)`, and truncation divides by the window at
`R/objective.R:112` (`ll <- ll - log(Fub - Flb)`), so a LEFT-TRUNCATED
survival model, which is what delayed entry is and which is routine,
meets the identical representability problem from the other side.
Closing that needs a log CDF as well, or a windowed log-difference slot.
The proposal here is the right first step and should not be sold as the
last one.

**And it is not the only way to close the class.** What
`frmtmb.spline` actually needs is for a floored likelihood to be
*visible*, and an `lccdf` slot achieves that by removing the floor. A
POST-FIT FAMILY HOOK would achieve it a different way, by letting a
family inspect its own fit and refuse. Core has no such hook today: the
family carries `post$mean_fn`, `post$var_fn` and `post$dev_fn`, none of
which core calls when a fit finishes, and `logLik()` reads
`object$opt$objective` directly (`R/methods-fit.R:233-240`). So
`logLik()` and `AIC()` on a `royston_parmar()` fit report whatever the
optimizer reached, and **no extension can gate them**. That is the one
limit `frmtmb.spline` documents rather than fixes; it refuses at every
entry point it owns and cannot reach those two. Either seam closes it,
and the hook would close it for every family rather than for survival
alone.

### Two limits recorded rather than proposed

Neither is a request. They are written down so that the next reader does
not spend the afternoon rediscovering them.

* **`logLik()` and `AIC()` cannot be gated from an extension.** See
  immediately above. `frmtmb.spline` refuses through `rp_floored()` and
  through all three curve functions, and a user who fits and reads
  `AIC()` without calling either still gets a floored number with no
  signal. Tracked against Part 1c.
* **A mapped random-effect block has no supported route to test.**
  `frmtmb_control()` takes no `map` argument, so nothing in an extension
  can construct a fit with a mapped `b` block and the behaviour of the
  cached joint covariance there is unmeasured. The nearest reachable
  analogue, a distributional parameter held fixed
  (`bf(y ~ s(x, k = 8), sigma = 0.5)`, which sets `betad_fixed_idx` and
  takes the same index-remapping path in `sp_coef_pos()`), works at a
  measured `cov_rel_error` of 1.55e-15 over 8 `predict()` calls. If core
  ever exposes `map`, that is the case to measure.

---

## Part 2. A penalized coefficient block a nonlinear body can consume

This is the real proposal, and this lane must NOT implement it.

### The shape of the problem

Every penalized smooth in core today ends up MULTIPLIED BY `Z`. The
frame builds `mgcv::smoothCon()` then `mgcv::smooth2random()`, splits
the basis into a null space that joins `X` and a wiggly part that
becomes a random-effect block with one variance per penalty, and the
objective adds `Z b` to a linear predictor. That is the whole
integration, and it is exactly right for `s(x)` and exactly useless for
the three consumers below, all of which need the spline's VALUE inside
an expression rather than its contribution added to a predictor.

### What core would export

One term, usable in a nonlinear body and in a family:

```r
ps(expr, k = 10, degree = 3, penalty = 2, by = NULL, id = NULL)
```

`ps()` in a nonlinear body declares a penalized coefficient block and
evaluates to the spline's VALUE at `expr`, on the tape, as an ordinary
quantity the body may then multiply, exponentiate or feed to a solver.
It is not added to anything.

At frame time `ps()` contributes:

| piece | where it goes | why |
|---|---|---|
| knot vector, degree | frozen on the frame | a basis must be the same on `newdata`; this is `poly()`'s frozen-basis rule |
| null-space columns | `beta`, as today | unpenalized directions are fixed effects |
| penalized columns | `b`, one block, one `theta` | the smoothing parameter IS `sigma^2 / var`, which is what makes core reproduce mgcv's `method = "ML"` |
| a tape closure | the per-call evaluation frame `ev` | this is the new part; NOT `nl_env`, see below |

The closure is what does not exist today. It must evaluate the basis at
`expr` ON THE TAPE, because in consumer (i) `expr` is itself a function
of parameters.

**And it has to live in `ev`, not in `nl_env`.** An earlier draft of this
table put it in `nl_env`, which is wrong for core as it stands.
`nl_env` is fixed at parse time (`R/parse.R:1357`,
`environment(nf) %||% env`), before any coefficient exists, so a closure
living there would capture a stale environment and could never see the
current draw of the penalized block. The body is evaluated at
`R/objective.R:286-293`:

    ev <- c(dparv[[lp[["resp"]]]][c(lp[["nl_pars"]], lp[["nl_dpar_refs"]])],
            lp[["data_list"]])
    eta <- tryCatch(eval(lp[["nl_body"]], ev,
                    ad_overload_env(lp[["nl_env"]], lp[["nl_body"]])), ...)

`ev` is rebuilt per objective call and already holds advectors, so a
`ps` closure placed there is found by ordinary R lookup, receives its
argument already evaluated as an advector, and returns one. **No second
evaluation path is needed.** One word in the table, and it is the
difference between a design that works and one that captures a stale
environment.

**The parse is the harder half, and it is where the work is.**
`nl_dpar()` (`R/parse.R:1246-1252`) learns everything it knows about a
body from one line:

    vars <- all.vars(body)          # R/parse.R:1247

and derives `nl_pars` and `datavars` from it. It never inspects call
structure, so core has no place today that would notice a `ps()` call at
all, let alone allocate its knot vector, its null-space columns and its
`theta`. That is a new parse pass over the body, and it is the same line
Part 3 below is about: a collector that cannot see past `all.vars()`
cannot see `ps()` either. The two are one fix at the same site.

### The basis on the tape is not a problem, and here is the measurement

A B-spline basis is a divided difference of truncated powers:

```
B_{i,k}(x) = (t_{i+k} - t_i) * sum_{j=i}^{i+k} (t_j - x)_+^{k-1} /
                                prod_{m != j} (t_j - t_m)
```

and the truncated power is spelled branch free as

```r
pos <- function(e) 0.5 * (e + abs(e))
```

`abs()` works on an RTMB advector. `pmax()` does not: RTMB refuses it
with "Comparison is generally unsafe for AD types", and RTMB exports no
`CondExp`. So this construction is not a convenience, it is the only one
available.

**Measured** (cubic, `ord = 4`, 8 interior knots on [0, 1.5], 37
evaluation points, against `splines::splineDesign`):

    value, no AD                      max abs error   6.82e-14
    value, no AD                      max rel error   3.64e-11
    value, x as an AD parameter       max abs error   4.64e-14
    d/dx, reverse mode, against
      splineDesign(derivs = 1)        max abs error   1.52e-13

The relative figure is larger than the absolute one only because basis
values near a boundary knot are themselves near zero. **No tape
comparison is needed anywhere**, and the derivative the tape produces is
`splineDesign`'s own derivative to 1.5e-13. The recursive Cox-de Boor
form, which does need comparisons, is not required.

One caveat the measurement also shows: the divided-difference form needs
DISTINCT knots. Repeated boundary knots, which `splineDesign` normally
uses, make the divided difference degenerate. The fix is the usual one,
place distinct boundary knots outside the data range, and `ps()` should
do that rather than expose it.

### The three consumers

#### (i) A warped population curve in an `nl = TRUE` body

The shape-invariant model of D'Alessandro, Thoresen and Sorensen (2026,
arXiv 2603.11728). One population curve `f`, and each subject sees it
shifted and scaled in both axes:

```r
frm(bf(y ~ a * ps(b * (t - c), k = 12) + d,
       a + b + c + d ~ 1 + (1 | subject), nl = TRUE),
    data = dat)
```

`ps()`'s argument contains `b` and `c`, which are parameters, so the
basis is evaluated at a point the tape controls. This is the consumer
that forces the whole design: nothing about `Z b` can express it,
because there is no fixed `Z`.

(The reference implementation, github.com/matteodales/snmmTMB, is
unlicensed. It was not read and nothing from it is here.)

**What the rest of core must say.**

* `predict(newdata = )`: works, and needs nothing new. The basis is
  frozen on the frame and the body is re-evaluated at the new `t`, which
  is what the nonlinear path already does.
* `predict(se.fit = TRUE)`: **refused**, and it is already refused,
  at TWO sites rather than one. The response-scale refusal is
  `R/predict.R:1558` ("se.fit is not supported on the response scale for
  a nonlinear predictor yet"); the link-scale one is
  **`R/predict.R:1179`** ("se.fit is not supported for the nonlinear
  predictor yet"), and that second one is the one a nonlinear
  `frm_lp_basis()` would have to lift, because the link scale is where a
  curve is read. Both must stay refused until it grows a nonlinear
  branch, and the reason belongs in the message: `A` is `d eta / d c`
  and for a nonlinear body that is a Jacobian, not a design matrix. It
  is computable, by taping the body's derivative with respect to the
  coefficients, and it is a bigger job than Part 1.
* `frm_curve()`: refused today by its own linearity probe, correctly.
  With a nonlinear `A` it would work unchanged, since everything after
  `A` and `V` is the same arithmetic.
* `simulate()`: works. The body is evaluated at drawn coefficients like
  any other nonlinear body.
* `importance`: the importance correction re-evaluates the objective at
  drawn random effects (`R/importance.R`). A `ps()` block's coefficients
  are random effects with one variance, so they are drawn like any
  other block. The correction is valid and should be marked `works`
  ONLY after it is run: the warping parameters `b` and `c` are also
  random effects here, and the Laplace approximation is being asked
  about a curve that moves.
* Compat: `ps` is a new `special`, and the rows that matter are
  `ps x nl` (`works`, and the reason `ps` exists), `ps x se.fit`
  (`refused`, with the nonlinear-predictor reason), `ps x predict`
  (`works`), `ps x REML` (`untested`), `ps x quadrature` (`refused` or
  `untested`: marginalizing a `ps` block by Gauss-Kronrod means
  integrating over a curve, and the node count for a 12-coefficient
  block is not a thing anyone should reach for by accident).

#### (ii) Shape-constrained smooths through exponentiated coefficients

A monotone increasing curve is a spline whose coefficient DIFFERENCES
are positive, which is Pya and Wood's construction: the free parameters
are unconstrained, the coefficients are their exponentials cumulated,
and the basis is untouched.

```r
frm(bf(y ~ ps(x, k = 12, shape = "increasing")), data = dat)
```

The coefficients now enter through `exp()` and `cumsum()`, so `eta` is
NOT linear in them and `Z b` is not what the objective computes. The
penalized block is the same block; only the map from block to
coefficients changes. That is why the two consumers belong to one
proposal: a `ps()` block should carry an optional coefficient map, with
the identity as the default, and the `shape` argument is that map.

**What the rest of core must say.**

* `predict()`: works, once the coefficient map is applied where the
  basis meets the coefficients rather than at fit time only.
* `predict(se.fit = TRUE)`: the delta method needs `d eta / d(free
  parameters)`, which is the basis times the Jacobian of the map. That
  Jacobian is closed form for the exponential-cumsum map and is a small
  amount of code. So this consumer CAN have standard errors where (i)
  cannot, and the answer should not be "refused" just because (i) is.
* `frm_curve()`: works through the same Jacobian.
* `VarCorr()`: reports the block's variance as usual, and it should keep
  meaning the same thing: the penalty is still on the FREE parameters,
  so the variance is still `sigma^2 / lambda`.
* `simulate()`: works.
* Compat: `ps(shape =) x se.fit` is `works` where `ps() x se.fit` inside
  `nl` is `refused`, and the table should carry both rows rather than
  one blurred one. A single `ps` row would be a lie in one direction or
  the other.

#### (iii) Spline-valued time-varying inputs for `frmtmb.ode`

`extensions/frmtmb.ode/R/ode.R:512-549` builds a time-varying input as a
PIECEWISE-CONSTANT function of time: the input's value is looked up in a
step function at each solver step. That is the right first
implementation and it is wrong for two things a pharmacokinetic or
motor-control model normally has. A step input has no derivative, so a
stiff solver's Jacobian is discontinuous at every knot; and a smoothly
varying input estimated from data is a smooth, not a staircase.

```r
frm_ode(..., input = ps(time, k = 20))
```

The ODE right-hand side evaluates the input at the solver's OWN time
points, which are chosen adaptively and are not data. So the basis must
be evaluable at an arbitrary tape-valued time, which is exactly the
capability (i) needs, and is why this consumer is in the same document.

**What the rest of core must say.**

* `predict()` and `simulate()`: work, because the solver is re-run.
* `se.fit`: refused for the same reason as (i), and the ODE extension
  already refuses it.
* The ODE extension's own compat rows for `input` gain a `ps` variant.
  The piecewise-constant one stays: it is cheaper, it is the right
  choice for an infusion that really is a step, and removing it would be
  a regression.
* One thing the ODE lane must decide and this document cannot: whether
  a `ps()` input's coefficients are estimated jointly with the ODE
  parameters (they are random effects, so yes) and what that does to the
  solver's step count. Measure before shipping.

### Cost, honestly

Part 1 is a return statement plus its documentation. Part 2 is a new
term in the parser, a frozen basis on the frame, a coefficient map, a
tape closure, and a nonlinear branch in the delta method. It should not
be one commit and it should not be one lane.

---

## Part 3. A parser defect found on the way (a CORE item)

This one is not about splines and does not wait on anything in Part 1 or
Part 2. It is a defect in core's nonlinear-body variable collector, it
has a one-line site and a small fix, and it should be filed and fixed on
its own. It is repeated here only because Part 2 runs into the same line
from the other direction: a collector that cannot see past `all.vars()`
cannot see a `ps()` call either.

### What it is

A call in FUNCTION POSITION inside a nonlinear body hides its own
arguments from the variable collector.

```r
all.vars(quote(f(x)(y)))              # "y"          -- x is gone
all.vars(quote(a * curry(tv)(zv)))    # "a" "zv"     -- tv is gone
all.names(quote(f(x)(y)))             # "f" "x" "y"  -- all.names sees it
```

`all.vars()` treats the function-position subtree as a function name and
drops it whole, arguments included. Core collects a nonlinear body's
data variables with `all.vars()`, so a variable that appears ONLY inside
a function-position call is never collected, never put in the body's
evaluation data, and the body then fails on it.

### Reproduced on 5dfdd84

```r
curry  <- function(u) function(v) u * v
curry2 <- function(u, v) u * v
d <- data.frame(tv = runif(60), zv = runif(60))
d$y <- 2 * d$tv * d$zv + rnorm(60, 0, 0.05)

frm(bf(y ~ a * curry2(tv, zv), a ~ 1, nl = TRUE), gaussian(), d)
#> fitted, a = 2.006

frm(bf(y ~ a * curry(tv)(zv), a ~ 1, nl = TRUE), gaussian(), d)
#> Error: The nonlinear formula body could not be evaluated:
#>        object 'tv' not found

frm(bf(y ~ a * curry(tv)(zv) + 0 * tv, a ~ 1, nl = TRUE), gaussian(), d)
#> fitted, a = 2.006
```

The third line is the diagnosis: naming `tv` anywhere `all.vars()` can
see it makes the identical model fit, to the same coefficient.

### Why it matters here

A curried spline is the natural spelling for a basis that is built once
and then evaluated: `ps(t, k = 12)(theta)` reads better than
`ps_eval(t, k = 12, theta)` and is what someone writing consumer (i)
would try first. The defect is not hypothetical for this proposal, it is
in its path.

### The fix

Replace `all.vars(body)` with a walker that descends into a
function-position call. Minimally:

```r
body_vars <- function(e) {
  if (!is.call(e)) return(all.vars(e))
  v <- character(0)
  if (is.call(e[[1L]])) v <- body_vars(e[[1L]])
  for (i in seq_along(e)[-1L]) v <- c(v, body_vars(e[[i]]))
  unique(v)
}
```

`body_vars(quote(a * curry(tv)(zv)))` gives `c("a", "tv", "zv")`, and
`body_vars()` agrees with `all.vars()` on every body that has no call in
function position, which is every body core is tested on today. The
change belongs wherever the nonlinear body's data variables are
collected, and it wants one test asserting the three lines above.

The refusal message is worth a second look while that is done. "object
'tv' not found" is what R said, not what happened; "`tv` is used in the
model formula but was not collected from `data`" would be nearer, and
the current wording sends a reader to look for a typo in their data
frame rather than at the shape of their expression.
