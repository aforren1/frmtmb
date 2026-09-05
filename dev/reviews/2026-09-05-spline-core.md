# Review: SPLINE-CORE lane

Reviewer's worktree read: `C:/Users/adf44/source/r/frmtmb-wt-spline-core`,
branch `wt-spline-core`, base 2210aa1 (0.51.0). Main at 316a28b, clean and
untouched by this review.

Private library
`scratchpad/rsc-lib`: frmtmb 0.52.0 (worktree), frmtmb.spline 0.2.0,
frmtmb.ddm 0.2.0, frmtmb.ode 0.1.0, frmtmb.sample 0.1.0, flexsurv 2.3.2,
gratia 0.11.2, brokenstick 2.7.0. R 4.6.1.

Everything below is reproduced by the reviewer unless it says otherwise.

---

## 1. The identity, ported independently from the paper

The paper is real and is what the lane cites: arXiv:2603.11728,
D'Alessandro, Thoresen and Sorensen, *A Semiparametric Nonlinear Mixed
Effects Model with Penalized Splines Using Automatic Differentiation*.
Fetched the HTML, extracted the LaTeX from the MathML annotations, and
worked from equations (3), (4) and section 2.3.1. `snmmTMB` was not read
and the lane's own reference port was read only after mine ran.

My port differs from the lane's in what it takes from the fit. The lane's
vignette reference reads `pt$U0`, `pt$Us` and `pt$knots` off the fitted
object and reassembles `dnorm()` terms. Mine rebuilds the knot vector from
`ps()`'s documented rule, rebuilds the second-difference penalty `S`, the
sum-to-zero constraint and the Wood (2004) eigensplit from scratch, and
writes equation (4) in the paper's own form, penalty included as
`lambda * theta' S theta` rather than as `||omega||^2`.

Three things fall out, all independent of the lane:

| check | result |
|---|---|
| my knot vector against `pt$knots` | max abs diff **0** |
| internal knots strictly inside the padded range | **11**, the paper's number for the SMOCC basis |
| rank of the constrained penalty | 13 range + 1 null = 14 = 15 - 1, matching `n_pen` 13 and `n_fixed` 1 |
| `theta' S theta - ||omega||^2`, my `S` and my constraint | **4.135e-12** |

The last row is the one that matters: it says the lane's eigensplit really
is a valid mixed-model representation of the second-difference penalty
under a sum-to-zero constraint, checked against a penalty matrix the lane
did not supply.

The density itself, `frmtmb:::build_objective(frame)(estimates)` against my
port of equation (4):

| point | frmtmb | my port | abs | rel |
|---|---:|---:|---:|---:|
| at the optimum | 3653.3851349042 | 3653.3851349043 | **9.4133e-11** | 2.5766e-14 |
| at a perturbed vector (mine, not the lane's) | 3807.2227442180 | 3807.2227442184 | 3.6380e-10 | 9.5555e-14 |

The optimum figure reproduces the lane's claim to every digit it printed
(9.413e-11, 2.577e-14). CONFIRMED.

### Table 1 and the beta0 offset

My fit of the lane's spelling, against the published table and its
intervals:

| parameter | frmtmb | paper | 95% interval | inside |
|---|---:|---:|---|:--:|
| beta0 | 68.7543 | 68.20 | 67.60-68.70 | **no** |
| beta1 | 1.8058 | 1.80 | 0.92-2.68 | yes |
| beta2 | 0.0002 | 0.00 | -0.01-0.01 | yes |
| beta3 | 1.0500 | 1.00 | 0.73-1.27 | yes |
| sd_b1 | 2.8685 | 2.86 | 2.56-3.17 | yes |
| sd_b2 | 3.3241 | 3.28 | 2.89-3.68 | yes |
| sigma | 1.0567 | 1.05 | 1.01-1.09 | yes |

Six of seven, as claimed. On the seventh the lane says the offset is a
centering convention. I can now put a number on that rather than an
argument:

```
mean over the data of the fitted curve   -0.556096
beta0 + mean(fitted curve)                68.1982
paper's beta0                             68.2
beta0(frmtmb) - 68.2                       0.554333
```

Re-centering the lane's fit to the convention "the fitted curve averages to
zero over the data" lands the intercept on the paper's 68.2 to **0.002 cm**.
That is not an argument, it is an identity, and it settles the bullet.

On my own reading of the paper's convention: the paper says only
"subject to sum-to-zero constraint" (section 4, the SMOCC basis) and, in
the discussion, "a sum-to-zero constraint to ensure identifiability of an
intercept term". It never says *on what*. The lane's finding doc asserts
the paper constrains the fitted basis; the paper's text does not say so.
But the standard meaning of "sum-to-zero constraint" in exactly this
literature (Wood, whose 2004 paper the reparameterization is taken from)
is the constraint `sum_i f(x_i) = 0` on the FITTED VALUES, and the 0.002 cm
residual above is decisive. The lane's explanation is right; its claim to
be reading it off the paper's text is stronger than the text supports. A
one-clause softening is on the punch list, nothing more.

### One difference from the paper the lane does not record

Section 2.3.1 of the paper places knots by a rule `ps()` does not
implement. For a spline argument whose domain is not bounded a priori the
paper rescales the argument to `[0, 1]` using bounds
`min(t) - 3*sigma_b`, `max(t) + 3*sigma_b` that are smooth functions of
the CURRENT variance estimates, so the paper's knots move with the fit.
`ps()` freezes knots on a padded data range with every nonlinear parameter
at zero, and documents that it does. Both give 11 internal knots and a
cubic basis; they are not the same basis, and `ps()` additionally uses
distinct boundary knots outside the range where the standard construction
repeats them.

This does not touch the identity check, which is against the lane's own
model. It does mean "recovers the paper's Table 1" is empirical agreement
between two related models rather than a reimplementation of one. The
agreement of six parameters plus a beta0 offset that resolves to 0.002 cm
is strong evidence they are equivalent in practice. Worth a sentence in
`?ps` or the vignette; the current text implies the knot rule is the
paper's.

---

## 2. `ps()` beyond the identity

### Gradient of the objective through the taped basis

`frmtmb:::build_objective(frame)` flattened to a vector, taped with
`RTMB::MakeTape()`, against `numDeriv::grad()` at two points on a
`y ~ lev + ps(t + shift, k = 10, pad = 0.3)` fit, p = 43.

| point | max abs AD - numDeriv | max rel | worst component |
|---|---:|---:|---|
| optimum | 4.62e-06 | 4.62e-06 | a `b` entry whose gradient is 4.6e-06, i.e. zero |
| perturbed | 1.69e-05 | **2.06e-08** | a `b` entry whose gradient is 819.27 |

The absolute figures are numDeriv, not the tape. Two controls say so.
Tightening numDeriv's step (`eps = d = 1e-6, r = 6`) makes it WORSE by two
orders (max |AD - numDeriv| 9.45e-03), which is step-size sensitivity in
the difference and not in the derivative. And a hand-placed central
difference at `h = 1e-5` on the worst component gives 819.26970029 against
the tape's 819.26970114, a disagreement of **8.5e-10** on a number of
magnitude 819, relative 1e-12. The taped divided-difference basis
differentiates correctly. PASS.

### A warp of exactly zero against a plain smooth

Does NOT reproduce `s()`'s logLik to 1e-8, and cannot. 200 unique
covariate values, `y = 2 + sin(2 pi t) + N(0, 0.15)`:

| model | logLik |
|---|---:|
| `bf(y ~ lev + ps(t, k = 12, pad = 0), lev ~ 1, nl = TRUE)` | 80.7565339765 |
| `y ~ s(t, bs = "ps", k = 12)` | 82.0567625064 |
| `mgcv::gam(y ~ s(t, bs = "ps", k = 12), method = "ML")` | 100.4151551725 |

The reason is the basis, and it is not hidden. `ps()` lays `k + ord`
UNIFORM knots with three outside each end, so with `pad = 0` its knots
still span [-0.321, 1.326] on data spanning [0.008, 0.997]; `s(bs = "ps")`
puts its knots over the data range and uses the repeated-boundary-knot
convention. Different bases are different models with different marginal
likelihoods, and there is no public way to hand `s()` the knot vector
`ps()` builds. The models agree where it is meaningful to compare them:

```
fitted values, ps(pad = 0) vs s(bs = "ps"):  max abs 1.230e-03   cor 0.999999774
```

The 18.4-nat gap between frmtmb's own `s()` and `mgcv::gam(method = "ML")`
on identical formula and data is a pre-existing cross-package difference in
what `logLik()` counts, not something this lane touched.

Verdict: the claim "a zero warp reproduces `s()`" is NOT establishable to
1e-8, the lane never makes it, and the fitted-curve agreement to 1.2e-03 is
what the shared basis space supports. No defect.

### predict at newdata, inside and outside the range

Inside: correct, and the frozen knots hold. A sub-grid agrees with the
matching rows of a full grid to **0**, and `predict(newdata = the original
data)` agrees with `predict(fit)` to **0**.

Outside: correct arithmetic, SILENT presentation, and this is the one
usability defect I found in `ps()`. On a fit whose `knot_range` is
[-0.331, 1.328] and whose full knot span is [-1.042, 2.039]:

| t | predict | what it is |
|---:|---:|---|
| 0.50 | 1.96172 | the curve |
| 1.35 | 2.98550 | **a partial sum**: past `knot_range`, the basis is no longer a partition of unity |
| 1.55 | 2.83324 | the same, decaying |
| 1.75 | 2.05797 | the same, nearly decayed |
| 3.00 | 1.64484 | exactly `lev`: the curve reads as zero |

No warning at any row. `ps_coverage_warning()` (R/ps.R:521) reports
coverage at FIT END on the FITTED data only; nothing reports it at
`predict(newdata =)`. A user extrapolating a growth curve past the padded
range gets a curve that bends smoothly to the intercept and is told
nothing. See the punch list.

### simulate()

Works, and the block is drawn: three simulated columns have SDs 0.698,
0.683, 0.691 against the data's 0.692.

`simulate(fit, newdata = nd)` returns the FITTED data's row count, not
`nd`'s. That is NOT this lane: `simulate.frmtmb_fit` has no `newdata`
argument at 2210aa1 or at 0.52.0 (R/predict.R:2607, signature
`(object, nsim, seed, re.form, censored, ...)`), so `newdata` lands in
`...` and is dropped on every shape. Reproduced identically on
`y ~ x + (1 | g)`, `y ~ s(x)` and a plain nonlinear body. Out of scope for
this lane; worth a separate issue.

### The importance correction

Refused, exactly as the compat row says, with the inherited reason:

```
`importance` cannot correct a nonlinear predictor (y.mu). A nonlinear body
mixes parameter values with raw data columns, and the corrected objective
evaluates the predictor once per draw ... Use importance = 0
```

But `?ps` says the opposite. R/ps.R:96-99 (and `man/ps.Rd:101`) reads

> The importance correction ([frm()]'s `importance =`) reweights the same
> Laplace integral and draws a `ps()` block like any other random-effect
> block, so it applies. Read its diagnostics...

while `R/compat.R` registers `r("ps()", "importance", "refused", ...)` and
the findings doc says refused. The help page tells the user to do something
that always errors. Punch list.

### Refusals by name

All fire, all name the term, all say why:

| asked | refused |
|---|---|
| `REML = TRUE` | yes, names `ps(t, k = 8)` |
| `quadrature = TRUE` | yes, names the term |
| `frmtmb_control(profile = TRUE)` | yes, names the term |
| `mvbf()` | yes, "not supported in a multivariate model yet" |
| `ps(by = id)` | yes, and says `by =` and `id =` are not silently ignored |
| `ps(k = 51)` | yes, with the measured reason |
| `ps(ps(t)))` | yes |
| `ps()` called as a function | yes |

---

## 3. `predict(se.fit = TRUE)` as a consumer of `frm_lp_basis()`

Ten shapes, fifteen comparisons, the lane's 0.52.0 against main's core
(0.51.0 built from `git archive 316a28b` into a separate library), same
seeds and same data on both sides.

| shape | max abs se change | max rel | max abs fit change |
|---|---:|---:|---:|
| gaussian GLMM, `re.form = NA` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| gaussian GLMM, `re.form = NULL` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| gaussian GLMM at newdata | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| distributional, mu | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| distributional, `dpar = "sigma"` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| ordinal (`cumulative()`) | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| smooth, `re.form = NA` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| smooth, `re.form = NULL` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| gp, in sample | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| gp at newdata | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| rr, `re.form = NA` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| rr, `re.form = NULL` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| mixture | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| mvbf, `resp = "y"` | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| car (`escar`) | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| nonlinear | refused on both, byte-identical message |

**Nothing moved.** Not "below 1e-12": bit-identical, on every shape
including the two the rewrite most plausibly could have disturbed (rr,
whose loadings live in `theta`, and gp, whose kriging variance the new
code returns separately in `extra_var`). The rewrite of
`predict(se.fit = TRUE)` onto `frm_lp_basis()` is behavior-preserving.

This is a stronger result than the lane's own 2.776e-17, which was measured
against the pre-rewrite code inside the lane rather than against main.
CONFIRMED and exceeded.

---

## 4. `lccdf` and the fit-end hook

### The five sites

All present, all where the findings doc says (line numbers as they now
stand):

| # | site | where |
|---|---|---|
| 1 | `frmtmb_family(lccdf =)` and its validation | R/families.R:292, :322-327 |
| 2 | `frmtmb_ad_overload()` wrapping beside `lcdf`'s | R/families.R:327 |
| 3 | `fam_lccdf()` arity shim and `has_lccdf()` | R/families.R:4148, :4160 |
| 4 | the use site in the censored term | R/objective.R:106-120 |
| 5 | the frame guard | R/frame.R:1301-1315 |
| 6 | `post$fit_check`, the hook the proposal did not enumerate | R/fit-end.R:18, called from R/fit.R:999 |

The five built-in families that declare it, read off the installed
package: `fam_gaussian`, `fam_lognormal`, `fam_exponential`,
`fam_weibull`, `fam_cox`. `fam_poisson` has `lcdf` and is type
"discrete"; `fam_inverse_gaussian` has `lcdf` and no `lccdf`;
`fam_Gamma` has neither. Exactly the claim.

### The gaussian table

Reproduced against `stats::pnorm(lower.tail = FALSE, log.p = TRUE)` as
truth, with `log(1 - F)` formed from the family's own `lcdf`:

| z | true | log(1 - F) | lccdf |
|---:|---:|---:|---:|
| 1.0 | -1.84102 | -1.84102 | -1.84102 |
| 5.0 | -15.06500 | -15.06500 | -15.06500 |
| 8.0 | -35.01344 | **-34.94504** | -35.01344 |
| 8.3 | -37.49422 | **-Inf** | -37.49422 |
| 20 | -203.91716 | -Inf | -203.91716 |
| 37 | -689.03059 | -Inf | -689.03059 |
| 100 | -5005.52421 | -Inf | -5005.52421 |
| 500 | -125007.13355 | -Inf | -125007.13355 |

Digit for digit the lane's table. The eight-digit claim at z = 5 also
holds: true -15.0649983939887, log(1 - F) -15.0649983938340, lccdf
-15.0649983939887.

Weibull at shape = 3: `lccdf` gives -45.57267 at q = 4 where log(1 - F)
is already -Inf, and -364.58135 at q = 8. The lane's figures.

The gradient in eta, `lccdf` taped with RTMB against a central difference
of the truth, and log(1 - F) central-differenced:

| z | true | log(1 - F) | lccdf (AD) |
|---:|---:|---:|---:|
| 5 | 5.18650 | 5.18663 | 5.18650 |
| 8 | 8.12137 | **0** | 8.12137 |
| 10 | 10.0981 | **NaN** | 10.0981 |
| 20 | 20.0498 | NaN | 20.0498 |
| 40 | 40.0250 | NaN | 40.0250 |

My log(1 - F) column differs from the lane's in exactly how it dies (0
and NaN where they reported 7.58 and Inf) because mine is a finite
difference of the dead value and theirs is the objective's own AD
gradient. Same conclusion reached twice: the old form's derivative is
unusable from z = 8 and the new one is exact.

### A family without `lccdf` still takes the old path

Two custom families, identical except for the slot:

```
lcdf only  logLik -210.2879925060
+ lccdf    logLik -210.2879925060
difference 0.000e+00 ; max coef difference 2.220e-16
```

The slot is inert where the tail is representable. That is the
regression-safety result and it is exact.

The guard's refusals all fire by name: an `lccdf`-only family fits right
censoring (same logLik to the digit) and is refused for LEFT censoring
and for `trunc(ub =)`; a family with neither is refused.

### What the slot buys, measured on a fit

A gaussian-shaped custom family, one row censored at 40, about 77 SD
above the fitted mean:

```
lcdf only : ERROR: The optimizer failed on this model
            (NA/NaN gradient evaluation)
+ lccdf   : logLik -448.111264  intercept 1.404001  sigma 3.326955
            at those coefficients: exact         -448.111264
                                   log(1-F) form -Inf
            |exact - reported| 2.842e-13
```

### The fit-end hook

Probed with a custom family whose `fit_check` records what it is given.

**Called with** `(fit, resp)`: the finished `frmtmb_fit` and the response
NAME as a string. Runs once per response, from `fit_end_checks()`
(R/fit-end.R:18), which `fit_assembled()` calls at R/fit.R:999 after
`check_convergence()`.

**What a family can do in it**: it sees `fit$estimates`, `fit$frame` and
`fit$opt`; `logLik(fit)` works inside it and so does `vcov(fit)`
(returned a 3 x 3). It does NOT see a populated `fit$cache$sdr` -
`names(fit$cache)` was empty in my probe - so a check that wants standard
errors pays for them itself. A warning raised in the hook reaches the
caller and the fit is still returned.

**It can also stop, and there is no guard.** `fit_end_checks()` calls
`fc(fit, resp)` with no `tryCatch`. A family whose `fit_check` throws,
deliberately or through a bug in the extension, destroys the whole fit
after `frm()` has done all the work. `royston_parmar()` is careful (it
wraps its own `rp_floored()` call in `try()` and only warns), but nothing
in core makes that the rule. See the punch list.

### `rp_floored()` is a diagnostic now

Confirmed in the source and on a fit. `rp_floored()`
(`extensions/frmtmb.spline/R/rp-check.R:194-199`) returns early for
`action = "report"` and then stops ONLY on `mono_rows`; the censored
count is computed, returned in `n_censored_floored`, and never refuses.

Reproduced on the review's 600-subject design and on a second seed of my
own, with an exact log-likelihood I wrote from the hazard-scale
definition (log h - H for an event, -H for a censor) rather than from the
lane's helper:

| seed | reported logLik | my exact | abs diff | max -log S | rp_floored() |
|---|---:|---:|---:|---:|---|
| 20260905 (the lane's) | -575.5379429 | -575.5379429 | **4.206e-12** | 55.7302 | did not refuse |
| 777 (mine) | -588.4725370 | -588.4725370 | **3.070e-12** | 55.0779 | did not refuse |

The lane claims 1.93e-12 on the first; I get 4.2e-12 with an
independently written reference, which is the same result. The max -log S
of 55.73 is the lane's 55.73. On the lane's seed the fit does warn "false
convergence (8)", the honest note their findings make; on my seed it
converges quietly, so that warning is a property of that design rather
than of the change.

I could not force a non-monotone fit in six designs, so I did not
independently exercise the monotonicity refusal. The lane's
`test-rp-floored.R` covers it and passes (section 8).

---

## 5. The phantom beta

Reproduced on main's core (0.51.0, built from `git archive 316a28b` into
a separate library) and on the lane's, same data, same seed. Witness: two
nonlinear parameters that are purely random effects,
`bf(y ~ (1 + a) * sin(...), a ~ 0 + (1 | id), p ~ 0 + (1 | id), nl = TRUE)`.

The brief says 0.50.0, where the lane first found it. I used 0.51.0
instead, for two reasons: it is the version the lane branched from, so it
is the comparison that decides whether the fix is a regression against
what the trunk has TODAY; and building it from `git archive` leaves the
main checkout untouched, where using the 0.50.0 in the user library would
have meant reasoning from a build I did not control. The defect is
present in both - `R/frame.R`'s two coefficient-naming lines are
unchanged between 0.50.0 and 0.51.0 - and 0.51.0 is the stronger test.

| | main 0.51.0 | lane 0.52.0 |
|---|---|---|
| `beta` names | **`a_`, `p_`** | none |
| `beta` length | 2 | 0 |
| `beta` values | 0, 0 | - |
| **logLik** | **-155.700739254** | **-155.700739254** |
| `theta` | -0.7067283, -1.2583186 | -0.7067283, -1.2583186 |
| `betad` | -1.2230762 | -1.2230762 |
| `vcov()` | 3 x 3, **contains NaN** | 1 x 1, no NaN |
| warning | "Some standard errors are not finite ... probably overparameterized" | none |
| Hessian, 3 smallest eigenvalues | **-7.816e-14, -9.663e-25**, 43.97 | 43.97, 49.04, 504.62 |

Everything the lane claims, and the eigenvalue row is the proof: the two
singular directions on main are exactly the two phantoms, and on the lane
the smallest eigenvalue IS main's third. The fit does not move; only the
two coordinates that never entered the likelihood are gone.

### The BEHAVIOR CHANGE wording

The failure mode for a user with a hand-written start vector is LOUD, not
silent, which is the important part:

```
0.51.0:  start$beta length 3 -> fits
         start$beta length 1 -> ERROR: start$beta must have length 3
0.52.0:  start$beta length 3 -> ERROR: start$beta must have length 1
         start$beta length 1 -> fits
```

Same logLik either way (-155.69131306). No silent misalignment, and the
error names the length it wants.

But NEWS does not say it. `NEWS.md:93-102` files this under `FIX`, not
under the `BEHAVIOR CHANGE` label the same file gives the 0.51.0 `mo()`
change, and its closing sentence -

> Affected models fitted to the same optimum; only the inference was lost.

- describes what was WRONG BEFORE, not what CHANGES NOW. A reader
scanning NEWS for breaking changes will not find this one, even though
the lane had to edit two of its own tests for it and the findings doc
states it plainly. Punch list.

---

## 6. Merges

### The textual merge

`git merge-file -p ours base theirs` for every file this lane shares with
a sibling, base 2210aa1, ours = spline-core, theirs = the sibling's
worktree:

| sibling | file | conflicts |
|---|---|---:|
| esicar | R/compat.R | **0** |
| esicar | R/frame.R | **0** |
| esicar | R/objective.R | **0** |
| esicar | NEWS.md | 1 |
| ce | R/families.R | **0** |
| ce | R/predict.R | **0** |
| ce | NEWS.md | 1 |
| priors | NAMESPACE | **0** |
| priors | NEWS.md | 1 |

Every `R/` file and NAMESPACE merge clean. The only conflicts are NEWS.md,
where all three lanes prepend a section to the same place; that is a
one-minute hand merge with no semantics in it.

`R/methods-fit.R` is named in the brief as a ce file. This lane does not
touch it, so it is not shared and there is nothing to merge.

Hunk ranges say why the R files are clean:

* **frame.R** - esicar at base 2090, 2180-2188, 2364-2371; spline-core at
  1298-1309, 1473-1487, 2021-2031, 2264-2269. Closest approach 69 lines.
* **objective.R** - esicar at 194-199 and 249-256; spline-core at 93-103
  and 284-290. No approach closer than 90 lines.
* **predict.R** - the two lanes interleave inside
  `predict.frmtmb_fit()` and still do not collide. See below.

### Does ce's `dpar_report_hook()` / `dpars_natural()` sit inside the region this lane rewrote?

**Adjacent on both sides, and disjoint.** In base coordinates:

* ce inserts the `hook <-` / `hook_val <-` definitions at base 1220-1229,
  immediately after `n <- ed[["n"]]` and BEFORE the `if (!se.fit)` early
  return.
* spline-core's rewrite runs from `has_rr <- ...` to `se_eta <-
  sqrt(var_eta)`, base 1226-1245 - i.e. it begins after that early
  return.
* ce's second change, the `out <- if (!is.null(hook)) ...` block at base
  1244-1251, sits immediately AFTER spline-core's region and CONSUMES
  `se_eta`.

So ce brackets the rewritten block without touching a line of it. That is
why the merge is clean, and the merge is also semantically right:
spline-core changes only HOW `se_eta` is computed, and I measured in
section 3 that the value is bit-identical. ce's hook therefore receives
exactly the `se_eta` it was written against.

`dpars_natural()` itself is new ce code at base 740-746 and 754-767, far
from anything this lane touches.

### Does ce's response-scale delta method go through `frm_lp_basis()` or around it?

**Half through, half around, and the half that goes around is the half
that should.**

* The `se_eta` factor now comes from `lp_basis_out()`, the same helper
  `frm_lp_basis()` returns from. ce's `abs(hook_val("deriv")) * se_eta`
  is therefore a delta method built on top of the seam.
* `hook_val("value")` and `hook_val("deriv")` go around it: they call
  `dpars_natural()`, which reads back through `predict()` recursively.
  That is correct as far as it goes, because `frm_lp_basis()` returns
  `d eta / d coef` for ONE predictor and a mixture's softmax reporting
  transform is a function of SEVERAL predictors.

There is a real gap here, and it is ce's rather than this lane's: ce's
own comment says it applies "the same one-predictor rule the link inverse
gets below", which is exactly the approximation a softmax over several
predictors does not satisfy. After the merge the machinery to do it
properly - a Jacobian over all the contributing predictors, with the
joint covariance at the right rows - EXISTS and is unused. Worth an issue
against ce after consolidation; not a merge blocker and not this lane's
to fix.

### Does esicar's `has_expand` gate a `ps()` block correctly?

**Yes, and this lane anticipated it.**

esicar sets `has_expand = has_rr || has_esicar` (frame.R:2371) and
derives the real gate from the blocks in `frame_needs_expand()`
(covstruct.R:1760). A `ps()` block registers with `covstruct = "smooth"`,
so it is neither `rr` nor esicar and never turns the gate on by itself -
which is right, because a smooth block's coefficient space equals its
parameter space.

When the gate IS on because something else in the model turns it on,
`expand_b()` (covstruct.R:1781) copies every other block through its
`else` branch, `cvec[bk$c_idx] <- b[bk$b_idx]`, so a smooth block passes
through unchanged.

And spline-core reads the right vector: frame.R:2345-2351 sets
`pt$c_idx <- bk$c_idx` with the comment "coefficient space, not parameter
space: an rr block elsewhere in the same model makes the two differ, and
the closure reads the expanded vector", and `ps_env()` (R/ps.R:251)
indexes `bvec[pt$c_idx]`, where `bvec` is exactly what esicar's
`needs_expand` branch produces.

### Proved, not argued: the four-lane tree

I assembled the merged tree (base 2210aa1 + spline-core's working tree +
the merged shared files + esicar's `R/covstruct.R` + ce's
`R/conditional-effects.R` and `R/methods-fit.R` + priors' `R/priors.R`
and `R/simulate-new.R` + the merged NAMESPACE) and installed it. It
builds clean. On it:

| check | result |
|---|---|
| `ps()` + `rr(sp + 0 | id, d = 2)` in one model | fits, logLik 8.86026080, `has_rr TRUE`, **`c_idx != b_idx`** (n_c 156 vs 81 b entries) |
| `ps()` + `car(W, type = "esicar")` in one model | fits, logLik -14.14576494, **`has_expand TRUE`**, `c_idx == b_idx` |
| ce's mixture `theta1` on the response scale, with `se.fit` | 0.333333, **in (0, 1)**, se 0.03849 computed from the rewritten `se_eta` |
| the ten-shape `se.fit` sweep, merged tree vs spline-core alone | **0.000e+00 on all fifteen**, same refusal message |

The first row is the one that matters: it is the case spline-core's
comment anticipated, the expansion is genuinely active (coefficient space
is nearly twice parameter space), and the `ps()` block indexes correctly
through it.

### Merge order

The consolidation merges esicar, priors and ce BEFORE this lane. That
order is the right one and the risk is low, for three reasons I measured
rather than assumed.
First, the merge is textually free. Every `R/` file this lane shares with
esicar or ce merges with zero conflicts, and the closest two independent
hunks come within 69 lines of each other in `frame.R` and interleave
without touching in `predict.R`. Only NEWS.md conflicts, three times,
which is a hand merge with no semantics in it.

Second, the two places the lanes genuinely meet are both safe in the
order given, and I built the merged tree to prove it rather than reading
it off. ce's reporting hook brackets the `se.fit` block this lane
rewrote without touching a line of it, and it consumes `se_eta`, whose
value I measured to be bit-identical before and after the rewrite on
fifteen comparisons across ten model shapes - so ce lands on the same
numbers whether it goes first or second. esicar's `has_expand` correctly
leaves a `ps()` block alone (a smooth block is neither `rr` nor esicar)
and correctly carries it through `expand_b()`'s pass-through branch when
something else in the model turns the gate on; this lane already indexes
`c_idx` rather than `b_idx` for exactly that case, and on the merged tree
a model with both a `ps()` term and an `rr` block fits with coefficient
space nearly twice parameter space.

Third, the order esicar-priors-ce-then-spline is the one that puts the
least at risk, and it is worth keeping for a reason beyond convenience:
this lane is the only one that bumps the version and the only one whose
extension (`frmtmb.spline` 0.2.0, `Depends: frmtmb (>= 0.52.0)`) has a
floor that the bump has to satisfy, so having it land last means the
floor is never briefly wrong in the trunk. The one thing to do rather
than assume is to re-run the ten-shape `se.fit` sweep after ce lands and
before this lane goes on top: ce's `dpar_report_hook()` is the only
sibling change that reads a quantity this lane recomputes, and a
zero-difference sweep is a thirty-second confirmation that the seam
still hands it what it expects. Nothing I measured suggests it will not.
---

## 7. Versions

| lane | DESCRIPTION Version | changed vs 2210aa1 |
|---|---|---|
| **spline-core** | **0.52.0** | yes |
| esicar | 0.51.0 | no |
| ce | 0.51.0 | no |
| priors | 0.51.0 | no |
| rdm-gng | 0.51.0 | no |
| main (316a28b) | 0.51.0 | - |

No collision: spline-core is the only lane that bumps it, and no sibling
touches DESCRIPTION at all.

`frmtmb.spline` 0.2.0 declares `Depends: frmtmb (>= 0.52.0)`, which
matches the bump exactly, and drops `Matrix` from Imports as claimed. The
other extensions' floors are untouched and all below 0.52.0 (ddm >=
0.49.0, ode >= 0.46.0, latent >= 0.48.0, sample declares none), so none
is broken by the bump.

One nit: core's `Suggests` gains `brokenstick` inserted between `ape` and
`bayesplot`, out of the alphabetical order the rest of the list keeps.

## Part 3, the parser, verified separately

Not on the brief's numbered list, but it is the lane's Part 3 and it is
cheap to check, so I did.

`nl_body_vars()` against `all.vars()` on twelve expression shapes, run on
the installed 0.52.0:

| expression | `all.vars()` | `nl_body_vars()` |
|---|---|---|
| `f(x)(y)` | y | **x,y** |
| `a * curry(tv)(zv)` | a,zv | **a,tv,zv** |
| `(g(h(p))(q))(r)` | r | **p,q,r** |
| `d$b` | d,b | d,b |
| `d@b` | d,b | d,b |
| `stats::rnorm(k)` | k | k |
| `stats:::rnorm(k)` | k | k |
| `function(u) u + w` | u,w | u,w |
| `m[, 1]` | m | m |
| `m[i, ]` | m,i | m,i |
| `a + b[[j]]` | a,b,j | a,b,j |
| `lst$el$deep` | lst,el,deep | lst,el,deep |

The walker differs from `all.vars()` in the function-position rows and
NOWHERE else, which is exactly the claim "no body that parsed before
collects a different set". The third row is a case the lane's own write-up
does not mention: nested function position, two levels deep, is also
collected.

The three spellings, main 0.51.0 against lane 0.52.0, same data:

| spelling | main 0.51.0 | lane 0.52.0 |
|---|---|---|
| `a * curry2(tv, zv)` | a = 1.9900537 | a = 1.9900537 |
| `a * curry(tv)(zv)` | **Error: object 'tv' not found** | **a = 1.9900537** |
| `a * curry(tv)(zv) + 0 * tv` | a = 1.9900537 | a = 1.9900537 |

Identical to seven figures across all three, which is the lane's claim
(their number is 2.005982 on their data; mine is 1.9900537 on mine, and
the three-way agreement is the point).

`m[, 1]` in a body fits on both versions (a = 2.9928165 each), so the
empty-argument fault was in the new walker and is fixed, with no
behavior change against the base.

## `frm_lp_basis()`, the three rows the extension adoption rests on

| claim | measured |
|---|---|
| rr carries theta columns at `re.form = NULL` | A is 150x62 with **11 theta columns** (`theta.1`..`theta.11`); at `re.form = NA` there are **0**, correctly |
| the reconstructed se equals `predict(se.fit)` | **0.000e+00** at both `re.form` settings |
| gp returns kriging variance separately in `extra_var` | 1.2e-06 at seen positions, **0.184 at an unseen one**; reconstructed se matches `predict()` to **0.000e+00**, and dropping `extra_var` would put it out by 0.1045 |
| the nonlinear branch's `A` is a Jacobian | 7 x 2, `coef_names` correct, `eta` matches `predict()` to **0**, and `A` matches a central difference of `predict()` to **3.474e-11** |

## The compat table

`ps()` is registered as `kind = "special"` (R/compat.R:535), which makes
it inherit the specials default rule, so `frm_compat()` now carries 68
`works`, 67 `conditional`, 72 `untested` and 8 `refused` rows mentioning
`ps()` or `frm_lp_basis`, of which the lane registered 11 explicitly.

I suspected the derived rows over-claimed and tested it. They do not. The
statuses are right: `student`, `weibull`, `lognormal` and `tweedie`
crossed with `ps()` all fit, and the three that refused
(`Gamma`, `zero_inflated_poisson`, `von_mises`) refused on my response's
shape rather than structurally.

Two notes fall out of that check, neither this lane's fault:

* the derived NOTE on a `<family> x ps()` row reads "Predictor specials
  build design columns before the family sees them", which is the one
  thing `ps()` does not do. Low priority, but it is a sentence attached
  to 68 rows.
* the pre-existing `nl x <family>` rows say **refused** for `student`,
  `weibull`, `lognormal` and `tweedie`, and all four FIT. That is stale
  on main too (I checked 0.51.0), so it is not this lane's, but the
  consolidator should know the table disagrees with itself on
  neighbouring rows.
---

## 8. Runs

Scope taken from `git -C ... diff --name-only 2210aa1` plus untracked: 27
tracked files changed, 11 untracked added (one of which is this review).
No sibling-owned file is touched - `R/priors.R`, `R/sampling-api.R`,
`R/conditional-effects.R`, `R/methods-fit.R` and `R/covstruct.R` are all
absent from the list - and nothing under `extensions/frmtmb.ode`,
`frmtmb.ddm`, `frmtmb.latent` or `frmtmb.sample` is touched. The scope
claim in the findings doc is accurate.

### The full core suite, one file per process, audited by name

113 processes, one per `test-*.R`, run in the package namespace
(`test_file(f, package = "frmtmb")`), `NOT_CRAN = true`:

```
expected 113   ran 113   unique 113   missing 0   extra 0   duplicates 0
files 113 | failed 0 | error 0 | warning 0 | skipped 87 | passed 6255
```

That reproduces the lane's headline figure exactly, to the passed count.

The files the contract audits by name, each its own process:

| file | passed | skipped | failed | error |
|---|---:|---:|---:|---:|
| test-nl.R | 14 | 0 | 0 | 0 |
| test-nlf.R | 74 | 0 | 0 | 0 |
| test-nl-lexical.R | 31 | 0 | 0 | 0 |
| test-nl-rtmb-scope.R | 103 | 0 | 0 | 0 |
| test-nl-body-vars.R (new) | 23 | 0 | 0 | 0 |
| test-predict-newdata.R | 12 | 0 | 0 | 0 |
| test-predict-lp-basis.R (new) | 30 | 0 | 0 | 0 |
| test-cens-trunc.R | 35 | 0 | 0 | 0 |
| test-cens-lccdf.R (new) | 27 | 0 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| test-bracket-access.R | 8 | 0 | 0 | 0 |
| test-compat.R | 265 | 0 | 0 | 0 |
| test-ps.R (new) | 43 | 0 | 0 | 0 |
| test-tmb-examples.R | 28 | 0 | 0 | 0 |

Every count matches the lane's table digit for digit.

A note on my own process, since it produced a wrong answer before a right
one: my first sweep ran the tests in the GLOBAL environment rather than
the package namespace, so every file calling an internal function
reported errors. That was my runner, not the lane. Fixed by passing
`package = "frmtmb"`, re-verified on `test-autoscale.R` (0 failed, 46
passed), and the whole sweep re-run from scratch to a fresh file. The
table above is that run.

### The extension suites, one process each

| package | files | failed | error | warning | skipped | passed |
|---|---:|---:|---:|---:|---:|---:|
| frmtmb.spline 0.2.0 | 7 | 0 | 0 | 0 | 0 | **265** |
| frmtmb.ddm 0.2.0 | 15 | 0 | 0 | 0 | 0 | **926** |
| frmtmb.ode 0.1.0 | 4 | 0 | 0 | 0 | 0 | **211** |
| frmtmb.sample 0.1.0 | 10 | 0 | 0 | 0 | 2 | 886 |

spline 265 and ode 211 are the lane's numbers exactly. ddm is 926 where
the lane reported 923 passed and 1 skipped: three assertions that skipped
for them ran for me, which is a Suggests package present in my library
and not theirs, not a difference in the code. Both are clean.
`frmtmb.sample` the lane did not report; it is clean too.
### as-cran

`R CMD build` then `R CMD check --as-cran --no-manual` against `rsc-lib`,
with `_R_CHECK_CRAN_INCOMING_=false` and pandoc from
`RSTUDIO_PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"`
on PATH.

| package | version checked | Status | ERROR | WARNING | NOTE |
|---|---|---|---:|---:|---:|
| frmtmb | 0.52.0 | **OK** | 0 | 0 | 0 |
| frmtmb.spline | 0.2.0 | **OK** | 0 | 0 | 0 |

Not one `WARNING`, `NOTE` or `ERROR` string appears anywhere in either
log. The stages that cost something:

| stage | core | spline |
|---|---|---|
| examples | [40s] OK | - |
| examples with `--run-donttest` | [73s] OK | - |
| tests | **[14m] OK** | **[30s] OK** |
| re-building vignette outputs | **[256s] OK** | [16s] OK |

The suites as they ran INSIDE the checks:

```
frmtmb        [ FAIL 0 | WARN 0 | SKIP 92 | PASS 6239 ]
frmtmb.spline [ FAIL 0 | WARN 0 | SKIP  1 | PASS  261 ]
```

Both differ slightly from my standalone runs (6255/87 and 265/0) because
a handful of tests gate on environment variables and on Suggests
packages that resolve differently inside a check's isolated library.
Nothing failed in any of the four runs.

The core suite took 14m inside the check where the lane reports 375s.
That is my machine, not the code: another lane's `roxygen2::roxygenise()`
and several R processes were running throughout. The RESULT is what
matters and it matches.

One scoping note, so the reader is not misled. Both tarballs were built
BEFORE the one edit I made to the worktree (see Edits below), so these
two OK statuses certify the lane's tree exactly as delivered. My edit
changes a roxygen comment block and the matching prose in a generated
`.Rd`; `R/ps.R` still parses (17 top-level expressions) and
`tools::parse_Rd()` still parses `man/ps.Rd` (36 sections), and neither
touches code that runs.

### The third pre-existing defect: `eval_dpars(b = NULL)`

Reproduced on main's 0.51.0 and confirmed fixed, with a correctness check
rather than only an absence of error:

| | main 0.51.0 | lane 0.52.0 |
|---|---|---|
| `eval_dpars(fit)` | mu[1:3] = 1.40350, 2.21388, 2.88821 | identical |
| `eval_dpars(fit, b = NULL)` | **Error: evaluation nested too deeply: infinite recursion** | mu[1:3] = 1.38248, 2.19286, 2.86719 |
| does it drop the RE contribution? | unreachable | `max|mu - predict(re.form = NA)|` = **0.000e+00** |

So it not only stops recursing, it returns exactly what
`?frmtmb-extension-api` documents `b = NULL` to mean, and the default
path is unchanged. Confirmed.
---

## Punch list

Ordered by what I would fix first. Nothing here is a correctness defect
in the numbers; the first two are documentation that contradicts measured
behavior, and the third is a robustness hole a future extension will find.

### 1. `?ps` told users the importance correction works. FIXED BY ME.

**R/ps.R:96-99** and **man/ps.Rd:101-104** said

> The importance correction ... draws a `ps()` block like any other
> random-effect block, so it applies.

It is always refused. `R/compat.R:1123` says `refused`, the findings doc
says refused, and the measured message is "`importance` cannot correct a
nonlinear predictor (y.mu)". The help page instructed the user to do
something that errors every time.

Edited in both files to match the compat row's reason. See Edits below.

### 2. `predict(newdata =)` past the knot span is silent

**R/ps.R:521** (`ps_coverage_warning()`), **R/predict.R:1215**.

The fit-end report covers the FITTED data. `predict()` at a new point
past `knot_range` returns a decaying partial sum, and past the outer knot
returns the body evaluated with the curve at exactly zero, with no
warning at any row (measured: 2.98550 and 2.83324 at t = 1.35 and 1.55,
then exactly `lev` at t = 3, on a fit whose `knot_range` ends at 1.328).

Either run the same coverage check over `newdata` and warn, or say so in
`?ps`'s "What the rest of the package does with it", which currently
only promises that `predict()` "works". The `pad =` argument and the
fit-end warning both exist because this cliff is real; `predict()` is the
one door it is not guarded at.

### 3. `fit_end_checks()` runs a family's hook with no guard

**R/fit-end.R:20-24**.

```r
fc <- resp$family[["post"]][["fit_check"]]
if (!is.function(fc)) next
fc(fit, resp$resp_name)
```

No `tryCatch`. A `fit_check` that throws - deliberately, or through a bug
in a third-party extension - destroys a completed fit after `frm()` has
done all the work. I confirmed it: a family whose hook calls `stop()`
makes `frm()` raise that error and return nothing.

`ps_coverage_warning()` on the line above has the same exposure through
its unguarded `eval_dpars(fit)` call.

`royston_parmar()` is careful about this (it wraps `rp_floored()` in
`try()` and only warns), but that is the extension being disciplined, not
core enforcing it. Wrap the hook so that a failing check degrades to a
warning naming the family. A fit that finished is worth more than a check
that did not.

### 4. NEWS does not say the phantom fix changes `start$beta`

**NEWS.md:93-102**.

Filed as `FIX`, not as the `BEHAVIOR CHANGE` the same file gives the
0.51.0 `mo()` change. The closing sentence, "Affected models fitted to
the same optimum; only the inference was lost", describes the old bug
rather than the new break. A model with a purely-random nonlinear
parameter loses one `beta` entry per such parameter, so an existing
script with a hand-written `start$beta` now errors (`start$beta must have
length 1`). The lane's findings doc states this plainly and the lane
edited two of its own tests for it; NEWS is what users read.

Add one sentence and the label.

### 5. The beta0 explanation attributes a convention the paper does not state

**dev/spline-core-findings.md** ("the paper constrains the fitted basis")
and **vignettes/case-studies.Rmd:1699-1704**.

The paper says only "subject to sum-to-zero constraint" and, in the
discussion, "a sum-to-zero constraint to ensure identifiability of an
intercept term". It never says on what. The conclusion is right - I
measured `beta0 + mean(fitted curve) = 68.1982` against the published
68.2 - so soften the attribution, not the claim. "The paper does not say
which sum-to-zero constraint it uses; re-centering ours to the fitted
values reproduces its intercept to 0.002 cm" is both weaker and true.

### 6. `ps()`'s knot rule is not the paper's, and the text implies it is

**R/ps.R:58-69** ("The knots are frozen") and
**vignettes/case-studies.Rmd:1673-1680**.

The paper's section 2.3.1 rescales the spline argument to [0, 1] using
bounds `min(t) - 3*sigma_b` and `max(t) + 3*sigma_b` that move smoothly
with the CURRENT variance estimates. `ps()` freezes knots on a padded
data range with the nonlinear parameters at zero. Both give a cubic basis
with 11 internal knots and the fits agree, but they are not the same
construction. One sentence in `?ps` saying so would stop a reader taking
`ps()` for a reimplementation.

### 7. The derived compat note is wrong for `ps()`

**R/compat.R:535**. Registering `ps` under `kind = "special"` attaches the
specials default note, "Predictor specials build design columns before
the family sees them", to 68 `<family> x ps()` rows. `ps()` builds no
design columns; that is the whole reason it exists. The statuses are
right (I tested four of them); only the sentence is.

### 8. `use_re` is duplicated verbatim

**R/predict.R:1085-1087** and **R/predict.R:2952-2954** are the same three
lines. The lane factored the three shared MESSAGE templates into
`check_re_form()`, `stop_unknown_response()` and `stop_newdata_missing()`
for exactly this reason; this is the fourth thing the two functions share
and the only one left inline. Being identical is why section 3 measured
zero, so this is tidiness, not risk.

### 9. A third `n_predict` assertion was left describing deleted machinery

**extensions/frmtmb.spline/tests/testthat/test-curve.R:48-66**, the test
named "a grouping block costs one probe, not one call per level".

The lane tightened two of the three `n_predict` assertions when the
probe went away - 11 to 1 at line 43, 32 to 1 at line 291 - and left
this one at `expect_lt(attr(cv, "check")$n_predict, 30L)`, above a
six-line comment that explains

> They are skipped in blocks of 24, and the measured count is 26 rather
> than the 9 the smooth alone would need: one chunk straddles the
> boundary between the group block and the smooth block...

There are no chunks any more. `n_predict` is the constant `1L` at
`R/curve-cov.R:84`, so the assertion passes vacuously and the comment
documents machinery this lane deleted. The test's title no longer
describes what it tests either: there is no probe. Tighten it to
`expect_equal(..., 1L)` and delete the comment.

(The findings doc says "both are pinned at 1", which is true of the two
it names; this is a third test it does not mention.)

### 10. `Suggests` is out of alphabetical order

**DESCRIPTION:40**. `brokenstick` sits between `ape` and `bayesplot`.

## Out of scope, but the consolidator should know

* **`simulate(newdata =)` is silently ignored on every model shape.**
  `simulate.frmtmb_fit` (R/predict.R:2607) has no `newdata` argument at
  2210aa1 or at 0.52.0, so it lands in `...`. A 5-row `newdata` returns
  80 rows on `y ~ x + (1 | g)`, on `y ~ s(x)`, on a plain nonlinear body
  and on a `ps()` model alike. Pre-existing; deserves its own issue.
* **`pi` in a nonlinear body is treated as a data column.**
  `bf(y ~ a * sin(2 * pi * t), ...)` fails with "variable lengths differ
  (found for 'pi')" on 0.51.0 and 0.52.0 alike, because `all.vars()` (and
  now `nl_body_vars()`) collects `pi` as a variable. Pre-existing, but
  this lane now owns the function that decides what a body's data
  variables are, so it is the natural place to fix it later.
* **`nl x <family>` compat rows are stale.** They say `refused` for
  `student`, `weibull`, `lognormal` and `tweedie`; all four fit. Stale on
  main too.
* **ce's reporting delta method is one-predictor.** After the merge,
  `frm_lp_basis()` provides exactly the multi-predictor Jacobian a
  softmax reporting transform needs, and ce's hook does not use it. An
  issue against ce, not against this lane.

## Edits I made to the worktree

One, and it is documentation only:

1. **`R/ps.R:96-100`** - replaced the roxygen paragraph claiming the
   importance correction "applies" with one saying it is refused and why,
   matching `R/compat.R:1123-1124`.
2. **`man/ps.Rd:101-105`** - the same replacement in the generated Rd, so
   the two stay consistent until roxygen is next run. `tools::parse_Rd()`
   parses the result (36 sections).

No other file in the worktree was touched. `git status` in the worktree
is otherwise exactly as I found it, and the main checkout is untouched
and still at 316a28b with a clean tree.
---

## Verdict: GO WITH FIXES

Every substantive claim in the lane's write-up reproduces, and several
reproduce more strongly than claimed. Nothing I measured contradicts a
number in `dev/spline-core-findings.md`. The headline results, all
re-derived rather than re-read:

* the SMOCC identity at **9.4133e-11** absolute, **2.5766e-14** relative,
  against a port of the paper's equation (4) that I wrote from the
  arXiv source with my own knot vector, my own second-difference penalty
  and my own sum-to-zero constraint - and the Wood (2004) eigensplit
  independently validated at `theta' S theta - ||omega||^2 = 4.135e-12`;
* `predict(se.fit = TRUE)` **bit-identical to main** (0.000e+00, not
  merely below 1e-12) on fifteen comparisons across ten model shapes,
  including the two that could most plausibly have moved, `rr` and `gp`;
* the gaussian `lccdf` table digit for digit, and the fit that the old
  `log(1 - F)` path could not complete at all now exact to 2.842e-13;
* the phantom-beta fix reproduced on both versions with the outer
  Hessian's eigenvalues as the witness, and the fit not moving to any
  printed digit;
* 113/113 core test files clean at **6255 passed**, and spline 265 and
  ode 211 exactly as reported.

The lane also found and fixed three defects that predate it and were
none of its business, and reported honestly on the things it did not do
(consumer (ii), `by =`, left truncation, simultaneous coverage) and on
the one place its change makes an optimizer complain where it used to
converge quietly.

The fixes are the reason this is not a plain GO, and none of them is a
correctness defect in a number:

* one help page contradicted the compat table and the measured behavior
  (**fixed by me**);
* `predict()` past a `ps()` knot span is silent about a cliff that the
  fit-end check exists to warn about;
* `fit_end_checks()` gives a family's hook no guard, so a throwing
  extension check destroys a completed fit - a new seam, worth closing
  before extensions start using it;
* NEWS files a breaking change to `start$beta` as a `FIX`;
* two documentation attributions to the paper are stronger than the
  paper's text.

All five are small and local. None of them is a reason to hold the lane
out of the trunk, and none of them interacts with a sibling.

### Merge order

The consolidation merges esicar, priors and ce before this lane, and
that order is right. Every `R/` file this lane shares with a sibling
merges with **zero conflicts**; only NEWS.md conflicts, three times, and
that is a hand merge with no semantics in it. The two places the lanes
genuinely meet are both safe and I built the four-lane tree to show it
rather than argue it: ce's `dpar_report_hook()` brackets the `se.fit`
block this lane rewrote without touching a line of it and consumes
`se_eta`, whose value I measured to be bit-identical before and after
the rewrite, so ce lands on the same numbers either way; and esicar's
`has_expand` correctly leaves a `ps()` block alone and correctly carries
it through `expand_b()`'s pass-through branch when something else turns
the gate on, which this lane already anticipated by indexing `c_idx`
rather than `b_idx`. On the merged tree a model with both a `ps()` term
and an `rr` block fits with coefficient space nearly twice parameter
space, a `ps()` + esicar model fits with `has_expand TRUE`, ce's mixture
weight comes back inside (0, 1) with a standard error built on the
rewritten `se_eta`, and the ten-shape sweep is 0.000e+00 against this
lane alone. Keeping spline-core last is also worth doing for a reason
beyond convenience: it is the only lane that bumps the version, and the
only one whose extension carries a floor (`frmtmb.spline` 0.2.0 needs
`frmtmb (>= 0.52.0)`) that the bump has to satisfy, so landing it last
means that floor is never briefly wrong in the trunk. The one thing to
do rather than assume is to re-run the ten-shape `se.fit` sweep after ce
lands and before this lane goes on top; it takes thirty seconds and it
is the only quantity a sibling reads that this lane recomputes.

---

# Punch re-check, 2026-09-05

Second pass over the same worktree after the lane's punch round. Main has
moved to **efbcbd4** (ce merged at d34d309, esicar at 61c479e); the lane
is still based at 2210aa1. `rsc-lib` rebuilt from the updated worktree;
current main installed separately into `rsc-main2-lib`; a third library
`rsc-mtree2-lib` holds the lane three-way merged ONTO efbcbd4, which is
what would actually land.

My own edit from the first pass (the `?ps` importance paragraph) survived
the punch round intact, now at `R/ps.R:111` and `man/ps.Rd:116`, with no
occurrence of the old "so it applies" text in either file.

## The three behavior items, reproduced with my own models

### (1) `ps_span_warning()` - CONFIRMED

`R/ps.R:553`, armed through a new `check` argument to `ps_env()`
(`R/ps.R:261`) that the closure reads at `R/ps.R:277`, set only at
`R/predict.R:1258-1259` as `check = !is.null(newdata)`. The objective's
`ps_env()` call (`R/objective.R:315`) and `eval_dpars()`'s
(`R/predict.R:587`) both take the default `FALSE`, so nothing arms it on
the tape.

On my own fit, `ps(t + sh, k = 10, pad = 0.3)`, knot span
[-0.31464, 1.315]:

| door | warnings | wanted |
|---|---:|---|
| `predict(newdata =)` inside the span | **0** | 0 |
| `predict(newdata =)` with 4 of 5 rows outside | **1** | >= 1 |
| in-sample `predict()` | **0** | 0, the fit-end report already spoke |
| `fitted()` | **0** | 0 |
| refit (does the tape see it?) | **0** | 0 |
| `simulate()` | **0** | 0 |

The message carries everything claimed:

```
4 of 5 predicted values of t + sh lie outside the frozen knot span of
ps(t + sh, k = 10, pad = 0.3) [-0.3146, 1.315]. The basis is a partial
sum there and exactly zero past the outer knot, so these are not
extrapolations of the fitted curve: the curve decays and the prediction
bends to whatever the rest of the body gives. Refit with a larger pad =
to cover the range you predict on
```

Row count `4 of 5` present, both span ends present, and it says what the
number IS rather than only that it is suspect. Far-field predictions at
t = 20, 50, 100 are 2.230507338 each against `fixef(fit)$lev` of
2.2305073: agreement is **exactly 0**, not merely under 1e-8.

`is.numeric()` at `R/ps.R:554` is what keeps it off the tape, and the
refit row above is the check that it works.

### (2) `fit_end_checks()` guards both checks - CONFIRMED, and better than claimed

`R/fit-end.R:18-44`. Both the coverage report and the family hook are
wrapped, each in its own `tryCatch`, the family one INSIDE the per-response
loop.

A family hook that throws:

```
fit returned         : TRUE (class frmtmb_fit)
logLik with hook     : -40.1279799371
logLik without hook  : -40.1279799371
|difference|         : 0        (claimed < 1e-10)
coefficients         : identical at tolerance 0
warnings             : 1, naming 'probeF' and post$fit_check
```

My first attempt to break the coverage report did not actually break it -
setting `knot_range` to a string makes `<` coerce rather than raise - so I
forced a real error instead by removing `fit$estimates`, which makes
`coef_b()` raise:

```
ps_coverage_warning() alone        : THROWS ("requires numeric/complex ...")
fit_end_checks() on the same object: returns without error
warning                            : names the coverage report
```

Two things I checked that the lane did not claim, both good:

* a hook that WARNS still delivers its warning to the caller - the guard
  catches errors only, which is right;
* on an `mvbf()` fit where response A's hook throws, response B's hook
  **still runs**. The `tryCatch` is inside the loop, so one family's
  failure does not silence the next one's check. That is the placement
  that matters and it is correct.

`post$fit_check` is documented at `man/frmtmb_family.Rd:86-93`, and the
documentation states the new contract explicitly: "Warn from it rather
than stopping; a hook that throws is caught, reported as a warning naming
the family, and the fit is returned regardless."

### (3) NEWS relabelled - CONFIRMED, and more thorough than asked

`NEWS.md:93-111`. The bullet is now "BEHAVIOR CHANGE, and a fix", and a
second bullet, "What changes for existing code", names the exact failure
(`start$beta must have length <n>`), the remedy (drop the element for the
parameter with no fixed design), where to look (`par_template(fit)`), and
states that fits do not move. It also records that two of the package's
own tests carried such a vector. That is more than the one sentence I
asked for and it is the right more.

## Items 5 to 10

| # | item | state |
|---|---|---|
| 5 | beta0 attribution | **fixed**, and the arithmetic is on the page |
| 6 | knot rule not the paper's | **fixed**, `R/ps.R:72` |
| 7 | derived compat note | **fixed**, 36 of 36 rows |
| 8 | `use_re` duplicated | **fixed**, `re_form_keeps()` |
| 9 | stale `n_predict` assertion | **fixed**, test renamed too |
| 10 | `Suggests` order | **fixed** |

Detail on the three worth checking rather than eyeballing:

**5.** The vignette no longer attributes a convention to the paper. It now
says the paper "does not say on what, so the right thing to do is measure
rather than attribute", and runs the arithmetic in a live chunk. The
findings doc carries the numbers (`dev/spline-core-findings.md:304`,
`:848`) and, to the lane's credit, the honest caveat I had not thought to
ask for: on an even 0-to-130-week grid the same re-centring gives 76.61,
because the children are not evenly spread over age, so "re-centred to
the fitted values" is a statement about the observed design and not a
universal one. That is a better answer than the one I asked for.

**7.** Measured on the installed package rather than read: 36 family x
`ps()` rows, all `works`, **0 still carrying** the old
"build design columns before the family sees them" note, and all 36 on a
single corrected note beginning "ps() builds no design column: it
evaluates a penalized spline on the tape and hands its VALUE to ...".
36 of 36, exactly as reported.

**8.** `re_form_keeps()` at `R/predict.R:62-66`, called at `:1101` and
`:2970`. Logically equivalent to the old three-liner (the old second line
was redundant), and the sweep below confirms it empirically at
0.000e+00.

**9.** The test is renamed "a grouping block costs nothing, because there
is no probe", the stale chunk-of-24 commentary is replaced by an accurate
historical note, and the assertion is now `expect_equal(..., 1L)`.

## The measurement the merge order needed

Current main **efbcbd4** carries ce's `dpar_report_hook()` and
`dpars_natural()`, which bracket the `predict(se.fit)` region this lane
rewrote. I installed three trees and ran the same ten-shape sweep on each:

* `rsc-main2-lib` - main at efbcbd4 (ce + esicar merged), 0.51.0;
* `rsc-mtree2-lib` - the lane three-way merged ONTO efbcbd4, which is what
  would land, 0.52.0;
* `rsc-lib` - the lane alone at base 2210aa1, 0.52.0.

| shape | merged vs current main | lane-alone vs current main |
|---|---|---|
| gaussian GLMM `re.form = NA` | 0.000e+00 | 0.000e+00 |
| gaussian GLMM `re.form = NULL` | 0.000e+00 | 0.000e+00 |
| gaussian GLMM at newdata | 0.000e+00 | 0.000e+00 |
| distributional, mu | 0.000e+00 | 0.000e+00 |
| distributional `dpar = "sigma"` | 0.000e+00 | 0.000e+00 |
| ordinal | 0.000e+00 | 0.000e+00 |
| smooth `re.form = NA` | 0.000e+00 | 0.000e+00 |
| smooth `re.form = NULL` | 0.000e+00 | 0.000e+00 |
| gp in sample | 0.000e+00 | 0.000e+00 |
| gp at newdata | 0.000e+00 | 0.000e+00 |
| rr `re.form = NA` | 0.000e+00 | 0.000e+00 |
| rr `re.form = NULL` | 0.000e+00 | 0.000e+00 |
| mixture, mu | 0.000e+00 | 0.000e+00 |
| mvbf | 0.000e+00 | 0.000e+00 |
| car (escar) | 0.000e+00 | 0.000e+00 |
| nonlinear | same refusal | same refusal |

**Shapes that moved beyond 1e-12: NONE**, in either column. And lane-alone
against merged is also identical everywhere, so the merge changes nothing
about the lane's own numbers.

### The shape that actually exercises ce's hook

The table above does not settle it on its own, because ce's hook fires
only for a mixture's `theta*` dpar on the RESPONSE scale, and the sweep's
mixture row is `mu`. So I ran that case directly:

| `predict(dpar = "theta1", type = "response")` | fit[1] | se[1] |
|---|---:|---:|
| current main efbcbd4 (has ce) | 0.33333333 | 0.03849003 |
| **lane merged onto efbcbd4** | **0.33333333** | **0.03849003** |
| lane alone (base 2210aa1, pre-ce) | -0.69314718 | 0.17320513 |

| comparison | max abs d(fit) | max abs d(se) |
|---|---:|---:|
| merged vs main | **0.000e+00** | **0.000e+00** |
| lane alone vs main | 1.026e+00 | 1.347e-01 |

That is the answer. ce's softmax reporting transform survives this lane's
rewrite EXACTLY, standard error included - and the standard error is the
part that matters, because it is `abs(hook_val("deriv")) * se_eta` and
`se_eta` is precisely the quantity the lane now computes through
`lp_basis_out()`. The lane alone differs by 1.03 only because it branched
before ce existed, not because it breaks anything.

On the link scale all three agree at 0.000e+00, as they should: ce's hook
does not touch it.

### Is the merge still conflict-free against current main?

Three-way merge, base 2210aa1, ours = the lane, theirs = efbcbd4:

| file | main changed it | conflicts |
|---|---|---:|
| **R/predict.R** | YES | **0** |
| **R/families.R** | YES | **0** |
| R/frame.R | YES | **0** |
| R/objective.R | YES | **0** |
| R/compat.R | YES | **0** |
| R/fit.R | no | 0 |
| R/parse.R | no | 0 |
| NAMESPACE | no | 0 |
| DESCRIPTION | no | 0 |
| NEWS.md | YES | 1 |

Still conflict-free on every `R/` file, including the two the brief names.
Only NEWS.md conflicts, as before. The merged tree installs clean and
reports version 0.52.0.

Why predict.R stays clean, in base coordinates: the lane's new
span-arming edit is at 1205-1217, inside the NONLINEAR early-return path,
and its se.fit block is 1226-1244; ce's hunks are 1190-1202, 1220-1228
and 1244-1250. The `@@` ranges touch because they include three lines of
context each, but the changed lines are disjoint, which is what
`git merge-file` and the bit-identical sweep both say.

## Runs

Every file one process, `NOT_CRAN = true`, package namespace.

| file | passed | failed | error | warning | skipped |
|---|---:|---:|---:|---:|---:|
| `test-ps.R` | **51** | 0 | 0 | 0 | 0 |
| `test-fit-end-hook.R` | **18** | 0 | 0 | 0 | 0 |
| `test-predict-newdata.R` | 12 | 0 | 0 | 0 | 0 |
| `test-predict-lp-basis.R` | 30 | 0 | 0 | 0 | 0 |
| `test-message-uniqueness.R` | 6 | 0 | 0 | 0 | 0 |

Those are the only two `test-predict*.R` files in the suite. `test-ps.R`
is 51 where it was 43 before the punch round, and `test-fit-end-hook.R`
at 18 is the lane's figure exactly.

Full core suite, one file per process, audited against `ls test-*.R`:

```
expected 114   ran 114   unique 114   missing 0   extra 0   duplicates 0
files 114 | failed 0 | error 0 | warning 0 | skipped 87 | passed 6281
```

**114/114, 6281 passed, 87 skipped** - the lane's numbers to the digit,
and 114 is 113 plus `test-fit-end-hook.R`.

`frmtmb.spline` 0.2.0, one process: 7 files, 0 failed, 0 error, 0 warning,
0 skipped, **266 passed**. The lane's figure exactly (266, up from 265
before the punch round).

### as-cran, re-run on the punched tree

| package | version | Status | ERROR | WARNING | NOTE |
|---|---|---|---:|---:|---:|
| frmtmb | 0.52.0 | **OK** | 0 | 0 | 0 |
| frmtmb.spline | 0.2.0 | **OK** | 0 | 0 | 0 |

Not one `WARNING`, `NOTE` or `ERROR` string in either log. Stages:
examples [40s] OK, `--run-donttest` [61s] OK, tests [11m] OK, vignette
re-build [184s] OK. The suites inside the checks:

```
frmtmb        [ FAIL 0 | WARN 0 | SKIP 92 | PASS 6265 ]
frmtmb.spline [ FAIL 0 | WARN 0 | SKIP  1 | PASS  262 ]
```

The vignette rebuild passing matters here, because item 5 put live
arithmetic on the page: the `smocc-beta0` chunk executes during the
rebuild, so "68.75433, 68.19824, 68.20000" is a number the check
reproduces rather than a number the lane typed.

## Punch list, final state

| # | item | state |
|---|---|---|
| 1 | `?ps` claimed the importance correction applies | **closed** (my first-pass edit, kept intact) |
| 2 | `predict()` silent past the knot span | **closed** - `ps_span_warning()`, verified on my own models |
| 3 | `fit_end_checks()` unguarded | **closed** - both checks wrapped, verified with a genuinely throwing hook and a genuinely throwing coverage report |
| 4 | NEWS filed a break as a FIX | **closed** - relabelled, plus a "what changes for existing code" bullet |
| 5 | beta0 attributed to the paper | **closed** - measured on the page, with a caveat I had not asked for |
| 6 | knot rule implied to be the paper's | **closed** - `R/ps.R:72` |
| 7 | derived compat note wrong for `ps()` | **closed** - 36 of 36 rows |
| 8 | `use_re` duplicated verbatim | **closed** - `re_form_keeps()` |
| 9 | stale `n_predict` assertion and comment | **closed** - test renamed, assertion tightened |
| 10 | `Suggests` out of order | **closed** |

Ten of ten. Two were answered better than they were asked: item 5 added
the even-grid caveat (76.61) that makes the re-centring claim honest about
its own scope, and item 3's guard turns out to be per-response, so one
family's failing hook does not silence the next one's.

### One residual, new, low severity

`frm_lp_basis()` does not arm the span check. Its two `ps_env()` calls
(`R/predict.R:3131`, `:3151`) take the default `check = FALSE`, so on an
identical past-span grid:

```
predict(newdata = past-span)                warnings: 1
frm_lp_basis(newdata = past-span)           warnings: 0
frmtmb.spline::frm_curve(newdata = past-span) warnings: 0
```

with the three returning the same values (max |d eta| = 0, far field
exactly `lev`). The numbers are right; only the diagnostic is missing.

It matters a little more than it looks, because `frm_curve()` is the
function a user calls to DRAW the curve with a band, which is exactly
when they reach past the fitted range - the lane's own vignette draws the
SMOCC band through this path. `predict()`, the door that now warns, is
the one a plotting user is least likely to use.

Not a gate. `frm_lp_basis()` is a seam that hands back pieces, and one
could argue the caller owns the diagnostic; but if the warning is worth
having at `predict()` it is worth having at `frm_curve()`. A follow-up
for whoever owns `frmtmb.spline`, not a reason to hold this lane.


### Addendum: main moved again during this re-check

The brief named efbcbd4 and I measured everything above against it. Main
finished this session at **ffe7f4b**, three commits further on:

```
ffe7f4b docs: the frmtmb.eam subsite path, with a redirect from the old one
a5619b7 frmtmb.ddm is renamed frmtmb.eam
9a1d3f3 Merge branch 'wt-rdm-gng'
eeab234 frmtmb.ddm 0.3.0: the racing diffusion model and the go/no-go diffusion model
```

None of it touches an `R/` file this lane touches, so the se.fit
measurement above stands unchanged: `R/predict.R` and `R/families.R` are
byte-identical at efbcbd4 and ffe7f4b.

It does add ONE shared file the earlier pass did not have to consider,
`vignettes/case-studies.Rmd`, which main edited for the rename and which
this lane extends with its SMOCC section. Re-running the three-way merge
against ffe7f4b:

| file | main changed it | conflicts |
|---|---|---:|
| R/predict.R, R/families.R, R/frame.R, R/objective.R, R/compat.R | YES | **0** each |
| R/fit.R, R/parse.R, NAMESPACE, DESCRIPTION | no | 0 each |
| **vignettes/case-studies.Rmd** | **YES** | **0** |
| NEWS.md | YES | 1 |

Still clean, and the result is right rather than merely conflict-free.
main's hunks in the vignette are at 1443-1455 and 1659-1665; the lane's
are at 19-24, 1604-1609 and 1643-1648. In the merged file:

```
frmtmb.eam mentions : 4
frmtmb.ddm mentions : 0
SMOCC mentions      : 46   (the lane's section 13, intact)
ps(age + shift, k = 15, pad = 0.25) present at line 1663
```

That is the outcome you want: the lane never names `frmtmb.ddm` itself -
its four occurrences are inherited from the base - so main's rename wins
those lines and the lane's new section is untouched beside it.

The NEWS.md conflict is one block (ours 3-119 against theirs 120-332),
both lanes prepending a section at the top of the file. Keep both, in
version order. No semantics.

So the verdict below is unchanged by main having moved, and the merge is
clean against the trunk as it actually stands rather than as the brief
described it.

### Second addendum: and again, to 32b70bd with priors merged

Main moved twice more while this re-check ran, finishing at **32b70bd**
(`Merge branch 'wt-priors'`). `R/predict.R` and `R/families.R` are still
byte-identical to efbcbd4, so every number measured above stands. priors
does add NAMESPACE to the shared set. Re-running the merge once more, now
against 32b70bd:

| file | main changed it | conflicts |
|---|---|---:|
| R/predict.R, R/families.R, R/frame.R, R/objective.R, R/compat.R | YES | **0** each |
| **NAMESPACE** | **YES** (priors) | **0** |
| vignettes/case-studies.Rmd | YES | **0** |
| R/fit.R, R/parse.R, DESCRIPTION | no | 0 each |
| tests/testthat/test-compat.R, test-tmb-examples.R | no | 0 each |
| NEWS.md | YES | 1 |

Still one conflict, still only NEWS.md. NAMESPACE takes the lane's three
exports (`frm_joint_cov`, `frm_lp_basis`, `ps`) alongside priors' without
a clash.

**One watch item, not a finding.** Main's working tree at 32b70bd carries
an UNCOMMITTED modification to `tests/testthat/test-compat.R` - somebody's
in-flight work, not mine and not this lane's. That is the one file where
this lane also edits (it adds `ps` to the specials vocabulary), so the
merge check above compares against the COMMITTED 32b70bd and cannot see
it. If that edit lands before this lane does, re-run the three-way merge
on that one file. Everything else is measured.

For the record on hygiene: I never wrote to the main checkout. My only
write outside my own review document, across both passes, is the `?ps`
documentation edit from the first pass, which the lane kept.
## Verdict on the punch round: GO

Every item I raised is closed, and closed by measurement rather than by
assertion: I reproduced the span warning and the throwing-hook case on my
own models, forced a real failure through `ps_coverage_warning()` after my
first attempt turned out not to break it, and read the compat note off the
installed package rather than the diff. Nothing regressed - 114/114 core
files, 6281 passed, and both as-cran checks OK with zero notes.

### Merging last: still right, and now measured against the real trunk

The first pass argued the merge order from the lane against a main that
predated ce and esicar. That argument is now unnecessary, because main
HAS them and I measured against it directly.

The three-way merge of the lane against efbcbd4 is conflict-free on every
`R/` file, `R/predict.R` and `R/families.R` included, both of which main
changed; only NEWS.md conflicts, once. The lane merged onto efbcbd4
installs clean and produces se.fit values **bit-identical to current
main** on all fifteen comparisons across ten shapes - nothing moved at
1e-12 or anywhere else. And on the one case that actually exercises ce's
`dpar_report_hook()`, a mixture's `theta1` on the response scale, the
merged tree reproduces current main to **0.000e+00 on both the fitted
value and the standard error**, where the lane alone differs by 1.03
simply because it branched before ce. The standard error is the
load-bearing half of that: ce computes it as
`abs(hook_val("deriv")) * se_eta`, and `se_eta` is precisely what this
lane re-routed through `lp_basis_out()`. It hands ce's hook exactly what
ce wrote against.

So spline-core should still go last, and the reason is now stronger than
convenience: it is the only lane that bumps the version, the only one
whose extension carries a floor the bump must satisfy
(`frmtmb.spline` 0.2.0 needs `frmtmb (>= 0.52.0)`), and the only one whose
interaction with an already-merged sibling has been measured end to end
on the tree that would actually land. The re-run I recommended after ce
landed is the one above; it is clean, so there is nothing left to hold
for.
