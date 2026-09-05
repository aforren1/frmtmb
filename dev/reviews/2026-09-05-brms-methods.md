# Review: the brms post-fit method tier (`wt-brms-methods`)

Reviewer's ledger. Every number below was reproduced on this machine
rather than read from `dev/brms-methods-tests.md`; where mine differs
from the lane's, both are shown and the difference is named.

Environment: Windows 11, R 4.6.1, brms 2.23.0, rstan 2.32.7,
StanHeaders 2.39.1, posterior 1.7.0. The worktree's core was installed
into a private library and every run used an explicit `--library=` /
`lib=`. The Stan programs were read from a private COPY of the lane's
`dev/stan-cache-methods`, so nothing in the lane's tree was written to.

Scope confirmed: `git diff --name-only 564e185` reports `.gitignore`
and `NEWS.md` only; the untracked set is
`dev/brms-methods-tests.md`, `dev/makevars-stan-methods`,
`tests/testthat/helper-brms-methods.R`,
`tests/testthat/test-brms-methods.R`. **Nothing under `R/` changed**, so
the claim that the full core suite is not required stands. Main is
clean at `5dfdd84` and was not touched.

---

## 1. The mechanism

### Route (b) reproduces

`brm(empty = TRUE)` scaffold + `rstan::sampling(algorithm =
"Fixed_param", init = list(pars))` on the cached `stanmodel`, assigned
into `$fit`, then `brms::rename_pars()`.

| shape | draws | max spread across draws | `posterior_epred[1, ]` vs `fitted()` | `log_lik[1, ]` vs frmtmb row density | `sum(log_lik)` vs `logLik` |
| --- | --- | --- | --- | --- | --- |
| r1 (`y ~ x + z, sigma ~ x`) | 10 x 9 | 0 | **0**, `identical()` TRUE | **0** | -239.538977525 both, diff 2.8e-14 |
| r16 (`y ~ x, zi ~ x`, ZIP) | 10 x 8 | 0 | **4.44e-16**, `identical()` FALSE | **8.88e-16** | -418.097315756 both, diff 3.4e-13 |

The mechanism holds. One correction to the ledger: the doc says the
epred identity is "**0** (exactly, not to a tolerance)" and that this
holds "on row 1, and then on every shape". It is exactly 0 on row 1;
on r16 it is 4.44e-16 and `identical()` is FALSE. The tier's own
`expect_exact_num()` at 1e-8 relative is the right assertion and passes;
the DOC overstates. See punch list P8.

### Route (a) is byte-identical

`brm(..., algorithm = "fixed_param", init = list(pars), chains = 1,
iter = 10, warmup = 0)` on r1, against route (b) on the same shape:

```
route (a) seconds: 171.9
posterior_epred identical (all draws): TRUE
log_lik          identical (all draws): TRUE
variables()      identical: TRUE      fixef() identical: TRUE
route (a) epred vs frmtmb fitted(), max abs: 0
```

The lane's 153.6 s for route (a) is 171.9 s here under load; the
conclusion is unaffected. Route (b) shape build measured 18.9 s (r1,
first shape in a cold process: frmtmb fit + scaffold + RDS read +
sampling) and 7.1 s (r16). The doc's "1.04 s per shape" is the
sampling-and-scaffold portion, not the whole shape; the argument for
route (b) rests on the compile, which is real, and is not weakened.

### `brms:::exclude_pars` is the only internal, and its absence is safe

`grep ':::'` over `test-brms-methods.R`, `helper-brms-methods.R` and
`helper-brms.R` returns nothing. The only namespace reach is
`helper-brms-methods.R:53-58`. `exclude_pars` is in brms's namespace and
is NOT exported, so the `exists()` guard is doing real work.

I ran the guard's NULL branch by hand (`pars = NULL, include = FALSE`):

- `rename_pars()` does not error.
- On r1 (no group-level effects) there is **no difference at all**: 9
  variables either way, same set.
- On rC0 (sleepstudy) the fit gains 80 variables: `r_1[i,j]`, `z_1`,
  `L_1`, `Cor_1` alongside the renamed `r_Subject` spellings. Exactly
  what the comment predicts.
- Critically, `ranef()` is `identical()` to the guarded path and
  `posterior_epred()` is `identical()` to it. **No assertion in the
  tier would change.**

So the guard degrades the object's tidiness and not the tier's results.
One coupling to note: `brms_dpars_of()` and `brms_dpar_is_scalar()`
read `brms::variables(shape$brmsfit)`, so they are the two helpers a
future brms that changes `exclude_pars()` could move. Neither is
affected today (the extra names are `r_/z_/L_/Cor_`, none of which is a
dpar name), but it is the seam to watch.

**Verdict on the mechanism: sound, and the claim is understated rather
than overstated.**

---

## 2. The 22 divergences, classified

`(D)` a frmtmb defect. `(P)` a paradigm difference that is right for a
frequentist fit. `(C)` a design choice needing a user decision.

The "pin" column answers the question that matters for merging: when a
`(D)` is FIXED, does a NAMED expectation flip (acceptable: the fix is
told to update the pin), or does the fix pass silently / break an
unrelated check (not acceptable)?

| # | divergence | class | evidence (measured here) | pin quality |
| --- | --- | --- | --- | --- |
| 1 | zero-inflated `conditional_effects()` plots `exp(eta)` | **D** | ce 0.701289576884 vs epred 0.608531126228 at grid point 1; `predict(type = "conditional")` returns the ce value to 0 | **acceptable** |
| 2 | mixture `theta1` response scale is the predictor | **D** | `range(predict(type="response", dpar="theta1"))` = [0.607, **1.112**]; 1.25% of rows exceed 1 | **acceptable** |
| 3 | `conditional_effects()` refuses a mixture with `theta ~ x` | C | error names the dpar it looked at: "No plottable predictors found for dpar 'mu1'" | acceptable |
| 4 | conditions on `signif(mean +/- sd, 3)` | **D** | held at -0.9360, 0.0267, 0.9900 vs exact -0.93621694, 0.02673023, 0.98967740; max abs 3.226e-4 | **acceptable, noisy** |
| 5 | `int_conditions=` accepted, ignored, warning blames `predict()` | **D** | grid and estimates `identical()` to the call without it; warning is "ignoring unknown arguments to predict(): int_conditions", `conditionCall` NULL | **acceptable** |
| 6 | ce frame lacks `y`, `z`, `cond__`, `effect1__` | **D** | `names(cf)` = `x, estimate__, se__, lower__, upper__` | **acceptable** |
| 7 | two-way grid row order reversed | C | both hold the same 300 points; neither order is wrong | acceptable |
| 8 | `re_formula = NULL` conditions on the first observed level | **D** | curve matches Subject **308** to 5.68e-14 and no other level; no `Subject` column in the frame | **acceptable** |
| 9 | ordinal ce: layout, ignored `categorical=`, effect key | **D** | `categorical = TRUE` returns an `identical()` frame **and emits no warning at all** | **acceptable** |
| 10 | `mo()` gets a 100-point grid | **D** | frmtmb 100 rows, brms 4; agree at the 4 defined levels | **acceptable** |
| 11 | nonlinear predictor refused a Wald band | P | no analytic delta-method SE exists; message names three routes, `method = "predict"` reproduces brms to 1e-16 | acceptable |
| 12 | `method=` vocabulary | C | `"epred"`/`"predict"` vs `posterior_*`; values agree where both resolve | acceptable |
| 13 | `predict()` is the conditional mean | P | equals `fitted()` to 0; `stats`/`lme4` convention, `simulate()` carries the draw | acceptable |
| 14 | `hypothesis()` returns the table, not a `brmshypothesis` | C | `frmtmb_hypothesis` data frame; brms idiom fails loudly | acceptable |
| 15 | `sum(log_lik())` conditional vs `logLik()` marginal | P | 75.72 nats on `(1\|q\|g)`; the Laplace correction is the frequentist quantity | acceptable |
| 16 | pearson residuals divide by different quantities | P | frmtmb `(y-mu)/sigma` to 0; brms divides by a Monte Carlo predictive SD and deprecates the type | acceptable |
| 17 | scalar dpar `posterior_linpred` natural vs link scale | P | both self-consistent; `posterior_epred(dpar=)` agrees exactly on both kinds | acceptable |
| 18 | `ranef()` keyed by block, `coef()` by factor | **D** | `names(ranef(fit))` = `"Days \| Subject"`, `names(coef(fit))` = `"Subject"`; `ranef(fit)$Subject` is NULL while `coef(fit)$Subject` is not | **acceptable** |
| 19 | `coef()` return type depends on the model | C | numeric / list / list of data frames across three shapes | acceptable |
| 20 | `transform = TRUE` probability vs count under `trials()` | P | `p * n` equals frmtmb exactly; different named quantities | acceptable |
| 21 | `se()`'s unused sigma reported 1 vs 0 | **D** | frmtmb returns `exp(0) = 1` for a parameter the density does not use | **acceptable** |
| 22 | multi-column linear predictors unreachable for `cs()` | C | categorical and mixture columns ARE reachable via `dpar=`; `cs()` has no spelling | acceptable |

**Tally: 10 (D), 6 (P), 6 (C).**

A 23rd asserted divergence is not in the brief's list of 22: finding 6,
`conditional_effects(dpar =)` enumerating different effects
(test line 894, 8 assertions). It is **(C)** and the doc already argues
frmtmb's side. The test file carries 21 `# DIVERGENCE` comment blocks,
two of which cover two divergences each ("TWO DIVERGENCES on one call"
at line 763, and the ordinal block covering three at line 922), plus two
explicit `# NOT a divergence to fix` blocks. The count is consistent;
the brief's 22 simply omits finding 6.

### The (D) items, pin by pin

Every one of the ten is pinned by asserting the CURRENT wrong output
against a NAMED expectation, so a fix produces a labeled failure that
points at the assertion to update. None of the ten is pinned in a way a
fix would silently satisfy. Detail where it is not obvious:

- **D1** (line 620): three named expectations flip -- `"frmtmb ce is
  exp(eta), r16"`, `"frmtmb ce is type=conditional, r16"`, and
  `expect_gt(max(cf$estimate__ / cb$estimate__ - 1), 0.15)`. The three
  that must SURVIVE a fix (`predict(response) is epred`, `ce
  method=predict is epred`, `frmtmb fitted() is epred`) are already
  written as agreement assertions. This is the model the other pins
  follow.
- **D4** (line 760): flips `expect_identical(zf, sort(signif(...)))`,
  `expect_false(all.equal(zb, zf))` and `expect_gt(gap, 1e-6)`. It also
  breaks `expect_false(anyNA(m))` at line 801, because the alignment key
  `kf` stops matching `kb` once `cf$z` is unrounded. That is a real
  failure for the right reason under an unhelpful label; acceptable, but
  the block needs one comment saying so. Punch list P5.
- **D10** (line 678): `expect_identical(nrow(cf$inc), 100L)` and
  `expect_gt(length(setdiff(cf$inc$inc, 0:3)), 90)` flip; the
  `match()`-aligned estimate comparison is the invariant and survives.
  Correctly separated.
- **D21** (line 245): the response-scale pin flips; the link-scale
  `expect_identical(fixef(...)$sigma[["(Intercept)"]], 0)` survives if
  the fix is at the reporting layer. Correctly separated.

### The one place a fix passes silently

Not in the assertions -- in the helper's exclusion lists.
`brms_ce_shapes()`, `brms_meanlink_shapes()`, `brms_linpred_shapes()`,
`brms_resid_shapes()` and `brms_dpars_of()` route the diverging shapes
AROUND the agreement loops. When D1 is fixed, `r16` and `rC16` must be
put back into `brms_ce_shapes()`
(`helper-brms-methods.R:577-581`); when D2 is fixed, `theta` must come
back into `brms_dpars_of()` (`helper-brms-methods.R:501`); when D21 is
fixed, `r14c`'s `sigma` (`helper-brms-methods.R:505-507`). **If the
list is not edited, nothing fails and the tier silently stops covering
the shape it was built for.** The comments above each list are good and
name the reason, but a comment is not an assertion. Punch list P1.

### Status per divergence: the premise is false

The brief asks whether the dev doc records Status per divergence "the
way `dev/brms-likelihood-tests.md` does". It does not, and neither does
that document: `grep -n Status dev/brms-likelihood-tests.md` returns
exactly one hit, line 3, the document-level status. There is no
per-divergence Status convention in this repository to follow.

What `dev/brms-methods-tests.md` has instead is a bolded
**"Which is right: ..."** verdict, on 9 of 22 findings, and a
**"What a porting user experiences"** note on 6. Thirteen findings carry
neither. Findings 11b, 12, 13, 14, 15, 16 and 18 in particular end with
prose that leans one way without saying which package is right or
whether frmtmb intends to change. That is the real gap, and it is worth
closing before merge because the (D)/(P)/(C) call is exactly what the
next person needs and exactly what is missing. Punch list P2.

---

## 3. Reproductions without Stan, from the installed core alone

All five reproduce. Numbers are mine.

**(1) zero-inflated `fitted()` vs `conditional_effects()$estimate__`.**

| x | ce `estimate__` | `predict(type = "response")` | ratio |
| --- | --- | --- | --- |
| -2.84476 (point 1) | 0.701289576884 | 0.608531126228 | **1.152430** |
| 0.00142 (point 50) | 1.837330 | 1.127389 | 1.629721 |
| 2.90569 (point 100) | 4.909230 | 1.334783 | **3.677924** |

`predict(type = "conditional")` returns the ce value to 0;
`(1 - zi) * conditional` returns the epred value to 0
(zi[1] = 0.132268400549). `conditional_effects(method = "predict")`
gives the epred value. So the lane's "15.2% high" is right **at the
first grid point only**: the error grows monotonically to **268% high**
at the last one. The NEWS bullet's "runs 15% above the mean of the data"
is wrong twice over -- it is not the mean of the data, and 15% is the
smallest error on the curve, not a typical one. Punch list P7.

**New, and beyond the lane's claim: the same defect is live on the
hurdle families, and needs no Stan to show.** On
`bf(y ~ x, hu ~ x) + hurdle_poisson()`, default `conditional_effects()`
against `predict(type = "response")`:

```
ce[1]  = 0.52647059213     epred[1] = 0.982645403146    ratio = 0.536
max ratio - 1 over the grid = 1.196
```

Here the default curve is 46% BELOW the expected response at one end
and 120% above at the other. The doc predicts this ("the mechanism
above predicts the same result for them") and stops. It is a two-line
frmtmb-only check. Punch list P3.

**(4) the `signif` rounding.** `ce_second_values()` at
`R/conditional-effects.R:87-94` is
`signif(mean(col) + c(-1, 0, 1) * sd(col), 3)`. On row 1:

```
exact   -0.9362169409490204   0.0267302282880733   0.9896773975251669
frmtmb  -0.9360               0.0267               0.9900
identical(got, sort(signif(exact, 3))) : TRUE
max abs difference: 3.226e-4
```

**(5) the `int_conditions` warning.** Text reproduces verbatim:
`ignoring unknown arguments to predict(): int_conditions`, class
`simpleWarning`, `conditionCall` **NULL** (so the warning cannot even be
traced to a call). Grid and estimates are `identical()` to the call
without the argument, and no held value is in `c(-1, 1)`.

**New, and it weakens the lane's own correction.** The doc's finding 3
says the argument "is not swallowed in silence" and treats the warning
as present. It is present on the gaussian shape only. On an ordinal
fit the same call produces **no warning at all**:

```
conditional_effects(fo, categorical = TRUE)      -> NO WARNING
conditional_effects(fo, int_conditions = ...)    -> NO WARNING
conditional_effects(fo, nosucharg = 1)           -> NO WARNING
```

The warning comes from the `method = "epred"`, `band = "wald"` branch
forwarding its dots to `predict()`; the ordinal/categorical branch does
not forward, so unknown arguments there are silently discarded. That
makes divergence 9's "`categorical=` is accepted and ignored" strictly
worse than the doc records, and it is a second, independent finding.
Punch list P4.

**(8) `re_formula = NULL`.** frmtmb conditions on **Subject 308**, the
first `rownames(ranef(fit)[[1]])`, matching to **5.68e-14** and matching
no other level. The frame does NOT say so:
`names(cf)` = `Days, estimate__, se__, lower__, upper__`, no `Subject`
column. The curve is 86.5 away from the `re_formula = NA` population
curve at its furthest point, so the choice is large and invisible.

**(2) the mixture theta link.**
`family(fit)$links` = `mu1=identity, sigma1=log, mu2=identity,
sigma2=log, theta1=identity`. `predict(type = "response", dpar =
"theta1")` is `identical()` to `type = "link"`;
`max|response - plogis(eta)| = 0.359462897571`.

**New:** the returned value is **not confined to [0, 1]**. On the
fitted data `range` is [0.606969187849, **1.111956599700**] and 1.25% of
rows exceed 1; on `x` in [-4, 4] it reaches 1.126. This is the evidence
that settles divergence 2 as (D) rather than a scale convention: a
mixing weight above 1 is not a value any convention makes right.

---

## 4. CI impact

**Does the methods tier run in `brms-likelihood.yaml` today? No.** The
run step names one file, `tests/testthat/test-brms-likelihood.R`, so
the methods tier is not executed anywhere in CI.

**Would the same gate work? Yes.** Both tiers gate on the identical
`skip_unless_brms_fit()` in `helper-brms.R`, which reads
`FRMTMB_BRMS_FIT_TESTS`, already set to `"true"` in that job's `env:`.
No new gate is needed for the tier to RUN.

**Is the on-disk cache shared within one job? Yes, by construction.**
The job sets `FRMTMB_STAN_CACHE: ${{ github.workspace }}/dev/stan-cache`,
and `helper-brms-methods.R` compiles only through `brms_stan_model()`
from `helper-brms.R`, which reads that variable. The keys are content
addressed on the Stan code plus the rstan version, so identical
programs collide by design.

**Measured, not assumed:** the 23 methods shapes need **23 distinct
programs**, all 23 are present in the lane's cache, and no cache entry
is orphaned. Against the sibling lane's `dev/stan-cache` (36 programs,
built by the log-density tier) exactly **one** is new:
`0cc4f6d974bc3704ca76b4bd619d0fdf`, the `rfac` shape `y ~ x * f` --
which is precisely the shape the lane added on purpose because no
log-density row has a factor covariate.

**So the wall time becomes:** cold, +1 Stan compile (55-70 s on the
box's own figures) on top of the existing 24, plus the methods tier's
own non-Stan work. Warm, the tier costs **121.8 s** measured here on a
loaded Windows box against the lane's 145.4 s; Linux is faster. The
90-minute timeout is not at risk either way.

**One real defect if the tier is bolted on unchanged.** The cache key
step hashes only

```
tests/testthat/test-brms-likelihood.R
tests/testthat/helper-brms.R
```

`actions/cache` saves at post-job **only when the primary key was not an
exact hit**. Add the methods tier without adding its two files to that
`md5sum()` list and the first run compiles the `rfac` program, restores
an exact key hit, saves nothing, and recompiles `rfac` on **every
subsequent run, forever**. The fix is one line: put
`test-brms-methods.R` and `helper-brms-methods.R` into the `h <-
tools::md5sum(c(...))` vector. Punch list P6.

**Does the methods tier need its own gate? No, and it should not have
one.** A second environment variable would let the two tiers drift out
of sync and would forfeit the shared cache, which is the entire economic
argument for route (b). What it needs is to be added to the existing
job's run step as a SECOND `test_file()` call in a SEPARATE process --
the lane is right that `test_file()` re-sources the helper and resets
the in-session model cache, and the existing step's `stopifnot` /
skip-detection wrapper should be applied to both. The job name and the
workflow comment ("Only this tier runs") both need updating.

---

## 5. `dev/makevars-stan-methods`

Two lanes have now independently created the same file:

- this lane, `dev/makevars-stan-methods` (untracked, added to
  `.gitignore`? **no** -- see below): `CXX17FLAGS = -O2 -Wall -std=gnu++17`
- the `wt-brms-rows` lane, `dev/stan-cache/makevars-cxx17.mk`:
  `CXX17FLAGS = -O2 -Wall -mfpmath=sse -msse2 -mstackrealign -std=gnu++17`

**What should survive consolidation: the documented recipe, not the
file.** The reasoning:

1. The repair is already IN the repository twice as documentation:
   `dev/brms-likelihood-tests.md` "Toolchain finding" and, in executable
   form, the `Force C++17 for rstan's generated models` step of
   `.github/workflows/brms-likelihood.yaml`. That workflow step is the
   canonical recipe and it works on the CI platform.
2. A checked-in `.mk` cannot be right for two platforms at once: the
   rows lane's variant carries `-mfpmath=sse -msse2 -mstackrealign`,
   which are Windows/i386-era flags with no meaning on the CI runner.
   Shipping either file invites someone to point `R_MAKEVARS_USER` at
   the wrong one.
3. It is environment configuration for a machine, not package content.
   `.gitignore`ing it is the right instinct; keeping two ignored copies
   in two worktrees is just duplication.

Recommendation: delete `dev/makevars-stan-methods` from the lane, and
add one paragraph to `dev/brms-methods-tests.md` (or a cross-reference
to the likelihood doc's toolchain finding) giving the one line a
developer must put in `~/.R/Makevars` or a file pointed at by
`R_MAKEVARS_USER`. If the lane prefers to keep a local file, it must at
least be ignored -- and it currently is NOT. Punch list P9.

**Concretely, the `.gitignore` line is wrong.** The diff adds
`dev/stan-cache-methods/`, which correctly ignores the 34 MB of compiled
programs. It does NOT ignore `dev/makevars-stan-methods`, so that file
is untracked-and-unignored and will be swept into the next `git add -A`.

---

## 6. Suites

Each in its own process, worktree loaded with `pkgload::load_all()`,
private library, the lane's cache copied read-only.

| file | gate | tests | assertions | failed | skipped | errors | wall |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `test-brms-methods.R` | `FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true` | 43 | **814** | 0 | 0 | 0 | 121.8 s |
| `test-message-uniqueness.R` | `NOT_CRAN=true` | 1 | **6** | 0 | 0 | 0 | 10.8 s |
| `test-bracket-access.R` | `NOT_CRAN=true` | 3 | **7** | 0 | 0 | 0 | 7.2 s |

43 tests and 814 assertions match the lane's claim exactly, including
the per-test distribution. Message-uniqueness 6 and bracket-access 7
match. `test-brms-likelihood.R` row files were out of scope and were not
run.

---

## 7. What the 814 assertions actually assert

Measured, not estimated: the file was re-run under a `ListReporter` and
every expectation was attributed to its source line, then each call site
was bucketed. All 814 are accounted for; none is unclassified.

| bucket | assertions | share |
| --- | --- | --- |
| exact cross-package identity at 1e-8 (`E`) | **549** | 67.4% |
| structural: shapes, names, classes, containers, refusals, mechanism preconditions (`S`) | **188** | 23.1% |
| divergence pins (`D`) | **77** | 9.5% |

`expect_exact_num()` contributes two expectations per call (a length
check and a relative-difference check), which is where the 204 static
`expect_*` calls become 814 at run time: the loops over 23 shapes and
their dpars do the rest.

Of the 549 exact identities, **543 have brms on one side**. The other
six are frmtmb-internal consistency identities (`coef()` is `fixef()$mu`
at line 451, the pearson denominator at 537, `predict(type =
"response")` is `fitted()` at 570). So the headline claim -- that this
tier compares frmtmb against brms exactly, at scale -- is honest.

Of the 77 divergence pins, **45 belong to the ten `(D)` items** and 32
to the `(P)` and `(C)` ones. So fixing all ten defects would require
editing 45 assertions, 5.5% of the file, all of them labeled and all of
them in blocks whose comments already say what they pin. That is the
number to put in front of whoever owns the fixes.

A property worth naming: the `(P)` divergences are pinned by EXACT
identities, not by inequalities. `p * n` equals frmtmb's answer
(line 206), `log(nat)` equals frmtmb's link-scale sigma (line 338),
`plogis(eta)` equals brms's theta (line 284), the brms `log_lik` sum
equals frmtmb's conditional sum (line 103). A paradigm difference that
can be written as an exact identity is a difference that is understood,
and the tier consistently reaches for one. That is the best thing in
this file.

---

## 8. Punch list

| # | file and line | what | severity |
| --- | --- | --- | --- |
| P1 | `tests/testthat/helper-brms-methods.R:501`, `505-507`, `534-537`, `543-545`, `551-553`, `561-565`, `577-581` | The exclusion lists are the ONLY place a fix can pass silently: fix D1 and `r16`/`rC16` must be put back into `brms_ce_shapes()` or the tier stops covering them, with nothing failing to say so. Name the finding each exclusion belongs to and state that fixing it requires re-adding the shape. | **medium** |
| P2 | `dev/brms-methods-tests.md`, all finding sections | 9 of 22 findings carry a bolded "Which is right"; 13 carry nothing. Add a one-line verdict -- defect, paradigm, or open choice -- to every one. The brief's premise that `dev/brms-likelihood-tests.md` already does this per divergence is false (one `Status:` line, at line 3), so this doc would be establishing the convention, not following it. | **medium** |
| P3 | `dev/brms-methods-tests.md:287-292`, plus a new brms-free test | Finding 1b stops at "the mechanism predicts the same result" for the hurdle families. Measured here: `hurdle_poisson`, `ce[1]/epred[1] = 0.536`, `max ratio - 1 = 1.196`. The defect is visible with frmtmb ALONE (`conditional_effects()` default vs `predict(type = "response")`), so it belongs in the ordinary suite, not behind the Stan gate. | **high** |
| P4 | `dev/brms-methods-tests.md:384-392` | Finding 3's own correction is shape-dependent. On an ordinal fit, `categorical =`, `int_conditions =` and any unknown argument produce **no warning at all**, because the categorical branch does not forward its dots to `predict()`. That makes divergence 9's ignored argument genuinely silent. | **medium** |
| P5 | `tests/testthat/test-brms-methods.R:797-801` | One comment: fixing the `signif` rounding also breaks the alignment key at line 801, so `expect_false(anyNA(m))` will fail alongside the three named pins. Expected, but it reads like an unrelated failure. | low |
| P6 | `.github/workflows/brms-likelihood.yaml:88-90`, and `1-2`, `36`, `113-115` | Before the methods tier joins this job, add `test-brms-methods.R` and `helper-brms-methods.R` to the `tools::md5sum()` vector. `actions/cache` saves only on a primary-key MISS, so without this the one new program (`rfac`) is compiled and discarded on every run forever. The run step, the job comment "Only this tier runs" and the workflow name also need updating, with the second `test_file()` in its own process. | **high, if wired in** |
| P7 | `NEWS.md:16-18` | "runs 15% above the mean of the data" is wrong twice: 15.2% is the error at the FIRST grid point (it reaches **268%** at the last), and the comparison is against the expected response, not the mean of the data. | **medium** |
| P8 | `dev/brms-methods-tests.md:88-92` | "equals `fitted(fit)` to **0** (exactly, not to a tolerance)" on "every shape" is true on r1 and false on r16 (4.44e-16, `identical()` FALSE). | low |
| P9 | `dev/makevars-stan-methods`, `.gitignore` | The file is untracked AND unignored: `.gitignore` gained `dev/stan-cache-methods/` only. Either delete it and document the recipe (recommended, see section 5) or ignore it. As it stands a `git add -A` commits a machine-local toolchain file. | **medium** |
| P10 | `tests/testthat/helper-brms-methods.R:48-58` | Record the measured fallback in the comment: without `exclude_pars`, `rename_pars()` does not error, the object gains `r_1`/`z_1`/`L_1`/`Cor_1` (80 extra names on sleepstudy), and `ranef()` and `posterior_epred()` stay `identical()`. Saves the next reader the measurement. | low |

No punch-list item is in the tier's assertion logic. Every one is
documentation, CI wiring, or a file-hygiene slip.

---

## 9. Verdict

**GO WITH FIXES.**

The mechanism is sound and I reproduced it end to end, including the
byte-identical route (a) comparison the brief asked for. The suite
reproduces exactly -- 43 tests, 814 assertions, 0 failures, 0 skips --
and so do message-uniqueness (6) and bracket-access (7). Nothing under
`R/` changed, confirmed by diff. Two-thirds of the assertions are exact
cross-package identities with brms on one side, which is a real and
unusual achievement for a package of this kind.

The fixes are P1 through P9. None touches test logic; P3, P6 and P7 are
the ones I would not merge without.

**The tier should merge BEFORE the (D) items are fixed, not after.**
Every one of the ten defects is pinned so that a fix flips a NAMED
expectation in a block whose comment already says what it pins, and the
45 assertions involved are labeled. Holding the tier back until `R/` is
repaired would mean the ten defects -- including a
`conditional_effects()` curve that runs 268% high on a zero-inflated fit
and 46% low on a hurdle fit -- stay undocumented and unpinned in the
meantime, which is the state this lane was created to end. Merge it,
fix P1 through P9, then open the ten `(D)` items as their own work with
this document as the specification.

### Edits made to the worktree

One file created: `dev/review-brms-methods.md` (this document). No other
file in `C:\Users\adf44\source\r\frmtmb-wt-brms-methods` was modified,
and nothing was committed. `C:\Users\adf44\source\r\frmtmb` was not
touched and remains clean at `5dfdd84`.

---

# Punch re-check, 2026-09-05

Verified against the diff and by running, not against the lane's prose.
Same rules as the first round: private library, the lane's cache read
from a private copy, nothing committed, main untouched.

Lane state: HEAD still `564e185`; the tracked diff is now
`.github/workflows/brms-likelihood.yaml`, `.gitignore`, `NEWS.md`;
untracked is the two test files plus the two dev docs.
`dev/makevars-stan-methods` is gone. Main is clean at `5dfdd84`.

## All ten applied

| item | verified how | result |
| --- | --- | --- |
| P1 | read the table, ran the guard, flipped a probe | **done**, see below |
| P2 | counted the Class verdicts and cross-checked every one against my own table | **done**, 26 findings, 10/8/8 |
| P3 | ran the hurdle block's data outside the tier | **done**, numbers confirmed |
| P4 | read `:830-867`, re-ran my own probe | **done**, and better than I had it |
| P5 | read `:932-936` | **done** |
| P6 | parsed the YAML, walked both run steps | **done**, one residual |
| P7 | read `NEWS.md:17-25` | **done** |
| P8 | read `dev/brms-methods-tests.md:88-94` | **done**, goes further than asked |
| P9 | file deleted, recipe at doc `:160-176` | **done** |
| P10 | read `helper-brms-methods.R:48-57` | **done** |

## P1: the table and the guard

`brms_exclusions()` at `helper-brms-methods.R:551` returns **32 rows**,
every one carrying `list`, `key`, `finding` and `class`:

```
                        C  D  P
  brms_ce_shapes        1  9  1
  brms_dpars_of         0  2  0
  brms_linpred_shapes   3  0  0
  brms_meanlink_shapes  0  0 10
  brms_resid_shapes     0  0  6
```

All **11 defect rows** land in the two lists `brms_exclusion_agrees()`
dispatches on (`:605`), so every one reaches a live probe. All four
derived lists now come from the table with no hard-coded exclusion left
anywhere (`helper:664, 679, 687, 693`), and the lengths are unchanged
from before the refactor: **ce 12, linpred 20, meanlink 13, resid 17**
against 23 registered shapes. Every key names a registered shape.

**The guard discriminates.** I copied the two helpers into a scratch
directory, injected one line into the COPY only,
`if (identical(shape$name, "r16")) return(TRUE)` at the head of
`brms_ce_agrees()`, and ran the guard block against it. The lane's own
file was left untouched, asserted in the same script. Result: exactly
one failure out of 45 expectations, and it is named:

```
1. Failure ('test-rbm-guard.R:40:7'): every deferred exclusion still has a live defect
finding 1b looks fixed: r16 now agrees with brms, so drop its row from
brms_exclusions() and let brms_ce_shapes() cover it again
```

Finding, key and the list to edit, all three. The other ten defect rows
still passed. This closes the silent-pass hazard that was the whole of
P1.

## P6: the workflow, parsed

`yaml::read_yaml()` parses it. Walked:

- workflow `name: brms-tiers` (`:46`), job id `brms-tiers` (`:49`), and
  **no job-level `name:`**.
- The FILE is still `brms-likelihood.yaml`, which is correct: the
  self-referencing `paths:` filters at `:30` and `:38` still match it.
  Renaming the file would have stopped the workflow triggering on its
  own edits.
- The cache key hashes **exactly four** files,
  `test-brms-likelihood.R`, `helper-brms.R`, `test-brms-methods.R`,
  `helper-brms-methods.R`; all four exist and nothing stale is in the
  list. I checked what else could add a program:
  `test-brms-agreement.R` is the only other file mentioning
  `make_stancode`, it has zero `brms_stan_model()` calls and is not run
  by this job, so four is the right set.
- `cache.path` is `dev/stan-cache` and `FRMTMB_STAN_CACHE` points at the
  same place.
- `any::MASS` present (`:85`).
- **Both** run steps: the same `stopifnot` asserting
  `FRMTMB_BRMS_FIT_TESTS` and `NOT_CRAN`, `load_all`, failure on
  `failed`/`error`, and `if (any(df$skipped)) stop(...)`. **The second
  step cannot pass on skips.** The methods step additionally requires
  `MASS`.

**Branch-protection consequence.** For a GitHub Actions job the required
status check is named by the job's `name:` if set, otherwise by the job
id. There is no job-level `name:`, and the job id changed from
`brms-likelihood` to `brms-tiers`. So a branch protection rule that
requires `brms-likelihood` **will never be satisfied again**, and the
failure mode is the confusing one: the check is never reported, so
protected pull requests sit at "Expected, waiting for status to be
reported" indefinitely rather than going red. The rule has to be renamed
in the same change that merges this. The in-file comment at `:44-45`
flags that it must be updated but not that the symptom is a hang.

## P3 and P4, re-measured

Hurdle, on the tier's own data (`test-brms-methods.R:783`):

```
ratio at grid point   1: 0.86193188603
ratio at grid point 100: 3.37134980095
crosses 1 along the curve: TRUE
ce == predict(type = "conditional"): 0
method = "predict": "method = 'predict' needs a family with a simulator"
```

0.862 and 3.371 exactly as claimed, the sign change confirmed, and the
zero-inflated workaround correctly shown to be unavailable. The block is
gated on `skip_on_cran()` alone, so it runs without Stan and without
brms, which is what P3 asked for.

P4 at `:830` is better than my own probe was. I used `tryCatch`, which
catches the first warning and lets the second escape; the lane found
that a Wald band calls `predict()` twice and used `capture_warnings()`
with `expect_gte(length(w), 1L)` plus an `all(grepl(...))` over the
whole vector. Three `expect_silent()` calls pin the ordinal path.

## Suites, one per process

| file | tests | assertions | failed | skipped | errors | wall |
| --- | --- | --- | --- | --- | --- | --- |
| `test-brms-methods.R` | **46** | **872** | 0 | 0 | 0 | 126.4 s |
| `test-message-uniqueness.R` | 1 | **6** | 0 | 0 | 0 | 5.2 s |
| `test-bracket-access.R` | 3 | **7** | 0 | 0 | 0 | 2.4 s |

46 / 872 / 0 matches the claim. The delta from the first round is +3
tests and +58 assertions: the guard (45), the hurdle block, and the
ordinal-no-dots block. The guard's 45 breaks down as 32 key checks, one
`expect_setequal` on the class vocabulary, one `expect_gt` on the live
count, and 11 probe outcomes.

## Residual items

| # | file:line | what | severity |
| --- | --- | --- | --- |
| R1 | `.github/workflows/brms-likelihood.yaml:101`, comment `:92-97` | **The brms version is not in the cache key.** rstan is; brms is not, and brms is what GENERATES the Stan code. The job installs `any::brms` unpinned and runs on a weekly cron. On any brms release every program's content-addressed key changes, all 25 recompile, and the `actions/cache` primary key is unchanged, so post-job sees an exact hit and saves nothing: 25 recompiles on every run until someone edits one of the four files. This is the P6 failure mode one variable over. One line: add `packageVersion("brms")` beside rstan's at `:101`. | **medium** |
| R2 | `tests/testthat/test-brms-methods.R:151-159` | `brms_exclusion_agrees()` returns `NA` for any list it has no dispatch for, and the loop is `if (isTRUE(agrees)) fail() else succeed()`, so **`NA` passes silently**. 19 of the 32 rows sit in the three non-dispatched lists. All 19 are P or C today so nothing is unprobed now, but reclassify one to D and the guard succeeds without ever probing it. One line: fail when a D row's probe returns `NA`. | low-medium |
| R3 | `helper-brms-methods.R:551` vs `dev/brms-methods-tests.md` | The class of a divergence is stated twice, in the table's `class` column and in the 26 Class lines, in two files, with nothing tying them together. They agree today; I checked all 26. | low |
| R4 | workflow `:144`, `:174` | Neither run step asserts a test COUNT (`nrow(df)`), so a file that produced zero tests would pass. testthat errors on a parse failure, so this is defense in depth rather than a live hole. | low |
| R5 | workflow `:174` | The second run step has no `if: always()`, so a failure in the identity tier aborts the job and the methods tier's status is never reported. Fail-fast is defensible now that they share a job; worth one deliberate line either way. | low |

## Updated verdict

**GO.**

All ten punch items are applied and verified against the diff. P1 in
particular is done properly: one table, every defect row probe-backed,
the derived lists provably unchanged, and a guard I demonstrated will
fail with a named, actionable message the moment a defect is repaired.
P6 is correct on every point I weighted it on, four files, no stale
entries, both steps gated, the second unable to pass on skips, and the
lane made the right call leaving the file name alone.

The five residuals are all small and none blocks the merge. **R1 is the
one to fix before the workflow lands**, because it silently reintroduces
the cost P6 was written to remove, on a weekly schedule, and nobody
would notice it except as a slow job.

My first-round conclusion is unchanged: merge the tier before the ten
`(D)` items are repaired. The guard makes that recommendation stronger
rather than weaker, because the exclusion lists can no longer go stale
without a test saying so.

### Edits made to the worktree in this round

None. This document is the only file I have written in
`C:\Users\adf44\source\r\frmtmb-wt-brms-methods`, and the probe-flip
experiment ran against copies in the scratchpad (`rbm-scratch-tt/`),
with the lane's own `tests/testthat/helper-brms-methods.R` asserted
unmodified. Nothing was committed; main remains clean at `5dfdd84`.
