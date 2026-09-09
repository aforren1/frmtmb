# The per-subject non-decision-time bound

Item 1.0a of `dev/extension-gaps-plan.md`, in frmtmb.eam. Written on
the `wt-ndt` worktree off 92e9330. R 4.6.1, Windows 11, 16 logical
cores. frmtmb 0.55.1 and frmtmb.eam both built from this worktree into
a private library at `C:/Users/adf44/source/r/ndt-lib`; the reference
build of 92e9330 at `C:/Users/adf44/source/r/reflib-r2` supplies every
BEFORE number.

Every script named below is in `dev/ndt-scripts/` and carries its own
seed in its header. Its log is beside it.

## The claim, and whether it holds

**A per-subject non-decision-time bound recovers the eam scale
design.** It does, and it beats the target the plan set.

| | Phase 0, one global bound | this change, `ndt_group(s)` | the plan's target |
|---|---|---|---|
| convergence code | 1 | 0 | 0 |
| maximum absolute gradient | 1.25e11 | 9.87e-04 | |
| positive definite Hessian | no | YES | yes |
| standard errors that are `NaN` | 7 of 7 | 0 of 7 | |
| log-likelihood | -7148.81 | **-7003.01** | at least -7027.4 |
| `sd(ndt)` natural, truth 0.02629 | 0.0058 | **0.02543** | recovered |
| population `ndt`, truth 0.25 | 0.2261, se `NaN` | 0.2469, se 0.00196 | |
| `ndt` z against the truth | (none) | 1.60 | |
| condition effect, truth 0.9 | 0.891, no interval | 0.911 (0.858, 0.964) | covers |
| `sd(mu)`, truth 0.35 | 0.334 | 0.336 | |
| `sd(bs)` on the log link, truth 0.20 | 0.173 | 0.143 | |
| whole `frm(se = TRUE)` call | 269 s | 108 s | |
| subjects below their OWN fastest response | (n/a) | 30 of 30 | |
| tightest such margin | (n/a) | 27.3 ms | |
| `diagnose()` | reports nothing | nothing to report | |

Construction: `extensions/frmtmb.eam/tests/testthat/test-scale.R`, row
`eam`, seed 20260908, run by `dev/ndt-scripts/ndt-scale-run.ps1` under
`FRMTMB_SCALE_TESTS=true NOT_CRAN=true`, one row per fresh process.
The records are in `dev/ndt-scripts/ndt-scale.tsv`.

The log-likelihood is 145.8 units better than the global bound and
24.4 units better than the bound-lifted arm, at the same parameter
count.

**Run twice, in independent processes, on two builds of this change**
(the second after the starting-value and mixture fixes below): the log
likelihood, every estimate and every margin agree to all six digits
recorded. `-7003.01`, `ndt` 0.246867 with se 0.0019639,
`sd(ndt)` 0.0254253, tightest margin 27.3347 ms. The wall clocks differ,
108.0 s and 151.9 s, because the second ran beside two other R
processes; read them as upper bounds.

**The instrument checks out.** The same runs reproduce the
`eam-unbounded` row from Phase 0 to every digit `dev/scale-findings.md`
records: log-likelihood -7027.37, population `ndt` 0.2437 with se
0.0045, `sd(ndt)` natural 0.0225, condition effect 0.909 (0.857,
0.962). That row also measures the plan's "22 ms to spare" claim
directly: the tightest subject at that arm's optimum sits **21.9 ms**
below its own fastest response. The plan's arithmetic was right.

## Punch round 1, and what it changed about the design

`dev/reviews/2026-09-09-ndt.md` raised five blockers. All five shared
one root: a shipped parameter's natural scale had changed and the
callers of that scale had not been swept. The fix is not five patches.
It is a smaller change than the one first shipped.

**The scalar bound went back into the LINK.** A per-ROW bound cannot be
a link parameter, and that is still true; a per-DATA-SET bound can, and
was, for the whole of 0.6.0. So there are now two parameterizations and
only the opt-in one moved:

- no `ndt_group()`: one bound, in the link, `ndt` is a TIME. Byte for
  byte the pre-change path.
- with `ndt_group()`: a per-row bound, a plain logit, `ndt` is a
  FRACTION, and the density multiplies.

Because `log(ndt / (ub - ndt))` and `log(f / (1 - f))` are the same
number, moving the bound back into the link changes no linear
predictor, no log-likelihood and no coefficient. What it restores is
every reading of `ndt` on its RESPONSE scale. Measured on both builds,
`dev/ndt-scripts/ndt-round1.R`, seed 4242, 400 rows, `min(rt)` =
0.3145218:

| | 0.6.0 | this branch |
|---|---|---|
| `predict(dpar = "ndt", type = "response")` | 0.27528808 | 0.27528808 |
| log-likelihood | -205.767667983 | -205.767667983 |
| `prior(normal(0.30, 0.01), class = "ndt")` -> `ndt` | 0.28073215 | 0.28073215 |
| `prior(normal(0.10, 0.01), class = "ndt")` -> `ndt` | 0.16063194 | 0.16063194 |
| `prior(normal(0, 0.1), class = "Intercept", dpar = "ndt")` -> eta | 0.44280466 | 0.44280466 |
| `bf(ndt = 0.2)` with `max_ndt = 0.25`, log-likelihood | -237.065325 | -237.065325 |
| the same, `fitted()[1]` | 0.793986 | 0.793986 |
| `bf(ndt = 0.2)` with a bare `wiener()` | refused | refused |

That is **B2 and B3 closed**, not mitigated: the numbers the review
measured at 0.296152 against 0.106955, and -255.9 against -326.3, are
now equal on both builds. It also closes N3, because
`predict(dpar = "ndt", type = "response")` means seconds on `gddm()`
and on every other family that has one bound.

The evidence that the scalar path is untouched is the diff itself. Nine
test files were edited in the first round to follow `ndt` onto the
fraction; eight of them are now reverted to their 92e9330 state, and so
is the helper that converted times to fractions for them. What remains
in `tests/` is `test-scale.R` and the new `test-ndt-bound.R`.

### Did I take the reviewer's `post$dpar_response` hook (F2)?

**No, and the reason is a measurement rather than a preference.** The
hook makes `predict(dpar = "ndt", type = "response")` report seconds
for a scalar bound. Putting the bound back in the link does that AND
three things the hook cannot:

- a `prior(class = "ndt")` is placed through `coef_placement()` with
  the dpar's own `link$linkinv` and `link$mu_eta`
  (`frmtmb/R/priors.R:1898-1925`), not through the reporting hook, so
  the hook leaves B2 exactly where it was;
- a `bf(ndt = )` constant runs `fam$links$ndt$linkfun(constant)` at
  PARSE time (`frmtmb/R/parse.R:1264-1268`), so the hook leaves B3
  where it was too, including the case the review flagged as worse,
  where a previously loud refusal had gone quiet;
- `st`'s bound is twice `ndt`'s and has no hook of its own.

The hook is also the strictly more delicate of the two: the review's
own notes on it record a `deriv` written against the wrong scale
reporting a standard error of 0.2595 where 0.0062 was right, and
`ndt_time()` double-scaling to 0.098872 against a truth of 0.295562.
The link carries its own `mu_eta`, so the delta-method standard error
is right by construction and there is nothing to keep in step.

So F2's finding is accepted in full and its remedy is superseded by a
larger one. The one thing the hook would still buy is seconds from
`predict()` under `ndt_group()`, where the bound is per row; that stays
a fraction, `ndt_time()` is the route, and `?wiener` says so.

### B1, and the case where the label is absent

The floor table is now keyed on the group's LABEL. The coercion hashes
the label with two polynomial rolling hashes packed into one double
below 2^49, so a factor, a character vector, a logical and integer
codes that name the same groups all give the same fit, and no code
depends on which other labels were present. The character and logical
refusals are gone with the defect that motivated them.

`dev/ndt-scripts/ndt-round1.R`, seed 808, the review's own design: the
same two rows `s2` and `s3` through four framings.

| newdata | `ndt_time()` | fitted mean |
|---|---|---|
| subset, levels as fitted | 0.286143, 0.285354 | 0.695233, 0.694445 |
| ... then `droplevels()` | 0.286143, 0.285354 | 0.695233, 0.694445 |
| ... then `relevel(ref = "s3")` | 0.286143, 0.285354 | 0.695233, 0.694445 |
| hand-built, levels reversed | 0.286143, 0.285354 | 0.695233, 0.694445 |

The review's `droplevels()` row read 0.279511, 0.286143. And the
character and integer spellings of the same grouping now give the same
log-likelihood, -449.3260018, as the factor.

**The lookup cannot fail open.** `floors[label]` returns `NA` on a
miss, and an `NA` bound reaching the density would be the same defect
one layer down, so the case is constructed rather than reasoned about.
A group label the fit never saw, on all three paths that reach the
lookup:

| path | result |
|---|---|
| `ndt_time(newdata =)` | refused |
| `predict(newdata =, type = "response")` | refused |
| `frm()` on the new rows with the fitted family | refused |

and a missing group, and a newdata frame with the column gone, are each
refused by name. `tests/testthat/test-ndt-bound.R` pins all of it.

The review notes `rlddm()` will inherit this coercion. It now inherits
the fixed one.

### N0, the kept bound

All four families now KEEP a settled bound. `wiener()` and
`wiener_gng()` rebuild themselves from their config at finalize, so the
install's own guard could not see the carried bound and is handed it
explicitly through `ddm_ndt_keep()`. Refinalizing without the row that
SET the bound, `dev/ndt-scripts/ndt-round1.R`, seed 4242:

| family | all rows | refinalized without the fastest row | kept |
|---|---|---|---|
| `wiener` | 0.3145218 | 0.3145218 | yes |
| `lba` | 0.3145218 | 0.3145218 | yes |
| `rdm` | 0.3145218 | 0.3145218 | yes |
| `wiener_gng` | 0.3145218 | 0.3145218 | yes |

The review measured `wiener()` at 0.334522 against 0.345458 before
this. The claim in this file about `influence()` is now true rather
than aspirational, and the test asserts that the dropped row is the one
that set the bound, so a re-derivation would be visible.

### A defect the punch round found, and fixed

**An `ndt_group()` no family reads was accepted and ignored.**
`mixture()` takes the union of its components' allow-lists and never
runs their `family_finalize`, so
`mixture(wiener(max_ndt = 0.4), lognormal())` with `ndt_group(g)`
assembled, fitted, and scored every row against the one scalar bound
with the grouping unread. That is the failure this package already
fixed once, when `rt | dec(u) + vint(1 - u)` fitted with the second
column unread.

A per-group bound leaves a table behind, so its absence is exactly the
condition to refuse on. `ddm_check_ndt_group_read()` is registered
through `frmtmb_register_frame_check()` and refuses at frame assembly,
naming the response and saying what a mixture can do instead. It cannot
false-alarm on a grouped model, because there the table exists.

### The rest of the punch list

| item | what changed |
|---|---|
| B4 | the NEWS "to the last bit" sentence is replaced by what was measured: the objective and its gradient agree with 0.6.0's to within 2.81 ulp at a FIXED parameter vector, and an optimizer on a flat ridge can still land elsewhere, with the `rdm` case quoted at a coefficient moving 0.068 while the log-likelihoods differ by 2.7e-9 |
| B5 | the `sd(ndt)` caveat is now a NEWS bullet and a `?wiener` section, with the floors-only estimator's 0.02748, the 100-trial identity, and the 8.44 log-likelihood units and 11.90 to 7.67 ms the component does buy |
| N1 | both absolute tolerances are out of `test-scale.R`. The 20 ms one is a `scale_z()` against the fit's own standard error, which is what that helper exists for. The `-7027.4` one is REMOVED rather than rewritten: it is a log-likelihood carried from another arm at another seed, the tier already RECORDS `logLik` in its own output, and no self-referential form of it exists that does not cost a second 12,000-row fit. The test comment says that |
| N2 | `ndt_time(newdata =)` evaluates the grouping itself, so an absent column reached R's `object 's' not found`. It now names the term and says to supply it |
| N4 | a group of one trial is still accepted, because any trial-count threshold would fire on correct models. What changed is that the counts are recorded, at `family(fit)$ndt_bound$sizes`, and `?wiener` states the bound's quality against the trial count |
| N5 | the `NA` branch of the coercion says where it fires: `frm()`'s `na.action` drops the row first, and `ndt_time(newdata =)` is the path that reaches it |
| N6 | the over-length lines in this file are all markdown tables, which `dev/scale-findings.md` also carries; left as they are |
| F1 | the `gddm()` exclusion is restated on SCOPE in `frm_compat("gddm")`, in `?gddm` and in NEWS: `?gddm` already states the contract, the family could run the same conservative `tapply` check it runs for `vreal()`, and what a per-group bound would cost there is one Fokker-Planck solve per group |


## What changed

`ndt_group()` is a new addition term: the grouping the bound is taken
per. With it the bound is the row's own group's fastest response and
`ndt` is a FRACTION of that bound, on a plain logit, which the density
multiplies back out. Without it the bound is one number, the whole data
set's fastest response or `max_ndt`, it stays in the link as it has
since 0.5.0, and `ndt` is a TIME. `max_ndt` and `ndt_group()` together
are refused, and so is an `ndt_group()` no family reads.

Four families take the grouping: `wiener()`, `lba()`, `rdm()`,
`wiener_gng()`. `gddm()` does not, on scope; see below.

### Files

| file | what |
|---|---|
| `R/ddm-shared.R` | the whole mechanism: `ddm_scaled_logit()` moved here and shared, the label hash and the `ndt_group()` coercion, the bound derivation and its refusals, the per-row scaler, the slot wrappers, the starting values, `ddm_ndt_install()` with its three states, `ddm_ndt_keep()`, and the exported `ndt_time()` |
| `R/wiener-family.R` | `ddm_scaled_logit()` moved out to be shared; `ddm_finalize()` takes the family and the aterms, keeps a settled bound, and hands the install the per-group table when there is one; the help's non-decision-time sections rewritten |
| `R/lba.R` | `lba_finalize()` takes `aterms` and defers to the shared install; help |
| `R/rdm.R` | the same, through `ddm_ndt_finalize()`; the wrap reaches `lcdf` and `lccdf` |
| `R/wiener-gng.R` | the same, with the bound over GO rows only; `gng_go_which()` split out so `ndt_group()` is subset by the same rows |
| `R/gddm.R` | help only: it keeps its scaled logit, and the exclusion is restated on SCOPE |
| `R/zzz.R` | registers the `ndt_group` addition term, its five compatibility rows, and the frame check that refuses a grouping no family read |
| `NEWS.md` | the bullets, under an unreleased heading |
| `tests/testthat/test-ndt-bound.R` | NEW, 79 assertions |
| `tests/testthat/test-scale.R` | the eam row carries `ndt_group(s)` and the recovery assertions; the sv row's condition assertion is replaced by a measurement; both absolute tolerances are gone |
| every other test file | UNCHANGED from 92e9330. Eight were edited in the first round to follow `ndt` onto the fraction and all eight are reverted, which is the evidence that a model without `ndt_group()` did not move |
| `R/frmtmb.eam-package.R` | imports `frmtmb_register_frame_check` |
| `man/*.Rd`, `NAMESPACE` | roxygenised |

### The design, and the routes rejected

The plan named the shape: `ndt` becomes a fraction of the group's floor
inside `lpdf`. What it did not settle is how the grouping arrives and
how every other slot keeps working. Four routes were considered.

1. **A per-row bound inside the LINK.** Rejected, and it is why the
   change is breaking at all. `linkinv(eta)` is handed a vector whose
   length is the rows being predicted, so a captured per-row bound
   recycles silently against newdata of another length, and
   `linkfun()` is called once on a SCALAR starting value
   (`?frmtmb_family`, slot call order step 4) where a vector bound has
   no meaning.

2. **`vint()` as the grouping.** Rejected for `wiener()` by its own
   declaration: `exclusive_aterms = list(c("dec", "vint1"))` refuses a
   model supplying both, and `dec()` is how a Wiener model carries its
   boundary. The grouping would have had to be `vint2` under a
   `vint(upper, group)` spelling for some families and `vint1` for
   others, which is the positional ambiguity `gddm()`'s own comment
   already calls out.

3. **A registered addition term, `ndt_group()`.** Taken.
   `frmtmb_register_aterm()` is the seam core provides, the term joins
   the compatibility vocabulary automatically, and
   `aterms_for_newdata()` re-evaluates it on newdata with the
   registered coercion, so a prediction carries its own grouping.

4. **Reporting the time through `post$dpar_response`.** The hook takes
   `value(dpars, dpar)` (`frmtmb/R/predict.R:1365`) and receives no
   `aterms`, so it cannot reach a per-row bound. That much is right and
   is why `ndt_time()` exists. What the first round got wrong is the
   conclusion it drew for the SCALAR bound, where the answer is not a
   hook at all but leaving the bound in the link; see the punch round
   above.

The per-row bound reaches the density as DATA. `aterm_data()` puts
`ndt_floor` beside the addition-term values, so the tape gathers a
vector instead of doing a lookup per evaluation. It is read out of the
table captured when the family was finalized rather than recomputed
from whatever `y` it is handed, which matters because
`assemble_frame()` runs on new rows for `influence()`, `frm_simulate()`
and the prior-predictive path: a leave-one-out refit that dropped a
group's fastest trial would otherwise have silently changed the model.
`ddm_ndt_install()` is idempotent for the same reason, and the guard is
the recorded bound itself.

Under a grouping, every slot that consumes a non-decision time is
WRAPPED rather than rewritten: `lpdf`, `lcdf`, `lccdf`, `sim` and
`post$mean_fn` are given `dpars` with `ndt` (and `st`) already
multiplied out, so each density below is still written in times and
knows nothing about the fraction. Without a grouping no wrapper is
installed at all and the slots are the ones 92e9330 shipped.
The wrappers preserve the formal COUNT, because frmtmb reads
`length(formals())` to decide whether a slot takes the family-level
extras (`frmtmb/R/families.R:1010`, `frmtmb/R/objective.R:38`); a `...`
formal would have been read as one argument. They also drop the
reserved `.eta_ndt` entry with the value it belongs to, so a robust
dpar accessor cannot recover a linear predictor that no longer maps to
the value beside it.

### `max_ndt` under the new scheme

Unchanged in meaning: one absolute upper bound, in the units of the
response, applied to every row, with `ndt` a time measured against it
exactly as before. Its
refusal is unchanged too: above the fastest response it is refused
outside a mixture, and `allow_unreachable = TRUE` still lifts that for
the mixture case. What is new is that `max_ndt` with `ndt_group()` is
refused, because the two set the same bound to different things.

## The behavioral failure, seen

The lane rule is that a test pinning a bug is worthless unless the bug
was seen. Two forms are recorded.

**The weak form.** `tests/testthat/test-ndt-bound.R` run against the
reference build of 92e9330 (frmtmb.eam 0.6.0), with the failure cap
lifted and `package = "frmtmb.eam"` so the internals resolve:
`FAIL 14 | WARN 0 | SKIP 0 | PASS 12`. Of the 14, eight are
`Addition term ndt_group() is not supported`, two are
`could not find function "ndt_time"`, and four are the bound and link
assertions failing on their values.

The first time this was run the reporter said
`FAIL 10 | WARN 0 | SKIP 0 | PASS 0` with "Maximum number of failures
exceeded". Ten was the CAP and zero was the namespace, not the file:
without `package =` every test touching an internal errors before it
can pass. Both numbers were wrong and the lane's own standing rule is
what caught them.

**The strong form**, a difference in the answer rather than in the
vocabulary: `dev/ndt-scripts/ndt-before-after.R`, seed 909. Two groups
of 250 trials, true non-decision times 0.18 and 0.40, so the slow
group's truth is above the FAST group's fastest response and a single
bound cannot express it.

| | global bound | per-group bound |
|---|---|---|
| global `min(rt)` | 0.21423751 | the same data |
| the two groups' own `min(rt)` | 0.21423751, 0.43508186 | |
| fitted `ndt` per group | 0.13517323, **0.21423751** | 0.17635632, 0.40275205 |
| slow group's relative error | **46.4 percent** | 0.69 percent |
| fast group's relative error | | 2.02 percent |
| log-likelihood | -158.327 | -48.218 |
| convergence code | 0 | 0 |
| maximum absolute gradient | 2.75e-08 | 6.41e-06 |

The slow group's fitted non-decision time under the global bound is
0.21423751 against a bound of 0.21423751: it is AT the wall to within
1e-6, and the fit reports convergence. 110.1 log-likelihood units, at
the same parameter count.

That run also shows the OTHER face of the Phase 0 standard-error
finding. On the LINK scale the pinned coefficient's Wald interval is
(-16724, 16770), a standard error of 8544 on a coefficient of 23.0; on
the RESPONSE scale the delta method multiplies by `mu_eta`, which has
vanished, and Phase 0 measured 7.2e-06. The same pathology reads as an
exploded standard error in one coordinate and a collapsed one in the
other. That is worth carrying into the `diagnose()` seam the plan has
filed: a check that tests either alone is fooled by one of them.

**The control that makes the strong form a before/after.** With no
`ndt_group()` the new parameterization is the old one written
differently, so `dev/ndt-scripts/ndt-smoke.R` runs the same six models
on both builds, seed 4242:

| quantity | reference 0.6.0 | this worktree |
|---|---|---|
| plain wiener log-likelihood | -249.5516581329 | -249.5516581329 |
| `ndt.(Intercept)` on the link | 1.649938274 | 1.649938274 |
| `mu.condb` | 0.9366881774 | 0.9366881774 |
| `fitted()[1]` | 0.7721772658 | 0.7721772658 |
| `variability = "sv"` log-likelihood | -249.3745150425 | -249.3745150425 |
| `rdm` log-likelihood | 50.0016052927 | 50.0016053032 |
| `lba` log-likelihood | 87.1928304152 | 87.1928304152 |
| `wiener_gng` log-likelihood | -116.0024236055 | -116.0024236055 |
| `predict(dpar = "ndt", type = "response")` | 0.2307934893 | 0.8388827078 |

Every number meant to be identical is identical to the digits printed,
including the `ndt` coefficient itself.

**Ten printed digits are not a bit, and the review was right to say
so.** Under `identical()` and a ulp count over fourteen models, nothing
is bitwise except `gddm` and `simulate()`; the worst cases are an `rdm`
coefficient at 0.0681 and an `rdm` vcov at 3.94e6. The cause is settled
in the change's favour: at a FIXED parameter vector the objective and
its gradient agree to 2.81 ulp or better, so the density is faithful
and what diverges is an optimizer path on a flat ridge, with the `rdm`
log-likelihoods differing by 2.7e-9 while a coefficient differs by
0.068. That is what NEWS now says. The response-scale `ndt` row of the
table above is from the first round's design and no longer moves at
all: the scalar bound is back in the link, so
`predict(dpar = "ndt", type = "response")` returns 0.2307934893 on both
builds.

## Two defects this change introduced, found by the suite and fixed

Both were silent-wrong-answer shaped, which is why they are recorded
rather than merely fixed.

**A bare `wiener()` inside `mixture()` would have read the FRACTION as
a time.** `mixture()` does not finalize its components: it builds one
family out of theirs and never calls their `family_finalize`. Through
0.6.0 that case failed loudly, because the bound lived in the link and
the link refused until it had one. With the bound in the density,
nothing refused: `dpars$ndt` arrived in `(0, 1)` and was used as
seconds. The family now carries a bound from CONSTRUCTION, either the
one `max_ndt` names or a refusing placeholder
(`ddm_ndt_preinstall()`), so the old loud failure is back and says what
to do. `ndt_group()` inside a mixture is refused for the same reason:
no per-group table would ever be built.

**A fixed starting fraction of 0.5 is not the old starting value.**
`0.5 * min(y)` under a bound of `min(y)` is a fraction of one half, but
under `max_ndt` ABOVE the fastest response, which is exactly the
mixture case, half the bound starts the optimizer at a non-decision
time no observed row can reach. Measured on `test-defects.R`'s mixture,
seed 3, `dev/ndt-scripts/ndt-mixture-debug.R`:

| | reference 0.6.0 | this change, start 0.5 | this change, start fixed |
|---|---|---|---|
| starting `ndt` | 0.0907 | 0.2000 | 0.0907 |
| log-likelihood | -82.776 | -79.424 | -82.776 |
| fitted `ndt` | 0.3906 | 2.68e-07 | 0.3906 |
| fitted drift | 0.248 | 51.7 | 0.248 |
| fitted boundary | 2.553 | 50.8 | 2.553 |
| finite residuals | yes | NO | yes |

The wrong start found a HIGHER log-likelihood at a degenerate point,
which is the failure mode a starting value cannot be checked by
convergence alone. `ddm_ndt_install()` now sets the start against
`min(y)` and the bound together, so it is the old starting TIME
whatever the bound is.

## What would have changed my mind, measured

### Is the -7027.4 target reachable under a per-subject bound?

Yes, and it is exceeded by 24.4 units. The premise was checkable before
any fitting: `dev/ndt-scripts/ndt-floors.R`, seed 20260908, at the
tier's own 30 x 400 design.

| quantity | value |
|---|---|
| global `min(rt)` | 0.226234796 |
| true `ndt` range | 0.1996 to 0.2947 |
| subjects whose true `ndt` is above the GLOBAL bound | 20 of 30 |
| subjects whose true `ndt` is above their OWN bound | **0 of 30** |
| smallest margin, own floor less true `ndt` | **26.61 ms** |
| true `ndt` as a fraction of its own floor | 0.747 to 0.888, median 0.849 |

The truth is inside the per-subject-bounded parameter space with
26.6 ms to spare at its tightest, and the fit lands there with 27.3 ms.

### Is a per-subject floor a biased ceiling at 400 trials?

**Yes, by 47 ms on average, and the fit does not treat it as the
truth.** `dev/ndt-scripts/ndt-floors.R`, 20 seeds (20260908 + 0:19)
times 30 subjects, 600 subject-draws:

| quantity | value |
|---|---|
| floor less true `ndt`, ms | min 8.61, q25 31.2, median 42.5, q75 58.2, max 141.8 |
| mean overshoot | 46.05 ms, sd 20.52 |
| true `ndt` / own floor | mean 0.8478, sd 0.0570 |
| subjects above the GLOBAL floor | 443 of 600 |

and it shrinks with trials, seed 20260908:

| trials per subject | mean overshoot | sd | max |
|---|---|---|---|
| 50 | 77.6 ms | 26.3 | 162.1 |
| 100 | 63.7 ms | 25.1 | 133.7 |
| 200 | 59.0 ms | 18.2 | 105.9 |
| 400 | 47.1 ms | 13.8 | 81.3 |
| 800 | 43.3 ms | 12.5 | 76.9 |

The floor is an upper bound, not an estimate, and the fitted fraction
absorbs the overshoot: at the eam optimum the population fraction is
0.8497 against a true median ratio of 0.849, and the recovered
population `ndt` is 0.2469 against 0.25.

### The argument for the change is what happens as data accumulates

This is the strongest evidence for the per-group bound and the review
produced it, not this lane. `dev/rev-ndt-floors.R`, seed 20260908, the
tier's own design and truths, 30 subjects at 100, 200 and 400 trials,
three estimators on the SAME data: `re` is what ships,
`ndt ~ 1 + (1 | s)` with `ndt_group(s)`; `nore` drops the random effect
and keeps the grouping; `global` is 0.6.0's single bound.

| trials | arm | logLik | `sd(ndt)` | per-subject RMSE | corr with truth |
|---|---|---|---|---|---|
| 100 | `re` | -1887.866 | 0.03198 | 20.91 ms | 0.7770 |
| 100 | `nore` | -1887.866 | 0.03198 | 20.91 ms | 0.7770 |
| 100 | `global` | -1926.759 | 0.00659 | 27.23 ms | 0.587 |
| 200 | `re` | -3704.555 | 0.02882 | 14.08 ms | 0.8970 |
| 200 | `nore` | -3704.934 | 0.02957 | 15.68 ms | 0.8716 |
| 200 | `global` | -3755.676 | 0.01241 | 20.51 ms | 0.721 |
| 400 | `re` | -7003.008 | 0.02543 | 7.67 ms | 0.9607 |
| 400 | `nore` | -7011.446 | 0.02748 | 11.90 ms | 0.9092 |
| 400 | `global` | -7149.635 | 0.00578 | 30.37 ms | 0.533 |

Truth: `sd(ndt)` 0.02629.

**The per-group arm's per-subject error halves as the data grows,
20.91 to 14.08 to 7.67 ms, and the global arm's does not, 27.23 to
20.51 to 30.37 ms.** The global bound does not converge on the truth
with more data, because more data lowers the global minimum and
tightens the ceiling on every subject at once. Every other column moves
the same way: the log-likelihood gap widens from 38.9 to 146.6 units
and the correlation with the truth goes 0.59, 0.72, 0.53 on the global
arm against 0.78, 0.90, 0.96 on the per-group one.

This is now the first table in the NEWS bullet and a `?wiener` section,
because it is the argument a user needs and the log-likelihood row
alone is not it.

Two things it also settles that this lane had left open.

**The floor's overshoot does not propagate into the estimate.** The
population bias is +6.59, +6.40 and +2.69 ms at 100, 200 and 400
trials, that is 1 to 2.6 percent of a 250 ms truth, and it SHRINKS with
trials, while the floor's own overshoot as a ceiling is 63.7, 59.0 and
47.1 ms. The fitted fraction absorbs it.

**And `nore` is the estimator the plan's criterion has to beat.** At
100 trials the two arms are IDENTICAL in every digit printed: the
variance component collapses and the shipped model IS the floors-only
model. At 400 trials the component is worth 8.44 log-likelihood units
and cuts the per-subject error by 36 percent, which is a real
contribution and a larger one than the 0.2 ms of 8.0 ms this lane's own
15 x 250 study found. But `sd(ndt)`, the statistic the plan asks for,
comes back at 0.02748 from the floors-only estimator against 0.02543
from the full model and a truth of 0.02629, so both "recover" it.

### But the recovered spread is the FLOORS, not the random effect

This the plan does not mention and a reader of the table above would
get wrong. `ndt_i = fraction_i * floor_i`, so the fitted per-subject
non-decision times vary between subjects even when the random effect
is exactly zero, purely because the floors do.

`dev/ndt-scripts/ndt-floor-artifact.R`, 15 subjects x 250 trials, 8
seeds (771 + 0:7), two designs. `sd_floor_only` is the spread a COMMON
fraction times the floors would give, that is the spread with the
random effect switched off.

**Design A**, `ndt` varies (log spread 0.12) and the boundary does not,
which is the scale design's shape:

| quantity | mean over 8 seeds | sd |
|---|---|---|
| true `sd(ndt)` | 0.02745 | 0.00421 |
| fitted `sd(ndt)` | 0.02534 | 0.00541 |
| the floors alone would give | 0.02526 | 0.00553 |
| `sd(floor)` | 0.02981 | 0.00638 |
| correlation, fitted with true | 0.9565 | 0.0195 |
| correlation, FLOOR with true | 0.9514 | 0.0300 |
| bias, ms | +0.46 | 2.05 |
| RMSE per subject, ms | 8.03 | 1.09 |
| RMSE with the random effect off, ms | 8.25 | 1.09 |

The fitted spread and the floors-alone spread agree to three decimals,
and the random effect buys 0.2 ms of the 8.0 ms per-subject RMSE. The
recovery in the headline table is real, and most of it is the floor
carrying the information rather than the random effect estimating it.
That is not a defect: the floor IS information about that subject's
non-decision time, which is the whole argument for the change. It does
mean "`sd(ndt)` is recovered" should not be read as "the `ndt` variance
component is recovered".

This lane's 15 x 250 study understated the component's contribution.
The table above, at the tier's own 30 x 400, is the number to use: 8.44
log-likelihood units and 11.90 ms to 7.67 ms, not 0.2 ms of 8.0.

**Design B**, the falsification: `ndt` CONSTANT at 0.25 for every
subject and the log boundary spread 0.25 instead, so the floors vary
and the truth does not.

| quantity | mean over 8 seeds | sd |
|---|---|---|
| true `sd(ndt)` | 0 | 0 |
| fitted `sd(ndt)` | 0.00640 | 0.00094 |
| the floors alone would give | 0.00640 | 0.00094 |
| `sd(floor)` | 0.00796 | 0.00133 |
| bias, ms | -23.58 | 6.06 |
| RMSE per subject, ms | 24.44 | 5.85 |

The random effect collapses to zero, correctly, and the fit STILL
reports a between-subject spread of 6.4 ms because the parameterization
cannot report zero spread when the floors differ. The population
non-decision time also comes back 23.6 ms low, 9.4 percent.

**The control that settles whether the bound caused it.** A limitation
only counts against the change if the thing it replaced did not have
it, so design B was run again with BOTH parameterizations on the same
data, `dev/ndt-scripts/ndt-artifact-control.R`, same 8 seeds:

| quantity, mean over 8 seeds | per-group bound | global bound |
|---|---|---|
| population `ndt` bias, ms | -23.58 (sd 6.06) | -23.87 (sd 6.05) |
| reported `sd(ndt)`, ms, truth 0 | **6.40** (sd 0.94) | **0.00** |
| per-subject RMSE, ms | 24.44 (sd 5.85) | 23.87 (sd 6.05) |
| log-likelihood | -2060.08 | -2060.85 |
| `sd(floor)`, ms | 7.96 (sd 1.33) | the same data |

The 23.6 ms bias is NOT the per-group bound's: the global bound has the
same one, -23.9 ms, so it belongs to the design (a wide boundary
spread with a superfluous `ndt` random effect) and not to the
parameterization. The two differ by 0.29 ms.

What IS the per-group bound's is the phantom spread: 6.4 ms where the
global bound reports exactly zero, because `ndt_i = fraction_i *
floor_i` cannot say "the same for everybody" when the floors differ.
The price of that on per-subject accuracy is 0.57 ms of RMSE out of 24,
and the per-group fit is nonetheless the better model on 6 of the 8
seeds (mean +0.78 log-likelihood units, per-seed differences 2.14,
0.78, -2.57, 1.96, 2.29, 0.95, -1.31, 1.97).

**So the trade is: a 6.4 ms spread reported where there is none, and
0.6 ms of per-subject RMSE, on the design built to hurt the change,
against the whole of the eam row.** A user who believes their subjects
share one non-decision time should leave `ndt_group()` off, which is
the default and is unchanged from 0.6.0.


### Did an existing test depend on the global bound being right?

No test asserted that the global bound was the right constraint. Three
asserted the MECHANISM and had to move, and five more read `ndt` as a
time:

- `test-family.R`'s "the scaled-logit ndt link round-trips and keeps
  the support" asserted `linkinv(eta) < max_ndt` and that a family
  built before `frm()` refuses its link with "bound is not set yet".
  Both are statements about where the bound lived. The link now
  round-trips on the fraction, and a family built before `frm()` has a
  usable link.
- `test-brms-parity.R` asserted `links$ndt$name == "scaled_logit"` and
  `linkinv(40) == max_ndt`. Rewritten, with the parity claim
  sharpened: `ndt`'s linear predictor is no longer comparable with
  brms's at all, where before it was a different transform of the same
  quantity.
- `test-defects.R` asserted the same "bound is not set yet" refusal.
- `test-variability.R` asserted `st`'s scaled logit onto
  `(0, 2 * max_ndt)`.
- `test-lba.R`, `test-rdm-gng.R`, `test-simulate-density.R` and
  `test-extension-api.R` read `family(fit)$links$ndt$linkinv(eta)` as a
  time or handed a density parameters in times. Those now go through
  `ndt_time()` or through the new `helper-ndt.R`.

Everything else that reads `ndt` reads it as
`min(rt) / (1 + exp(-eta))`, which is unchanged arithmetic, so
`test-defects.R`'s two mixture tests and `test-family.R`'s two
recovery tests pass without an edit.

### One thing that DID move, and it is not the bound

The `eam-sv` scale row's condition effect no longer covers its truth.
That row used to assert coverage and now does not, so this is the
measurement that replaced the assertion. `dev/ndt-scripts/ndt-sv-arm.R`,
the tier's own seed 20260908, the same 12,000 rows fitted both ways:

| arm | log-likelihood | `mu.condb`, truth 0.9 | covers | `sv`, truth 0.4 | `sd(ndt)`, truth 0.0263 | population `ndt` bias |
|---|---|---|---|---|---|---|
| `sv`, per-group bound | -6926.65 | 0.8223 (0.7621, 0.8825) | no | 0.3105 | 0.0282 | +3.8 ms |
| `sv`, global bound | -7129.88 | 0.8565 (0.7936, 0.9195) | yes | 0.5589 | 0.0077 | -22.7 ms |
| plain, per-group bound | -6928.70 | 0.7982 (0.7449, 0.8514) | no | | | |

BOTH parameterizations underestimate the condition effect on this
draw, by 0.078 and 0.044, and the global arm covers only because its
interval is 8 percent wider and its point estimate 0.034 higher.
Neither recovers `sv`, which comes back at 0.31 and 0.56 against 0.4.
What moved is the drift-against-variability trade-off: with `ndt`
pinned 22.7 ms low the between-trial drift spread is absorbed by `sv`,
and with `ndt` right it is not. The non-decision time itself recovers
under the bound and does not under the global one, which is what this
lane changed.

One seed, one arm. The plan puts recovery at Phase 2's 60 replicates
and this is exactly the row that needs them. The tier now asserts
convergence, a finite interval and the `ndt` recovery, and RECORDS the
condition effect with the numbers above in the file.

## What was NOT done, and why

### `gddm()` does not get the per-subject bound

Not a scope decision. It could not use one. `gd_densities()`
(`R/gddm.R:339-341`) solves the Fokker-Planck equation once per
CONDITION and reads every distributional parameter at the FIRST ROW of
that condition:

    i <- d$first[j]
    pj <- lapply(ctl$dpars, function(nm) gd_at(dpars[[nm]], i))

so a per-row `ndt` is already collapsed to one value per condition
before any bound applies, and `gd_shift()` takes a scalar shift. A
per-row bound would change nothing that reaches the density.

`gddm()` therefore keeps its scaled logit, its `ndt` stays a TIME, and
`ndt_group()` is refused BY NAME through `accepts_aterms`, with a
`frm_compat("gddm")` row that says why. `ndt_time()` refuses a `gddm()`
fit rather than returning a number it has no bound for.

**A defect found and not fixed.** A `gddm()` model whose `ndt` (or any
other dpar) varies WITHIN a condition is silently fitted with the first
row's value for the whole condition. That is the same class of failure
this lane was sent to remove, in another family, and it is not about
the bound. It needs either a refusal at frame assembly when a dpar's
linear predictor is not constant within a condition, or a solve per
distinct parameter vector. It is not fixed here because it is a
different mechanism from item 1.0a and changing the solver is more than
this lane should do unreviewed.

### `frmtmb.learn`'s `rlddm()`

Item 1.0b, another lane. Nothing here blocks it: `ndt_group` is
registered by frmtmb.eam's `.onLoad()`, so a package depending on
frmtmb.eam has the term already, and `ddm_ndt_spec()`,
`ddm_ndt_scaler()`, `ddm_ndt_install()` and `ddm_ndt_preinstall()` are
the four functions it needs. They are `@noRd`, so reaching them means
either a second export seam in `R/extension-api.R` or a copy. The
export seam is the better answer and this lane did not make it, because
its shape depends on how `rlddm()` carries its own per-trial state.

### `utils` is used without being declared

`R/lba.R:453` and `R/rdm.R:401` call `utils::head()` and `DESCRIPTION`
lists only `RTMB` and `stats` under `Imports`. Pre-existing, and
`R CMD check --as-cran` does NOT flag it here, so it is a latent thing
rather than a live one; new code in this change avoids `utils::`
anyway. Worth a one-line `DESCRIPTION` fix by whoever owns the next
release.

## Which version bump this needs

Not chosen here. What the bump has to cover: a shipped distributional
parameter's meaning changed, a post-fit call that returned seconds now
returns a fraction, a new exported function, a new addition term, and
one family deliberately left on the old scheme. The NEWS bullets are
under `# frmtmb.eam (development version)`.

## What the plan should say instead

`dev/extension-gaps-plan.md` is not edited here. Item 1.0a's row should
read:

- the deliverable landed and beat its own target: -7003.01 against a
  target of -7027.4 and a baseline of -7148.8, with a positive
  definite Hessian and `sd(ndt)` at 0.0254 against 0.0263;
- the acceptance criterion "recovers `sd(ndt)`" is met in appearance:
  an estimator with the random effect on `ndt` switched off returns
  0.02748 against a truth of 0.02629 where the shipped model returns
  0.02543, and at 100 trials per subject the two are identical in every
  digit. The statistics that DO separate the model from its floors are
  the per-subject error and the log-likelihood, and those are what
  Phase 2 should assert;
- the row should carry the trial-count table: the per-group arm's
  per-subject error halves with more data, 20.91 to 14.08 to 7.67 ms,
  and the global bound's does not, 27.23 to 20.51 to 30.37 ms. That is
  the argument for the item and it belongs in the plan;
- item 2.1 is unblocked and should carry a second question: at design
  B above the parameterization reports a 6.4 ms spread that is not
  there and a population `ndt` 23.6 ms low. Sixty replicates would say
  whether that matters at the designs the field produces;
- the `eam-sv` row's condition effect is no longer covered at the
  tier's seed, under EITHER bound. That is a Phase 2 question and the
  tier records it rather than asserting it;
- item 1.0b can take `ndt_group()` as it stands, but needs an export
  seam from frmtmb.eam for the four internal functions;
- a new item, small: `gddm()` reads every dpar at the first row of its
  condition, so a dpar that varies within a condition is silently
  ignored.

## Tests and the check

### Punch round 1

The files this round's changes reach, one per process, `NOT_CRAN=true`,
failure cap lifted, `package = "frmtmb.eam"`
(`dev/ndt-scripts/ndt-reached.ps1`, log `ndt-reached-log.txt`). The
BASELINE column is `dev/suite-baseline.tsv`, joined per file, which is
what that file exists for.

| file | baseline | this branch | |
|---|---|---|---|
| test-ndt-bound.R | (new) | 79 | |
| test-family.R | 51 | 51 | back to baseline |
| test-defects.R | 59 | 59 | back to baseline |
| test-brms-parity.R | 13 | 13 | back to baseline |
| test-variability.R | 140 | 140 | back to baseline |
| test-lba.R | 109 | 109 | back to baseline |
| test-rdm-gng.R | 221 | 221 | back to baseline |
| test-simulate-density.R | 69 | 69 | unchanged |
| test-extension-api.R | 14 | 14 | back to baseline |
| test-gddm-family.R | 77 | 77 | unchanged |
| test-bracket-access.R | 1 | 1 | unchanged |
| test-message-uniqueness.R | 4 | 4 | unchanged |
| test-surface.R | 50 | 50 | unchanged |
| test-sampling.R | 97 | 94, SKIP 1 | the private library, not a regression; the review reran it at 97 |

`FAIL 0` on every one. The six rows marked "back to baseline" are the
first round's edits reverted, and their being back at the baseline
COUNT as well as at the baseline content is the check that the scalar
path really is the pre-change one.

The eam scale row, rerun on the reworked build under
`FRMTMB_SCALE_TESTS=true`: `FAIL 0 | WARN 0 | SKIP 2 | PASS 7`, one
more than before because the population `ndt` is now asserted as a
`scale_z()`. It reproduces the headline to every digit for the THIRD
time, on a third build: log-likelihood -7003.01, `ndt` 0.246867 with se
0.0019639, `sd(ndt)` 0.0254253, tightest margin 27.3347 ms, 30 of 30
subjects below their own floor, condition effect 0.910877 (0.858162,
0.963592), positive definite Hessian, maximum gradient 9.87e-04.
Record at `dev/ndt-scripts/ndt-scale-r1.tsv`.

The whole 22-file suite was NOT rerun this round, per the round's cost
rule: eight of the nine files the first round edited are byte-identical
to 92e9330 again, and the review holds current numbers for the rest.

`R CMD check --as-cran` on the built tarball, with pandoc and TinyTeX
on PATH and `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`: **Status: 1 NOTE**,
and the NOTE is the expected one, "Skipping checking math rendering:
package 'V8' unavailable". `checking tests ... [18m] OK`, which is the
whole suite in one process. Log at
`dev/ndt-scripts/ndt-check-log.txt`.

### The first round

The whole suite, one file per process (`dev/ndt-scripts/ndt-suite.ps1`,
log at `ndt-suite-log2.txt`): **22 files, FAIL 0, 1454 passing
assertions**, every process exit code 0. That run predates the punch
round and its per-file counts for the reverted files are the fraction
era's, not this branch's; the table above supersedes it.

The three scale rows, one per fresh process: `eam` FAIL 0 PASS 6,
`eam-unbounded` FAIL 0 PASS 1 (WARN 1, its pre-existing gradient
warning), `eam-sv` FAIL 0 PASS 5 (WARN 1, the same).


## The scripts, and what each one is for

Every one takes its seed from its own header and writes its log beside
it in `dev/ndt-scripts/`.

| script | what it establishes | log |
|---|---|---|
| `ndt-floors.R` | the per-subject floors against the truths, and the floor's bias as a ceiling at 50 to 800 trials | (stdout, quoted above) |
| `ndt-smoke.R` | the bit-identical control: no `ndt_group()` is the 0.6.0 parameterization | (stdout, quoted above) |
| `ndt-before-after.R` | the behavioral failure, both builds, two groups 250 trials | (stdout, quoted above) |
| `ndt-group-check.R` | the per-group plumbing: floors, refusals, newdata, `simulate()` | (stdout) |
| `ndt-mixture-debug.R` | the starting-value defect this change introduced, both builds | (stdout, quoted above) |
| `ndt-sv-arm.R` | the `eam-sv` condition effect under both bounds | `ndt-sv-arm-log.txt` |
| `ndt-floor-artifact.R` | designs A and B: how much of the recovered spread is the floors | `ndt-artifact-log.txt` |
| `ndt-artifact-control.R` | design B under BOTH bounds, which is what exonerates the change | `ndt-artifact-control-log.txt` |
| `ndt-round1.R` | the punch round: the scalar path against 0.6.0 (B2, B3), the level orderings and the absent label (B1), the kept bound (N0), and N2, N4, N5 | `ndt-round1-log.txt` |
| `ndt-scale-run.ps1` | the three eam scale rows, one per fresh process | `ndt-scale-log2.txt`, `ndt-scale.tsv` |
| `ndt-suite.ps1` | the whole suite, one file per process | `ndt-suite-log2.txt` |
