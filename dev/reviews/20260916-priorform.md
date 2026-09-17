# Review: lane wt-priorform (item 2.6c), round 1

Adversarial pass on the worktree at base aa9227e. Builds: BASE
`C:/Users/adf44/source/r/rellib-r3` (read only), LANE
`C:/Users/adf44/source/r/priorform-lib` (read only; not rebuilt), brms
2.23.0 from the user library. Every script below is under `dev/`, named
`priorform-rev-*`, and writes its log or table beside itself.

Before measuring, the lane library was checked against the worktree
source function by function (`priorform-rev-libmatch.R`): 973 top-level
functions in `R/*.R` deparse identically to the installed namespace,
0 differ, 0 missing.

The user decided two policies during this review, and they change how
two items are read:

1. The duplicated group-level refusal stays exactly as brms has it,
   including the animal model `(1 | gr(id, cov = A)) + (1 | id)`. That
   design is not ranked here. Cases where brms ACCEPTS and the lane
   refuses are still reported (F3).
2. Duplicate prior specifications on one slot will be refused as brms
   refuses them. Section 3 lists the blast radius for the punch round.

## Verdict

Mergeable after punch. No BLOCKER. Two MAJOR findings: F1 is a new
public function that gives a wrong answer, and F2 is a silent wrong
model that the lane's own fix pattern already covers.

## 1. Findings

### F1 MAJOR: `validate_prior()` shows a prior that `frm()` does not apply, and its table refits a different model

`fill_prior_table()` writes the rows in specification order, then fills
the rows below a class row only where the row is still at its default.
`resolve_priorlist()` applies later-wins per parameter. The two agree
only when the specifications go from less to more specific. When a
coefficient-specific or group-specific specification comes FIRST, the
table shows the earlier one on that row, and `frm()` applies the later
class-wide one.

Construction: `priorform-rev-roundtrip.R`, seed 20260916, lane build,
17 cases. `frm(prior = pl)` against `frm(prior = validate_prior(pl, ...))`:

| case | fit identical | resolved entries identical | logLik difference |
|---|---|---|---|
| `set_prior("normal(3, 0.01)", coef = "x") + set_prior("normal(0, 0.01)", class = "b")` | FALSE | FALSE | 123 |
| `set_prior("normal(0, 5)", class = "sd", group = "g") + set_prior("normal(0, 0.05)", class = "sd")`, beta with `(1 | g)` in mu and phi | FALSE | FALSE | 1.88 |
| the same two pairs in the other order | TRUE | TRUE | 0 |
| 13 other cases (distributional, dpar-only sd, bounds then density, density then bounds, mv, nonlinear with bounds, smooth, gp theta, cumulative, sratio with `cs()`, lkj, sigma dpar) | TRUE for 12 | TRUE for 13 | 0 |

For the first case the table prints `normal(3, 0.01)` with source
`user` on row `b x`, while `resolve_prior_input()` on the original
prior gives `beta2 = normal(0, 0.01)`. A user who reads the table is
told the wrong prior, and a user who passes the table back gets a fit
123 log-likelihood units away.

The lane's test "a validate_prior() table fits the model its prior
fits" uses only the class-then-coefficient order, so it cannot see
this.

Note for the punch round: brms applies the more specific specification
whatever the order. The duplicate-prior decision refuses same-slot
pairs only, so `coef = "x"` then class `b` stays legal and this
disagreement survives it. Either make the table follow the resolver
(later wins), or make the resolver follow brms (the more specific
wins). The second one changes `frm()` fits, and the user must decide it.

The remaining non-identical case is `mv`: resolved entries are
identical and the logLik differs by 2.3e-13. The table's row order
changes the order of the prior terms that are added together. It is a
NIT, but it means `identical()` is too strict for a round-trip test on
a multivariate model.

### F2 MAJOR (pre-existing, not introduced): `x:s(z)`, `x * s(z)` and `I(s(z))` fit a different model without a message

brms refuses all three with "The term ... is invalid", which is the
same rule the lane adopted for `cs()`. Base and lane build the same
wrong model. Construction: `priorform-rev-smoothint.R`, seed 20260916,
n = 200, `y = 1 + x sin(z) + e`:

| formula | brms | base and lane (identical) |
|---|---|---|
| `y ~ x:s(z)` | refuse | X = (Intercept), x; 0 smooths; logLik -94.309 |
| `y ~ x * s(z)` | refuse | the same model |
| `y ~ I(s(z))` | refuse | X = (Intercept) only; logLik -212.868 |
| `y ~ x:gp(z)`, `x * gp(z)` | refuse | refused, "could not find function gp" (a poor message, not a wrong model) |
| `y ~ s(z, by = x)` | accept | smooth fitted, logLik -37.577 |

The smooth and the interaction are dropped without a message. This is
the same defect class as the `x * cs(g)` drop that the lane fixed.
`cs_term_refusal()` extends to `s`, `t2`, `te`, `ti` and `gp` with a
change to the name test only. The lane did not introduce this defect,
so it does not block the merge. It is ranked MAJOR because it gives a
wrong model without a message.

### F3 MINOR: the duplicate refusal fires where brms accepts

`priorform-rev-dup.R brms|ref|lane`, seed 20260916, 54 designs.
Output: `priorform-rev-dup-{brms,ref,lane}.tsv`. The worker's generated
set combined DISTINCT pool terms. For that reason, it had no exact
twins and no reversed interaction, and its "0 false alarms" does not
cover these designs:

| design | brms | base | lane |
|---|---|---|---|
| `(1 | g) + (1 | g)` in mu (`mv_dup_one`), in an nlpar (`nl_dup_in_a`) | accept | accept | refuse |
| `(1 | gr(id, cov = A)) + (1 | gr(id, cov = A))` | accept | accept | refuse |
| `(1 | mm(g1, g2)) + (1 | mm(g1, g2))` | accept | accept | refuse |
| `(1 | g:h) + (1 | h:g)` | accept | accept | refuse |

What each package builds (`priorform-rev-twins.R`):

- brms collapses an exact twin to ONE block. `standata` has `M_1`
  only for `(1 | g) + (1 | g)` and for `(x | g) + (x | g)`. The
  cause is that `terms()` removes duplicate labels.
- Base frmtmb fitted TWO blocks for each twin, and the logLik equals
  the single-term logLik (-184.773190 both times). The split was not
  identified.
- brms keeps `(1 | g:h) + (1 | h:g)` as two blocks (`M_1 M_2`), the
  same unidentified split. The lane's sort-the-factors key refuses it.

The lane refusal is kinder than either fit. However, the policy is
"exactly as brms has it". So either collapse exact twins, as brms
does, and compare interaction groups as strings, or record both as
divergences. The worker's findings must also stop calling the rate 0
for these shapes.

Also recorded, not ranked, on frmtmb's own glmmTMB structures, where
brms has no verdict: the lane refuses `ar1 + diag`, `rr + diag`,
`us + diag`, `cs + diag`, `homcs + diag`, `equalto + (0 + f | g)`,
`equalto + diag`, `exp + diag` and `gr(prec =) + (1 | id)` on the same
factor and coefficients. Base accepts all of them. Identification on
base (`priorform-rev-ident.R`, 150 groups, 4 seeds, `pdHess` and max
theta SE):

- `gr(cov) + plain`: pdHess in 4 of 4, max SE 0.29 to 0.32. Identified.
- `equalto + us`: pdHess in 4 of 4, max SE 27 to 161. Weak.
- `ar1 + diag`: 2 of 4. `rr(d=1) + diag`: 3 of 4, SE above 2e4.
  `exp + diag`: 1 of 4. `us + diag` (a control that is not
  identified): 1 of 4.

Only the known-matrix design is clearly identified at this size, and
brms's rule covers that design. No vignette, Rd example or test in core
or the extensions writes any of these designs (grep of vignettes; every
Rd example was run, section 2).

### F4 MINOR: the `mvbind(...)` `nl` pass-through fix has no guard

The lane fixed `bf(mvbind(y1, y2) ~ a, a ~ 1, nl = TRUE)`, which built
linear per-response formulas on base (`priorform-rev-bfedge.R`: base
prints `y1 ~ a`, lane prints `y1 ~ a (nonlinear)`). Putting the old
call back in the lane's `bf()` (`priorform-rev-revert.R bfnl`) leaves
all 795 assertions of the six affected files passing. Add an assertion
on `$responses[[1]]$nl` or on the fitted parameter set.

### F5 MINOR: ledger P15 and P27 hide a sampling-route gap

The ledger records both as divergences because "the fit route is
flat", and gives brms's defaults as "frmtmb.sample's". On
`route = "sample"` with frmtmb.sample attached
(`priorform-rev-route.R`):

- P15 rescor: `(flat)`, where brms gives `lkj(1)`.
- P27 intercept with `offset(off)`, where off = 10: `student_t(3, 2, 2.5)`,
  where brms gives `student_t(3, -8, 2.5)`. The sampling default does
  not subtract the offset.

Neither is this lane's code. The reason text is still wrong, and two
default-prior defects in frmtmb.sample should be filed. The other 19
"cannot transfer" rows and 6 divergences were read against brms's test
source and the lane's code, and each reason held.

### F6 MINOR: a guard lost its purpose

`test-review-v29.R`, "blocks sharing a term label are all reported, not
the first twice", guarded two blocks whose VarCorr labels are the same.
With `(1 | id2)` the labels differ, so the test no longer covers the
repeated-label code paths that `R/methods-fit.R` still carries
(comments near lines 841 and 957). A shared label can still be built
across predictors, for example `(1 | g)` in `mu` and in `sigma`. Build
the guard on that design so that it covers what it was written to
cover.

### F7 NIT

- The refusal message loses `gr()`: `(1 | gr(g, dist = "student")) +
  (1 | g)` reads "appears in both (1 | g) and (1 | g)". The label is
  taken after `gr()` is stripped.
- `cs(0 + time || g)` goes to the category-specific branch and fails
  with "cs() needs an sratio, cratio, or acat family". The check tests
  `"|"` and not `"||"`. Base does the same.
- `print(empty_prior())` prints an empty table header. brms prints
  nothing. This is documented.
- `as.brmsprior(data.frame(prior = c("normal(0,1)", ""), lb = c(NA, "0")))`:
  brms marks the bounds-only row `(vectorized)`, and the lane marks it
  `user` with `(flat)`.
- `bf(y ~ x, phi = sigma ~ z)` silently uses `sigma`. brms does the
  same, so this is no divergence. It is only noted.
- Instrument note: the gp fit in `priorform-rev-bitwise.R` A differed
  at the 1e-11 level in the first run, when two R processes ran at the
  same time. The run on base repeated identically (`ref2`), the lane
  rerun equals base (`lane2`), and three fits in one process are
  identical. The cause is load, and it is a gp reproducibility note
  outside this lane.

## 2. Claims that held, and how they were attacked

- **No silent change in fits.** `priorform-rev-bitwise.R` part A fits
  28 formulas through `bf()`, `lf()`, `nlf()`, `mvbf()`, `mvbind`, `nl`
  flags (explicit FALSE, shared `a + b ~`), sratio and acat with
  `cs()`, cumulative with RE, beta with RE in mu and phi, and with and
  without sd priors. It also fits zi, smooth, gp, `cs()` and `ar1()`
  covariances, `||`, `gr(cov)` on a copy, `mm()`, theta and
  natural-sigma priors. On a quiet run, 28 of 28 are `identical()`
  between base and lane in logLik, optimizer vector, objective,
  sdreport `cov.fixed` and estimates (compare file
  `priorform-rev-bitwise-A-compare.txt`, with the reruns above). Part B
  runs every Rd example of core, coupling, eam, latent, learn, ode and
  spline with `set.seed(1)`: 114 fits, 114 identical. The only example
  errors are the 6 on base that call functions the lane adds.
- **`cs()` covariance structure is not touched by the new `cs()` term
  logic.** `cs(0 + time | g)` alone, in an ordinal model, twice on two
  factors, beside `cs(x)`, and `cs(x) * cs(0 + time | g)`: the verdict
  is the same on base and lane, and `cs_cov` fits bitwise identically.
- **Refusals match brms on the designs brms refuses.** 10 of 10 brms
  refusals in `priorform-rev-dup.R` are refused by the lane, including
  `gr(dist = "student")`, `||` plus an intercept, a factor slope plus
  an intercept, `|ID|` keys and a sigma predictor. Same group in two
  dpars, two responses, two nlpars, `|ID|` across predictors, `mm()`
  beside a member, `car()`, `s(g, bs = "re")` and `gp()` beside
  `(1 | g)`: all accepted by all three. `gr(by = )` is refused by
  frmtmb on base for another reason (not supported). No miss.
- **The animal model result.** Measured again independently with
  `stancode()`, not `default_prior()`: brms refuses
  `(1 | gr(id, cov = A)) + (1 | id)` and `(x | gr(id, cov = A)) + (0 + x | id)`.
- **Resolution did not change when sd rows gained `dpar`.** The table
  change does not touch `resolve_priorlist()`. `beta_phi_re_prior`
  (class sd over the mu and phi blocks) and the group-level variant fit
  bitwise identically on base and lane. The sampling route defaults
  read the class-wide key, as before.
- **`set_prior()` recycling against brms.** `priorform-rev-print.R`,
  39 calls diffed against brms: lengths 2 and 4 recycle the same way.
  Lengths 1, 2, 3 and zero-length arguments are refused in both, with
  different wording. `$prior`, `$class`, `$lb`, `$ub`, `$source` and
  `$coef` agree. Every single-row print line agrees with brms character
  for character, except two rows frmtmb refuses by design (sd with
  `coef`, and dpar with nlpar). Tables differ only by brms's `tag`
  column.
- **Guards fail closed** (`priorform-rev-revert.R`, a base function
  put back in a lane process, six test files, control 795 pass 0
  fail). Reverting the duplicate refusal fails 7 assertions, the `cs()`
  parse fails 4, the nested refusal 6, `bf()` 4 plus 3 errors, `lf()`
  2, `as_priorlist()` 1, `set_prior()` 8, and `block_dpar()` 1. Print:
  the first revert passed, but the instrument was open because the S3
  table was changed and the namespace binding was not.
  `priorform-rev-revert2.R` replaces the binding and fails 13
  assertions in 5 tests. These are behavioral failures, not "could not
  find function". The exception is F4.
- **The `bf()` edges.** `bf(bf(y ~ x), nl = TRUE)` skips the
  `nl = TRUE` check in `bf()`, but `frm()` refuses it later with a clear
  message. `bf(~ ~x)` is refused as brms refuses it. `y ~ x + ~z` and
  `sigma ~ (~z)` are accepted in all three, which is the top-level-only
  rule.
- **Stray warning.** The block tests `coef` and `group` index
  resolution, and does not test fitted values. The new generator still
  gives two slopes and two correlated blocks, which the block needs.
  `expect_no_warning()` is pinned at seed 77, where the correlations
  are -0.06 and -0.45 and far from the boundary, so the pin has little
  risk.

## 3. Blast radius of refusing duplicate prior specifications

brms's slot is (class, coef, group, resp, dpar, nlpar). Runtime
measurement: `priorform-rev-priortrace-run1.R` traces
`resolve_priorlist()` on the lane build and logs any priorlist with two
specifications on one slot. It ran on every test file that mentions a
prior: 31 in core, and 9 in the extensions (coupling, eam, learn,
spline 1 each; sample 5), one process each, all green
(`priorform-rev-priortrace-log.txt`). It also ran on all 162 Rd examples
(`priorform-rev-priortrace-rd.R`). The instrument was checked: it fires
on the known bounds-plus-bounds test. Hits
(`priorform-rev-priortrace.tsv`):

| where | slot | kinds | path |
|---|---|---|---|
| core `test-priors-bounds-grcov.R` 227-246, "a later bounds-only prior tightens an earlier one" | `b`, nlpar `guess` | bounds, bounds | `frm()` and `resolve_prior_input()` |
| frmtmb.sample `test-sample-direct.R` | `Intercept` | density, density (identical `student_t(3, -0.3, 3)`) | `frm_sample()` then `resolve_prior_input()` |
| frmtmb.sample `test-sampling-ported.R` | `Intercept` | density, density (identical `student_t(3, 2.4, 2.5)`) | `frm_sample()` then `resolve_prior_input()` |
| Rd examples, all packages | none | | |

The frmtmb.sample hits come from `frm_sample()`'s own stacking, not from
a user writing the pair twice (the source call was not pinned; a direct
`sample_resolve_priors()` on a plain MAP fit does NOT duplicate,
`priorform-rev-sampledup.R`). So the refusal must sit at the user
boundary (`frm(prior =)`, `frm_sample(prior =)`, `validate_prior()`),
not in `resolve_priorlist()`. A refusal there breaks those two sample
tests.

Code and documentation that state or depend on later-wins or on
bounds-only tightening (grep):

- `R/priors.R`: roxygen at 59-62 ("later specifications override
  earlier ones"), 108-113 ("a later bounds-only specification
  tightens"), 138-143 (the bounds section), 1462 (`validate_prior()`),
  1503 (`fill_prior_table()` order rule), 1885 (`entry_bounds()`, the
  tightening mechanism), 1947 and 1962 (`resolve_priorlist()` header
  and `claim()`). Rendered in `man/set_prior.Rd` 101, 159, 189-190 and
  `man/validate_prior.Rd` 57.
- `extensions/frmtmb.sample/R/sample.R`: `drop_superseded()` 1163-1180
  (keeps a bounds-only spec on purpose), `merge_prior_inputs()` 1231,
  `prior_stack()` 1277-1300, and the `sample_resolve_priors()`
  precedence text 1303-1320. These merge the fit's own prior, the
  call's prior and the defaults, and stay legitimate internal later-wins
  after a user-level refusal.
- Tests that assert same-slot later-wins: `test-priors-bounds-grcov.R`
  227-246 (must become a refusal test or a single spec with both
  bounds).
- Tests that use cross-slot pairs brms ACCEPTS, which must keep
  passing: `test-setprior.R` 40-46 (class `b`, then `coef = "x"`),
  `test-prior-compat.R` 885-888 (class `b` density, then `coef = "x"`
  bounds), `test-lkj.R` 292 (class `cor`, then class `theta`, which brms
  does not have), `test-brms-formula-priors.R` "a validate_prior() table
  fits the model its prior fits".
- No vignette writes a same-slot pair (grep of every `set_prior(` and
  `prior =` in core and extension vignettes; `frmtmb.Rmd` 391-394 uses
  four different classes).
- F1 interacts with this. A same-slot refusal does not make
  `fill_prior_table()` correct for coef-then-class order.

Not ranked, as briefed. Measured on brms 2.23.0 by the coordinator: brms
refuses same-slot identical rows, two densities, and a density followed
by bounds only. It accepts class `b` plus `coef = "x"`.

## 4. Scripts

| script | what | output |
|---|---|---|
| `priorform-rev-libmatch.R` | lane library against source | stdout |
| `priorform-rev-dup.R` | 54 designs, brms, base, lane | `priorform-rev-dup-*.tsv`, `-log.txt` |
| `priorform-rev-twins.R` | what brms and base build for twins | stdout |
| `priorform-rev-ident.R` | identification of frmtmb-only refused designs on base | `-log.txt` |
| `priorform-rev-bitwise.R`, `-bitwise2.R`, `-bitwise-compare.R` | bitwise fits, part A and Rd examples | `priorform-rev-bitwise-*.rds`, `-compare.txt` |
| `priorform-rev-gp.R`, `-gp2.R` | gp determinism control | stdout |
| `priorform-rev-roundtrip.R` | `validate_prior()` round trip, 17 cases | `-log.txt` |
| `priorform-rev-print.R` | print and recycling against brms | `priorform-rev-print-{brms,lane}.txt` |
| `priorform-rev-smoothint.R` | `x:s(z)` and related | `-log.txt` |
| `priorform-rev-bfedge.R` | `bf()` edge spellings | `-log.txt` |
| `priorform-rev-route.R` | ledger P15 and P27 on the sample route | `-log.txt` |
| `priorform-rev-revert-extract.R`, `-revert.R`, `-revert2.R` | fail-closed checks | `-revert-log.txt`, `-revert2-log.txt` |
| `priorform-rev-priortrace-run1.R`, `.sh`, `-rd.R` | duplicate-prior blast radius | `priorform-rev-priortrace.tsv`, `-log.txt` |
| `priorform-rev-sampledup.R` | frm_sample prior stack without sampling | `-log.txt` |

## Recheck, round 1

Scope: the round-1 changes (section 2P of `dev/priorform-findings.md`)
plus a regression check. Scripts are `dev/priorform-rev2-*`. The lane
library was checked against the worktree source again
(`priorform-rev-libmatch.R`): 978 functions identical, 0 differ, 0
missing.

### Verdict

Mergeable after a short punch. Nothing blocks. The new findings are
all MINOR. Every round-0 MAJOR is fixed, and the fixes are guarded.

### Remaining findings

**R1 MINOR: the duplicate-prior key refuses two different slots.**
`prior_slot_label()` compares the brms print name, which joins class,
group, resp, dpar, nlpar and coef with underscores. So
`set_prior(class = "b", dpar = "sigma")` and
`set_prior(class = "b", coef = "sigma")` both read `b_sigma` and are
refused. brms 2.23.0 accepts the pair (`priorform-rev2-slots.R`, first
row). A covariate named `sigma` is not hypothetical here:
`vignettes/inputs.Rmd` fits `bf(y ~ sigma + (1 | g))`. The same collision
can happen between a dpar or nlpar and any coefficient whose name starts
with that parameter and an underscore. Fix: compare the tuple (class,
bare coef, group, resp, dpar, nlpar), not the label. Keep the label for
the message.

**R2 MINOR: F1 is still visible in one `validate_prior()` row.** With
`set_prior("normal(0, 5)", class = "sd", dpar = "phi")` and a
class-wide `set_prior("normal(0, 1)", class = "sd")`, on
`bf(yb ~ x + (1 | g), phi ~ (1 | g))`, the fit gives the phi block
`normal(0, 5)` in both written orders. That matches brms. The table
still prints `normal(0, 1)` (vectorized) on row `sd / g / phi`
(`priorform-rev2-sdphi.R`). The written dpar-specific class row is
appended at the end of the table, and the vectorized pass looks only
upward. Only the display is wrong. Passing the table back fits
identically, because vectorized rows are skipped.

**R3 MINOR, a miss that base shares:** `(1 | mm(g1, g2)) + (1 | mm(g1,
g2, weights = cbind(wt1, wt2)))` is refused by brms ("Duplicated
group-level effects", group `mmg1g2`). Base and lane accept it, because
the mm group key is the whole label with its arguments
(`priorform-rev2-special.R`, `tw_mm_args`).

**R4 MINOR, shared by base:** `y ~ x + s(z) + s(z)` builds two smooths
with the same column name `s(z).fx1` on base and lane. logLik is
-252.0751169332, where the single smooth gives -252.0751169328. brms
collapses the pair to one term, as `terms()` does. The lane now
collapses `s(z) * s(z)` but not `s(z) + s(z)`, so the twin rule is
consistent for bar terms and not for smooths.

**Carried, not ranked, shared by base, filed for later:**

- A class `sd` prior with a group and no dpar reaches the `phi` block
  as well. brms applies it to `mu` only: `sd group g` then `sd` gives
  brms `sd/g/phi = student_t(3, 0, 2.5)` and the lane `normal(0, 5)`.
- In a multivariate model, `sd group = g` with no resp is accepted and
  broadcast. brms refuses it as matching no parameter.

The resolver did not change in round 1 for either behavior.

**NIT:** on `frm_sample()`, a MAP fit's `coef = "x"` density plus a
bounds-only `coef = "x"` spec on the call keeps the density and adds the
bound (`priorform-rev2-entry.R`). brms's `update(prior =)` puts the new
row in place of the old one (`brms:::update.brmsfit`, lines 87-90).

### What held

- **Entry points.** A same-slot pair is refused by `frm()`,
  `update(prior =)`, `validate_prior()`, `frm_simulate()`,
  `par_template()` and `frm_sample(prior =)`.
  - `frm_bootstrap()`, `frm_allfit()` and `refit()` take no prior and
    reuse the fit's own, which was already checked.
  - `frm_multiple()` passes its dots to `frm()`.
- **Slot keys against brms** (`priorform-rev2-slots.R`, 11 cases). The
  lane agrees with brms on 10 of 11; R1 is the exception.
  - Accepted by both: `sd` plus `sd group`; `Intercept` plus
    `Intercept dpar = sigma`; nlpar `b` plus nlpar `coef`; the
    multivariate `resp` and `coef` combinations.
  - Refused by both: a repeated `sd`; `coef = ""` against a missing
    coef; `nlpar coef = "Intercept"` twice; a density plus bounds on
    one slot.
- **Precedence against brms.** Validated tables, with brms's inherited
  rows filled from their parents:
  - The lane matches brms in both written orders for `b` coef against
    class, `cor` group against class, `b dpar = sigma` against `b`,
    nlpar `coef` against nlpar, and `Intercept dpar = sigma` against
    `Intercept`.
  - `sd group = g dpar = phi` against `sd group = g` also matches brms.
  - `sd` with `group` plus `coef` could not be compared, because frmtmb
    refuses coef on sd by design.
  - The analogous `update()` case matches brms: brms keeps the old user
    rows beside the new ones, so a stored `coef = "x"` prior still beats
    a new class `b` prior. The lane does the same
    (`priorform-rev2-entry.R`: `beta2 = normal(0, 1)` survives a call
    `class = "b"` `normal(0, 5)`).
- **The stated fit changes.** `dev/priorform-punch-order.R` reproduces
  every number: coef-then-class goes from -261.336784 to -409.935633,
  sd group-then-class from 211.855093 to 215.282086, and both
  class-first orders are `identical()`. A rerun `.rds` is `identical()`
  to the recorded one. My 17 round-trip cases now give identical fits
  and entries for 14. The two bounds-plus-density cases are refused by
  policy. `mv` has identical entries and a 2.27e-13 logLik difference.
- **The special refusal, no false alarm on 22 designs of mine**
  (`priorform-rev2-special.R`, brms, base and lane).
  - Accepted by brms and fitted by the lane identically to base:
    `s(z, by = x)`, `s(z, by = fg)`, `t2(x, z)`, `mo(mz) * w`,
    `w:mo(mz)`, `s(z) + offset(off)`, `s()` in an nlpar, `s()` in
    `sigma`, `poly(x, 2) + s(z)`, `s(z, g, bs = "fs")`, `mmc` inside a
    bar, and `s(z) + (1 | g)`.
  - Refused by both: `sigma ~ w * s(z)` and `(x + s(z))^2`.
  - The round also fixed three silent drops on base. `y ~ x +
    s(z) * s(z)`, `y ~ x + (s(z))` and `y ~ (x + s(z))` fitted `y ~ x`
    on base (-272.236). On the lane they fit the smooth, and the fit is
    `identical()` to `y ~ x + s(z)`.
  - `gp(x, by = fg)` and `me()` are refused for other reasons on both
    builds.
- **Twins and nesting, 13 designs.** The lane agrees with brms on 12 of
  13; R3 is the exception.
  - Two blocks: `(1 | g:h) + (1 | h:g)`.
  - Refused: `(1 | g:h) + (1 | g/h)`, `(x || g) + (x | g)`,
    `(1 + x | g) + (x + 1 | g)`, and `|p|` with `|q|` keys.
  - Collapsed: `(1 | g/h)` twice, `(x || g)` twice, `(1 | p | g)` twice,
    twins in an nlpar, and whitespace-only variants.
  - Every collapsed twin has an `identical()` optimizer result to the
    single term (`priorform-rev2-twinll.R`, 6 pairs).
  - Two `gr()` twins could not be compared, because base and lane
    already refuse `gr()` without cov and with `dist` plus cov.
- **Guards fail closed** (`priorform-rev2-mutants.R`, six test files,
  control 730 pass 0 fail).
  - Specificity order: identity fails 4 assertions, reversed rank 4,
    group weight removed 4, coef weight removed 3.
  - Slot check: no-op fails 5 plus 1 error; paren-insensitive key 1;
    bounds-only specs exempt 2.
- **No regression from round 0.** Part A: 28 of 28 fits are `identical()`
  to base in logLik, parameters, objective, `cov.fixed` and estimates,
  gp included. Part B, every Rd example of core and six extensions: 114
  of 114 identical (`priorform-rev2-bitwise-compare.txt`).
