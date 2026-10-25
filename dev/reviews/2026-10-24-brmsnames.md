# Review of lane `brmsnames` (2026-10-24)

Worktree `C:\Users\adf44\source\r\frmtmb-wt-brmsnames`, base 18b536fc
(frmtmb 0.74.0), nothing committed. Lane library
`C:/Users/adf44/source/r/wt-brmsnames-lib`; base library
`C:/Users/adf44/source/r/rellib-r12`; brms 2.23.0. Record:
`dev/brmsnames-1024-findings.md`.

Reviewer scripts are `dev/review-brmsnames-*.R` and
`dev/review-brmsnames-run.sh`; their logs are in
`dev/review-brmsnames-log/` and their saved names in
`dev/review-brmsnames-out/` and `dev/review-brmsnames-suite/` (all
gitignored). Every script takes the arm (`base` or `lane`) as its first
argument and states its seed in its header.

## Verdict: FIX

The three code changes are correct, and I could not break them. The
fixes this review asks for are in the record and the tests: the change
also fixes two silent wrong answers (`drop1()`/`step()` and
`re_formula`) that NEWS does not mention and no test pins (finding 1),
and two NEWS sentences say more than the code does (findings 2 and 3).
None needs a code change in the lane's three hunks. After findings 1
to 4 are done, SHIP.

## What I checked and found correct

- **The record.** `dev/brmsnames-findings.md` (the older lane's
  tracked record) is unchanged: `git status` and
  `git diff 18b536fc --stat` do not list it. Every lane script is
  `dev/brmsnames-1024-*`.
- **The installed build is the source.** `deparse()` of
  `unsub_specials`, `parse_linpred`, `prior_coef_brms` and
  `brms_fixef_rows` in the lane library is identical to the worktree's
  `R/` (`review-brmsnames-unsub.R`, log `unsub-lane.txt`).
- **The lane's tests, lane library, one file per process,
  `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`**
  (`review-brmsnames-run.sh lane tests1`). Every `lib:` line names
  `wt-brmsnames-lib`.

  | file | pass | fail | err | skip | warn |
  |---|---|---|---|---|---|
  | test-brms-names-specials.R | 177 | 0 | 0 | 0 | 0 |
  | test-v15.R | 41 | 0 | 0 | 0 | 0 |
  | test-setprior.R | 27 | 0 | 0 | 0 | 0 |
  | test-brms-priors.R | 103 | 0 | 0 | 0 | 0 |
  | test-priors-autocor-classes.R | 62 | 0 | 0 | 0 | 0 |
  | test-priors-bounds-grcov.R | 50 | 0 | 0 | 0 | 0 |
  | test-prior-sd-coef.R | 37 | 0 | 0 | 0 | 0 |
  | test-get-prior-route.R | 80 | 0 | 0 | 0 | 0 |
  | test-compat.R | 763 | 0 | 0 | 0 | 0 |
  | test-message-uniqueness.R | 6 | 0 | 0 | 0 | 0 |
  | test-brms-names.R | 138 | 0 | 0 | 0 | 0 |
  | test-brms-shapes.R | 82 | 0 | 0 | 0 | 0 |
  | test-brms-methods.R (hypothesis) | 1006 | 0 | 0 | 0 | 0 |
  | sample test-draws-names-specials.R | 4 | 0 | 0 | 0 | 0 |
  | sample test-draws-sd-names.R | 12 | 0 | 0 | 0 | 0 |

  There is no `test-hypothesis*.R` or `test-priors.R` file in core;
  `hypothesis()` is tested in `test-brms-methods.R` and others.
- **The name table, rebuilt for terms the lane did not cover.** brms
  2.23.0's names from `standata()`, `default_prior()` and a real
  brmsfit built as the lane builds one (`brm(empty = TRUE)`, a stanfit
  of brms's own code and data, 1 chain, 40 iterations, seed 1,
  `rename_pars()`), compiled with `R_MAKEVARS_USER` set
  (`review-brmsnames-brms.R`, data seed 2410 in
  `review-brmsnames-cases.R`). frmtmb by `review-brmsnames-frm.R`.
  "SAME" means identical in order. `set_prior` runs every coef row of
  brms's table through frmtmb's `set_prior()` and `validate_prior()`;
  `hyp` runs `hypothesis()` on every `b_` and `sd_` name of brms's
  `variables()`.

  | # | term | fixef / vcov, lane | base fixef | set_prior not ok, lane (base) | hyp not ok, lane (base) |
  |---|---|---|---|---|---|
  | 1 | `log1p(xp)` | SAME | SAME | 0 (1) | 0 (0) |
  | 2 | `abs(x)` | SAME | SAME | 0 (1) | 0 (0) |
  | 3 | `exp(x):z` | SAME | `.frm_expx:z` | 0 (1) | 0 (1) |
  | 4 | `I(exp(x))` | SAME | `I.frm_expx` | 0 (1) | 0 (1) |
  | 5 | `x + (exp(x) \| g)` | SAME | SAME | 0 (1) | 0 (1) |
  | 6 | `(1 + exp(x) \| gr(g, by = fg))` | SAME | `.frm_expx` | 0 (2) | 0 (3) |
  | 7 | `mvbf`, `exp(x)` in both | SAME | `y_.frm_expx` ... | 0 (2) | 0 (2) |
  | 8 | nl, `sigma ~ exp(z)` | SAME | order and `.frm_` | 0 (1) | 0 (1) |
  | 9 | categorical, `mu2 ~ exp(x)` | SAME | `mu2_.frm_expx` | 0 (1) | 0 (1) |
  | 14 | nl, `a ~ 1 + exp(z) + (1 \| g)`, `b ~ 0 + z` | SAME | `a_.frm_expz` | 0 (1) | 0 (1) |
  | 13 | `mvbf`, nl `yn` and linear `y` | order SAME, names DIFF | order and `.frm_` | 0 (2) | 3 (4) |

  Case 13's remaining difference is not the lane's (finding 6). On case
  8 brms gives `sigma_Intercept a_Intercept b_Intercept sigma_expz`
  and the lane gives the same; base gave `a_Intercept b_Intercept
  sigma_Intercept sigma_.frm_expz`. `variables()` differs from brms on
  every case only by brms's centered `Intercept` draw, which frmtmb's
  ML `variables()` does not list (not this lane).
- **`brms_rename()` collisions.** Three pairs that rename to one brms
  name, `log(xp)` beside a column `logxp`, `I(x - 1)` beside `IxM1`,
  and `poly(x1, 2)` beside `poly(x, 12)` (`polyx121` twice): brms
  2.23.0 stops in `standata()` with "Internal renaming led to
  duplicated names ... Occured for: 'log(xp)', 'logxp'", and frmtmb
  stops with the same message on both arms (cases 10 to 12), from the
  `brms_rename(cn, check_dup = TRUE)` the frame already calls
  (`R/frame.R:3978`). So no two columns of one design reach the prior
  table under one name. `y ~ mo(k) + mok` with a data column `mok` is
  also refused by frmtmb on both arms, where brms lists one `b mok`
  row (`review-brmsnames-mok.R`); frmtmb's refusal is the safer of
  the two and is not this lane's.
- **`unsub_specials()` does not undo what the alias is for.** One
  formula with user functions named `ar1()`, `toep()` and `diag()`
  (defined in the formula's own environment), `exp(x)`, and the
  structures `exp(pos + 0 | g)` and `ar1(0 + ft | g2)`
  (`review-brmsnames-unsub.R`, seed 77): all four fixed terms are kept
  (`Intercept ar1x toepx diagx expx`), both blocks keep their
  structures, `logLik` is -101.0629 (df 10) on both arms, and the
  population prediction equals `X %*% fixef()` to 0. `predict()` and
  `fitted()` on newdata, `re_formula = NA`, `update(. ~ . - toep(x))`,
  `update(newdata =)`, `drop1()`, `conditional_effects()`,
  `emmeans()`, `avg_slopes()`, `frm_sample()` (draws
  `b_ar1x b_toepx b_diagx b_expx`), `default_prior()` and
  `set_prior(coef = "ar1x")`, and `saveRDS()`/`readRDS()` then
  `fitted()` on newdata all work on the lane. Also `diag(1 + t | g)`
  and `toep(0 + ft | g2)` beside `exp(x)`, and `diag(1 + exp(x) | g)`.
  Nested and odd spellings agree with brms's `standata()` names on the
  lane: `base::exp(x)`, `exp(exp(x / 3))`, `exp(x)^2`, `log(exp(x))`
  (`review-brmsnames-nsexp.R`); base gave `.frm_exp.frm_expxD3`,
  `.frm_expx`, `log.frm_expx` for the last three.
- **The prior rename.** The old design-column spelling still moves the
  right coefficient for class `b`: with `normal(3, 0.001)` on
  `poly(x, 2)1`, `I(z^2)`, `log(xp)` and `factor(k)2`, and on brms's
  `polyx21`, `IzE2`, `logxp`, `factork2`, the one coefficient at 3 is
  always the named one, `validate_prior()` lists the user row under
  brms's name with no appended duplicate (9 rows each time),
  `prior_string()` takes both, and `frm_sample(prior =)` takes both
  (draws mean of `b_polyx21` and `b_IzE2` near 3)
  (`review-brmsnames-prior.R`, seed 31). On base, `polyx21` and `IzE2`
  were refused by `validate_prior()`, `frm()` and `frm_sample()`. No
  code reads a prior table row by position: the consumers are
  `fill_prior_table()` (by key) and the printed tables; no extension
  calls `prior_table()`, and no test indexes a prior table by number.
- **The nonlinear order.** Every `fixef()`/`vcov()` consumer in core
  and in all eight extensions reads by name (`multiple.R` pools by
  `rn`; `get_vcov()`, `get_coef()`, `influence.R` and frmtmb.latent's
  `hmm_scale_vcov()` use `interop_vcov()`/`vcov_estimated()`, which
  are in internal order and unchanged); a grep for a numeric index on
  `fixef()`, `vcov()` or `summary()$fixed` finds only a roxygen line.
  The lane's own runs of all nine suites in both tiers are clean.
- **The suite-wide diff.** `dev/brmsnames-1024-suitediff.R frmtmb`
  rerun on the lane's records reproduces the record: 274 files, 3753
  fits, 359 changed names, 53 distinct, none unaccounted. I re-recorded
  three files on both arms with a copy of the recorder that writes to
  my own directory (`review-brmsnames-suitenames.R`; `test-v15.R`,
  `test-emmeans.R`, `test-edgecases.R`): 27, 26 and 68 fits, and the
  names are `identical()` to the lane's records on both arms
  (`review-brmsnames-sncmp.R`). The diff does not record `fixef()`
  order, so it says nothing about change 3; the lane says so.

## Findings

### 1. MAJOR: two silent wrong answers fixed, and neither recorded nor pinned

The restored call fixes more than names. Neither fix is in NEWS, the
findings or a test.

- `drop1()` on a term that calls `exp()` (or any function named like a
  structure) reported `Df 0` and the full model's AIC, so the term
  looked free to drop at no cost; `step()` acted on that table.
  `y ~ exp(x) + z`, seed 3 (`review-brmsnames-drop1.R`):

  ```text
  base                          lane
  <none>          607.84        <none>    607.84
  .frm_exp(x)  0  607.84        exp(x)  1 625.69
  z            1  605.84        z       1 605.84
  step(): y ~ exp(x) + z        step(): y ~ exp(x)
  ```

  The AIC of `frm(y ~ z)` refitted is 625.689, the lane's row. Base's
  `step()` kept `z` although dropping it lowers the AIC to 605.84.
- `fitted()`/`predict(re_formula = ~ (1 + exp(x) | g))` stopped on base
  with "the re_formula term (1 + exp(x) | g) matches no group-level term
  of this fit. ... The fit's group-level terms are:
  (1 + .frm_exp(x) | g)". On the lane it works, and the prediction of
  `y ~ exp(x) + (1 + exp(x) | g)` with that `re_formula` equals the
  full prediction and equals, to every printed digit, the fit of the
  same model on a precomputed column `xx = exp(x)`
  (`review-brmsnames-reform.R`, seed 77: 1.226339906 1.745194277
  1.149327590 in all three).

Required: a NEWS bullet for each (a bug fix, with the 0.74.0 behavior
and the numbers), and two tests in `test-brms-names-specials.R`, seen
to fail on rellib-r12 (they do; logs `drop1-base.txt`,
`reform-base.txt`): `drop1()` of `y ~ exp(x) + z` has `Df` 1 for the
`exp(x)` row and its AIC equals `AIC()` of the refit without it (a
ratio, not a tolerance); and `fitted(re_formula = ~ (1 + exp(x) | g))`
is `identical()` (or equal at the ratio the suite uses) to `fitted()`
with `re_formula = NULL`.

### 2. MINOR: the prior table sorts as brms only within a class

The second bullet says "The table rows now sort by brms's name, as
brms sorts them." Within one class, the lane's order equals brms's.
Across classes it does not: frmtmb lists `Intercept` before `b` and
`sd` before `cor`, where brms lists `b`, `cor`, `Intercept`, `sd`,
`sigma`. That holds on all nine of my cases with those classes
(`review-brmsnames-order.R`), and it was so on 0.74.0, so it is not a
regression. Required: say "within a class" in NEWS (or sort the classes
as brms does, which is outside this lane).

### 3. MINOR: the design-column coef works for class `b` only

The roxygen of `set_prior(coef =)` and NEWS say the design column name
"also works", with no class named. For class `sd`,
`set_prior(class = "sd", coef = "I(x^2)", group = "g")` and
`coef = "poly(x, 2)2"` are refused ("No random-effect SDs match ...
group "g": "Intercept", "IxE2""), while `coef = "exp(x)"` is accepted,
because only the parentheses-dropped spelling matches there. On a
`cs(I(x^2))` term `class = "b", coef = "I(x^2)"` is refused
(`review-brmsnames-sdprior.R`, seed 9). Both arms refuse these, so it
is not a regression, but the new sentence promises it. Required: name
class `b` in the roxygen and NEWS sentence, or let the sd resolver
take the design spelling as `prior_coef_match()` does for `b`.

### 4. MINOR: no test pins the old spelling that NEWS promises

The only test in any suite that writes a design-column coef with
parentheses is `coef = "s(x).fx1"` (`test-brms-parity-defects.R:484`).
Nothing pins `coef = "I(x^2)"` or `"poly(x, 2)1"` after the table
moved to brms's names. Required: a test that a tight prior written
with the design spelling and with brms's spelling moves the same,
named coefficient, as `review-brmsnames-prior.R` does.

### 5. MINOR (pre-existing, to file): `fixef()` of draws, brms order

NEWS bullet 3 is about the fit. brms's `fixef()` is the draws one, and
frmtmb.sample's `fixef(draws)` keeps the template order on both arms
(`review-brmsnames-nldraws.R`, seed 12):

```text
bf(y ~ a + b * x, a ~ exp(z / 2), b ~ 1 + z, sigma ~ z, nl = TRUE)
  fit  : sigma_Intercept a_Intercept a_expzD2 b_Intercept b_z sigma_z
  draws: a_Intercept a_expzD2 b_Intercept b_z sigma_Intercept sigma_z
bf(y ~ x, sigma ~ z)
  fit  : Intercept sigma_Intercept x sigma_z
  draws: Intercept x sigma_Intercept sigma_z
```

The older record files the order of `variables()` and
`posterior_summary()` (its "Found and NOT fixed" item 1), not of the
`fixef()` matrix, which is a summary and could follow
`brms_fixef_rows()` without touching the draws. Required: file it in
the findings' "Defects found and not fixed", and say in NEWS bullet 3
that it is the fit's methods.

### 6. MINOR (pre-existing, to file): mv nonlinear prefix order

brms's prefix order is `dpar`, `resp`, `nlpar` (`vars_prefix()`,
`R/brmsterms.R:662`), so `mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + exp(z),
b ~ 1, nl = TRUE), bf(y ~ x + exp(z)))` has `yn_a_Intercept` and
`b_yn_a_Intercept`. frmtmb, on both arms, names it `a_yn_Intercept` and
`b_a_yn_Intercept`, and `hypothesis()` refuses 3 of brms's names
(case 13). The lane's order fix holds on this model. Required: file it
(it is a naming defect outside the placeholder scope).

### 7. MAJOR (pre-existing, confirmed): the mo() coding defect

Confirmed with brms and a fit (`review-brmsnames-mo.R`, seed 5, n =
300). brms 2.23.0 `standata()`: `mo(k0)` with 0..3, `mo(k)` with 1..4
and `mo(k2)` with 2..5 all give `Jmo = 3` and `Xmo_1` in 0..3. frmtmb
codes `D = max(x)` (`R/frame.R:3809-3810`), since 0.15.0 (65fed080,
2026-08-31, the commit that added mo()).

| term | frmtmb logLik (df) | Intercept (SE) | b (SE) |
|---|---|---|---|
| `mo(k0)`, 0..3 (brms's model) | -427.181869 (5) | 2.0834 (0.1138) | 1.5159 (0.0547) |
| `mo(k)`, 1..4 | -427.181869 (6) | 1.8520 (NaN) | 1.1948 (NaN) |
| `mo(k2)`, 2..5 | -427.181869 (7) | 1.9612 (NaN) | 0.9340 (NaN) |

Severity, by output:

- Point predictions are right: `fitted()` at k = 1..4 agrees across
  the three to 1e-6.
- ML coefficients are wrong: `Intercept`, `b` and the simplex are not
  identified, so `fixef()` reports values that are neither brms's nor
  unique, with NaN standard errors. A warning names the cause.
- ML prediction intervals are lost: 3 of 4 predictions have NaN
  standard errors (warning).
- `frm_sample()` samples another model than brms, with only ESS
  warnings. Intercept 0.912 (sd 1.330) against 2.098 (0.117), `b_mok`
  1.426 (0.336) against 1.510 (0.053), and a fourth cell
  `simo_mok1[1]` with posterior mean 0.176 on a category no row is
  in (2 chains, 1000 iterations, seed 11).
- frmtmb refuses negative integers, which brms takes after it
  subtracts the minimum.

Integer scales that start at 1 (Likert items) are the common input for
`mo()`, so this is a wrong Bayesian model on common data. It is out of
this lane's scope, as the lane says. Required: a lane in the next round
that subtracts the minimum, stores it for newdata (brms's
`attr(x, "min")`), and takes negative integers.

### 8. NIT: dead alias machinery and stale text

After `unsub_specials()` no formula refers to `.frm_<name>`, so these
are now dead or stale:

- the `assign()` of `.frm_<name>` into `prot` and `shared`
  (`R/parse.R:2022-2023`);
- the `shared` documentation (`R/parse.R:1736-1737`) and the
  `shared_env` comment (`R/parse.R:2502-2504`);
- `check_cs_in_bar()`'s `.frm_cs` branch and its comment
  (`R/frame.R:694-701`);
- `dev/rtmb-pitfalls.md` item 9, which describes the alias as the
  protection without saying it is undone after `splitForm()`.

Required: update pitfalls item 9 and the two comments; remove the dead
binding and branch, or leave them with a comment that says why.

### 9. NIT: the frmtmb.sample floor is still open

The new sample test fails on frmtmb 0.74.0 (the lane's seen-fail log)
and `extensions/frmtmb.sample/DESCRIPTION` still has
`frmtmb (>= 0.74.0)`. The lane names this. Required at consolidation:
raise the floor to the release that carries the change, or skip the
test on an older core.

## Overlaps with other lanes

The lane touches `R/parse.R`, `R/priors.R`, `R/brms-shapes.R`,
`NEWS.md`, `man/set_prior.Rd` and `extensions/frmtmb.sample/NEWS.md`,
all named in the round's overlap list. `R/parse.R` is the likely
merge conflict with sarfcor (`fcor()`, the structure seam); the hunks
are local to `parse_linpred()` after `splitForm()`.

## Short check of fix round 1

Verdict: SHIP. Finding R1 (minor) is a one-line change; make it before
merge or file it. Finding R2 is a nit.

The fix round is recorded in "Fix round 1" of
`dev/brmsnames-1024-findings.md`. It now also changes `R/brms-names.R`,
`R/frame.R`, `dev/rtmb-pitfalls.md`,
`extensions/frmtmb.sample/R/methods-draws.R` and
`tests/testthat/test-simplex-mvbf.R`. New reviewer scripts:
`dev/review-brmsnames-fr1.R`, `-fr1vp.R`, `-fr1cmp.R` and
`-fr1draws.R`, data seed 21, plus case 15 in
`review-brmsnames-cases.R`. Logs are in `dev/review-brmsnames-log/`.

### What I checked

- **The installed build matches the source.** In the lane library,
  `deparse()` of the following functions is identical to the
  worktree: `unsub_specials`, `parse_linpred`, `parse_one_response`,
  `prior_table`, `prior_rows_sd_order`, `prior_coef_brms`,
  `prior_coef_match_brms`, `prior_spec_slot`, `sd_spec_reach`,
  `sd_coef_shared`, `resolve_priorlist`, `brms_lp_prefix`,
  `brms_fixef_rows`, `check_cs_in_bar`, and frmtmb.sample's
  `fixef.frmtmb_draws`.
- **Tests.** I ran these on the lane library, one file per process,
  `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`
  (`review-brmsnames-run.sh lane fr1`). Every `lib:` line names
  `wt-brmsnames-lib`. Every file has fail 0, err 0, skip 0, warn 0:

  | file | pass |
  |---|---|
  | test-brms-names-specials.R | 253 |
  | test-simplex-mvbf.R | 26 |
  | test-setprior.R | 27 |
  | test-brms-priors.R | 103 |
  | test-brms-formula-priors.R | 210 |
  | test-prior-compat.R | 199 |
  | test-prior-location-dpar.R | 55 |
  | test-prior-mv-resp.R | 36 |
  | test-prior-sd-coef.R | 37 |
  | test-prior-sd-scope.R | 25 |
  | test-prior-support.R | 150 |
  | test-priors-autocor-classes.R | 62 |
  | test-priors-bounds-grcov.R | 50 |
  | test-get-prior-route.R | 80 |
  | test-update-prior-drop.R | 29 |
  | test-nlf.R | 77 |
  | test-compat.R | 763 |
  | test-message-uniqueness.R | 6 |
  | test-cs-bar-brms.R | 9 |
  | sample test-draws-names-specials.R | 9 |
  | sample test-draws-sd-names.R | 12 |

- **F2, the row order.**
  - The `lncor` rows: brms 2.23.0's `default_prior()` of the two
    `logistic_normal()` responses of `test-simplex-mvbf.R`
    (`smv_data(2)`) lists `V` then `Y`. The lane lists `V Y` and
    0.74.0 listed `Y V`. On the rows the two packages share, the
    lane's whole table is in brms's order.
  - Order against brms on 8 models (`review-brmsnames-fr1.R`): one
    exp model with a random slope, `sigma ~ z`, a 2-response `mvbf`
    with random effects, a nonlinear model with a random effect,
    `s() + mo() + (1 | g)`, `cs()` with `acat`, `ar()`, and `student`.
    `default_prior()` is in brms's order on the shared rows for all 8
    on the lane and for none on 0.74.0.
  - The same for `validate_prior()` with priors each model takes, on
    8 cases (`review-brmsnames-fr1vp.R`; the nl case has two nlpar
    priors, the mv case has resp priors). All 8 agree with brms on
    the lane and none does on 0.74.0.
  - Only the order changed. On 22 tables I compared the whole
    rows, `prior`, `class`, `coef`, `group`, `resp`, `dpar`, `nlpar`,
    `lb`, `ub` and `source`, between the arms
    (`review-brmsnames-fr1cmp.R`; for 0.74.0 the coef is put in
    brms's spelling with the alias dropped). The sorted rows are
    identical in all 22 and the order differs in all 22.
  - The 12 cases of my name table (1 to 9, 13 to 15) are also in
    brms's order on the shared rows (`review-brmsnames-order.R`).
  - The suite diff rerun reproduces the record: core 4743 changed
    names, of which 4380 are "prior (row order)", and frmtmb.sample
    1175, all of them "prior (row order)". No name is unaccounted.
    That classifier compares the coef multiset only; the full-row
    comparison above is the independent check.
  - Nothing reads prior table rows by position. A grep of `R/`,
    `extensions/*/R` and every test directory for a numeric row index
    on a prior column finds only `pr$prior[2]` and `pr$class[1]` in
    `test-brms-formula-priors.R` and in the `$<-.frmtmb_priorlist`
    example. Both edit a user's `set_prior()` list, which is not
    sorted. No `head()` of a prior table appears anywhere.
- **F3, one slot per coefficient** (`review-brmsnames-fr1.R`). These
  pairs are refused as duplicates on the lane: `I(x^2)` + `IxE2`,
  `poly(z, 2)1` + `polyz21`, `log(xp)` + `logxp` (class `b`),
  `I(x^2)` + `IxE2` (class `sd`, group `g`), and `(Intercept)` +
  `Intercept` on `0 + Intercept`. Six different coefficients in mixed
  spellings are accepted, and a `normal(m, 0.001)` prior on each puts
  each estimate at its own m: `logxp 1; IxE2 2; polyz21 3; polyz22 4;
  factork2 5; factork3 6`. No two coefficients of one block can share
  a slot. The frame refuses a fixed or group-level design whose
  columns rename to one name: `(1 + log(xp) + logxp | g)` gives
  "Duplicated group-level effects are not allowed", and brms refuses
  it too.
- **F6, a model with two responses, both nonlinear, with sigma and a
  random effect** (case 15: `mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + z +
  (1 | g), b ~ 1, sigma ~ z, nl = TRUE), bf(y ~ c + d * x, c ~ exp(z),
  d ~ 1 + (1 | g), sigma ~ x, nl = TRUE))`). I compared against a
  compiled brmsfit (`review-brmsnames-brms.R 15`).
  - `fixef()` and `vcov()` are identical to brms in order and names:
    `sigma_yn_Intercept sigma_y_Intercept yn_a_Intercept yn_a_z
    yn_b_Intercept sigma_yn_z y_c_Intercept y_c_expz y_d_Intercept
    sigma_y_x`.
  - `variables()` differs only by brms's centered `Intercept_sigma_*`
    draws.
  - The prior table is in brms's order on the shared rows, every
    brms coef row is accepted, and `hypothesis()` takes every brms
    `b_` and `sd_` name.
  - Case 13 (finding 6) now matches brms everywhere.
  - On 0.74.0, case 15 does not fit: "could not find function
    ".frm_exp"" (see R2).
- **F5, `fixef()` of the draws** (`review-brmsnames-fr1draws.R`,
  `frm_sample()` with 1 chain, 200 iterations, seed 3). On cases 8,
  13 and 15 the row names are identical to brms's `fixef()` and to
  the fit's. `summary = FALSE` has its columns in the same order and
  keeps `nchains`. Each reordered column still holds its own `b_`
  draws (4 of 4, 6 of 6 and 10 of 10 compared with `identical()`).
  `review-brmsnames-nldraws.R`, which showed finding 5, now gives the
  fit's order for both models.
- **F8, the alias machinery removed.** I reran
  `review-brmsnames-unsub.R` on the lane. It uses user functions
  `ar1()`, `toep()` and `diag()` from the formula's own environment,
  `exp(x)`, `exp(pos + 0 | g)` and `ar1(0 + ft | g2)`. The output is
  unchanged from the first review: `logLik` -101.0629 (df 10); the
  population prediction equals `X %*% fixef()` to 0; and the
  newdata predictions are 1.879339 2.485717 3.557874 5.014105. These
  also all work: `update()`, `update(newdata =)`, `drop1()` (`ar1(x)`
  `Df 1`), `conditional_effects()`, `emmeans()`, `avg_slopes()`,
  `default_prior()`, `set_prior(coef = "ar1x")`, `frm_sample()`
  (`b_ar1x b_toepx b_diagx b_expx`), `saveRDS()`/`readRDS()`, and the
  `diag()`/`toep()`/`diag(1 + exp(x) | g)` fits. `test-cs-bar-brms.R`,
  which tests the removed `.frm_cs` path, passes (9).

### R1. MINOR: a smooth's two spellings are still two slots

NEWS now says the two spellings of a coefficient "name one slot, so
writing both is refused as a duplicate". That is false for the
unpenalized column of a smooth. On `y ~ s(x) + z`,
`c(set_prior("normal(3, 0.001)", coef = "s(x).fx1"),
set_prior("normal(-3, 0.001)", coef = "sx_1"))` is accepted.
`validate_prior()` lists one user row (`sx_1`, `normal(-3, 0.001)`),
and the fit puts `sx_1` at -3. So the first prior is dropped with no
word. 0.74.0 does the same, so it is not a regression.

The cause: `prior_spec_slot()` keys by `brms_rename()`, which makes
`s(x).fx1` into `sx.fx1` and not `sx_1`. Required:
key by `prior_coef_brms()`, which already maps `.fx<k>` to brms's
`_<k>`, and add the pair to the duplicate test; or file it and narrow
the NEWS sentence.

### R2. NIT: one more 0.74.0 crash fixed, not in NEWS

0.74.0 cannot fit case 15: "could not find function ".frm_exp"". The
`exp()` in the second response's nonlinear formula did not find the
alias. The lane fits it. Add it to the bug-fix bullets or to the
findings.
