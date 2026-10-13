# Lane moparity: brms features still refused or diverging

Round opened 2026-10-07, base 48c1e6c9 (frmtmb 0.69.0, frmtmb.sample
0.17.0), worktree `frmtmb-wt-moparity`, private library
`wt-moparity-lib`. "Base" is `rellib-r7`. Every number below is from a
log in `dev/moparity-log/` (gitignored) written by the script named
beside it; the test runs are in `dev/release/moparity-*-files/`.

## Punch round 2

Review "Re-check after punch round 1" in
`dev/reviews/2026-10-08-moparity.md`: NOT MERGEABLE on B2; r1, r2,
r3. "Seen to fail" is the round-1 lane build.

**B2. An empty row of a prior table takes its slot.**
`update_merge_prior()` (R/confint.R) reads the slot of every row of a
`brmsprior` or `frmtmb_prior_rows` table before the translation, which
skips empty rows (new helper `prior_table_slots()`), and drops the
stored specification of each of those slots as well as of the
translated rows'. A prior with no rows (`NULL`, `empty_prior()`, a
zero-row table) keeps the stored prior. A legacy named list stored and
a table given is refused as a mix of spellings. Comment, `?update`
and NEWS now say what brms does (update.R lines 150 to 155). The
reviewer's script on the fixed build
(`dev/moparity-rev2-update.R`, `dev/moparity-log/p2-rev2-update.txt`),
against brms 2.23.0 (`dev/moparity-rev2-log/brms-update.txt`):

| case | brms rows, b_x | lane before | lane now |
|---|---|---|---|
| frmtmb `default_prior()` table, sigma edited | b flat, b_x 0.6133 | b x kept, 0.0520 | sigma only, 0.6239 |
| brms `default_prior()` table, sigma edited | Intercept the table's default, sigma | (as above) | Intercept `student_t(3, 0.9, 2.5)`, sigma, 0.6239 |
| one row `b x`, prior `""` | b x flat, Intercept kept, 0.625 | b x kept, 0.0511 | Intercept kept, 0.6239 |
| `brms::empty_prior()` | (kept) | kept | kept, 0.0511 |

(brms's b_x is a posterior mean, frmtmb's an ML estimate under the
same prior rows.) Duplicate new rows and a new row for a removed term
stay refused, the accepted divergence. Test "a prior table's rows take
their slots, empty rows too" replaces "... keeps a flat one", and the
brmsprior test adds an empty `brms::prior_string("", class = "b", coef
= "x")` row. Seen to fail: tier `b2pre`, 26 pass, 3 fail; now `p2a`,
29 pass.

**r1. A correlated `rp_cor_case()`.** The test's data now has a true
random-effect correlation of 0.6 (`u[2, ] <- 0.6 * u[1, ] + 0.8 *
u[2, ]`), iter still 1000, gate |z| < 5. Over sampler seeds 1 to 24
(`dev/moparity-rp-gate.R` with `RP_DESIGN=30,240,0.8,0.6,0.6`,
`dev/moparity-log/rp-gate-cor06-sum.txt`): base max |z| 3.38, lane
3.60, 0 of 48 runs at 5 or more; smallest centered ESS 45. Power: the
reviewer's transposed factor (`b = L' z`) fails the gate on 12 of 12
seeds (`dev/moparity-rev2-log/gate-newcor-1-1.txt`) and on the test's
seed 33 at |z| 11.53, which passes unmutated at 2.76
(`dev/moparity-log/p2-gate-mut-33.txt`, `p2-gate-ok-33.txt`). Test
file: lane 371 pass (`p2rp`); base passes the gate, its 3 other
failures are the lane's gp non-centering (`p2rpbase`).

**r3. A column or nonlinear parameter named like a special.** Every
special-name test in R/parse.R that used `all.names()` now asks for a
CALL (`calls_function()`, new `called_functions()`): the `mo()`
interaction, group-level and nested checks, `mi()`, `mmc()`, the
crossed-autocorrelation refusal and its term name, and the
autocorrelation check of a nonlinear body. `dev/moparity-r3-names.R`
tries 32 special names (mo, mi, me, mmc, mm, cs, gp, s, t2, te, ar,
ma, arma, cosy, unstr, car, sar, fcor, cens, trunc, se, rate, trials,
dec, thres, sm, rr, us, diag, ou, toep, cv) as a column in `y ~ v *
z`, `y ~ v:z`, `y ~ v + (1 + v || g)`, `y ~ v + (1 + v | g)` and as a
nonlinear parameter: 23 failures on 0.69.0 and on the round-1 build
(mo and mi in all four formulas; ar, ma, arma, cosy, unstr crossed and
as a nonlinear parameter), 0 now (`dev/moparity-log/r3-names-base.txt`,
`r3-names-pre.txt`, `r3-names-post.txt`). The other names already
passed: their checks look at the call head. New test file
`test-special-names.R` (same fits as with a plain name, to 1e-10
relative logLik; and the real specials still refused where they
should be). Seen to fail: tier `r3pre`, 2 errors (the third failure
there was a fault of the test, fixed).

**r2** (every `||` term as brms reads it) is filed in
`dev/test-backlog.md`, section "Filed by lane moparity after punch
round 1".

**Runs after the round** (lane library, one file per process, gates
on). `test-update-prior-drop.R` 29, `test-special-names.R` 37 (tier
`p2a`); frmtmb.sample `test-reparam.R` 371 (`p2rp`); all clean. Core
and frmtmb.sample suites, tier `p2full`: 278 of 278 files, frmtmb 227
files 18722 pass, frmtmb.sample 51 files 2803 pass; 0 fail, 0 err, 0
warn; the one skip is frmtmb.sample's `test-scale.R`. R CMD check
--as-cran after the last source edit, run one after the other: frmtmb
`Status: 1 NOTE` (V8), frmtmb.sample `Status: OK`. No process of this
lane is left running.

## Punch round 1

Review `dev/reviews/2026-10-08-moparity.md`: NOT MERGEABLE on B1;
minors m1 to m10. All addressed below; every "seen to fail" run is the
lane build installed before the fix (`wt-moparity-lib`), tier logs
`dev/release/moparity-<tier>.log`.

**B1. `update(prior = )` with a brmsprior or a prior table.**
`update_merge_prior()` (R/confint.R) now reads both sides through
`as_priorlist()` before it merges, so a `brmsprior`, a
`default_prior()`/`get_prior()` table and `set_prior()` merge the same
way. `NULL`, and a table with no prior in it, keep the stored prior
(wrong for a table with empty rows, which brms reads as resetting
their slots; corrected in punch round 2, B2). A
named list of prior objects (the legacy spelling) merges by name with a
stored named list; a mix of that spelling and specifications is
refused, because one prior cannot hold both. The reviewer's script on
the fixed build (`dev/moparity-rev-update.R`,
`dev/moparity-log/p1-rev-update.txt`): `brms::prior(exponential(2),
class = sigma)` now gives b_x 0.081058 with all three rows, equal to
the `set_prior()` spelling, where it gave 0.45576 with the sigma row
only. Tests: three in `test-update-prior-drop.R` (brms object, also on
a fit made with a brmsprior; a `default_prior()` table and a flat one;
the named list and the refused mix). Seen to fail: tier `b1pre`, 18
pass, 8 fail; after: `b1post`, 26 pass.

**m1.** frmtmb.sample NEWS now states the logged numbers
(`dev/moparity-log/dir-sample-cmp.txt`): Dirichlet comparison, means
within 1.71 MCSE and sds within 2.51; `mo(x) * f2` (two levels) and
`mo(x) + (mo(x) | g)`, within 2.12.

**m2.** The `chol_L` comment in R/covstruct.R cites the logged base
run (`dev/moparity-log/gp-cmp.txt`: 1286 of 2000 transitions at the
maximum tree depth, R-hat 1.52, bulk ESS 7). The "R-hat 1.56 to 1.89 at
600 iterations" came from lane surface's backlog entry, with no log
here; it is removed. `dev/moparity-gp-init.R` was rerun with its log
kept (`dev/moparity-log/gp-init.txt`): chain 1 of `last.par.best` has
326 divergences and mean sigma 0.2558, chains 2 to 4 have 24 to 38,
the `random` chains 9 to 27, as the item 7 text says.

**m3. The correlated `test-reparam.R` gate.** Measured with
`dev/moparity-rp-gate.R` (per seed: z of every mean and sd, in units
of the two chains' own MCSEs from `posterior::mcse_mean()` and
`mcse_sd()`), summary `dev/moparity-rp-gate-sum.R`:

| data | arm | runs | max \|z\| | runs with max \|z\| >= 5 | min ESS centered |
|---|---|---|---|---|---|
| test's (12 groups, sds 0.7, 0.4), iter 1000 and 2000 | base | 47 (+1 NA) | 13.30 | 10 | 1 |
| same | lane | 48 | 5.43 | 2 | 2 |
| 30 groups, 240 rows, sds 0.8, 0.6, iter 1000 | base | 24 | 4.16 | 0 | 82 |
| same | lane | 24 | 4.25 | 0 | 82 |

(`dev/moparity-log/rp-gate-orig-sum.txt`,
`dev/moparity-log/rp-gate-30_240-sum.txt`.) On the test's data no
gate on the centered chain's draws is sound: a centered chain of 1000
sticks on some seeds (ESS 1 to 2; seed 18 on base gives an NA), so its
MC error cannot be estimated from the chain. On 30 groups the centered
chain mixes (ESS at least 82 over 48 runs) and the MCSE z is close to
calibrated (median per-run max 2.79 over 132 statistics). The test now
uses that data, at 1000 iterations again (0.69.0's count, half the
lane's 2000), and gates |z| < 5 on every mean and sd, the group-level
values included; a wiring bug that moves a mean by one posterior sd is
about 10.8 z units at the median ESS (3.8 at the smallest). Pilots
`dev/moparity-log/rp-gate-pilot-*.txt`: 12 groups with 240 rows and 20
groups with 400 rows still gave centered ESS 10 and 22. The test's seed
33: lane 371 pass, base the gate passes too (tiers `p1s2`, `p1sbase`;
base's other 3 failures are the lane's gp non-centering, expected).
posterior's "ESS has been capped" warning (base, one statistic) is
allowed by name: a capped ESS only enlarges the MCSE.

**m4.** `default_prior_notes(fit, user_pl)` (frmtmb.sample R/sample.R)
leaves out every simplex that a Dirichlet row of the run's prior
covers (resolved entries, as `mo_simplex_nlp()` reads them), for
`simo` and `sbhaz`; the note now reads "on each mo() simplex without a
simo row of its own (...)". Tests: `test-mo-simplex-draws.R` (two
simplexes, a row on one, then on both: no note) and `test-cox-draws.R`
(a `sbhaz` row removes the `sbhaz` note). Seen to fail: tier `m4pre`
(2 fail, 1 error).

**m5.** `?frm` "Monotonic effects", `?set_prior` class `"simo"` and
NEWS say it in plain terms: the MAP mode under `dirichlet(a)` is the
mode of `dirichlet(a + 1/2)` on the weights, pulled further toward the
uniform simplex; `dirichlet(1)` is not the absence of a prior; only no
`simo` row applies no density. `?frm` quotes the review's sizes
(`dirichlet-conv.txt`: 60 rows, 0.30 0.38 0.32 against ML 0.05 0.73
0.22; 400 rows, at most 0.003).

**m7. `(mo(x) * f || g)`.** `expand_double_verts()` (R/parse.R) now
makes a `||` term with `mo()` ONE `diag()` block over the whole term,
which is brms's reading of `||`: `(mo(x) * f2 || g)` has `(Intercept),
mox, f2b, mox:f2b`, brms's columns (`dev/moparity-rev-dblbar.R`,
`dev/moparity-log/p1-dblbar.txt`). A `||` term without `mo()` keeps
lme4's expansion (`(f2 || g)` still has `f2a, f2b`). The NEWS item and
`?frm` say so. Test in `test-mo-factor-re.R`: the six `sd_g__` names,
one `diag` block, and the logLik of the indicator spelling within 1e-8
relative. The gated brms test of `(1 + mo(x) || g)` (now one block
rather than two) still agrees with brms's program. Seen to fail: tier
`m78pre` (the refusal "mox:fa").

**m8.** The quadrature refusal for group-level `mo()` is its own
message ("quadrature = TRUE does not support group-level mo() terms
..."), checked after the dim-1 test, so `(0 + mo(x) | g)` names the
real reason. Test in `test-mo-factor-re.R`, with the absent case
(`mo(x) + (1 | g)` fits with quadrature). Seen to fail: tier `m78pre`.

**m9.** `test-cox-chart.R`: a weight at the face is below 1e-8 of the
largest weight, the MAP baseline's smallest weight above 1e-3 of its
largest. `test-cox-draws.R`: the same two as ratios to the largest
weight. `test-dirichlet-prior.R`: the Jacobian gap is relative to the
summed magnitude of the density's terms (no `max(1, .)`), the
missing-Jacobian inverse must exceed 1e4 times the largest gap, and the
MAP-versus-ML distance must exceed 1e3 times the distance between two
searches for the MAP. `test-reparam.R` gp gradient: the finite
difference's own error is measured (steps 1e-5 and 5e-6, plus the
rounding of fn) and the gradient must be within 10 times it.

**m6, m10** are filed in `dev/test-backlog.md` ("Filed by lane
moparity after punch round 1").

**For consolidation (record only).** R/fit.R `mo_search()` tail
conflicts with lane optrap: keep this lane's step 3
(`cox_face_leave()`) before optrap's tail, and take optrap's
settle-always tail (`obj_settle()`/`obj_state_at()`) in place of this
lane's `probed` restore; optrap's step 2 keeps `if (is.na(tm$b))
next`, so a cox term stays out of it. optrap also edits the
`mo_search()` roxygen where this lane adds item 3. Delete the
reviewer's library `C:/Users/adf44/source/r/moparity-rev-lib` (brms
2.22.0 only) afterwards. frmtmb.sample's floor must move to this
round's frmtmb (`simplex_frame_terms()` is new); bump minor for both.

**Runs after the round** (lane library, one file per process, gates
on). Affected files: `test-update-prior-drop.R` 26, `test-mo-factor-re.R`
44, `test-dirichlet-prior.R` 24, `test-cox-chart.R` 17 (tiers `p1a`,
`p1b`); frmtmb.sample `test-reparam.R` 371, `test-mo-simplex-draws.R`
40, `test-cox-draws.R` 12 (`p1s`, `p1s2`); all 0 fail, 0 err, 0 skip, 0
warn. Core and frmtmb.sample suites, tier `p1full`: 277 of 277 files,
frmtmb 226 files 18682 pass, frmtmb.sample 51 files 2803 pass; 0 fail,
0 err, 0 warn; the one skip is frmtmb.sample's `test-scale.R` (scale
tier). R CMD check --as-cran (`dev/moparity-check/`): frmtmb `Status: 1
NOTE` (V8 unavailable for the HTML manual's math), frmtmb.sample
`Status: OK`. A first frmtmb check run beside frmtmb.sample's failed in
`R CMD build`, which copied a file the other check held open; the rerun
alone is the one reported.

## Summary

| item | outcome |
|---|---|
| 1. `mo(x) * f`, `(mo(x) \| g)` | implemented; log densities agree with brms's compiled programs (1e-6 relative) |
| 2. `simo` (and `sbhaz`) Dirichlet priors | implemented in both charts; MAP and sampler validated |
| 3. `update(prior = )` | merges as brms does; new rows still checked |
| 4. `cox()` baseline | moved to the `mo()` chart, a face search, the sampler's softmax + `dirichlet(1)` |
| 5. multi-scale floors | `sd_floored()` on nine structures; bitwise unchanged above log sd -308.45 |
| 6. `parnames()` | decided from the dispatching frame |
| 7. exact `gp()` sampling | non-centered, as brms's `zgp` |
| 8. `as_tmbstan()` | contract kept (raw tape); the message names only density-free simplexes |

Version bump: minor for both packages (new features and breaking
changes). frmtmb.sample's floor is the frmtmb of this round: it calls
`simplex_frame_terms()`, new, and reads Dirichlet entries
`resolve_prior_input()` now returns.

## Item 1: `mo()` crossed with a factor, and in group-level terms

### What brms 2.23.0 does (`dev/moparity-brms-code.R`, `brms-code.txt`)

- `y ~ mo(x) * f` with a three-level `f`: brms stops building the model,
  "replacement has 2 rows, data has 3" in `frame_sp()`.
- `y ~ mo(x):f` and `y ~ f + mo(x):f`: brms makes three simplexes, three
  `bsp`, and Stan data `Csp_1`, `Csp_2`, `Csp_3`, but its Stan code reads
  `(bsp[1]) * mo(simo_1, ...) * Csp_1[n] + (bsp[2]) * mo(simo_2, ...) *
  Csp_1[n] + (bsp[3]) * mo(simo_3, ...) * Csp_1[n]` and declares only
  `Csp_1`. The cause is `frame_sp()`: `attr(mm, "covars")` has one entry
  per formula TERM and `Ic` is recycled over the model-matrix COLUMNS.
  brms#1828 is this bug, acknowledged by the owner ("You are right, this
  is indeed a bug in brms"), milestone 2.24.0; brms 2.22.0 fitted it.
- `y ~ mo(x) * f2` with two levels: correct (one covariate column).
- `y ~ mo(x) + (mo(x) | g)`: `(bsp[1] + r_1_2[J_1[n]]) * mo(simo_1, .)`,
  so the group-level coefficient varies the population coefficient on
  the SAME simplex. `y ~ (mo(x) | g)` with no population term: brms
  computes `r_1_2` and never uses it; its own check ("Special
  group-level terms require corresponding population-level terms")
  runs only when the formula has some population-level special term.

brms is clearly wrong on the factor case and on the silent drop, so
frmtmb fits the model brms's Stan data describe (brms 2.22.0's) and
refuses the drop in brms's own words.

### Base (`dev/moparity-refusals.R`, `refusals-base.txt`)

```text
== y ~ mo(x) * f
ERROR: mo() interactions support numeric multipliers only: f is factor; expand it to numeric indicator columns first
== y ~ mo(x) + (mo(x) | g)
ERROR: mo() is only supported as a standalone additive term or a two-way ':'/'*' interaction: (mo(x) | g)
```

### The change

- `R/parse.R`: an interaction term with a `mo()` call is expanded by
  `terms()` (`mo_term_entries()`): each label with a `mo()` factor is a
  special term, the other labels go to the fixed part, and each label
  keeps its factors in the order `terms()` of the whole formula gives
  them. Three-way terms (`mo(x) * z * w`) are now read; `mo(x):mo(w)`
  is still refused. A bar term with `mo()` on its left has each
  call replaced by a placeholder symbol `.frm_mo<k>` (`mo_bar_rewrite()`).
- `R/frame.R`: `mo_sp_design()` builds brms's sp model matrix: one
  formula over every `mo()` term of the predictor, each `mo()` call a
  column of ones, so `model.matrix()`'s own rule decides contrasts
  versus full indicators, and its column names are brms's (it reads the
  special terms alone, as brms's `sp_model_matrix()` does: `z * mo(x)`
  is `mox:z`, `z + mo(x):z` is `z:mox`). One simplex and one
  coefficient per column;
  `lp$mo_mm` keeps the terms, levels and contrasts for new data. For the
  group-level term, the placeholder is a column of ones in the model
  frame, the component's columns are renamed as brms names them (`mox`,
  `mox:z`), and `mo_re_link()` maps each to the population column it
  varies (`lp$Z_mo`), refusing a term with none.
- `R/objective.R`: `mo_re_zterm()` adds `(Z_cols b) * (m - 1)` to `Z b`,
  `m` the term's transform (`mo_transform()`), so the tape of a model
  without such a term is unchanged.
- `R/predict.R`: `lp_Z()`, the design at the estimates with the
  group-level columns scaled, for every numeric reader of `lp$Z` as a
  design (in-sample eta, the delta-method Jacobian, `se-check.R`'s
  separation certificate); new data rebuild the placeholder and scale
  it at the new rows' categories; factor multipliers on new data come
  from `mo_sp_newdata()`.
- `quadrature = TRUE` and `importance` refuse a group-level `mo()` term
  (their designs are formed off the tape).

### Validation

Against brms's compiled programs at frmtmb's estimates
(`tests/testthat/test-mo-factor-re.R`, gated; `brms_lp_check()` of
`helper-brms.R`: log density equal to 1e-6 relative after the flat
Dirichlet constants, gradient below 1e-3, the joint density and the
`z` gradient for random effects):

- `mo(x) * f` (frmtmb, factor) against brms's `fb + fc + mo(x) +
  mo(x):fb + mo(x):fc` (numeric indicators, where brms 2.23.0 is right),
  constant `3 lgamma(3)`.
- `mo(x):f` against `mo(x):fa + mo(x):fb + mo(x):fc`.
- `mo(x) + (mo(x) | g)`, joint, constant `lgamma(3)`.
- `mo(x) * f + (1 + mo(x) || g)` against the indicator spelling, joint.

All pass (`moparity-t2`, `moparity-t3`: 39 of 39). Names against brms's
`default_prior()` (`dev/moparity-brms-names.R`, `brms-names.txt`):
`z + mo(x):z` is `z:mox` in brms (frmtmb 0.69.0 said `mox:z`, now
`z:mox`), `f2 + mo(x):f2` is `f2a:mox`, `f2b:mox`, `mo(x) * z * w` is
`mox, mox:z, mox:w, mox:z:w`; the test asserts the simo coefs equal.
A logical multiplier is `moincome:xlTRUE` in brms and now here
(`dev/moparity-brms-lgl.R`).

frm_sample() against brms 2.23.0, 4 chains of 2000, seed 1
(`dev/moparity-dir-sample.R`, compared by `dev/moparity-dir-sample-cmp.R`,
`dir-sample-cmp.txt`): on `y ~ mo(x) + (mo(x) | g)` the largest |z| of
a mean or sd difference over its joint Monte Carlo error is 1.12
(`bsp_mox`, three weights, `sd_g__Intercept`, `sd_g__mox`), on
`y ~ mo(x) * f2` 2.11; 0 divergences on both sides.

`conditional_effects()` on `mo(x) * f` equals `fitted(newdata)` at its
grid exactly (max difference 0), and runs on `(mo(x) | g)` with and
without `re_formula = NULL` (`dev/moparity-ce.R`, `ce.txt`).

Seen to fail on rellib-r7 (`dev/release/moparity-base2-files/`):
`test-mo-factor-re.R` 9 errors ("numeric multipliers only", "only
supported as a standalone additive term", the three-way refusal);
`test-open-issues.R`'s updated test, 1 error.

## Item 2: the Dirichlet on a simplex

### brms 2.23.0 (`dev/moparity-brms-simo.R`, `brms-simo.txt`)

`dirichlet(c(5, 1, 1))` on `coef = mox1` sets `con_simo_1`; a scalar is
recycled (`dirichlet(2)` gives `2 2 2`); a class-wide row is refused
("do not correspond to any model parameter"); a wrong length is refused
("Invalid Dirichlet prior. Expected input of length 3."); any other
density is refused ("Currently 'dirichlet' is the only valid prior for
simplex parameters"); bounds are refused. `default_prior()` lists one
`dirichlet(1)` simo row per simplex (`mox:z1`, `mox1`, and `mox1` with
`dpar = sigma`). frmtmb follows each rule, with brms's words.

### The densities (`R/priors.R`, `dirichlet_logdens()`)

- Softmax chart (the sampler's): `log p(z) = lgamma(sum a) - sum
  lgamma(a) + sum a_k log w_k` (the Jacobian of `z -> w_2..w_D` is
  `prod w`).
- The fit's stereographic chart: the conformal factor `(2 / (1 +
  |z|^2))^(D-1)`, the orthant map `dw = 2^(D-1) prod |u_k| dS`, and the
  `2^D` sheets give `log p(z) = C + sum (2 a_k - 1) log |u_k| - (D - 1)
  log(1 + |z|^2) + (D - 2) log 2`.

Checks (`dev/moparity-dirichlet-jac.R`, `dirichlet-jac.txt`): against
finite-difference Jacobians at 20 points per D for D = 3, 4, 5 in both
charts, worst relative gap 6.08e-10; the density integrates to 0.999999
(stereographic) and 1.000000 (softmax) over the chart for D = 3.

A concentration below 1/2 is unbounded at a face in the fit's chart (its
factor is `|u_k|^(2 a_k - 1)`), so a MAP fit has no mode and `frm()`
refuses it (`check_map_dirichlet()`); `frm_sample()` takes any.

### MAP (`dev/moparity-dirichlet-map.R`, `dirichlet-map.txt`)

`y ~ mo(inc) + z`, `dirichlet(c(2, 1.5, 3))`:

```text
MAP objective -fn -173.444657110993 | hand log posterior -173.444657110993 | gap -2.84e-14
max |gradient| at the fit: 5.032834e-05
optim of the hand function, best of 10 starts: value -173.444657110987 (fit -173.444657110993, fit - best -5.83e-12), max |param - fit| 2.14e-07
ML weights: 4.83298e-01 5.16702e-01 2.64498e-13
```

### The sampler (`dev/moparity-dir-sample.R`, `dir-sample-cmp.txt`)

`dirichlet(c(5, 1, 1))` against brms 2.23.0 under the same prior, 4
chains of 2000: weak simplex (`y = 0.05 x + noise`) largest |z| 1.42,
strong (`ls ~ mo(income)`) 2.50; 0 divergences on both sides. The
default `dirichlet(1)` steps aside for a simplex with a row of its own
(`mo_simplex_nlp(fit, entries)`), and the row's entry is read in the
softmax chart (`mo_sample_entries()`).

Class `"sbhaz"` (a `cox()` baseline, brms's class of the same name and
default) takes the same density with no coef.

Seen to fail on rellib-r7: `test-dirichlet-prior.R` 6 errors and 1
failure ("Cannot parse prior 'dirichlet(c(5, 1, 1))'", "Unsupported
prior distribution 'dirichlet'", no simo rows in `get_prior()`);
frmtmb.sample `test-mo-simplex-draws.R` 3 errors.

## Item 3: `update(prior = )`

brms 2.23.0 (`dev/moparity-brms-update.R`, `brms-update.txt`): with
stored rows `normal(0, 0.05)` on `b x` and `normal(0, 10)` on
`Intercept`, `update(prior = exponential(2) on sigma)` carries all three
(b_x mean 0.0484); a new `b x` row replaces the stored one (b_x 0.6168);
a row on `coef = nosuch` is dropped without a word; `update(~ . - x,
prior = <sigma>)` drops the `b x` row. brms binds the new rows to the
stored USER rows and drops duplicates on (class, coef, group, resp,
dpar, nlpar).

frmtmb (`update_merge_prior()`, `R/confint.R`): the same merge on
`prior_spec_slot()`; the stored rows are marked droppable, the new rows
`update_strict` and checked as `frm()` checks them, so `nosuch` is an
error (a deliberate divergence: brms drops a row the user just wrote,
which `brm()` itself refuses). `prior = NULL` keeps the stored prior,
as in brms. Seen to fail on rellib-r7 (`moparity-t1base`):
`test-update-prior-drop.R` 8 failures (the stored rows lost).

## Item 4: `cox()`'s baseline simplex

`R/families.R`: the baseline reads `mo_simplex(raw, chart)`, the chart
from the response's `bhaz_chart` addition term (NULL: the fit's chart;
`mo_chart_frame()` sets "softmax" on the sampler's copy).
`R/fit.R`: `mo_search_terms()` lists the baseline with no coefficient,
and `mo_search()`'s step 3 (`cox_face_leave()`) tries one objective
value per weight at a face; if the likelihood rises from it, the model
is refitted from the simplex mixed a tenth of the way to the
barycenter. The first lane build had no step 3 and seed 22 stopped
0.0215 below the maximum at a saddle with three SEs lost; base did not
(its softmax was elsewhere).

The study (`dev/moparity-cox-study.R`, 200 seeds, right censoring,
n = 300; the reference maximum is a polish of the fit by the
multiplicative EM update on the weights with the intercept folded into
their scale, alternated with Newton on the slope; `cox-sum.txt`):

```text
base: seeds 200 (distinct 200) | gap > 1e-6 1, > 1e-4 1, max 0.000186 | codes 0=195 1=5 | warned 6 | weight < 1e-6 at the fit 128, at the polished max 129 | evaluations total 21581 median 109
lane: seeds 200 (distinct 200) | gap > 1e-6 0, > 1e-4 0, max 1.8e-08 | codes 0=200 | warned 0 | weight < 1e-6 at the fit 129, at the polished max 129 | evaluations total 22726 median 108
paired lane - base logLik: > 1e-6 1, < -1e-6 0, max gain 0.000186, max loss -1.47e-09
base seeds with a warning: 40 129 144 162 168 196
base seeds short by more than 1e-6: 162
```

The base's five code-1 fits are the plateau ("singular convergence
(7)", weights 6e-27 to 2e-54); seed 162 lost an SE and stopped 1.86e-4
short. Cost: 1.05 times the evaluations.

The sampler (`dev/moparity-cox-sample.R`, `cox-sample-cmp.txt`), seeds
2 and 5, 4 chains of 2000: base sampled the coordinates with no density
(seed 2: 239 divergences, `sbhaz_raw_2` mean -3.3e11, R-hat 3.07); the
lane samples the softmax under `dirichlet(1)` and reports `sbhaz[k]`,
0 divergences, largest |z| against brms 2.23.0 2.09 over b_x and the
five weights on both seeds.

Seen to fail on rellib-r7: core `test-cox-chart.R` (3 errors on the
first spelling; see the rerun below), frmtmb.sample `test-cox-draws.R`
2 failures and 2 errors (R-hat 2.31, `sbhaz_raw` names, 270 transitions
over the tree depth).

## Item 5: floors for blocks with several scales

Sweep over every structure at log sd -1137.64, each scale alone and all
at once, b off zero, on RTMB's tape (`dev/moparity-floor-sweep.R`,
`floor-sweep-base.txt`, `floor-sweep-lane.txt`): base had 42 `+Inf` or
`NaN` values on the tape, over `us` with 2 and 3 coefficients, `us_t`
with 2, `cs`, `toep`, `ar1`, `hetar1`, `gr(cov =)` and `gr(prec =)` with
two coefficients (and the one-scale cases below); the lane has 7, all
in structures this item does not touch.

`sd_floored(log_sd) = exp(log_sd) + 1e-150` now gives those nine their
standard deviations (`us`, `us_t` d >= 2, `cs`, `toep`, `ar1`,
`hetar1`, `us_sigma(floored = TRUE)` for `gr_cov` and `gr_prec`). The
threshold (`dev/moparity-floor-threshold.R`): `exp(x) + 1e-150 ==
exp(x)` for every log x of -308.450495 or more (x 1.101e-134). The
test asserts, on the tape, `identical()` objective AND gradient against
the same functions with `exp()` at 42 values of every structure's log sds
on [-308.4, 4] and at its start, and that the old density is `NaN` or
`-Inf` there (`test-covstruct-floor.R`, 940 of 940; 81 failures on
rellib-r7).

A whole model's tape, one process, the lane's tape against one built
with `exp()` in place of the floor (`dev/moparity-floor-ident.R`,
`floor-ident.txt`; `y ~ x + (x | g)`, the fixture of
`test-reparam.R`'s correlated case): fn and gr identical at 50 of 50
random points, joint and Laplace; at log sd -1137.64 the two differ
(the control: -3834.218 against -Inf).

A SAMPLER is another matter, and by design. On that fixture the
centered chain of `frm_sample(bf(y ~ x + (x | g)), ...)` differs from
base's from early in warmup, although the fit, the inits and the first
warmup iterations are identical (`dev/moparity-rp-cmp.R`,
`dev/moparity-rp-cmp2.R`). With the floor taken out of the lane's
registry in the same process the chain is base's draw for draw
(`dev/moparity-rp-cmp3.R`, `rp-cmp3.txt`: "lane without the floor vs
base: identical TRUE"), so a warmup proposal went below log sd -308.45,
where the floor is meant to act. That moved `test-reparam.R`'s
"the two routes agree on the correlated model" gate, an sd comparison of
two 500-draw chains, from 0.03 inside its bound (base: 0.476 against
0.506) to 0.03 outside it (0.533). At 2000 iterations both builds pass
by more than half the bound (0.071 and 0.138 against 0.299;
`dev/moparity-rp-iter.R`, `rp-iter-*.txt`), so the test's case now
samples 2000.

The fixtures. `dev/moparity-audit-run.R` wraps `RTMB::MakeADFun()` to
record the smallest theta every objective, gradient and Hessian call
evaluates, over every test file of the eight packages, one process per
file, the gates on (`dev/moparity-audit.sh`, `dev/moparity-audit-sum.R`,
`audit-sum.txt`):

```text
files with a record: 368 | distinct 368
             pkg files with_theta  calls    min_theta files_failing
          frmtmb   226        148 229640 -4318.120009             0
 frmtmb.coupling    11          3  11546   -14.937899             0
      frmtmb.eam    29          2    121    -2.512445             0
   frmtmb.latent    10          2   4137    -2.244591             0
    frmtmb.learn    15          6   2655 -2492.816399             0
      frmtmb.ode    11          4    974    -2.993922             0
   frmtmb.sample    51         51  10075 -1448.009280             0
   frmtmb.spline    15          9  24466  -160.337223             0
files below -308.45: 15
```

Those 15 files were run again recording only the positions the floor
changes, the log sds of the nine structures, found through the blocks
`build_objective()` keeps in its closure (`dev/moparity-audit-run2.R`,
`dev/moparity-audit2.sh`, `dev/moparity-audit2-out/all.txt`): in 14 the
smallest is -64.69 or more (their deep values are other structures'
thetas: `gp()` length scales, one-scale blocks). In one,
`test-predfix.R`, it is -4318.12, on the `(1 + x | g)` block
(`dev/moparity-audit-run3.R`). Its fits, in one process with the floor
and with `exp()` in its place (`dev/moparity-predfix-floor.R`,
`predfix-floor.txt`): the four gaussian random-slope fits at scales 1
to 1e-6 identical (estimates and logLik), the separated bernoulli
`yb ~ z + x + (1 + x | g)` of seed 533 identical by default and under
`autoscale = FALSE`; under `autoscale = TRUE` it was an error ("NA/NaN
gradient evaluation" in the pre-fit) and is now a fit at logLik 3.4e-13
with both sds at 0 and two warnings (the pre-fit's false convergence,
the separation). The test accepts either. So "bitwise unchanged above
the floor" holds at every point every optimizer of every test
evaluated, and one fixture's path went below it, which the floor turned
from an error into a boundary fit.

Not fixed, and why:

- One-scale `dnorm()` blocks (`us` and `us_t` with one coefficient,
  `diag`, `homdiag` and `smooth`) give `NaN` on the tape (`dnorm(b, 0,
  0)`), not `+Inf`: a sampler rejects a `NaN` and does not stick, and
  these are the most common blocks, so their tape was left alone.
- `hsgp`: its spectral sds legitimately fall below 1e-134 at long
  length scales and high frequencies, so a floor there would change
  fits.
- `gr(prec =)` with one coefficient, `car()` and `spde()` are
  parameterized by a PRECISION `exp(-2 theta)`, which overflows rather
  than underflows; no floor of that form is bitwise invisible.

## Item 6: `parnames()`

`dev/moparity-parnames.R` (a fresh process per case, 26 cases: four
load orders by four call routes, a fit and draws; `parnames-base.txt`,
`parnames-lane.txt`). Base: 24 of 26 give one warning; frmtmb's own
generic captured before brms loads and called after (`saved`) gives 0
on a fit, and on draws the call errors. Lane: 25 of 26 give one; the
draws `saved` case errors as on base ("no applicable method for
'parnames' applied to an object of class \"frmtmb_draws\"").

The method now asks the frame that dispatched it (`via_brms_generic()`,
in both packages; `posterior_samples()` used the same idea). The first
spelling asked `isNamespaceLoaded("brms")` and gave 2 warnings when
brms's generic was captured and brms then unloaded (the `unloaded`
case); it now recognizes brms's generic by the namespace the function
was defined in, which a captured function keeps, so brms is never
loaded to find out. Tests: core and frmtmb.sample
`test-generic-collision.R` (the `saved`, `brms_ns`, `brms_first` and
`unloaded` cases), seen to fail on rellib-r7 on the `saved` case (NWARN
0).

Not fixed, found: the draws `saved` case is a dispatch hole, not a
warning count. When brms is loaded before frmtmb.sample, frmtmb.sample's
`S3method(parnames, frmtmb_draws)` is registered on the generic the
`parnames` binding resolves to, brms's, so frmtmb's own generic, if
someone holds it, has no method for draws. Contrived; filed.

## Item 7: exact `gp()` sampling

`R/covstruct.R`: `chol_L` for `gp`, `exp(theta[1]) t(chol(gp_corr()))`
on the tape (RTMB's `chol()`), and `ncp_mode()` admits it as it admits
`hsgp` (its length scales enter the factor; class `lscale` priors
them). `dev/moparity-gp-ncp.R`, lane surface's `y ~ gp(x)` (60 points,
data seed 5) and a `gp(x, by = f)` on the same `x`, both routes, 4
chains of 1000, seed 1; `dev/moparity-gp-cmp.R` (`gp-cmp.txt`) sets
each against brms 2.23.0 on the same data. `|z|` is a mean or sd
difference over the joint Monte Carlo error, over b_Intercept, sigma,
sdgp and lscale:

```text
== base gp formula: divergences 0, at max tree depth 1286, seconds 227 | max R-hat 1.517, min bulk ESS 7 | max |z| mean 1.14 sd 2.99
== base gp fit: divergences 0, at max tree depth 1103, seconds 227 | max R-hat 2.243, min bulk ESS 5 | max |z| mean 0.93 sd 2.00
== base gpby formula: divergences 1, at max tree depth 1958, seconds 278 | max R-hat 1.325, min bulk ESS 11 | max |z| mean 2.35 sd 2.05
== base gpby fit: divergences 0, at max tree depth 1977, seconds 325 | max R-hat 1.475, min bulk ESS 8 | max |z| mean 2.35 sd 1.80
== lane gp formula: divergences 70, at max tree depth 0, seconds 81 | max R-hat 1.016, min bulk ESS 312 | max |z| mean 1.94 sd 1.41
== lane gp fit: divergences 365, at max tree depth 0, seconds 83 | max R-hat 1.210, min bulk ESS 13 | max |z| mean 1.43 sd 2.23
== lane gpby formula: divergences 88, at max tree depth 0, seconds 73 | max R-hat 1.034, min bulk ESS 162 | max |z| mean 2.02 sd 1.33
== lane gpby fit: divergences 39, at max tree depth 0, seconds 77 | max R-hat 1.006, min bulk ESS 375 | max |z| mean 0.89 sd 1.24
== gp brms: divergences 36, seconds 95 | max R-hat 1.015, min bulk ESS 413
== gpby brms: divergences 106, seconds 82 | max R-hat 1.063, min bulk ESS 49
```

The centered base agrees with brms only because its Monte Carlo errors
are wide (bulk ESS 5 to 11). Non-centered, three of four runs mix as
brms's do, with divergences of brms's order (brms itself has 36 and
106), and the time falls by a factor of 2.8 to 4.2. The gp `fit` route
is the exception, explained below.

Found on the way, not fixed: on the fit route the chain started exactly
at the ML mode (chain 1 of `init = "last.par.best"`) sticks (326 of 500
divergent, sigma 0.256 against 0.29), while the jittered chains and
every `init = "random"` chain mix (9 to 38 divergences each;
`dev/moparity-gp-init.R`, log `dev/moparity-log/gp-init.txt`, written at
punch round 1). The ML mode of this gp is at a small sigma the
posterior does not favor; mode initialization is a policy question for
the sampler lane.

## Item 8: `as_tmbstan()`

Kept: `as_tmbstan()` hands over the fit's own tape, by its documented
contract, and brms has no counterpart that argues otherwise. With item
2, a MAP fit with a `simo` row HAS a proper density on that simplex in
the fit's chart (the Dirichlet times the chart's Jacobian integrates to
one, above); the folds still keep a chain on one of the `2^D`
equivalent sheets. The message now names only the simplexes with no
density on the tape (`mo_fit_priored()`), and `?as_tmbstan` says so.

## Tests and suites

New test files: core `test-mo-factor-re.R`, `test-dirichlet-prior.R`,
`test-cox-chart.R`; frmtmb.sample `test-cox-draws.R`. Tests added to
core `test-covstruct-floor.R` (2), `test-update-prior-drop.R` (3),
`test-generic-collision.R` (1); frmtmb.sample `test-mo-simplex-draws.R`
(3), `test-reparam.R` (1), `test-generic-collision.R` (1). Changed
because the behavior changed on purpose: core `test-open-issues.R`
(factor and character multipliers fit; a logical one is `xlTRUE`),
`test-prior-compat.R` (class simo carries over; a class-wide row and a
non-Dirichlet density are refused in brms's words), frmtmb.sample
`test-reparam.R` (gp is eligible; the correlated case samples 2000
iterations, see item 5).

Seen to fail on rellib-r7 (`dev/release/moparity-t1base-files/`,
`-base2-files/`, `-base3-files/`):

```text
test-update-prior-drop.R   pass=10 fail=8 err=0
test-covstruct-floor.R     pass=859 fail=81 err=0
test-generic-collision.R   pass=69 fail=1 err=0   (saved: NWARN 0)
test-mo-factor-re.R        pass=0 fail=0 err=9
test-dirichlet-prior.R     pass=1 fail=1 err=6
test-cox-chart.R           pass=4 fail=4 err=2    (seed 162 gap 1.86e-4 and an SE warning; seed 40 code 1, "singular convergence (7)")
test-open-issues.R         pass=39 fail=0 err=1
test-prior-compat.R        pass=186 fail=0 err=2
frmtmb.sample test-mo-simplex-draws.R   pass=27 fail=0 err=3
frmtmb.sample test-cox-draws.R          pass=0 fail=2 err=2 warn=4 (R-hat 2.31, sbhaz_raw names)
frmtmb.sample test-reparam.R            pass=292 fail=2 err=1
frmtmb.sample test-generic-collision.R  pass=64 fail=3 err=0 (first spelling; see item 6)
```

All eight suites, one file per process, `NOT_CRAN`,
`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS` and `FRMTMB_FUZZ`
on, the lane's library first (every log's `lib:` line names
`wt-moparity-lib`), on the final build (`moparity-full2`,
`dev/moparity-suite-sum.R`, `full2-sum.txt`):

```text
tier moparity-full2 : 368 files, 0 without a RESULT line
             pkg files  pass fail err skip warn
          frmtmb   226 18669    0   0    0    0
 frmtmb.coupling    11   542    0   0    5    0
      frmtmb.eam    29  1746    0   0    3    0
   frmtmb.latent    10   360    0   0    2    0
    frmtmb.learn    15   501    0   0    2    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    51  2725    0   0    1    0
   frmtmb.spline    15   590    0   0    1    0
```

Every skip is an extension's `test-scale.R` (the scale tier,
`FRMTMB_SCALE_TESTS`, not run). The ported brms suite shows no stale
ledger row: no `dev/brmsport-verdicts.tsv` row flips. The first full
run (`moparity-full1`, before the last fixes) had 3 core failures and
1 error, all fixed: duplicate message templates (`stop_mo_one_variable()`,
`stop_mo_crossed()`), `z * mo(inc)` named `z:moinc` where brms says
`moinc:z` (the names now come from `model.matrix()` on the special
terms alone, as brms's do), and `test-prior-compat.R`'s stale
expectation of the simo refusal; and frmtmb.sample's `test-reparam.R`
gate (item 5).

OpenBLAS 0.3.26 as BLAS and LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`, the copy in `dev/ciharden-out/Rob0.3.26-lapack`;
`dev/moparity-run-par.sh` with `RSCRIPT`), the 39 files the change
touches or reaches (`dev/moparity-log/ob-jobs.txt`; `moparity-ob1`) and
the 12 changed after it (`moparity-ob2`), `dev/moparity-suite-sum.R`:

```text
tier moparity-ob1 : 39 files, 0 without a RESULT line
        frmtmb    31 3376    0   0    0    0
 frmtmb.sample     8  587    1   0    0    0   (test-reparam.R's correlated gate, before the 2000-iteration change)
tier moparity-ob2 : 12 files, 0 without a RESULT line
        frmtmb     8  476    0   0    0    0
 frmtmb.sample     4  415    0   0    0    0
```

R CMD check --as-cran (`sh dev/moparity-check.sh <package>`, built and
checked in `dev/moparity-check/`, `R_LIBS` the lane's library, rellib-r7,
the user library): frmtmb **Status: 1 NOTE** (the HTML manual's
skipped math rendering, V8 unavailable); frmtmb.sample **Status: OK**.
After the core tarball was built, NEWS.md and one roxygen paragraph of
`frm()` (the names of `mo()` interaction columns) and the comment above
`check_special_mult()` were reworded; no
code changed, and the installed library has the final source.

## What was not done, and defects found

- `mi()` multipliers stay numeric-only: brms's sp model matrix would
  code a factor there too, but the lane's item is `mo()`.
- A group-level `mo()` term in `mm()`, `gr(cov =)`, `gr(prec =)`,
  `gr(by =)`, `gr(dist = "student")` or a structured covariance is
  refused by name (`assemble_frame()`); each builds its design in a way
  the column scaling has not been built into.
- brms's `variables()` lists `simo_` for a fit; frmtmb's
  `variables(fit)` lists `b_Intercept bsp_mox sigma` on 0.69.0
  (`dev/moparity-variables.R`) and
  here. Not changed (the draws carry `simo_`).
- frmtmb.sample's draws keep a block's raw `theta_k` columns; the sd of
  a block is reached through other surfaces. Not changed.
- The fit route of `frm_sample()` on an exact `gp()`: the chain started
  at the ML mode sticks (item 7).
- `parnames()` on draws through frmtmb's own generic held from before
  brms loaded, with frmtmb.sample loaded after brms: "no applicable
  method" (item 6), on 0.69.0 too.
- The floors left out (item 5): the one-scale `dnorm()` blocks give
  `NaN` on the tape at an underflowing sd (rejected by a sampler, not
  stuck), `hsgp` legitimately reaches spectral sds below 1e-134, and
  the precision blocks (`gr(prec =)` with one coefficient, `car`,
  `spde`) overflow rather than underflow.
- `update(prior = )` refuses a NEW row that matches nothing, where
  brms drops it silently (item 3). A documented divergence.
- A `dirichlet()` on a MAP fit is the density of the chart coordinates,
  as every frmtmb MAP prior is, so `dirichlet(1)` there pulls toward
  the uniform simplex and is not the absence of a prior. Documented in
  `?frm` and NEWS.

Backlog items of "Filed at the 0.69.0 release" this lane closes:
`update(prior =)` merging; the `simo` row refused; `as_tmbstan()` on
`mo()` fits (decided: the contract stays, the message is narrowed);
`mo(x) * f` and group-level `mo()`; cox()'s softmax baseline; a dense
block with several scales; `parnames()`'s generic; exact `gp()`
sampling. Not edited in `dev/test-backlog.md`, which other lanes of
the round edit too; consolidation marks them.

## Files

Core: `R/objective.R` (`mo_chart_frame()`, `build_objective()`; new
`simplex_frame_terms()`, `mo_transform()`, `mo_re_zterm()`),
`R/fit.R` (`fit_assembled()`, the `frm()` roxygen "Monotonic
effects", `frmtmb_control()` roxygen `mo_search`,
`mo_search_terms()`, `mo_search()`; new `cox_face_leave()`),
`R/priors.R` (the `set_prior()` roxygen, `set_prior_one()`,
`frmtmb_prior_classes`, `brms_direct_prior_classes`,
`brms_prior_class_refusal()`, `parse_prior_dist()`,
`prior_dist_string()`, `prior_table()`, `prior_drop_unmatched()`,
`resolve_priorlist()`, `prior_base_logdens()`, `resolve_priors()`;
new `prior_dirichlet()`, `parse_dirichlet_arg()`,
`check_map_dirichlet()`, `dirichlet_dist()`, `dirichlet_logdens()`,
`simo_targets()`, `sbhaz_targets()`), `R/families.R` (`fam_cox()`, the
`cox()` roxygen in the families page), `R/covstruct.R` (`us_sigma()`;
the registry `nll` of `us`, `us_t`, `ar1`, `cs`, `hetar1`, `toep`,
`gr_cov`, `gr_prec`; `ncp_mode()`; new `sd_floored()`, the `gp`
`chol_L`; the non-centering comment block), `R/frame.R`
(`dpar_frame_rhs()`, `assemble_frame()`; new `mo_sp_design()`,
`mo_sp_newdata()`, `mo_ph_unwrite()`, `mo_ph_rename()`,
`frame_has_mo_re()`, `mo_re_index()`, `mo_re_link()`), `R/parse.R`
(`parse_linpred()`; new `stop_mo_one_variable()`, `stop_mo_crossed()`,
`colon_factors()`, `mo_term_entries()`, `mo_bar_rewrite()`),
`R/predict.R` (`pred_design()`, `eval_dpars()`, `lp_eta_design()`,
`lp_delta_A()`; new `lp_Z()`), `R/confint.R` (`update.frmtmb_fit()`
and its roxygen, `parnames.frmtmb_fit()`; new `update_merge_prior()`,
`via_brms_generic()` in place of `parnames_owner_warned()`),
`R/methods-fit.R` (`summary_mo_frame()`), `R/se-check.R`
(`sep_certificate()`), `R/importance.R` (`check_importance_scope()`),
`R/compat.R` (the `mo()` note), `R/sampling-api.R` (exports and
documentation), `NAMESPACE`, `NEWS.md`, `man/` (frm, frmtmb-families,
frmtmb-priors, frmtmb-sampling-api, frmtmb_control, set_prior,
update.frmtmb_fit).

frmtmb.sample: `R/sample.R` (`prior_augmented_obj()`,
`ncp_objective()`, `mo_sample_frame()`, `mo_sample_pars()`,
`mo_simplex_nlp()`, `mo_draws_to_chart()`, `draws_mo_cols()`,
`ncp_reasons`, the default-prior notes, `frm_sample()` (roxygen, the
`mo_soft` test), `check_laplace()`, `as_tmbstan()` and its roxygen;
new `mo_fit_priored()`, `mo_sample_entries()`), `R/methods-draws.R`
(`posterior_samples.frmtmb_draws()`, `parnames.frmtmb_draws()`,
`via_brms_generic()` in place of `ps_via_brms_generic()`), `NEWS.md`,
`man/as_tmbstan.Rd`, `man/frm_sample.Rd`. `R/draws-brms.R` is not
edited: its `draws_natural_cols()` reads `draws_mo_cols()`, which now
lists the `cox()` baseline too.

Tests: as listed above. Scripts: `dev/moparity-*.R`, `dev/moparity-*.sh`
(including `dev/moparity-run-par.sh`, a copy of
`dev/release/run-par.sh` that takes `RSCRIPT` and finds the root from
`dev/`).
