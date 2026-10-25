# Lane bhaz: brms's spelling of the Cox baseline hazard

Base 18b536fc (frmtmb 0.74.0 plus two CI commits), worktree
`frmtmb-wt-bhaz`, branch `wt-bhaz`, library
`C:/Users/adf44/source/r/wt-bhaz-lib` (core and frmtmb.survival from
this worktree) in front of `rellib-r12`. Nothing is committed. Logs
are under `dev/bhaz-log/`, `dev/bhaz-bitwise-out/` and
`dev/release/bhaz-*` (all gitignored); every number below was pasted
from them.

## 1. What brms 2.23.0 does

Scripts: `dev/bhaz-brms-probe.R` (log `dev/bhaz-log/brms-probe.txt`)
and `dev/bhaz-brms-probe2.R` (`brms-probe2.txt`). Data set A: seed 11,
n = 200, `x ~ N(0, 1)`, then a 3-level `g`, then event
times `rexp(exp(-0.5 + 0.7 x))`,
censoring `rexp(0.2)`; times run from 0.0044252 to 14.42168. brms
2.23.0, splines2 0.5.4. No compile: `standata()`, `stancode()`,
`default_prior()`.

What the source says, and the runs confirm:

- `resp_bhaz(gr = NA, df = 5, ...)` (`R/formula-ad.R`). `df` and `...`
  are "flags": brms evaluates them in the formula's environment and
  passes them to `splines2::mSpline()` (`Zbhaz`) and
  `splines2::iSpline()` (`Zcbhaz`) through `bhaz_basis_matrix()`
  (`R/data-response.R`). `intercept` defaults to TRUE there; `degree`
  is splines2's default, 3. So `bhaz(degree =, intercept =)` work in
  brms through `...`, though only `gr` and `df` are documented.
- The basis: M-splines of degree 3 with an intercept, `df` basis
  functions, `df - degree - intercept` internal knots on equally spaced
  quantiles of the response (`bhaz(df = 6)`: knots 0.51273, 1.35622;
  the default `df = 5`: one knot, 0.85633, the median).
- `Boundary.knots` default (`bhaz_basis_matrix()`): `c(max(min(y) -
  d / 50, 0), max(y) + d / 50)` with `d = max(y) - min(y)`; here
  `c(0, 14.710025)`.
- `bhaz(knots =)` replaces the quantile knots, and `df` is then
  ignored without a message: `bhaz(df = 8, knots = c(0.5, 1.1))` gives
  `Kbhaz` 6.
- `bhaz(Boundary.knots = c(0.5, 10))`, inside the data: splines2
  warns "Some 'x' values beyond boundary knots may cause
  ill-conditioned basis functions" twice and extrapolates. 63 of 200
  rows of `Zbhaz` have a negative entry, but no row is negative for
  every simplex (the review's `dev/review-bhaz-bnd.R`, 0 of 63). The
  stronger fact: 62 rows lie below the lower knot and 1 above, and at
  equal weights the cumulative baseline hazard there reaches -1.082, so
  S(t) > 1; for some weights the hazard is negative.
- A knot outside the boundary: "Internal knots must be set inside
  boundary." `bhaz(df = 3)`: "The specified df was too small."
  `bhaz(df = 4.5)`: `Kbhaz` 4 (`as_one_integer()` truncates).
  `bhaz(foo = 1)`: accepted, passed to splines2, ignored.
  `bhaz(df = 6) + bhaz(df = 7)`: "Each addition argument may only be
  defined once."
- A variable in a flag (`bhaz(knots = kn)`) must be in `data` or
  `data2`, or brms stops with "The following variables can neither be
  found in 'data' nor in 'data2'", although brms then evaluates the
  flag in the formula's environment. Literal vectors work.
- `bhaz(gr = g)`: one simplex per level of `g`. `standata()` adds
  `ngrbhaz` 3, `Jgrbhaz` (table 55 69 76) and a 3 x K `con_sbhaz`; the
  Stan code declares `array[ngrbhaz] simplex[Kbhaz] sbhaz` and computes
  `bhaz[n] = Zbhaz[n] * sbhaz[Jgrbhaz[n]]`; `default_prior()` adds one
  `sbhaz` row per group (`a`, `b`, `c`) under the global
  `dirichlet(1)`. The basis itself is shared by the groups. An unnamed
  first argument is `gr`: `bhaz(6)` is one group named "6".
- `bhaz()` with `gaussian()`: "Argument 'bhaz' is not supported for
  family 'gaussian(identity)'."
- `brms::cox()` has one argument, `link`. `brms::cox(bhaz = list(df =
  6))` stops in brms itself with "unused argument (bhaz = list(df =
  6))", before `frm()` sees it. brms 2.22.0's NEWS: "Support stratified
  `cox` models via the new addition term `bhaz`" (#1489); before that
  the basis was a `cox()` argument. So the backlog's second item is
  brms's own error, and brms 2.23.0 has no `cox(bhaz =)` spelling.

### The two bases are the same

frmtmb.survival 0.1.0's `bhaz_spec()`, `mspline_design()` and
`ispline_design()` (rellib-r12) against brms's `standata()`, on the
200 times above:

| basis | max abs diff / max abs, M | same, I | knots equal | boundary equal |
|---|---|---|---|---|
| `cox()` vs brms default | 3.862e-16 | 5.551e-16 | TRUE | TRUE |
| `cox(df = 6)` vs `bhaz(df = 6)` | 2.337e-16 | 6.661e-16 | TRUE | TRUE |
| `cox(df = 5, degree = 2, intercept = FALSE)` vs the same `bhaz()` | 3.750e-16 | 4.441e-16 | TRUE | TRUE |

and, with brms's own knots fed to frmtmb's spec (`brms-probe2.txt`):
`knots = c(0.5, 1.1), Boundary.knots = c(0, 16)` 1.140e-16 / 4.441e-16;
`Boundary.knots = c(0, 16)` 3.862e-16 / 4.441e-16. With boundary knots
inside the data the rows inside agree (8.7e-17 / 7.8e-17 on 137 of 200
rows); outside them frmtmb's recursion gives zero rows where splines2
extrapolates.

The knot equality is exact (`all.equal(tolerance = 0)`), so the bases
differ only in the rounding of two implementations of one recursion.

## 2. What changed

### Core: `frmtmb_register_aterm(option =)` (a new argument of an export)

A registered addition term may now be an OPTION of the family rather
than data: `frmtmb_register_aterm(name, option = function(family,
call, env))`. The parser (`parse_response()`) keeps the term's call as
written, out of `aterms`, and `parse_one_response()` calls
`apply_aterm_options()` once the response's family is known, before
anything reads the family. The function gets the family, the call and
the formula's environment, and returns the family the fit uses. A
family whose `accepts_aterms` does not name the term is refused before
the function runs ("Addition term `bhaz()` is not supported for family
'gaussian'"), which is brms's rule for `bhaz()`. A duplicate is
refused by the existing duplicate check, and a function that returns
no family is refused. `update()`, `predict()` on new data and refits
parse again from the written family, so the function must be a pure
function of (family, call); `cox()`'s is (see below).

Why a new seam and not the old one: a registered term is per-row data,
evaluated against the data and coerced to numbers, with a fixed arity.
`bhaz(df = 6)` is a family argument, `bhaz()` has 0 to 5 arguments,
and `knots` is a vector that is no column. Guardrail 2 forbids
reaching into core, so the extension needed a public seam. It is
useful on its own: any family argument a brms user writes on the
response side has this shape.

New exports: none. `frmtmb_register_aterm()` gains `option = NULL`;
the default keeps every existing registration as it was.

`R/brms-vocabulary.R`: `bhaz = "frmtmb.survival"` in
`extension_aterm_homes` (its own line); `bhaz_family_pointer` and the
`bhaz` entry of `brms_aterm_refusals` are gone, and the table is an
empty named character vector. Without frmtmb.survival, `bhaz()` stops
with "The bhaz() addition term is in the extension package
frmtmb.survival, which is not loaded. Install it and run
library(frmtmb.survival), then use bhaz()".

### frmtmb.survival

- `.onLoad()` registers `bhaz` with `option = cox_bhaz_option`, before
  the compat rules, and `cox()` declares `bhaz` in `accepts_aterms`.
  A compat row `bhaz() x cox: works` says what the term does.
- `cox(link, df, degree, intercept, knots = NULL, Boundary.knots =
  NULL, bhaz = NULL)`. `knots` and `Boundary.knots` are new: brms sets
  them only through `bhaz()`, and the frmtmb spelling of every
  `bhaz()` needs them. `bhaz =` exists only to refuse brms's pre-2.22.0
  spelling by name.
- `cox_build()` validates the written arguments, each labeled with the
  place it was written (`cox(df =)` or `bhaz(df =)`), and records what
  `cox()` itself was given (`fam$cox_args`). `cox_bhaz_option()`
  matches the call against brms's signature `(gr = NA, df, degree,
  intercept, knots, Boundary.knots)`, evaluates each argument in the
  formula's environment, refuses an argument set in both places with
  two values, merges, and rebuilds the family from what `cox()` was
  given plus `bhaz()`. Applied twice it gives the same family.
- `bhaz_spec()` takes `knots` and `boundary`. With neither, the code
  path is the old one (bitwise check below).

### What frmtmb refuses that brms takes, and why

| spelling | brms 2.23.0 | frmtmb |
|---|---|---|
| `bhaz(gr = g)`, `bhaz(g)` | one baseline per level | refused: one baseline; fit each stratum, or add `g` to the formula |
| `bhaz(df = 8, knots = <2 knots>)` | ignores df | refused, names both |
| `bhaz(df = 6)` with `cox(df = 5)` | (no `cox(df =)`) | refused, names both |
| `bhaz(Boundary.knots =)` inside the data | extrapolates; cumulative hazard down to -1.08 | refused |
| `bhaz(foo = 1)` | ignored | refused |
| `bhaz(df = 4.5)` | truncated to 4 | refused (`sv_check_count()`) |
| `bhaz(int = FALSE)`, `inter =`, `i =`, `Bound =`, `B =` | dropped without a message (fix round 1) | refused, names the argument meant |
| `bhaz(periodic =)`, `per =`, `derivs =`, `integral =` | change the basis | refused |
| `bhaz(NA, 6, <knots>)` | third unnamed argument ignored | refused |
| `bhaz(6)` | one stratum named 6, df 5 | refused, suggests `bhaz(df = 6)` |
| `cox(bhaz = list(df = 6))` | "unused argument" | refused with the two spellings |

What frmtmb takes that brms does not: `bhaz(knots = kn)` with `kn` in
the formula's environment and not in `data2`.

Also changed: the refusal of a `df` below `degree + intercept` says
"degree (3) without an intercept" when there is none; the old message
said "degree + 1 ... with an intercept" there too. Where it is raised
was changed in the first version of this lane and changed back in fix
round 1 (section 9).

## 3. A brms model and its frmtmb spelling give the same fit

Script `dev/bhaz-brms-fit.R`, log `dev/bhaz-log/brms-fit.txt`.
Data set B: seed 11, n = 200, with no `g` drawn before the times, so
not data set A of section 1: times 0.00044 to 8.754, 146 events. brms
2.23.0 / rstan 2.32.7, 2 chains x 1000 iterations, seed 1,
`R_MAKEVARS_USER` set; each model compiled in
40-42 s. For each model: frmtmb's basis from its own frame
(`dry_run = "frame"`) against `standata()`; the two frmtmb spellings
against each other; and, at each of brms's 1000 posterior draws, brms's
`log_lik()` summed over rows against frmtmb's `cox()` `lpdf` (events)
and `lccdf` (right-censored rows) on frmtmb's own basis, with the
draw's `sbhaz` mapped to the softmax chart.

| model | basis M / I | frmtmb spellings | max rel. diff over 1000 draws | draws above frmtmb's MLE |
|---|---|---|---|---|
| `bhaz(df = 6)` | 3.866e-16 / 5.551e-16 | identical logLik, estimates, vcov | 1.494e-16 | 0 |
| `bhaz(knots = c(0.5, 1.1), Boundary.knots = c(0, 16))` | 1.113e-16 / 4.441e-16 | identical | 1.493e-16 | 0 |
| `bhaz(df = 5, degree = 2, intercept = FALSE)` | 3.102e-16 / 4.441e-16 | identical | 1.406e-16 | 0 |

frmtmb's MLE log likelihoods: -189.9092473950, -189.9111125281,
-201.7588985333; the highest brms draw in each: -190.1157643592,
-190.2726843239, -201.9319516707. `b_x`: frmtmb 0.7564 (SE 0.1044),
0.7565 (0.1044), 0.7653 (0.1045); brms posterior mean 0.7353 (sd
0.1121), 0.7186 (0.1073), 0.7451 (0.1025). The posterior means sit
below the MLE by about a fifth of a standard error, which the
`dirichlet(1)` prior and the short chains allow (Rhat max 1.017, 1.01,
1.054; the third run warned of low ESS). The log likelihood is the
same function in both packages to 1.5e-16 relative; that is the
claim, and it does not depend on the sampler.

## 4. Stratified baseline: not built

`bhaz(gr = )` is refused. It is not cheap on the existing seams:

1. The group is per-row DATA. An option term carries none, and a
   family sees data only through `aterms`; `gr` would need a second
   registered data term or an option term that also adds a column.
2. The baseline would be one simplex per group. `simplex_pars` declares
   one simplex per extra parameter, read with one `mo_simplex()`, and
   core's prior code refuses `class = "sbhaz"` with a `group`
   (`R/priors.R`, "class = \"sbhaz\" is the whole baseline simplex").
   brms names the draws `sbhaz[<group>, k]`, sets priors with `group =`,
   and frmtmb.sample's draw names and `dirichlet(1)` default, the
   `mo_search()` face search, and `cox_baseline()` all assume one
   simplex. One extra parameter per level (`sbhaz_a`, ...) would fit,
   but its prior class and draw names would not be brms's.
3. The number of simplexes is a property of the data, which
   `extra_pars(y, aterms)` could read only after (1).

So it touches the seam's data path, the simplex seam, priors and
frmtmb.sample. The refusal says what the model is and the two ways to
get close to it.

## 5. Unchanged fits, bitwise

Script `dev/bhaz-bitwise.R` (arm `lane` or `base`), compared by
`dev/bhaz-bitwise-cmp.R`; output `dev/bhaz-bitwise-out/`. Seed 33,
n = 400, Weibull times with a 20-level frailty. Nine fits: `cox()`,
`cox(df = 6)`, `cox(df = 9)`, `cox(df = 5, degree = 2, intercept =
FALSE)`, `cox(link = "identity")`, `(1 | g)` frailty, `trunc(lb =)`,
interval censoring, `family = "cox"`. Three processes per arm, every
lane run against every base run: 81 comparisons, 81 identical in
logLik, estimates, vcov and `cox_baseline()`, 0 differ.

## 6. Tests

### New and changed test files

- `extensions/frmtmb.survival/tests/testthat/test-cox-bhaz.R` (new, 9
  blocks, 93 assertions): each `bhaz()` spelling against its `cox()`
  spelling (`identical()` logLik, estimates, vcov); the absent cases
  (`bhaz()`, `bhaz(gr = NA)`, `bhaz(df = 5)` change nothing; a `cox()`
  argument `bhaz()` does not name is kept; the same value in both
  places is taken; boundary knots that cover the data, and exactly at
  its range, are taken); the formula environment; `update()` and new
  data; every refusal by name; `cox(bhaz = )`; the basis against brms
  2.23.0's `standata()` (16 ulps of the basis's own maximum; measured
  5.6e-16 relative, about 2.5 ulps).
- `tests/testthat/test-aterm-option.R` (new, 3 blocks, 14 assertions):
  the core seam with a test family (`tshifted`, `tshift()`), under
  `local_frmtmb_registries()`: the fit of the returned family, the call
  and environment the function sees, the term absent from `aterms` and
  from new data, the allow-list refusal before the function runs (call
  count 0, then 1 for a family with no allow-list), a function that
  returns no family, a duplicate, argument validation, and a data term
  beside an option term.
- `tests/testthat/test-guardrail1.R`: the two `bhaz` probes (gaussian,
  which brms refuses) became five with `family = 'cox'`, two of them
  refused by name once frmtmb.survival is loaded (`bhaz(df = 6)` on
  gaussian, `bhaz(gr = g)`); the block "the bhaz() refusal names what
  has no equivalent" is gone, and its assertions became three in "an
  unloaded extension's addition term is refused by name" (the home
  table, the refusal, and bhaz absent from `brms_aterm_refusals`). The
  file has 35 assertions before (HEAD's file on rellib-r12) and after
  (this file on the lane build). HEAD's file fails on the lane build,
  as it should (the extension arm now resolves bhaz).

### Seen to fail on rellib-r12 (`dev/bhaz-run-test.R ... base`)

- `test-cox-bhaz.R`: error in all 9 blocks. The behavioral form:
  "Addition term `bhaz()` is not supported: it sets brms's M-spline
  basis ..." from `check_unsupported_aterm()` (8 blocks), and "unused
  arguments (knots = c(0.5, 1.1), Boundary.knots = c(0, 16))" from
  `cox()` (`dev/bhaz-log/t-cox-bhaz-base.txt`).
- `test-guardrail1.R`: pass 29, fail 2, error 1. The extension arm
  with frmtmb.survival loaded refused all three `cox` probes of bhaz
  ("frmtmb.survival: bhaz: ERR Addition term `bhaz()` is not supported
  ..."), and `refused` held "custom", "sar", "fcor" against the
  declared "custom", "bhaz", "sar", "fcor" (`t-guardrail1-base.txt`).
- `test-aterm-option.R`: error in all 3 blocks, "unused argument
  (option = ...)", the weak (missing-argument) form; the behavioral
  form of the seam is the two files above.

### Runs on the final build

All with `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`, the lane library
(core and frmtmb.survival) in front of rellib-r12 (the other seven
extensions); every `lib:` line names the lane's frmtmb.

One file per process, all nine suites, 8 at a time
(`dev/release/run-par.sh bhaz-final`, log
`dev/release/bhaz-final.log`): 443 of 443 files ran, fail 0, error 0.

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 275 | 21368 | 0 | 0 | 15 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 32 | 1844 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 12 | 679 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 501 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 65 | 3294 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 7 | 391 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 15 | 498 | 0 | 0 | 0 | 0 |

The skips are other tiers' gates (`FRMTMB_DRMTMB_FIT_TESTS`,
`FRMTMB_FUZZ`, the scale tier). Named files: core `test-guardrail1.R`
35, `test-cox-home.R` 8, `test-compat.R` 763, `test-compat-register.R`
103, `test-compat-extensions.R` 2,
`test-aterm-option.R` 14, `test-parse.R` 44, `test-registry-helper.R`
7, `test-message-uniqueness.R` 6; frmtmb.survival `test-cox-bhaz.R` 93,
`test-cox-home.R` 12, `test-cox.R` 45, `test-saved-fits.R` 56;
frmtmb.sample `test-cox-draws.R` 12 and `test-compat-preflight.R`
86; all fail 0.

frmtmb.survival before (its 14 unchanged files on rellib-r12, tier
`bhaz-base-sv`): 405 pass, 0 fail. After: 405 + 93 = 498. Core: one
file added (14), `test-guardrail1.R` 35 before and after; the full
core suite was not run on rellib-r12.

A first full run (tier `bhaz-gated`) found one failure of this lane's,
fixed before the final run: core's `test-message-uniqueness.R` (fail
2) saw a second "Duplicated addition term" template, which the option
branch of `parse_response()` had copied; the option term now goes
through the existing duplicate check.

One session per suite, three orders (`dev/release/run-sessions.sh
bhaz-sess 8`, log `dev/release/bhaz-sess.log`): 27 of 27 sessions,
fail 0 and error 0 in every one. 157 LEAK lines, every one tagged
"[with load of ...]"; 0 untagged. Core: forward 21368, reverse 21518,
preload 21612 pass (the counts differ by order because some core tests
assert more with extensions loaded); every extension's three sessions
have the same count as its per-file run.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`): `test-cox-bhaz.R` 93, `test-cox.R` 45, survival
`test-cox-home.R` 12, `test-aterm-option.R` 14, `test-guardrail1.R` 35,
all fail 0 (`dev/bhaz-log/ob-*.txt`). The OpenBLAS copy of R was
deleted afterwards.

### R CMD check --as-cran

Built and checked in `dev/bhaz-check/<pkg>/` (`dev/bhaz-check.sh`),
`NOT_CRAN=true`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`.

- frmtmb.survival: `Status: 1 NOTE` (V8, math rendering of the HTML
  manual). Tests: `[ FAIL 0 | WARN 0 | SKIP 4 | PASS 459 ]`, which is
  0.74.0's PASS 366 and SKIP 4 plus the 93 new assertions.
- frmtmb: `Status: 1 NOTE` (V8, the same note). Tests:
  `[ FAIL 0 | WARN 0 | SKIP 228 | PASS 18851 ]` in 23 minutes; the
  skips are the gated tiers, which the check environment does not set.

## 7. What I did not do, and open items

- No stratified baseline (section 4).
- No `cox(bhaz =)` acceptance: brms 2.23.0 has no such argument; it
  is refused with the two spellings.
- `bhaz()` evaluates its arguments in the formula's environment and
  does not read `data2`. brms checks that a flag's variables are in
  `data` or `data2`, then evaluates the flag in the formula's
  environment, so a variable that is only in `data2` fails in brms
  too. Not pursued.
- The refusal of a `df` below `degree + intercept` in `bhaz_spec()` is
  now reachable only from a direct internal call (tests do that); it
  stays as a guard.
- Ported brms ledger (`dev/brmsport-verdicts.tsv`): no row mentions
  `cox` or `bhaz` (grep of every `dev/brmsport-*` file), so no row
  flips. brms's own blocks "standata handles cox models correctly"
  (`tests.standata.R`) and the cox block of `tests.stancode.R` are not
  in the ported ledger.
- `dev/test-backlog.md` ("brms spellings of `cox()`") is left for the
  consolidation: the first item is fixed; the second is brms's own
  error (section 1), and `cox(bhaz = )` is now refused by name.
- Overlaps with other lanes: `R/brms-vocabulary.R` (one new line in
  `extension_aterm_homes`, the bhaz refusal block removed), `R/parse.R`
  (`frmtmb_register_aterm()`, `parse_response()`,
  `parse_one_response()`, new `apply_aterm_options()`), `NEWS.md`,
  `man/frmtmb_register_aterm.Rd`, `tests/testthat/test-guardrail1.R`,
  `dev/core-boundary.md` (the sentence on the tables of refusals, which
  the sarfcor lane also edits, and the seams table), `R/compat.R`
  (fix round 1), `dev/upstream-bugs.md` (fix round 1, brms-36),
  frmtmb.survival (`R/cox-family.R`, `R/zzz.R`, package doc, NEWS,
  man). I did not touch frmtmb.survival's `test-saved-fits.R`.

## 8. Release needs

- Core: a minor bump (a new argument of an exported seam, and the
  `bhaz()` refusal text). frmtmb.survival: a minor bump (new `cox()`
  arguments and the `bhaz()` term), with its floor on core raised to
  the core release that has `frmtmb_register_aterm(option =)`:
  frmtmb.survival's `.onLoad()` calls it, so it does not load with
  frmtmb 0.74.0. Measured: this frmtmb.survival in front of
  rellib-r12's frmtmb stops with ".onLoad failed in loadNamespace() for
  'frmtmb.survival' ... unused argument (option = cox_bhaz_option)".
  `Depends: frmtmb (>= 0.74.0)` in its DESCRIPTION is left for the
  release to raise (no version numbers in a lane). The other
  extensions need no change and no new floor.

## 9. Fix round 1 (review `dev/reviews/2026-10-24-bhaz.md`, FIX)

1. **Partial argument names (major).** `cox_bhaz_option()` no longer
   uses `match.call()`. `cox_bhaz_args()` names the arguments as brms
   2.23.0 does, measured on data set A by `dev/bhaz-brms-partial.R`
   (log `dev/bhaz-log/brms-partial.txt`):
   - taken as brms takes them: `d = 6` (df, K 6), `g = g` (gr),
     `kn =`, `k =` (knots, K 6), `deg = 2` (degree 2), `bhaz(NA, 6)`
     (df 6), `warn.outside =`, `warn =`, `w =` (no change to the basis
     in brms, so taken and ignored);
   - refused with the argument meant, where brms drops the name
     without a message: `int = FALSE`, `inter = FALSE`, `i = FALSE`
     (brms keeps intercept TRUE; `int` and `i` reach splines2's
     `integral`, and `int = TRUE` changes the basis), `Boundary =`,
     `Bound =`, `B =` (brms keeps the default boundary 0, 14.71);
   - refused: `periodic`, `per`, `derivs`, `integral` (they change
     brms's basis; `per = TRUE` moved the knots to 0.3277 0.6411 1.0638
     2.0981), `de = 2` and `df = 6, d = 7` (brms: "argument 2 matches
     multiple formal arguments"), `Knots =` and `foo =` (dropped by
     brms), a third unnamed argument (`bhaz(NA, 6, c(0.5, 1.1))`: brms
     ignores it and fits K 6 on quantile knots).
   The defect was seen on the lane's first build in the review's
   `dev/review-bhaz-partial.R` (`int = FALSE` fit intercept FALSE, max
   abs dZbhaz 5.17; `Bound = c(0, 30)` fit boundary 0, 30, 0.457), and
   every spelling above is pinned in `test-cox-bhaz.R`, block "bhaz()
   reads argument names as brms 2.23.0 reads them".
2. **Wording.** The NEWS bullet, the `bhaz_spec()` comment, `?cox` and
   the compat note now say the cumulative baseline hazard can be
   negative (S(t) > 1), with the review's numbers (62 rows below the
   lower knot and 1 above, -1.082 at equal weights, no row negative for
   every simplex). `warn.outside` is taken (item 1). `bhaz(6)` says the
   first argument is `gr`, that brms reads one stratum named 6 with df
   5, and suggests `bhaz(df = 6)`.
3. **`?frmtmb_register_aterm` example.** Replaced with `knownsd()`: a
   family that declares the term, and an option that builds the family
   again. Outside `\dontrun{}`, so R CMD check runs it. Run with
   `example()` (`dev/bhaz-example.R`, log `dev/bhaz-log/example.txt`):
   logLik -103.5655 against the least-squares reference -103.5655. The
   section now says why the family must be rebuilt (the density is a
   closure over the constructor's arguments).
4. **frm_compat rows.** frmtmb.survival adds `bhaz() x royston_parmar:
   refused`, as `frm()` does. The generic rows come from four core rules
   shared by every addition term (`kind:aterm` against `kind:aterm`,
   `kind:covstruct`, `kind:special`, and `mvbf x kind:aterm`); their
   notes now say what a registered data term (`dec()`) and an option
   term (`bhaz()`) each do. Statuses are unchanged, so no row of any
   extension changes status. `test-aterm-option.R` backs the "works"
   of an option term with `weights()` and in `mvbf()` (identical
   logLik against the shifted family), and checks the two notes.
5. **The df check.** `cox()` no longer checks `df` against
   `degree + intercept`; the merged basis is checked in
   `cox_bhaz_option()` (labeled `bhaz()`), and `cox()` alone is checked
   at frame assembly, before the fit, as at 0.74.0 (`cox_stop_df()`,
   one template). `cox(df = 3) + bhaz(intercept = FALSE)` fits and is
   identical to `cox(df = 3, intercept = FALSE)`; `bhaz(degree = 2)`
   with `cox(df = 3)` likewise. The BREAKING note is gone: the only
   change left is the message's wording without an intercept.
6. Floor: left to the consolidation.
7. **Survival NEWS.** The refusals moved to "Where `bhaz()` differs
   from brms" under "New features"; `cox(bhaz =)` and the df message
   are "Other changes". The basis figure is "6.7e-16 of their size or
   less", the largest of the measured maxima (6.661e-16, section 1;
   the review's 6.7e-16).
8. **Findings slips.** `test-compat-preflight.R` is listed under
   frmtmb.sample; sections 1 and 3 name data sets A and B;
   `dev/core-boundary.md` and `R/compat.R` are in the overlap list.
9. **`dev/core-boundary.md`.** The history sentence is cut to
   "`bhaz()` is in `extension_aterm_homes`; frmtmb.survival refuses
   `bhaz(gr = )` by name."
10. **Not changed, by instruction:** without frmtmb.survival, `bhaz()`
    on any family says "install it and run library(frmtmb.survival),
    then use bhaz()", and with it loaded the same call on gaussian is
    refused. `dec()` has the same shape with frmtmb.eam.

Outside the lane: brms 2.23.0's interval-censored `cox()` scores the
interval with the cumulative hazard at `Y` for both ends, so the
interval's log probability is -Inf. Added to `dev/upstream-bugs.md` as
brms-36 with `dev/bhaz-brms-icens.R` (log `dev/bhaz-log/brms-icens.txt`:
row 2, Y 1.2816, rcens 1.9225, both ends -1.019042). Not posted. The
number 36 may collide with another lane's new entry at consolidation.

### Runs after the fixes

One file per process (`dev/release/run-par.sh bhaz-fix1`, log
`dev/release/bhaz-fix1.log`), 24 of 24 files, fail 0, error 0, warn 0:
core `test-aterm-option.R` 19, `test-guardrail1.R` 35, `test-compat.R`
763, `test-compat-register.R` 103, `test-compat-extensions.R` 2,
`test-message-uniqueness.R` 6; frmtmb.ode `test-compat.R` 15;
frmtmb.sample `test-cox-draws.R` 12, `test-compat-preflight.R` 86;
frmtmb.survival's 15 files 550 (405 unchanged + `test-cox-bhaz.R` 145).

One session, three orders (`run-sessions.sh bhaz-fix1-sess 6 frmtmb
frmtmb.survival`): frmtmb forward 21373, reverse 21523, preload 21617;
frmtmb.survival 550 in all three; fail 0 and error 0 in all six. 110
LEAK lines, 0 untagged. Core's counts are the first round's plus the 5
new assertions of `test-aterm-option.R`.

R CMD check --as-cran (`dev/bhaz-check.sh`): frmtmb `Status: 1 NOTE`,
frmtmb.survival `Status: 1 NOTE`, both the V8 note on the HTML manual;
examples OK (the new example runs), tests OK (core 24 min, survival
27 s). I deleted the `.Rcheck` directories before reading the testthat
summary lines, so this round's check has no PASS counts; the
one-file-per-process and session counts above are the record.

## 10. Fix round 2 (two nits after SHIP)

1. A full name and its shortened form together, `bhaz(degree = 2,
   deg = 3)` and `bhaz(knots = c(0.5, 1.1), kn = c(1, 2))`: brms keeps
   the full name and drops the short one (`dev/bhaz-brms-partial.R`,
   appended cases: degree 2; knots 0.5 1.1), while `bhaz(df = 6, df =
   7)` stops in brms ("formal argument "df" matched by multiple actual
   arguments"). Still refused; the message now says "brms keeps
   `degree =` and drops `deg =` without a message" for the first kind
   and "brms stops on this too" only for an exact name given twice.
   Both kinds are pinned in `test-cox-bhaz.R`.
2. The `knownsd()` example of `?frmtmb_register_aterm` registers inside
   a function that calls `frmtmb_local_registries()`, under
   `requireNamespace("withr")`. `dev/bhaz-example.R` (log
   `dev/bhaz-log/example.txt`): logLik -103.5655 against -103.5655;
   the aterm registry is `identical()` before and after the example,
   `knownsd` is not registered after it, and the compat features are
   identical.

Reruns on the reinstalled build: `test-cox-bhaz.R` 147 pass,
`test-aterm-option.R` 19, frmtmb.survival `test-message-uniqueness.R`
4; fail 0, error 0, warn 0. R CMD check not rerun, by instruction.
