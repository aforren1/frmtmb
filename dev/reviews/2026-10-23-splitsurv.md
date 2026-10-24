# Review: lane splitsurv (curve inference into core, frmtmb.survival, frmtmb.spline retired)

Reviewer, 2026-10-23. Worktree `frmtmb-wt-splitsurv`, base 19b162b9
(0.73.1), uncommitted. Lane library `wt-splitsurv-lib`, base library
`rellib-r11` (frmtmb 0.73.0, frmtmb.spline 0.10.2). Every number below
was produced by a command I ran; scripts are `dev/review-splitsurv-*.R`
and `.sh`, logs in `dev/review-splitsurv-out/` (gitignored). Nothing in
the lane's files was edited.

## Verdict

| move | verdict |
|---|---|
| (a) `frm_curve()`, `frm_curve_deriv()`, `frm_curve_feature()` into core, `post$curve_check` hook | SHIP |
| (b) `royston_parmar()`, `rp_floored()` into frmtmb.survival | SHIP (nits 6, 10) |
| (c) `cox()`, `cox_baseline()` into frmtmb.survival | FIX (findings 1, 2, 5). The guardrail 3 numbers are right, and by `dev/core-boundary.md`'s own rule they make (c) a weak candidate; see "Guardrail 3" |
| seams `simplex_pars =`, `aterm_data =` | FIX (findings 3, 4) |
| retirement of frmtmb.spline | SHIP (nit 10) |

No blocker and no major finding. Every equality claim reproduces
bitwise, every test file passes on the lane build, guardrail 1 holds
for every live spelling, and frmtmb.survival reads only exports.

## Findings

### 1. minor: the 21 lost assertions of `test-brms-families.R` are misdescribed, and 4 are not carried

The record says the 21 are "4 links built, 17 roster links refused"
and that `test-cox-links.R` "runs the same loops". The roster has 17
links, of which cox() accepts 4, so the refusal loop lost 13, not 17.
The other 4 come from a third test that the record does not name.

```
$ paste <(cut -f3,5 dev/splitsurv-before-out/frmtmb--test-brms-families.tsv) \
        <(cut -f3,5 dev/splitsurv-after-out/frmtmb--test-brms-families.tsv) | awk -F'\t' '$2!=$4'
"every (family, link) pair brms accepts constructs"               278  ...  274
"every roster link brms refuses for a family is refused by name"  574  ...  561
"link and link_<dpar> agree with links on every registry family"  268  ...  264
$ Rscript dev/review-splitsurv-roster.R
17
refused by core loop: 13  by survival loop: 13
```

`test-cox-links.R` reproduces the first two (4 built, 13 refused) and
the live brms link set. It does not reproduce the 4 assertions of
"link and link_<dpar> agree with links on every registry family" for
cox(): `$link` is one string, no `link*` field is stored on the list,
`$link` equals `links$mu$name`, and `$linkinv` equals
`links$mu$linkinv`.

Fix: correct the sentence in `dev/splitsurv-findings.md` (4 + 13 + 4),
and add those 4 assertions for `cox()` to `test-cox-links.R`, or say
why they are dropped.

### 2. minor: a saved cox() or royston_parmar() fit, reloaded where frmtmb.survival is not installed, fails with "could not find function"

I saved an RP fit, a non-monotone RP fit and a `cox(df = 6)` fit with
frmtmb.survival loaded. Then I read them in fresh sessions in two
cases: frmtmb.survival installed but not attached, and not installed.
For "not installed", the lane's frmtmb was copied alone into a scratch
library ahead of rellib-r11, which has no frmtmb.survival. I deleted
the copy afterwards.

```
$ Rscript dev/review-splitsurv-saved.R ns      # installed, not attached
survival ns loaded after readRDS: TRUE ; attached: FALSE
frm_curve(rp good)             OK rows=10 first=1.697905204
frm_curve(rp bad) [gate]       ERROR: rp_floored(): this fit's reported likelihood is not the model's. ...
update(rp, . ~ 1)              FIT logLik=-775.000388314777
$ Rscript dev/review-splitsurv-saved.R absent  # not installed
survival ns loaded after readRDS: FALSE ; attached: FALSE
frm_curve(rp good)             ERROR: could not find function "rp_floored"
frm_curve(rp bad) [gate]       ERROR: could not find function "rp_floored"
update(rp, . ~ 1)              ERROR: could not find function "sv_rp_valid_y"
summary(cox)                   ERROR: could not find function "cox_cbhaz_design"
```

When the package is installed, `readRDS()` loads its namespace, and the
hook and the refusal work. That arm is correct. When the package is not
installed, R replaces the namespace with the global environment and
every closure fails on its first internal call. For RP this is not new,
because frmtmb.spline behaved the same way. For cox() it is new: on the
base build a cox() fit needed nothing but core. Guardrail 1 is about
spellings, not saved objects, so this does not break the letter of the
rule. But the user meets an internal name and is not told which
package to install.

Fix: core knows the package that built each family
(`frm_family_package()`). In `curve_family_check()`, and in one shared
place that summary(), update() and the other methods pass through,
refuse with the wiener-style message if that package is neither
`"frmtmb"` nor loadable (`requireNamespace()`): "This fit's cox family
is from frmtmb.survival, which is not installed". Add a test for it.

### 3. minor: frmtmb.sample still names cox() by kind, so a third-party declared simplex gets an unannounced default prior

`extensions/frmtmb.sample/R/sample.R:1977-1983` announces the
`dirichlet(1)` default only when `identical(t$kind, "cox")`. A
third-party family that declares a simplex gets the same default on
the sampler (correct), but the announcement does not list it. I used
the family of `test-simplex-pars.R`; the true weights are
(0.6, 0.4, 0):

```
$ Rscript dev/review-splitsurv-tsample.R
frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
  Intercept          student_t(3, -1.6, 2.5)
  b                  (flat), as brms leaves slopes
variables: b_Intercept tw[1] tw[2] tw[3] lp__
draw means: 0.590505 0.404846 0.004649
```

The third weight is at 0 in the ML fit (3.3e-16). On the draws it is
0.0046, which is the `dirichlet(1)` default at work, and the
announcement does not name it. `set_prior("dirichlet(20)", class =
"tw")` reaches the sampler (draw means 0.306, 0.467, 0.226). The
record says sample's "R code is unchanged ... whose `kind` is still
`"cox"`". It does not say that this is the one remaining special case
of cox() by name in the ecosystem. NEWS says a declared simplex is
"drawn by frmtmb.sample ... under brms's dirichlet(1)". That is true,
but it is silent.

Fix: in `sample.R`, announce every non-`mo()` entry of
`simplex_frame_terms()` that has no user row, by its brms stem
(`sub("\\[1\\]$", "", t$names[1])`), as the `simo` note does. Then
`kind` is no longer read. Add the third-party case to a sample test.
That is a patch to frmtmb.sample, which this lane already touches.

`simplex_frame_terms()` in core reports `kind = <family name>`, for
example `"tsimplex"` in `test-simplex-pars.R`. That is generic and is
no remaining special case in core. In core `R/`, a grep for
`cox|royston|\brp_|sbhaz|bhaz|\bsp_[a-z]` finds only doc text, the
`extension_family_homes` entry, and `"sbhaz"` as brms's prior class
(set_prior-time checks, `frmtmb_prior_classes`, and the one allowed
name in `check_simplex_pars()`). That class name is brms's, so keeping
it is right.

### 4. minor: `simplex_pars =` is not checked against the extra parameters it names

`check_simplex_pars()` cannot see `extra_pars`'s output, because it is
a function, and nothing checks it at the frame either. A value that
names no extra parameter gives:

```
$ Rscript dev/review-splitsurv-typo.R
fit, correct declaration     FIT logLik -582.8445619
fit, typo in declaration     ERROR: not enough degrees of freedom to define contrasts
get_prior, typo              TABLE classes: Intercept,tw
prior on typo class          ERROR: not enough degrees of freedom to define contrasts
```

`family_simplex_terms()` computes `D = length(par_template[[zeta]]) + 1`,
which is 1 for a missing component. `mo_search_terms()` keeps it,
because `length(z) == D - 1 == 0`, and `get_prior()` lists a prior
class for a simplex that does not exist. This seam is meant for third
parties, so a typo must be refused by name.

Fix: where frame assembly evaluates `extra_pars(y, aterms)`
(`R/frame.R` near line 2984), refuse when a value of
`simplex_pars` is not one of its names, or when it has fewer than 1
coordinate. Name the family and the simplex. Add the guard-absent case
(the correct declaration) beside it in `test-simplex-pars.R`.

### 5. minor: the floor `frmtmb (>= 0.73.1)` is satisfied by a core that cannot run frmtmb.survival

The record calls it a placeholder. But 0.73.1 is the released base
(19b162b9), and its `frmtmb_family()` has neither `aterm_data =` nor
`simplex_pars =` (see `git diff 19b162b9 -- R/families.R`, signature
hunk). An install over 0.73.1 satisfies the Depends and then fails at
load, as the lane measured on 0.73.0. The next release's install script
(`dev/rel073-install.sh` lists `spline`) and its floor script also need
frmtmb.survival. The record names the floor script only.

Fix: write the real floor at consolidation, or `0.73.1.9000` now, so
that the placeholder cannot be satisfied by the released core. Add
frmtmb.survival in place of frmtmb.spline to the next release's install
and floor scripts.

### 6. nit: a compat note that the equality claim kept verbatim is now false

`extensions/frmtmb.survival/R/zzz.R`, row `royston_parmar x fitted`:
"core::cox() refuses for the same reason." cox() is no longer in core.
The lane changed frmtmb.learn's equivalent note but kept this one, so
that "note differs 0" holds for `royston_parmar`. Fix the text and
record one differing note.

### 7. nit (pre-existing, not this lane's): brms's `bhaz()` spellings for cox() have only generic refusals

```
$ Rscript dev/review-splitsurv-bhaz.R base|lane
Addition term `bhaz()` is not supported (supported: weights(), ...)   # both builds
$ Rscript dev/review-splitsurv-g1.R surv
frm(brms::cox(bhaz=df6))           ERROR: unused argument (bhaz = list(df = 6))   # all four arms
```

The behavior is the same on base and lane. Now that cox() has a home,
frmtmb.survival could register a by-name refusal or an equivalent for
`bhaz(df =, knots =, Boundary.knots =)`, with `gr =` refused, because
frmtmb has no stratified baseline. Add it to `dev/test-backlog.md`; it
does not block.

### 8. nit: the `dirichlet()` refusal moved later, which is no regression against brms

```
$ Rscript dev/review-splitsurv-prior.R base
frmtmb set_prior(class=sigma)            ERROR: dirichlet() is the density of a simplex; ...
$ Rscript dev/review-splitsurv-prior.R lane
frmtmb set_prior(class=sigma)            no error
frmtmb validate_prior(p1)                ERROR: dirichlet() is the density of a simplex; ... or a simplex a family declares, not class = "sigma"
frmtmb get_prior-free frm(p1)            ERROR: (same)
frmtmb set_prior(class=Intercept,dpar=sigma) ERROR: ... not class = "Intercept"
brms set_prior(class=sigma)              no error
brms validate_prior(b1)                  no error
brms stancode(b1)                        no error
brms set_prior(class=b)                  no error
brms stancode(b2)                        no error
```

brms 2.23.0 refuses neither spelling at `set_prior()`,
`validate_prior()` or `stancode()`. The lane still refuses at the first
point where the model is known, and `validate_prior()` refuses too, so
the timing is no worse than brms's. Two things are worth one line in
NEWS. First, `validate_prior()` also refuses. Second,
`class = "Intercept", dpar = "sigma"` is still refused at
`set_prior()`, while the same prior written as `class = "sigma"` is
refused at the fit.

The new dot refusal of `print()` is consistent with core: on a curve,
`print(cv, foo = 1)` and `print(cv, digits = 3)` both give
"print() has no argument `...`. It takes: x", as `print.frmtmb_fit`
does. frmtmb.spline ignored both. NEWS says so. It could name `digits`,
because users of data frames will try it.

### 9. nit: core docs that name cox() without its package

`R/fit.R:374` (rendered in `man/frm.Rd:418`) says "a `cox()` fit spends
slot 1 on its baseline hazard", and `R/confint.R:2441` says
`cox(df = 6)`. Neither says frmtmb.survival (rule 6 of the brief). Add
"(frmtmb.survival)" to the first, because it is user-facing.

### 10. nit: pkgdown redirects and comments

`pkgdown::check_pkgdown()` reports "No problems found" for core and for
frmtmb.survival (pandoc on PATH). Building the redirects into a scratch
destination (`pkgdown:::build_redirects()`) writes all 11
`frmtmb.spline/` pages, and `reference/cox_baseline.html` refreshes to
`https://aforren1.github.io/frmtmb/frmtmb.survival/reference/cox_baseline.html`.
Gaps:

- The old subsite's `articles/index.html` and `authors.html` have no
  redirect.
- `frmtmb.spline/news/index.html` goes to frmtmb.survival's NEWS, which
  does not hold frmtmb.spline's history. The record says that history is
  in git, which is true, but the redirect implies otherwise. Point it at
  core's NEWS section, or accept it and say so.
- The comment above `redirects:` says the block is "Generated by
  dev/docsci-gen-redirects.R". These 12 lines are hand-written by
  `dev/splitsurv-retire-pkgdown.R`.
- `pkgdown.yaml`'s comment "frmtmb.sample suggests frmtmb.eam,
  frmtmb.latent and frmtmb.ode" omits frmtmb.survival. The install
  order itself is right: survival is installed before sample.

## What I checked and found right

### Guardrail 3, recomputed

```
$ (line counts from git show 19b162b9:<path>; code = not blank, not a comment)
curve: 1846 / 805      rp: 1139 / 411
$ git log --oneline 4f5ea39f..19b162b9 | wc -l                       -> 37
$ git log --oneline 4f5ea39f..19b162b9 -- <spline curve*.R>          -> f5770ca0
$ git log --oneline 4f5ea39f..19b162b9 -- <royston-parmar.R rp-*.R>  -> (none)
$ git log --oneline 4f5ea39f..19b162b9 -G"cox|bhaz" -- R/
2edf7ab5 (4 lines naming cox)  12b92070 (8)  0a8d9a02 (6)  58b97a2c (77)
$ grep -i -c "cox\|bhaz" dev/brmsport-verdicts.tsv                 -> 0
```

The added lines of 2edf7ab5, 12b92070 and 0a8d9a02 are family-list
and doc entries, so "1 substantive commit" (58b97a2c) is right. The
dependency column is right: gratia is used only in
`test-curve-gratia.R` and `test-curve-difference.R`, each behind
`skip_if_not_installed("gratia")`, so its new core Suggests is
justified. splines2 is used only by survival's `test-cox.R`. No core
test or vignette now uses splines2. The test-time column is the lane's
loaded-machine sum. My one-per-process run puts the nine curve files at
38.9 s of elapsed time, which agrees with "about a minute".

Reading the numbers: (c) is 174 code lines, 0 dependencies of its own,
1 substantive commit in five rounds, and 0.09% of core's test time.
`dev/core-boundary.md` calls such a feature "a weak candidate". The
lane says the same, and says that the gain is the seam. I agree. Note
that the seams are useful without the move: `simplex_pars =` and
`aterm_data =` could ship in core with cox() left there. The user
decides.

### Guardrail 1, four sessions

```
$ Rscript dev/review-splitsurv-g1.R base|nosurv|ns|surv
                          base        nosurv (lane)                       ns (namespace only)  surv
frm(family = "cox")       FIT -193.92 ERROR: The cox family is in the     FIT -193.92          FIT -193.92
frm(family = "Cox")       FIT         extension package frmtmb.survival,  FIT                  FIT
brmsfamily("cox")         OK          which is not loaded. Install it     OK                   OK
frm(brmsfamily("cox"))    FIT         and run library(frmtmb.survival),   FIT                  FIT
frm(brms::cox())          FIT         then use family = cox()  (all five) FIT                  FIT
frm(brms::brmsfamily("cox")) FIT                                          FIT                  FIT
get_prior(family = "cox") OK          (same refusal)                      OK                   OK
frm(family = cox())       FIT         could not find function "cox"       could not find ...   FIT
```

The log likelihood is the same number (-193.917912936552) in every arm
that fits, base included. `cox()` and `cox_baseline()` without the
package are plain R "not found" errors, as `wiener()` is without
frmtmb.eam. That is expected. `frm_compat()` gives 7650 rows without
frmtmb.survival and 7955 with it.

### Equality claims, rerun

All three of the lane's dump scripts were run again into
`dev/review-splitsurv-out/`, and the lane's compare scripts were run
with only the directory changed:

```
curve-before.log: frmtmb 0.73.0 C:/Users/adf44/source/r/rellib-r11 / frm_curve from frmtmb.spline
curve-after.log:  frmtmb 0.73.1 C:/Users/adf44/source/r/wt-splitsurv-lib / frm_curve from frmtmb
fit-before.log:   frmtmb 0.73.0 ... rellib-r11 ; cox from frmtmb
fit-after.log:    frmtmb 0.73.1 ... wt-splitsurv-lib ; cox from frmtmb.survival
```

All 19 curve outputs are "bitwise" or have an identical error message.
All 11 fits are "identical: TRUE; warnings identical: TRUE". The
`frm_compat()` result is 7955 = 7955 pairs, the same pair set, status
differs 0, note differs 5 (core's reworded frm_curve rows), features
identical, and 0 unresolved. Every line is the same as in the record.
The formals of all seven moved functions are `identical()` between base
and lane (`dev/review-splitsurv-formals.R`).

### Guardrail 2

`grep -rn "frmtmb:::\|getFromNamespace\|asNamespace\|:::" extensions/frmtmb.survival/{R,tests,vignettes}`
finds only comments that name `brms:::`/`flexsurv:::` and
`asNamespace("brms")` in `test-cox-links.R`. Every name in
frmtmb.survival's `importFrom(frmtmb, ...)` and every `frmtmb::` call
is exported by the lane's core. Core's NAMESPACE gains `frm_curve`,
`frm_curve_deriv`, `frm_curve_feature` and the two print methods, and
loses `cox` and `cox_baseline`. No other export was added.

### Tests, run by me (lane library first, NOT_CRAN=true, one file per process)

`bash dev/review-splitsurv-par.sh dev/review-splitsurv-out/tests dev/review-splitsurv-jobs.lst 4`.
Every `lib:` line says `frmtmb 0.73.1 C:/Users/adf44/source/r/wt-splitsurv-lib`.

| package | file | tests | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| frmtmb | test-compat.R | 36 | 767 | 0 | 0 | 0 | 0 |
| frmtmb | test-compat-extensions.R | 1 | 2 | 0 | 0 | 0 | 0 |
| frmtmb | test-cox-home.R | 2 | 8 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-check.R | 2 | 8 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-deriv.R | 12 | 48 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-difference.R | 13 | 71 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-gratia.R | 5 | 13 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-illcond.R | 1 | 4 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-new-levels.R | 8 | 57 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve-scale.R | 1 | 0 | 0 | 0 | 1 | 0 |
| frmtmb | test-curve-span.R | 13 | 73 | 0 | 0 | 0 | 0 |
| frmtmb | test-curve.R | 20 | 85 | 0 | 0 | 0 | 0 |
| frmtmb | test-simplex-pars.R | 5 | 33 | 0 | 0 | 0 | 0 |
| frmtmb.sample | test-cox-draws.R | 3 | 12 | 0 | 0 | 0 | 0 |
| frmtmb.sample | test-simulators.R | 10 | 54 | 0 | 0 | 0 | 0 |
| frmtmb.survival | (13 files) | 64 | 345 | 0 | 0 | 0 | 0 |

The survival file counts are the same as the record, file by file.
`test-curve-scale.R` skips outside the scale tier, as frmtmb.spline's
`test-scale.R` did. The survival suite in one session (`run-session.R`
with `FRMTMB_RELLIB` set to the lane library) gives
`forward|reverse|preload files=13/13 pass=345 fail=0 err=0 skip=0 warn=0`,
with no LEAK line.

### Test accounting

The 16 frmtmb.spline test files all have a home (record table, checked
against the files). `test-surface.R`'s "print methods say what was
checked" is in core's `test-curve.R`, the same code without the unused
`sp_surface_fit()` call. The 8 cox() tests of `test-famgaps.R` (36
assertions, 100 to 64) are in survival's `test-cox.R` by name. The
mixture, simulate-row, mv and `update()` cases that core's tests lost
are in survival's `test-cox.R`, and core's versions now use
`weibull` / a test family. The `cox_face_leave()` test is now
`simplex_face_leave()` on the test family, because the internal cannot
be reached from the extension. The cox() face behavior is still
covered by survival's "a cox() baseline at a face reaches the maximum,
silently" and by the `cox_face40` fit equality. The only gap is
finding 1.

### Cases the lane did not test

- A polynomial curve, `y ~ poly(x, 2, raw = TRUE)`: `frm_curve()`
  equals `frm_linpred(se.fit = TRUE)` exactly (difference 0, SE ratio
  error 0). The derivative agrees with `b1 + 2 b2 x` to 1.6e-9 SEs. The
  maximum is 0.3321744, which is the closed form `-b1 / (2 b2)`
  (`dev/review-splitsurv-curve2.R`, `-feat.R`).
- `frm_curve()` on a cox() fit, `time | cens(cens) ~ s(x, k = 6)`, has
  no `curve_check`, so it is not refused. The centered log hazard ratio
  follows `sin(2 pi x)` within its SEs at all 9 grid points.
- `update()` of an RP fit works with the package attached, and also
  with only its namespace loaded (finding 2 covers the absent case).
- A third-party `simplex_pars` family on frmtmb.sample: findings 3 and
  4.

### Other checks

- R's NEWS reader returns no entries for "# frmtmb.survival
  (development version)" alone
  (`tools:::.build_news_db_from_package_NEWS_md()` gives `NULL`; with
  "# frmtmb.survival 0.1.0" it gives 2 entries). That confirms the
  NOTE. It does not matter: the release renames the heading. Core's
  own dev section is dropped by the same reader, silently, because
  older sections exist. I did not rerun R CMD check.
- YAML: `yaml::read_yaml()` parses all 14 workflows and both
  `_pkgdown.yml`. The matrix entry, the sample job's sibling install,
  `pkgdown.yaml`'s extension list, `run-sessions.sh`,
  `run-session.R`, `build-docs.R`, `run-check.ps1` and `run-suite.ps1`
  all have frmtmb.survival where frmtmb.spline was. `run-scale.ps1`
  runs core's `test-curve-scale.R`. All 7 `helper-scale.R` copies
  (core plus 6 extensions) have one sha256, and their content equals
  frmtmb.spline's at 19b162b9 with line endings ignored. The workflow's
  count of 7 is still right. One cost to note: survival's Suggests
  includes brms (for one `test-cox-links.R` test), so
  `check-frmtmb-survival.yaml` installs brms and rstan.
- The new package matches the other extensions: DESCRIPTION, NEWS,
  README, package help with a "What this package reads from frmtmb"
  section, `helper-warnings.R` identical to core's, `.Rbuildignore`,
  `_pkgdown.yml`, the check workflow modeled on the spline one, and the
  coverage matrix entry. flexsurv, rstpm2, splines2 and survival are in
  Suggests.

## Overlaps for consolidation

As the lane lists: `extension_family_homes` (one entry per line),
`R/compat.R`, `_pkgdown.yml`, `test-coverage.yaml`, `run-sessions.sh` /
`run-session.R`, `test-compat-extensions.R`, README, NEWS, and
frmtmb.sample (DESCRIPTION Suggests, `check-frmtmb-sample.yaml`, two
tests, and finding 3 if it is fixed there).

## Short check of fix round 1

Verdict: **SHIP.** F1, F2, F3, F4, F6, F9 and F10 are fixed as the
review asked. F5 (the floor and the next release's scripts) is left for
consolidation, as the record says. No new finding.

The lane build was installed at 12:10, after every changed source
(`find ... -newer`). Every `lib:` line names
`frmtmb 0.73.1 C:/Users/adf44/source/r/wt-splitsurv-lib`.

### F2: saved fits, read where frmtmb.survival is not installed

The fits were saved with frmtmb.survival attached
(`dev/review-splitsurv-saved2.R save`). They were then read in a
separate process whose library holds a copy of the lane's frmtmb ahead
of rellib-r11, which has no frmtmb.survival
(`... saved2.R absent`; `survival installed: FALSE`; the copy is
deleted). `REFUSED[ext]` means an error of class
`frmtmb_extension_missing` that says "The <family> family of this model
is from the extension package frmtmb.survival, which is not installed
...":

| method | cox() | royston_parmar() | user family (global env) |
|---|---|---|---|
| print, summary, vcov, confint | REFUSED[ext] | REFUSED[ext] | OK |
| logLik, AIC, fixef, family, nobs | OK | OK | OK |
| frm_linpred(type = "link"), and with newdata | OK | OK | OK |
| frm_linpred(newdata, se.fit = TRUE) | REFUSED[ext] | REFUSED[ext] | OK |
| predict(), predict(newdata =) | cox()'s own stated refusal (no simulator) | REFUSED[ext] | "no simulator yet" (the family has none) |
| fitted, residuals | REFUSED[ext] | REFUSED[ext] | OK |
| simulate | cox()'s own stated refusal | REFUSED[ext] | "no simulator yet" |
| conditional_effects(effects = "x") | REFUSED[ext] | REFUSED[ext] | OK |
| anova(f, f) | OK | OK | OK |
| drop1, update(. ~ 1) | REFUSED[ext] | REFUSED[ext] | OK |
| frm_curve | REFUSED[ext] | REFUSED[ext] | OK |

No method fails with "could not find function". The methods that still
answer for cox() and RP (logLik, AIC, fixef, family, nobs, the link
predictor, anova of two fits) need none of frmtmb.survival's functions.
cox()'s `predict()` and `simulate()` give the family's stored reason,
which is a correct answer.

The user family was built by `frmtmb_family()` at top level. It records
no package (`attr(ufam, "frmtmb_package")` is `NULL`, and
`frm_family_package()` is `"frmtmb"`), so no method refuses it.

A mixture of a core component and an extension component,
`mixture(lognormal(), royston_parmar(df = 1, knots = numeric(0), bknots = ))`,
records the component packages `frmtmb frmtmb.survival`
(`dev/review-splitsurv-mix.R`). Read without frmtmb.survival,
print, summary, fitted and update refuse with "The
mixture(lognormal, royston_parmar) family of this model is from the
extension package frmtmb.survival". logLik answers. A
`mixture(lognormal(), royston_parmar())` without fixed knots is refused
before any of this, at construction, by RP's own "no knots yet" message.

### F3, F4

- F3: the third-party family's run of `frm_sample()` now announces
  "tw: dirichlet(1), brms's default, on the simplex the tsimplex()
  family declares ...". With `set_prior("dirichlet(20)", class = "tw")`
  the note is not shown. The draws are unchanged at the same seed
  (means 0.590505 0.404846 0.004649). `sample.R` no longer reads
  `kind == "cox"`; `t$kind` is used only in the text.
- F4 (`dev/review-splitsurv-typo2.R`): a typo in the declaration is
  refused in `frm()`, `get_prior()` and a fit with a prior, with "The
  tsimplex family declares simplex_pars = c(tw = "tw_rwa"), and its
  extra_pars() gives no parameter called 'tw_rwa'. It gives 'tw_raw'".
  A simplex with no coordinate is refused by name in `frm()` and
  `get_prior()`. The correct declaration, the guard-absent case, fits
  (logLik -582.8445619, as before), and `get_prior()` lists `tw`. The
  reason for running the check at frame assembly is right: `extra_pars`
  is `function(y, aterms)`, and cox()'s length comes from the basis
  built from the data, so construction cannot see the names.

### F1, F6, F9, F10

- F1: the record says 4 + 13 + 4. `test-cox-links.R` has "link and
  linkinv agree with links on cox()" with the 4 assertions; the file now
  has 4 tests and 31 assertions.
- F6: the RP `fitted` note reads "cox(), in this package too, refuses
  for the same reason". The compare now shows note differs 6, status
  differs 0.
- F9: `R/fit.R:374` and `R/confint.R:2441` name frmtmb.survival.
- F10: the redirects now include `frmtmb.spline/articles/index.html`
  and `authors.html`, and `news/index.html` points to
  `frmtmb.survival/articles/spline-history.html`. The article's body is
  frmtmb.spline's NEWS at 19b162b9 with every heading two levels down
  (`diff` with `sed 's/^#/###/'` is empty). pkgdown lists it as
  `articles/spline-history.html`, and `.Rbuildignore` excludes
  `vignettes/articles`. `check_pkgdown()` reports no problems for core
  or frmtmb.survival. The comments above `redirects:` and in
  `pkgdown.yaml` are corrected.

### Equality, rerun on the current build

The after arm of each dump was rerun against my base dumps. All 19
curve outputs are bitwise or have an identical message, and all 11 fits
and `get_prior` are identical, warnings included. `frm_compat()` gives
7955 = 7955 pairs, status differs 0, note differs 6 (the 5 frm_curve
rewordings and F6), and 0 unresolved.

### Tests (one file per process, NOT_CRAN=true)

All pass, with fail 0, err 0 and warn 0:

- core: test-extension-home.R 23, test-simplex-pars.R 37,
  test-curve-check.R 8, test-cox-home.R 8, test-curve.R 85,
  test-compat.R 767, test-compat-extensions.R 2, test-arg-refusal.R 122,
  test-mixture-reference.R 12, test-ordinal-mixture.R 99,
  test-simulate-density.R 620, test-famgaps.R 64, test-update-pool.R 21,
  test-families-extreme.R 108;
- frmtmb.sample: test-simplex-family-notes.R 4, test-cox-draws.R 12,
  test-simulators.R 54;
- frmtmb.survival: all 14 files, including the new test-saved-fits.R
  (27).

I did not rerun the full core suite. The new checks sit in
`build_objective()`, `response_mean()` and `sdr_of()`, which every fit
reaches, so the consolidation's full run is the one to read.
