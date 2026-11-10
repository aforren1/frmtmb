# Lane compat: the survey's compat-row corrections

Round opened 2026-11-09, base e0cce906 (frmtmb 0.80.0), measured on
`rellib-r19` (before) and `wt-compat-lib` in front of `rellib-r19`
(after). Brief: `round-briefs/r081-compat.md`. Nothing is committed.

## 1. The brief's premise does not hold: the ten rows were applied

The brief says the ten proposals of `dev/equiv-data/compat-proposals.tsv`
"have not been applied since 0.78.0". They were applied in 0.79.0, by
three commits of 2026-11-04:

- 0357639e (refusals of models brms builds): `cens()` x `disc`, `cens()`
  x `equidistant`, `vint()` x `sum_to_zero`, `trials()` x `cs_pred()`
  refused; `von_mises` x `mi()` works (group `mi_families`).
- 56fc92d3 (nonlinear bodies under ordinal, Gamma and fcor): the nl x
  family rows, `cs_pred()` x `nl` refused, `fcor()` x `nl` conditional.
- 22eb806c (latent residual): `poisson` x `arma()` conditional.

NEWS 0.79.0 lists them ("`frm_compat()` rows, each checked by running
both packages: `cens()` x `disc`, ..."). `dev/compat-rows.R` on
rellib-r19 (log `dev/compat-rows-base.log`) prints, for 0.80.0:

| id | pair | proposed | table on 0.80.0 |
|---|---|---|---|
| p167 | poisson x arma() | works | conditional |
| p158 | fcor() x nl | works | conditional |
| p171 | zero_inflated_beta x nl | works | works |
| p169 | von_mises x mi() | works | works |
| p166 | nbinom1 x nl | works | works |
| p014 | cs_pred() x nl | refused | refused |
| p146 | trials() x cs_pred() | refused | refused |
| p010 | cens() x disc | refused | refused |
| p011 | cens() x equidistant | refused | refused |
| p059 | vint() x sum_to_zero | refused | refused |

p167 and p158 read "conditional", not the proposed "works". That is the
better row: the latent-residual note states the conditions (frm()
warns on ordinal rows and where a Laplace correction moves sderr or
ar[1]; it refuses `cov = FALSE` after the inverse link), and the
`fcor()` x `nl` note keeps a family whose term would be a latent
residual refused. The survey's specs meet both conditions (gaussian
for p158, `arma(cov = TRUE)` for p167), so under either status the
pair agrees with brms. I left both rows as they are.

### Why the record still says "compat row wrong" on 0.80.0

`dev/equiv-flow.R` scores a pair against `row$status`, which
`dev/equiv-measure.R` reads from `dev/equiv-data/sample.tsv`: the
status FROZEN when the sample was drawn on 0.77.0. Lines 234-245,
264-273, 287-290 and 476-477 turn that frozen status into "compat row
wrong". So a pair whose row was corrected after the draw reads "compat
row wrong" on every later build, and the 0.80.0 record
(`dev/equiv-data/results.tsv`, run with `dev/rel080-equiv-out`) still
lists all ten, with the proposals regenerated from it. That is why the
brief reads them as open.

Shown by construction: the stock instrument on rellib-r19
(`bash dev/compat-equiv-run.sh stock <rellib-r19>
dev/compat-equiv-base-stock-out`) gives "compat row wrong" for all
ten. `dev/compat-equiv-live.R` runs the same flow, with one change: it
scores against the build's own `frm_compat()` row
(`bash dev/compat-equiv-run.sh live ...`). On rellib-r19, with no
change of mine:

```
p010 refused     (frozen untested) refusal correct
p011 refused     (frozen untested) refusal correct
p014 refused     (frozen untested) refusal correct
p059 refused     (frozen untested) refusal correct
p146 refused     (frozen works)    refusal correct
p158 conditional (frozen refused)  equivalent  L1 const -5.68e-14 grad 1.88e-06
p166 works       (frozen refused)  cannot measure (brms has no nbinom1)
p167 conditional (frozen refused)  equivalent  L1 const -3.41e-13 grad 1.04e-08
p169 works       (frozen refused)  equivalent  L1 const  1.14e-13 grad 1.87e-05
p171 works       (frozen refused)  equivalent  L1 const  0        grad 1.96e-05
```

(extracted from `dev/compat-equiv-base-live-out/*.tsv`.) Every one has
left "compat row wrong" on 0.80.0 already.

Unfixed, the consolidation rerun (`RECORD=1`) reads the same frozen
statuses and lists these ten as "compat row wrong" again, whatever the
table says. Two ways out: score against the live row (what
`dev/compat-equiv-live.R` does), or refresh the `status` column of
`sample.tsv` from the current frame while keeping the draw. The stratum
(and so the weight) of a draw is a property of the frame at draw time,
so the first is the smaller change. The coordinator asked for it after
the first report; section 9 has the fix and its measurements.

### The refusals name the reason and the term

From the live run's logs (`dev/compat-equiv-base-live-out/<id>.log`):

- p010, p011 (`y | cens(c) ~ ...` on cumulative): "cens()/trunc() need
  a family with a CDF, and 'cumulative' supplies none;
  frm_compat("cens()") lists the families that do." brms: "Argument
  'cens' is not supported for family 'cumulative(logit)'."
- p059: "cumulative: the addition term `vint()` is not one this family
  reads, so writing it would change nothing about the fit. This family
  takes `thres()`, `weights()`." brms: "Argument 'vint' is not
  supported ..."
- p146: "Category specific effects are not supported for this family.
  cs() needs an ordinal family: cumulative(), sratio(), cratio(),
  acat() or hurdle_cumulative()". brms: the same first sentence.
- p014: "cs() is not supported in a nonlinear formula body: b1 * x +
  cs(z). brms 2.23.0 cannot compile it either." brms's program fails
  its compile (semantic error, line 74).

No refused pair fits without a word, and no refusal hides its cause.

## 2. The rule behind the four refused rows

The four refusals of 0.79.0 were one row per sampled pair. Their
neighbors disagreed with them. `dev/compat-ordinal-grid.R` on the
0.80.0 frame (`dev/compat-frame-base-out`) shows, for the core addition
terms against the features only an ordinal family has:

- `cens()` x `disc` and x `equidistant` refused, `cens()` x
  `sum_to_zero` untested;
- `vint()` x `sum_to_zero` refused, x `disc` and x `equidistant`
  untested;
- `trials()` x `cs_pred()` refused, and `cens()`, `trunc()`, `se()`,
  `mi()` and `rate()` x `cs_pred()` WORKS (the kind-level rule
  `kind:aterm` x `kind:special`);
- `thres()` x `cens()`, `trunc()`, `se()`, `trials()`, `mi()` refused
  (a loop of its own), `thres()` x `rate()` WORKS (`kind:aterm` x
  `kind:aterm`), x `vint()` and `vreal()` untested.

The reason is one fact. The five ordinal families declare the same
addition terms (`dev/compat-ordinal-accepts.R`, rellib-r19):

```
cumulative         index subset thres weights
hurdle_cumulative  index subset thres weights
sratio             index subset thres weights
cratio             index subset thres weights
acat               index subset thres weights
```

and `thres()`, `disc`, `equidistant`, `sum_to_zero` and `cs_pred()`
belong to them alone (with all eight extensions loaded, no family but
these five has a row other than "refused" with any of the five; awk
over `dev/compat-frame-base-out/compat.tsv`). So `trials()`, `cens()`,
`trunc()`, `se()`, `mi()`, `vint()`, `vreal()` and `rate()` can never
meet one of them.

### Measured, every pair the rule covers

`dev/compat-ordinal-measure.R` (seed 20261109, n = 80; log
`dev/compat-ordinal-out-base.log`, table
`dev/compat-ordinal-out/measure-base.tsv`) spells each of the 8 terms
with each of the 5 features on an ordinal family (cumulative; sratio
for `cs()`), and runs frm() and brms 2.23.0's `make_stancode()`.
`weights()` is the control:

```
control weights(): frm fits 5 of 5 ; brms builds 5 of 5
other terms: frm refuses 40 of 40 ; names the term 40 of 40 ; brms refuses 40 of 40
```

Every brms refusal is "Argument '<term>' is not supported for family
'<family>(logit)'". Every frm() refusal names the term: the
accepts_aterms message for `trials()`, `vint()`, `vreal()`, `rate()`,
the CDF message for `cens()` and `trunc()`, the se() message, and
"Argument 'mi' is not supported for family ..., as in brms" for `mi()`.
The other three ordinal families declare the same list, so the
spelling on cumulative and sratio stands for them.

### The change

`R/compat.R`:

- a group `ordinal_only = c("thres()", "disc", "equidistant",
  "sum_to_zero", "cs_pred()")`;
- one rule per term, `r(<term>, "group:ordinal_only", "refused", ...)`,
  for the eight terms above, in the "disc and the threshold
  structures" section;
- removed, because the rule now decides them: the rows `cens()` x
  `disc`, `cens()` x `equidistant`, `vint()` x `sum_to_zero` (which
  carried `override = TRUE` for no tie), `trials()` x `cs_pred()`, and
  the `thres()` loop over `cens()`, `trunc()`, `se()`, `trials()`,
  `mi()`.

The rule's signature is (name, group) = (3, 2). It beats the kind-level
"works" rules and the `vint()`/`vreal()` x `*` untested blankets (3, 0),
and no (3, 2) or (3, 3) rule disagrees with it (`frm_compat()` builds,
and `frmtmb_compat_validate()` in test-compat.R passes).

### Every row that moves

`dev/compat-dump.R` (runs `dev/equiv-dump-compat.R` on a library list)
before and after, diffed by `dev/compat-frame-diff.R` (log
`dev/compat-frame-diff.log`):

```
before: 12012 pairs, 7518 core
after:  12012 pairs, 7518 core
pairs only before: 0 ; only after: 0

status moves: 31 (core 31 )

untested -> refused    works -> refused
                 25                   6
  cens() x cs_pred()                 works     -> refused   core=TRUE
  cens() x sum_to_zero               untested  -> refused   core=TRUE
  cs_pred() x mi()                   works     -> refused   core=TRUE
  cs_pred() x rate()                 works     -> refused   core=TRUE
  cs_pred() x se()                   works     -> refused   core=TRUE
  cs_pred() x trunc()                works     -> refused   core=TRUE
  cs_pred() x vint()                 untested  -> refused   core=TRUE
  cs_pred() x vreal()                untested  -> refused   core=TRUE
  disc x mi()                        untested  -> refused   core=TRUE
  disc x rate()                      untested  -> refused   core=TRUE
  disc x se()                        untested  -> refused   core=TRUE
  disc x trials()                    untested  -> refused   core=TRUE
  disc x trunc()                     untested  -> refused   core=TRUE
  disc x vint()                      untested  -> refused   core=TRUE
  disc x vreal()                     untested  -> refused   core=TRUE
  equidistant x mi()                 untested  -> refused   core=TRUE
  equidistant x rate()               untested  -> refused   core=TRUE
  equidistant x se()                 untested  -> refused   core=TRUE
  equidistant x trials()             untested  -> refused   core=TRUE
  equidistant x trunc()              untested  -> refused   core=TRUE
  equidistant x vint()               untested  -> refused   core=TRUE
  equidistant x vreal()              untested  -> refused   core=TRUE
  mi() x sum_to_zero                 untested  -> refused   core=TRUE
  rate() x sum_to_zero               untested  -> refused   core=TRUE
  rate() x thres()                   works     -> refused   core=TRUE
  se() x sum_to_zero                 untested  -> refused   core=TRUE
  sum_to_zero x trials()             untested  -> refused   core=TRUE
  sum_to_zero x trunc()              untested  -> refused   core=TRUE
  sum_to_zero x vreal()              untested  -> refused   core=TRUE
  thres() x vint()                   untested  -> refused   core=TRUE
  thres() x vreal()                  untested  -> refused   core=TRUE

note-only changes: 9 (core 9 )
  cens() x disc                      refused
  cens() x equidistant               refused
  cens() x thres()                   refused
  cs_pred() x trials()               refused
  mi() x thres()                     refused
  se() x thres()                     refused
  sum_to_zero x vint()               refused
  thres() x trials()                 refused
  thres() x trunc()                  refused
```

31 + 9 = 40, the 40 pairs measured above; no other row moves, and no
extension-only pair moves. Core status counts (awk over the two dumps,
`in_core_only`): before conditional 2485, refused 975, untested 1199,
works 2859; after conditional 2485, refused 1006, untested 1174, works
2853. Of the 40, the survey's sample holds only the four refused pairs
(p010, p011, p059, p146), whose notes change and statuses do not.

No row moves to "works", and no fit changes: frm() refused all 40
before and after.

## 3. The nl x family rule

Settled in 0.79.0 (56fc92d3). The table no longer has an `nl` x
`kind:family` refusal; each core family has its own row, measured
(`nl_same`, `nl_ord`, `nl_multi`, `nl_start` notes, test-nl-family.R
and test-nl-brms.R). `dev/equiv-data/nl-sweep.tsv` (recorded on
0.80.0) reads every family "works" in the table but asym_laplace
("conditional"), with "nl honored" for every single-location family
(D1, the ordinal drop, is gone: cumulative -395.9068 linear and nl).
The six "nl differs from both" rows are the multi-location families
(categorical, multinomial, dirichlet_multinomial, dirichlet,
dirichlet2, logistic_normal) and asym_laplace: the first six because
brms writes one body into every location predictor, which is not the
linear model with per-category coefficients, and the table says so
(`nl_multi`, an L1 identity in test-nl-brms.R); asym_laplace's row is
conditional for its kink. The summary's line "rows the table says
refused: ;" is empty because no family row is refused. p171 and p166
sit on `nl_same` and the nbinom1 row; both are "works" and pinned in
test-nl-family.R. This lane moves no nl row.

## 4. Tests

- `tests/testthat/test-compat.R`, new test "an addition term no ordinal
  family takes refuses its features": the rule's term list equals the
  core addition terms minus what the ordinal families declare (read off
  `compat_core_family_accepts()`, so the two cannot drift); each of the
  8 x 5 pairs reads "refused"; the group holds the five features; each
  pair's frm() refusal names the term; and the absent case, `weights()`
  with each feature, fits and keeps "works" for `weights()` x `thres()`.
  The loops read the expected lists, not the rule, so a missing rule
  fails them.
  SEEN TO FAIL on rellib-r19 (`FRMTMB_RELLIB` unset, log
  `dev/release/compat-base1-files/frmtmb--test-compat.txt`):
  `pass=816 fail=33 err=0`: 31 status failures (25 "untested", 6
  "works", the same 31 pairs as the frame diff), the rule list and the
  group. The refusal checks pass on the base, as they should: frm()
  already refused. On the lane build: `pass=851 fail=0 err=0 skip=0
  warn=0`.
- Pins of the 0.79.0 rows that had no `frm_compat()` assertion beside
  the behavior: `poisson` x `arma()` "conditional" in test-latres.R
  (beside the poisson `arma(cov = TRUE)` fit), `fcor()` x `nl`
  "conditional" in test-nl-family.R (beside the fcor body fit), and
  each `mi_families` family x `mi()` "works" in test-mi-families.R
  (beside each family's mi() fit; von_mises is p169). These pass on
  both builds; they guard against a regression and do not pin a
  defect. p166, p171 and p014 were pinned already (test-nl-family.R).

## 5. What I did not do

- Change any of the ten survey rows: all ten were right on 0.80.0.
- Edit the record (`dev/equiv-data/`) or rerun it with `RECORD=1`. The
  instrument changed after the first report, on the coordinator's
  request (section 9); `dev/equiv-flow.R` did not.
- Move `weights()`, `subset()` and `index()` x `disc`, `equidistant`,
  `sum_to_zero` from "untested". The ordinal families declare them and
  the control fits (weights() 5 of 5), but a fit that runs is not a
  measurement of agreement with brms, and the brief asks for the
  surveyed rows.
- OpenBLAS pass: no numerical code changed (only the table and tests).

## 6. The lane build: remeasure and suites

### The ten pairs on the lane build

`bash dev/compat-equiv-run.sh live "<wt-compat-lib>;<rellib-r19>"
dev/compat-equiv-lane-live-out` (every log's `lib:` line reads
`wt-compat-lib/frmtmb 0.80.0`), extracted from the tsvs:

```
"p010" "refused" "refusal correct" "" NA NA
"p011" "refused" "refusal correct" "" NA NA
"p014" "refused" "refusal correct" "" NA NA
"p059" "refused" "refusal correct" "" NA NA
"p146" "refused" "refusal correct" "" NA NA
"p158" "conditional" "equivalent" "L1" -5.6843418860808e-14 1.87812502971951e-06
"p166" "works" "cannot measure" "" NA NA
"p167" "conditional" "equivalent" "L1" -3.41060513164848e-13 1.04238987441718e-08
"p169" "works" "equivalent" "L1" 1.13686837721616e-13 1.86907602606201e-05
"p171" "works" "equivalent" "L1" 0 1.96498519257141e-05
```

(columns: id, live status, verdict, level, L1 constant, L1 gradient.)
Against the record (`dev/equiv-data/results.tsv`, 0.80.0): all ten
were "compat row wrong"; on the live status 5 are "refusal correct", 4
"equivalent" at L1 and p166 "cannot measure" (brms has no nbinom1). The
same ten through the stock instrument (`... stock ...
dev/compat-equiv-lane-stock-out`) read "compat row wrong", all ten, as
on rellib-r19: that verdict is fixed by the frozen status (section 1).
The numbers are identical to the rellib-r19 live run: my change does
not touch these pairs' statuses (p010, p011, p059 and p146 change note
only).

`dev/compat-ordinal-measure.R` on the lane build (log
`dev/compat-ordinal-out-lane.log`) repeats the base reading:

```
control weights(): frm fits 5 of 5 ; brms builds 5 of 5
other terms: frm refuses 40 of 40 ; names the term 40 of 40 ; brms refuses 40 of 40
```

### Every test file, one per process

`bash dev/release/run-par.sh compat-all dev/compat-all.jobs 3` with
`NOT_CRAN=true FRMTMB_RELLIB=<wt-compat-lib>` (log
`dev/release/compat-all.log`, "COMPAT-ALL ran 509 of 509"); every
`lib:` line names `wt-compat-lib` for core; the extensions load from
rellib-r19 over it. Summed by awk from the RESULT lines:

```
frmtmb files=317 pass=23361 fail=0 err=0 skip=239 warn=0
frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=33 pass=1848 fail=0 err=0 skip=5 warn=0
frmtmb.latent files=12 pass=687 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=86 pass=3967 fail=0 err=0 skip=18 warn=0
frmtmb.spatial files=8 pass=478 fail=0 err=0 skip=6 warn=0
frmtmb.survival files=15 pass=552 fail=0 err=0 skip=0 warn=0
```

The skips are the gated tier's. Of the 55 files with a skip, three
read the compat table (grep for `frm_compat` or `compat_`); gated
(`FRMTMB_BRMS_FIT_TESTS=true`, `dev/release/compat-gated.log`):
test-fcor.R pass=57, frmtmb.sample test-sampling-ported.R pass=214,
frmtmb.spatial test-sar.R pass=126, each fail=0 err=0 skip=0 warn=0.
The known red test-se-tier.R passed in this run (pass=95 fail=0).

### Every suite in one session, three ways

`bash dev/release/run-sessions.sh compat-sess 3` with `NOT_CRAN=true
FRMTMB_RELLIB=<wt-compat-lib>` (log `dev/release/compat-sess.log`,
"COMPAT-SESS ran 27 of 27"), the SESSION lines:

```
SESSION frmtmb forward files=304/317 pass=23361 fail=0 err=0 skip=239 warn=0
SESSION frmtmb reverse files=304/317 pass=23511 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=304/317 pass=23604 fail=0 err=0 skip=241 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam reverse files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam preload files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward files=81/86 pass=3967 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.sample reverse files=81/86 pass=3967 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.sample preload files=81/86 pass=3967 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
```

Every LEAK line carries `[with load of ...]` (extension loads by
test-arg-refusal.R, test-cox-home.R, test-families-extreme.R and their
kin, and rxode2's and sf's environment variables); none is untagged.

### R CMD check

`bash dev/compat-check.sh` (core only; no extension changed), built and
checked in `dev/compat-check/` with `--as-cran`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` = wt-compat-lib, rellib-r19, user library. Log kept as
`dev/compat-check.log`; the tree is deleted.

```
* checking HTML version of manual ... [30s] NOTE
Skipping checking math rendering: package 'V8' unavailable
Status: 1 NOTE
```

The one NOTE is the expected V8 math-rendering note.

## 7. Release notes

- NEWS: one bullet under `# frmtmb (development version)`, `##
  Compatibility table`: rows that now read refused. No row now reads
  works, so there is no second bullet. Not BREAKING: no fit that ran
  before refuses now (frm() refused all 40 pairs on 0.80.0).
- Bump: patch. Floor: none (no new dependency, no extension change).

## 8. Files

Package: `R/compat.R`, `tests/testthat/test-compat.R`,
`tests/testthat/test-latres.R`, `tests/testthat/test-mi-families.R`,
`tests/testthat/test-nl-family.R`, `NEWS.md`.

dev (new): `compat-findings.md`, `compat-rows.R`,
`compat-ordinal-grid.R`, `compat-ordinal-accepts.R`,
`compat-ordinal-measure.R`, `compat-equiv-live.R`,
`compat-equiv-run.sh`, `compat-dump.R`, `compat-frame-diff.R`,
`compat-check.sh`, and the job lists `compat-all.jobs`,
`compat-changed.jobs`, `compat-gated.jobs`.

dev (changed, sections 9 and 10): `equiv-measure.R`, `equiv-summary.R`,
`equiv-specs.R` (p154's entry).
dev (new, section 9): `compat-equiv-measure-orig.R` (the instrument as
at e0cce906, kept for the comparison), `compat-noreg-run.sh`,
`compat-noreg-compare.R`, `compat-live-status.R`,
`compat-summary-check.R`, `compat-p154-probe.R`.

## 9. The instrument scores against the live row

Asked by the coordinator after the first report: keep `sample.tsv`'s
status as the stratum, score each verdict against the build's own
`frm_compat()` row, record both, and generate the proposals from the
live row.

### The change

- `dev/equiv-measure.R`: reads the live row with `frm_compat()`, passes
  `eq_measure()` a copy of the sample row whose `status` is the live
  one, and flattens with the ORIGINAL row, so the `status` column of
  the per-pair tsv (which `equiv-summary.R` reads as the stratum) stays
  the sample's. Two new columns, `stratum_status` and `live_status`, go
  into the tsv and so into results.tsv; the rds carries both too.
  `dev/equiv-flow.R` is not touched (lane translator changed it).
- `dev/equiv-summary.R`, one hunk of five lines after the proposals are
  selected (line 202), away from translator's two hunks (the
  weak-evidence rule near line 86 and its sentence near line 285):
  a proposal's `status_now` is the pair's `live_status` where the
  per-pair tsv has one. Strata, weights and rates still read `status`.

### 1. The ten pairs leave "compat row wrong"

The edited stock instrument on the lane build
(`dev/compat-diff-edit-out/<id>.log`, `lib:` wt-compat-lib):

```
verdict p171 : equivalent level L1 - every level agrees
verdict p169 : equivalent level L1 - every level agrees
verdict p167 : equivalent level L1 - every level agrees
verdict p166 : cannot measure level  - frmtmb fits; no brms analogue: brms has no nbinom1 family
verdict p158 : equivalent level L1 - every level agrees
verdict p146 : refusal correct level  - both refuse
verdict p059 : refusal correct level  - both refuse
verdict p014 : refusal correct level  - both refuse (brms at the compile)
verdict p011 : refusal correct level  - both refuse
verdict p010 : refusal correct level  - both refuse
```

The same as `dev/compat-equiv-live.R` gave (section 6).

### 2. No regression where the stratum equals the live row

25 pairs whose sample status equals the live status, across all four
strata (untested p001 p002 p007 p008 p015 p040 p051; conditional p062
p064 p065 p071 p073 p100; works p119 p126 p129 p130 p141 p152; refused
p153 p156 p157 p161 p162 p172), verdicts equivalent at L1, L1+L2 and
L4, refusal correct, cannot measure, and p002's measured defect. Each
was measured by the instrument as at e0cce906
(`dev/compat-equiv-measure-orig.R`) and by the edited one, on the lane
build, in separate processes (`bash dev/compat-noreg-run.sh ...
dev/compat-noreg <ids>`), and compared column by column, `seconds`
excepted (`dev/compat-noreg-compare.R`, log
`dev/compat-noreg-compare.log`):

```
columns only in B: stratum_status live_status
pairs with any difference outside seconds: 0 of 25
```

Every other column is identical as text, the L1 constants and
gradients included.

### Every pair whose live row differs from its stratum

`dev/compat-live-status.R` over all 172 pairs, no fit
(`dev/compat-live-out/base.log` on rellib-r19 and `lane.log` on the
lane build; the two lists are identical, 17 pairs, every pair has a
live row):

```
p010 cens()                 x disc             untested    -> refused     record: compat row wrong
p011 cens()                 x equidistant      untested    -> refused     record: compat row wrong
p013 cs_pred()              x frm_curve        untested    -> works       record: equivalent
p014 cs_pred()              x nl               untested    -> refused     record: compat row wrong
p043 s()                    x prior            untested    -> works       record: equivalent
p049 thres()                x nl               untested    -> works       record: equivalent
p059 vint()                 x sum_to_zero      untested    -> refused     record: compat row wrong
p070 |ID|                   x mmc()            conditional -> works       record: equivalent
p087 gr_cov                 x mmc()            conditional -> works       record: equivalent
p095 me()                   x prior            conditional -> works       record: equivalent
p146 trials()               x cs_pred()        works       -> refused     record: compat row wrong
p154 cens()                 x cosy()           refused     -> conditional record: refusal correct
p158 fcor()                 x nl               refused     -> conditional record: compat row wrong
p166 nbinom1                x nl               refused     -> works       record: compat row wrong
p167 poisson                x arma()           refused     -> conditional record: compat row wrong
p169 von_mises              x mi()             refused     -> works       record: compat row wrong
p171 zero_inflated_beta     x nl               refused     -> works       record: compat row wrong
```

All 17 measured with both instruments on the lane build
(`dev/compat-diff-compare.log`):

```
p010   compat row wrong           untested / refused DIFFERS: verdict reason proposal
p011   compat row wrong           untested / refused DIFFERS: verdict reason proposal
p013   equivalent                 untested / works identical
p014   compat row wrong           untested / refused DIFFERS: verdict reason proposal
p043   equivalent                 untested / works identical
p049   equivalent                 untested / works identical
p059   compat row wrong           untested / refused DIFFERS: verdict reason proposal
p070   equivalent                 conditional / works identical
p087   equivalent                 conditional / works identical
p095   equivalent                 conditional / works identical
p146   compat row wrong           works / refused  DIFFERS: verdict reason proposal
p154   refusal correct            refused / conditional DIFFERS: verdict reason proposal
p158   compat row wrong           refused / conditional DIFFERS: verdict reason proposal
p166   compat row wrong           refused / works  DIFFERS: verdict proposal
p167   compat row wrong           refused / conditional DIFFERS: verdict reason proposal
p169   compat row wrong           refused / works  DIFFERS: verdict reason proposal
p171   compat row wrong           refused / works  DIFFERS: verdict reason proposal
pairs with any difference outside seconds: 11 of 17
```

(the verdict column shows the original instrument's.) What it does:

- 6 pairs (p013, p043, p049, p070, p087, p095): the flow only asks
  whether a row is "refused", so untested or conditional to works
  changes nothing; identical.
- 10 pairs leave "compat row wrong" as in part 1.
- **p154 regresses**: "refusal correct" to "compat row wrong", with a
  proposal "refused". Its live row `cens()` x `cosy()` is conditional
  since the latent residual of 0.79.0: refused on gaussian and student,
  working on a family with no residual. Its spec is gaussian (drawn as a
  refused-stratum control), so it does not meet the condition, both
  packages refuse, and the flow, which cannot read a condition, calls
  the row wrong. The proposal is false: the row's refusal half is
  right, and its working half is what the spec never tried.

  Spelled on poisson, where the condition holds,
  `dev/compat-p154-probe.R` (seed 154, 16 of 200 rows right-censored
  at 4; log `dev/compat-p154-probe.log`) runs the same flow without
  editing the specs:

  ```
  verdict: equivalent | level: L1 | reason: every level agrees
  l1 const: -1.13687e-13 grad: 1.79e-11 points: 3 ok: TRUE extra: -1.70530e-13 -4.54747e-13
  ```

  So the conditional row is right. The first addendum report left the
  respelling to the coordinator; section 10 applies it.

### 3. The summary reads strata from the sample, proposals from the live row

`dev/compat-summary-check.R` rebuilds the 172 per-pair tsvs from the
record (measured verdicts, summary columns dropped), and
`dev/equiv-summary.R` runs on them without `--record`
(`dev/compat-sumA-out`). Its "Verdicts by stratum" through "agree,
among measured pairs only" block is the record's line for line, except
the corrected-weights line (89.0 against the record's 88.9; that line
reads the frame dump `dev/equiv-frame-out/compat.tsv`, absent in this
worktree, and is not touched by the change). Overlaying the 42 pairs
measured with the edited instrument (`dev/compat-sumB-out`), the
summary changes only in:

```
< | untested | 53 (46) | 0 (0) | 4 (4) | 5 (3) | 0 (0) | 8 (7) |
> | untested | 53 (46) | 4 (4) | 0 (0) | 5 (3) | 0 (0) | 8 (7) |
< | works | 36 (31) | 0 (0) | 1 (1) | 0 (0) | 0 (0) | 3 (3) |
< | refused | 0 (0) | 15 (15) | 5 (5) | 0 (0) | 0 (0) | 0 (0) |
> | works | 36 (31) | 1 (1) | 0 (0) | 0 (0) | 0 (0) | 3 (3) |
> | refused | 4 (4) | 14 (14) | 1 (1) | 0 (0) | 0 (0) | 1 (1) |
< | refused | 20 | 95.0 [76.4, 99.1] | 95.0 [76.4, 99.1] | 95.0 [76.4, 99.1] (20) | 95.0 [76.4, 99.1] |
> | refused | 20 | 95.0 [76.4, 99.1] | 95.0 [76.4, 99.1] | 100.0 [83.2, 100.0] (19) | 95.0 [76.4, 99.1] |
< - agree, among measured pairs only: 96.9
> - agree, among measured pairs only: 97.4
< Compat-row proposals: 10 (10 from sampled pairs, ...
> Compat-row proposals: 1 (1 from sampled pairs, ...
```

(and the cost line, which sums wall seconds). The strata are unchanged:
the `status` column of results.tsv counts 57 conditional, 20 refused,
60 untested and 35 works unique pairs, and the draws per stratum stay
70, 70, 40 and 20. The headline does not move (89.0), because the flow
counts a "compat row wrong" whose behavior agrees with brms as agree.
compat-proposals.tsv holds one row:

```
"p154"	"cens()"	"cosy()"	"conditional"	"refused"	"compat row wrong"
```

with `status_now` the live "conditional". It no longer proposes any
applied row; its one proposal is the p154 artifact above, gone after
section 10.

### Not done

- `RECORD=1`, and any edit of `dev/equiv-data/` or of the generated
  blocks of `dev/equiv-findings.md`.

## 10. p154 respelled on poisson

Decided by the coordinator: the pair is the feature pair `cens()` x
`cosy()`, and the family is a spelling, as p158's spelling was corrected
before. `dev/equiv-specs.R`, p154's own entry only (lines 969-979 here;
translator's nearest hunks start at its lines 891 and 1066 of the base,
so the files merge per hunk): `eq_dat("poisson", seed = 154)`, rows
above 4 right-censored at 4, `family = "poisson()"`. A comment in the
entry keeps the fact the old spelling measured: on gaussian and student
both packages refuse `cens()` beside `cosy()` (brms: "Invalid addition
arguments for this model"; frm(): "cosy(): cens() cannot be combined
with a residual correlation term", the record's p154 row).

The edited stock instrument on the lane build
(`dev/compat-p154-out/p154.log`, `lib:` wt-compat-lib):

```
pair p154 : cens() x cosy() [ stratum refused | live conditional ]
verdict p154 : equivalent level L1 - every level agrees
```

From its tsv (status, stratum_status, live_status, verdict, level,
l1_const, l1_grad, l1_points, l1_extra, proposal):

```
"refused" "refused" "conditional" "equivalent" "L1" -1.13686837721616e-13 1.78811410123103e-11 3 "-1.71e-13,-4.55e-13" ""
```

The stratum stays "refused".

The summary on the rebuilt record with the 43 edited pairs overlaid
(`dev/compat-sumC-out`) proposed nothing, and that exposed a latent
defect: `prop$source <- "sampled pair"` on a 0-row frame stops the
summary ("replacement has 1 row, data has 0"). Fixed in the same hunk
of `dev/equiv-summary.R` (`rep("sampled pair", nrow(prop))`). Then:

```
Compat-row proposals: 0 (0 from sampled pairs, 0 from the nl x family sweep; dev/compat-sumC-out/compat-proposals.tsv).
```

and compat-proposals.tsv is its header alone. Against the overlay of
section 9 (`dev/compat-sumB-out`) the summary changes only in the
refused stratum's verdict row (`| refused | 5 (5) | 14 (14) | 0 (0) |
0 (0) | 0 (0) | 1 (1) |`, p154 from compat row wrong to equivalent),
the levels (L1 72 unique pairs, 88 draws; none 38, 39), the proposal
count and the cost line. Every rate line is unchanged; the headline
stays 89.0.

## 11. Final verification (m1-m4)

Run 2026-11-10 on the lane build, after the review's four MINORs were
fixed.

### The m1-m4 changes

- m1: `R/compat.R` `compat_ordinal_only_rules()` derives the
  ordinal-only refusals for every registered addition term that no
  ordinal family declares, in place of the hand-written per-term rules.
- `R/parse.R` `registered_aterm_of()` and `aterm_spelling()`: a key is
  the term name plus an index 1..arity, no longer the key with its
  digits stripped.
- Tests: `tests/testthat/test-compat.R` and
  `tests/testthat/test-autocor.R`; `NEWS.md`.
- m2: `dev/equiv-flow.R`, both "both refuse" branches: a live row of
  refused or conditional gives "refusal correct".
- m3: `dev/equiv-measure.R` and the new `dev/equiv-live-row.R`: the
  live row is read with the extensions of the frame dump, in a child
  process.
- m4: `dev/equiv-summary.R` writes a generated line on stratum against
  live status; `dev/equiv-findings.md` section 2.6 has a paragraph on
  which status is which.

### 0. The installed core is the worktree's

`dev/compat-m-nscheck.R` (log `dev/compat-m-nscheck.log`) evaluates
the three functions from `R/compat.R` and `R/parse.R` and compares them
with the installed namespace:

```
lib: C:/Users/adf44/source/r/wt-compat-lib 0.80.0
compat_ordinal_only_rules    source:TRUE installed:TRUE identical:TRUE
registered_aterm_of          source:TRUE installed:TRUE identical:TRUE
aterm_spelling               source:TRUE installed:TRUE identical:TRUE
```

No reinstall was needed (installed 02:37:49, `R/parse.R` last changed
02:36:18).

### 1. No regression, 42 pairs, on the final instrument

`bash dev/compat-noreg-run.sh` (started 02:40:13) into
`dev/compat-noreg2-orig-out` and `dev/compat-noreg2-edit-out`. Its
first outputs (02:40:17) are newer than `dev/equiv-flow.R` (02:32:09),
`dev/equiv-live-row.R` (02:32:47) and `dev/equiv-measure.R`
(02:32:56). Every log's `lib:` line reads
`wt-compat-lib/frmtmb 0.80.0`.

The edited arm of p008 and p015 halted in the first batch, about 7
seconds into the run:

```
Error: no frm_compat() row for cens() x bounds
Error: no frm_compat() row for cs_pred() x prior
```

(kept as `dev/compat-m-failed/p0{08,15}-edit-first.log`). The child
`dev/equiv-live-row.R p008,p015` alone gives "untested" for both
(`dev/compat-m-live-p008.tsv`), and the two ids rerun with
`COMPAT_JOBS=1` measure normally ("stratum untested | live untested",
equivalent at L1, both arms). So a child failed at start under parallel
load; its message was lost. Instrument defect, not fixed here:
`dev/equiv-measure.R` reports any child failure (`rc != 0`, no file) as
"no frm_compat() row", which names the wrong cause. It should report
the child's exit status apart from a missing row.

`Rscript dev/compat-noreg-compare.R dev/compat-noreg2-orig-out
dev/compat-noreg2-edit-out` (log `dev/compat-noreg2-compare.log`):

```
pairs in A: 42 ; in B: 42 ; both: 42
columns only in B: stratum_status live_status
pairs with any difference outside seconds: 11 of 42
```

The 25 pairs whose stratum equals the live row are identical in every
column but the two new ones. Of the 17 that differ in status, 6 (p013,
p043, p049, p070, p087, p095) are identical, and 11 differ as section 9
says: p010, p011, p014, p059, p146 "compat row wrong" to "refusal
correct"; p158, p167, p169, p171 to "equivalent"; p166 to "cannot
measure"; p154 (respelled on poisson, section 10) "refusal correct"
(the original arm reads the old gaussian spec) to equivalent at L1,
`l1_const -1.13686837721616e-13`, `l1_grad 1.78811410123103e-11`,
`l1_points 3`, as in section 10.

### 2. The summary

`Rscript dev/compat-summary-check.R dev/compat-sumD-out
--live=dev/compat-live-out/lane-ext.tsv dev/compat-noreg2-edit-out`
("rebuilt 172 pairs; overlaid 42"), then `Rscript dev/equiv-summary.R
dev/compat-sumD-out` (log `dev/compat-sumD-summary.log`):

```
Strata are the compat status at the draw (sample.tsv); verdicts are scored against the measured build's own frm_compat() row. 17 of 172 pairs differ.
- **agree with brms (the headline): 89.0, Wald [84.5, 93.6]; bootstrap parametric [84.2, 93.2], nonparametric [84.2, 93.3]**
Compat-row proposals: 0 (0 from sampled pairs, 0 from the nl x family sweep; dev/compat-sumD-out/compat-proposals.tsv).
```

`compat-proposals.tsv` is its header alone. The verdict rows by
stratum and the rate lines are those of section 10.

### 3. Every test file, one per process

`NOT_CRAN=true FRMTMB_RELLIB=<wt-compat-lib> bash
dev/release/run-par.sh compat-m-all dev/compat-all.jobs 3` (log
`dev/release/compat-m-all.log`, "COMPAT-M-ALL ran 509 of 509"; the job
list equals every `test-*.R` of core and the eight extensions). All
509 `lib:` lines name `wt-compat-lib` for frmtmb; the extensions load
from rellib-r19. Summed by awk from the RESULT lines:

```
frmtmb files=317 pass=23381 fail=0 err=0 skip=239 warn=0
frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=33 pass=1848 fail=0 err=0 skip=5 warn=0
frmtmb.latent files=12 pass=687 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=86 pass=3967 fail=0 err=0 skip=18 warn=0
frmtmb.spatial files=8 pass=478 fail=0 err=0 skip=6 warn=0
frmtmb.survival files=15 pass=552 fail=0 err=0 skip=0 warn=0
```

No failure, so no rerun on rellib-r19 was needed. test-compat.R
pass=868, test-autocor.R pass=195, test-se-tier.R (the known red)
pass=95, each fail=0.

### 4. Every suite in one session, three ways

`NOT_CRAN=true FRMTMB_RELLIB=<wt-compat-lib> bash
dev/release/run-sessions.sh compat-m-sess 3` (log
`dev/release/compat-m-sess.log`, "COMPAT-M-SESS ran 27 of 27", 03:13
to 04:13, so the core sessions took under 1 hour). Every session's
`lib:` names `wt-compat-lib` for frmtmb.

```
SESSION frmtmb forward files=304/317 pass=23381 fail=0 err=0 skip=239 warn=0
SESSION frmtmb reverse files=304/317 pass=23546 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=304/317 pass=23649 fail=0 err=0 skip=241 warn=0
```

The 24 extension sessions read the same counts as section 6, each
fail=0 err=0 warn=0. DETAIL lines: 0. LEAK lines: 46, all 46 tagged
`[with load of ...]`; none untagged.

### 5. R CMD check

`bash dev/compat-check.sh` (core, `--as-cran`, 04:13 to 04:44; log
kept as `dev/compat-m-check.log`, the tree deleted):

```
* checking HTML version of manual ... [26s] NOTE
Skipping checking math rendering: package 'V8' unavailable
Status: 1 NOTE
```

The expected V8 note only.

### Review round 2: SHIP

`dev/reviews/2026-11-10-compat.md`, section "Round 2": **SHIP**; m1 to
m4 closed with evidence, and R/parse.R correct for every registered
term. One follow-up, filed in `dev/test-backlog.md` under "Filed in the
round of 2026-11-09 / From lane compat":

- **R1, MINOR, predates the lane.** Two registered addition terms can
  claim the same key, and `frmtmb_register_aterm()` does not refuse it.
  `x` (arity 3) beside `x2`: `y | x2(z) + x(u, v, w)` lets x's second
  argument overwrite x2(z) silently; the other order gives a misleading
  "Duplicated addition term" refusal. `a` (arity 12) beside `a1` does
  the same at a11 and a12. `reward1` beside reward, or `vint1` beside
  core vint(), is accepted. No shipped extension collides. Fix:
  `frmtmb_register_aterm()` refuses a name whose keys collide with
  another registered term's keys or a core key.

The reviewer's cost remark on m2: a conditional row that should be
refused, where both packages refuse, now passes as "refusal correct".
The reason text names it ("the refused half of a conditional row"),
but the summary does not count such pairs; a count line would make
them readable.

### Files of this verification

New: `dev/compat-m-nscheck.R`, its log, `dev/compat-m-failed/`,
`dev/compat-m-live-p008.tsv`, `dev/compat-noreg2-compare.log`,
`dev/compat-sumD-out/`, `dev/compat-sumD-check.log`,
`dev/compat-sumD-summary.log`, `dev/compat-m-*.driver.log`,
`dev/compat-m-check.log`, `dev/release/compat-m-all*`,
`dev/release/compat-m-sess*`. Changed: `dev/test-backlog.md` (R1
entry) and this file. Steps 0 to 5 changed no package code, test or
`dev/equiv-*` file; the subsection below changes `dev/equiv-measure.R`
and `dev/equiv-live-row.R` only.

### The live-row child failure: diagnosis and fix

Asked by the coordinator after this verification: the consolidation
reruns all 172 pairs with `RECORD=1` in parallel, so the p008/p015
failure of step 1 must not recur or be misreported.

**Diagnosis.** The two children died 90 ms apart (parent logs
02:40:21.98 and 02:40:22.07, 9 s into the first batch), and neither
printed a word. What the evidence rules out:

- A shared output file: the row file was `tempfile()` in each parent's
  own session folder, so no two children shared a path.
- A lost R error: a child's stdout and stderr reach the parent's log
  (`dev/compat-m-stderr-probe.log`: "child-err", "Error: boom", rc 1),
  so an R error, a load failure or a corrupt lazy-load message would
  have shown.
- A missing row: run alone, the child gives "untested" for both pairs
  (`dev/compat-m-live-p008.tsv`), and the core-only and all-loaded
  tables agree on all 172 pairs (`dev/compat-live-out/lane-ext.log`).
- A race between children, or on the library or the stan cache (the
  child does not touch the stan cache): `bash dev/compat-m-race.sh 3 6`
  (3 children at the same instant, 6 rounds, read by
  `dev/compat-m-race-read.R`, outputs `dev/compat-m-race1/`) gave 18 of
  18 with rc 0, all 8 extensions loaded and no NA row.
  `bash dev/compat-m-race-measure.sh 5 ... p001 p008 p015` (the first
  batch's three edit arms, parent and child, on the unfixed instrument,
  `dev/compat-m-race2/`) gave 15 of 15 verdicts and no failure.

So it did not reproduce in 33 parallel runs. The signature (an exit
with no output, two children at the same instant) is that of an end
from outside R: an injected self-kill (`tools::pskill()`) gives exactly
that, exit status 1 and no output. The cause stays unproven. The
children's command line, `...\x64\Rscript.exe dev/equiv-live-row.R p008
<tempfile>`, carried no worktree path, so a process sweep that matches
on a path could not tell whose they were.

**The fix**, in the two instrument files only:

- `dev/equiv-measure.R`, `live_row()`: the child writes to a per-pair,
  per-PID, per-attempt file in the parent's session folder
  (`live-<id>-<pid>-<attempt>.tsv`), and its stdout and stderr go to a
  log beside it. A child that exits non-zero, or exits 0 without a row
  file, is reported as "live-row child, attempt k failed: exit status
  <n>[, but no row file]; output: <last 8 lines | (none)>", and retried
  once. Two failures stop with "dev/equiv-live-row.R failed twice for
  <id>: ...". "no frm_compat() row" now means only that the child ran
  and the table has no unique row. The child starts with the script's
  full path, so its command line names the worktree.
- `dev/equiv-live-row.R`: writes `<out>.part` and renames it, so a
  child stopped mid-write leaves no file that reads as complete. An
  extension that is installed but fails to load now stops with its
  error; `requireNamespace(quietly = TRUE)` skipped it without a word,
  and the child read a smaller table.

**Injected failures.** `bash dev/compat-m-inject.sh p008 kill error
silent0 once none` runs the fixed `dev/equiv-measure.R` in a scratch
copy of the instrument whose `equiv-live-row.R` is
`dev/compat-m-inject-row.R` (log `dev/compat-m-inject.log`; the copy is
deleted):

```
== mode kill: exit status 1
live-row child, attempt 1 failed: exit status 1; output: (none)
live-row child, attempt 2 failed: exit status 1; output: (none)
  dev/equiv-live-row.R failed twice for p008: exit status 1; output: (none)
== mode error: exit status 1
live-row child, attempt 1 failed: exit status 1; output: Error: injected: lazy-load database is corrupt | Execution halted
...
== mode silent0: exit status 1
live-row child, attempt 1 failed: exit status 0, but no row file; output: (none)
...
== mode once: exit status 0
live-row child, attempt 1 failed: exit status 1; output: (none)
verdict p008 : equivalent level L1 - every level agrees
== mode none: exit status 0
verdict p008 : equivalent level L1 - every level agrees
```

**Parallel rerun on the fixed instrument.**
`dev/compat-m-race-measure.sh`, 3 rounds each: p008 p015 p154 and p146
p002 p166 at 3 jobs at once, and p008 p015 p001 p002 at 4 jobs at once
(logs `dev/compat-m-race3{a,b,c}.log`, outputs `dev/compat-m-race3/`).
30 of 30 measurements ended rc 0, and no log has a retry line. Each
round's tsvs against `dev/compat-noreg2-edit-out` (p008 and p015 there
are the single-job reruns), by `dev/compat-noreg-compare.R` (log
`dev/compat-m-race3-compare.log`): "pairs with any difference outside
seconds: 0" in all 9 rounds, the L1 constants included.

**Overlap.** Lane translator's worktree changes `dev/equiv-findings.md`,
`equiv-flow.R`, `equiv-lib.R`, `equiv-specs.R` and `equiv-summary.R`.
It does not change `dev/equiv-measure.R` and has no
`dev/equiv-live-row.R` (read-only `git status`), so the two files merge
without overlap.

New files: `dev/compat-m-race.sh`, `dev/compat-m-race-read.R`,
`dev/compat-m-race-measure.sh`, `dev/compat-m-inject.sh`,
`dev/compat-m-inject-row.R`, and their logs and output folders.
