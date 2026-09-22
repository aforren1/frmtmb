# Lane wt-brmsport: brms's own bin-1 suite, ported (item 2.6b)

Worktree `C:\Users\adf44\source\r\frmtmb-wt-brmsport`, base `e031e8c`
(frmtmb 0.59.0, frmtmb.sample 0.7.0). No file in `R/` was changed. The
lane adds a generated test tier, its harness, fixtures and a ledger.
Every run reads the private library `C:/Users/adf44/source/r/brmsport-lib`
(this worktree's two packages, installed by `dev/brmsport-install.R`
because the shared build `rellib-r3` predates commit `4e179c0`'s last
change to `R/`), then `pinlib` (StanHeaders 2.32.10, asserted by the
runner), then the user library (brms 2.23.0, testthat 3.3.2).

brms's source was re-fetched into `dev/brms-suite/` from the URL in
`dev/brms-suite-audit.md` section 1; its sha256 is
`b5f5bb56...992c37a5`, the recorded value, and
`dev/brmsport-blocks.R` stops if it is not.

This record is at punch round 2 (`dev/reviews/20260917-brmsport.md`,
"Recheck, round 1"). Section 10 lists each round-0 finding and what was
done; section 11 does the same for the recheck.

**THE TOTALS BELOW ARE THE STATE AT 0.60.0 AND ARE NO LONGER CURRENT.**
Lane `wt-shapes` fixed items 2.6d and 2.6f against them, and the ledger
was re-recorded, re-verdicted and regenerated: 192 passes became 231,
112 defects became 83, and the `pending 2.6d` class is gone. The
current figures, the 39 rows that moved and the five that changed
verdict without becoming passes are in `dev/shapes-findings.md`
section 6; `dev/brmsport-ledger.tsv` and
`dev/brmsport-log/ledger-summary.md` are regenerated in place. Rule 3
of section 1 and defect L5, L10, S1 and S4 of section 5 are DONE.

## 1. The user's rules and the divergence list

**The user decided three rules on 2026-09-17**, and the ledger applies
them:

1. Error message WORDING need not match brms. A refusal frmtmb makes
   for the same case, in its own words, is a PASS.
2. frmtmb objects must NOT carry brms's class names. An assertion of a
   brms class is a divergence, citing this decision.
3. The return shapes of `fitted`, `residuals`, `fixef`, `ngrps`, `vcov`
   and `summary` should match brms, and frmtmb.sample should offer
   `nsamples()` and `posterior_samples()`. Those rows are DEFECTS, item
   2.6f.

The divergence list was written before the fit-requiring files were
ported, from decisions that already existed. Punch round 1 corrected
its citations:

| # | decision | applies to | citation |
|---|---|---|---|
| D1 | `predict()` becomes brms's predictive summary; until then a `predict()` shape assertion is **pending 2.6d** | `predict has reasonable outputs` | user, 2026-09-16; plan row 2.6d |
| D2 | brms's errors are classed; an assertion needing the class is **pending 2.6e** | none: no bin-1 assertion passes `class =` | plan row 2.6e |
| D3 | `summary = FALSE` on `fixef()`, `ranef()`, `coef()`, `VarCorr()` of a FIT is refused: an ML fit has no draws | `ranef(fit2, summary = FALSE)` | `dev/brmsnames-findings.md` section B |
| D4 | `hypothesis()` on a fit: `Evid.Ratio` and `Post.Prob` are `NA`; `robust` and `scope` are refused | `hypothesis has reasonable ouputs` | `dev/brmsnames-findings.md`, "The hypothesis() object" |
| D5 | a sampled parameter with no brms counterpart keeps frmtmb's name; brms's derived draws columns are absent | rows 978, 984 | `dev/brmsnames-findings.md` A10 |
| D6 | WITHDRAWN as a decision. `inc_warmup` is a recorded GAP (`dev/brms-api-diff.md` (c), "Blocked, not small"), so row 107 is cannot transfer, absent | | |
| D7 | `conditional_effects()` on draws refuses `spaghetti`, `surface`, `select_points`, `too_far` | `conditional_effects` on draws | `dev/brmsnames-findings.md` "Found and NOT fixed" 9 |
| D8 | frmtmb has no `me()`, `gr(by =)`, `thres()`, and `arma()` needs `cov = TRUE`; an assertion reading a changed fixture term cannot transfer | fit3, fit5, fit6 terms | `dev/brms-suite-audit.md` section 5 |
| D9 | frmtmb objects do not carry brms's class names | the 18 `class` rows | **user decision, 2026-09-17, rule 2** (it cited 2.5b before, which is about generics, not class names) |
| D10 | Stan-side objects and switches have no counterpart | `testmode`, `save_pars()`, prior draws, `fit$fit@sim` | `dev/brms-api-diff.md` (c) |

## 2. What was built

| file | role |
|---|---|
| `tests/testthat/helper-brms-suite.R` | the harness (`brms_port()`, `brms_port_own()`, `brms_setup()`), the gate, the shims, the six fit fixtures |
| `extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R` | a generated copy of the above |
| `extensions/frmtmb.sample/tests/testthat/test-brms-suite-helper-copy.R` | a gated tier TEST that every object in the copy is identical to core's; hand-written, not generated |
| `extensions/frmtmb.sample/tests/testthat/helper-brms-suite-draws.R` | the draws fixtures; separate because core does not suggest frmtmb.sample |
| `tests/testthat/test-brms-suite-*.R` (10) | core and both tiers, GENERATED |
| `extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R` (5) | sample and both tiers, GENERATED, and the helper-copy test |
| `dev/brmsport-blocks.R` | locates the 95 bin-1 blocks by file and ordinal, stops unless 494 |
| `dev/brmsport-gen.R` | writes the 14 test files from brms's source and the verdicts |
| `dev/brmsport-record.sh`, `dev/brmsport-run.R` | one process per file, recording every assertion's outcome |
| `dev/brmsport-ledger.R` | records plus verdicts into `dev/brmsport-ledger.tsv` and the generated summary; stops on anything unclassified or contradicted |
| `dev/brmsport-verdicts-manual.tsv`, `dev/brmsport-verdicts-manual-fit.tsv` | verdict, class and reason of every assertion that does not hold |
| `dev/brmsport-verdicts-own.tsv` | the 27 rule-1 rows: frmtmb's own pattern for brms's case |
| `dev/brmsport-punch1-verdicts.R` | applies punch rounds 1 and 2's moves to the manual files, idempotently |
| `dev/brmsport-punch2-defects.R` | the two defects filed in punch round 2 (P1, P2) |
| `dev/brmsport-tier.sh` | runs the tier asserting the verdicts |
| `dev/brmsport-guards.R` | every harness rule firing and not firing, the copy check and the copy test, own-words specificity, the gate |
| `dev/brmsport-globcheck.ps1` | the gated glob on the tree and on an empty directory |
| `dev/brmsport-probe*.R`, `dev/brmsport-defects.R`, `dev/brmsport-fixtures.R`, `dev/brmsport-draws.R` | constructions |
| `dev/release/run-gated.ps1` | globs `test-brms-suite-*.R` in both packages and throws on an empty match |

**How an assertion is ported.** The generator copies brms's source text
of each top-level statement. An assertion becomes
`brms_port(id, verdict, reason, <brms's expect_* call>)`, or
`brms_port_own(id, pattern, note, <call>)` for a rule-1 row; any other
statement becomes `brms_setup(id, <statement>)`. The id is
`<topic>:<line>` in brms's file. The only rewrites are two text
substitutions in `tests.data-helpers.R` and the bindings at each file's
top: `brm` is `frm()`, `standata` is `frm(dry_run = "frame")` read
through brms's Stan data names, `update` drops brms's sampling
switches, `rename_pars` is the identity, `fit1` to `fit6` are fixtures.

**What "holds" means.** testthat would count it a pass AND none of
these, each detected automatically:

- the condition an `expect_error()`, `expect_warning()` or
  `expect_message()` caught is "could not find function", "object not
  found", "of mode 'function' was not found", "no package called", "not
  an exported object" or R's "is missing, with no default";
- the caught condition is frmtmb refusing an ARGUMENT ("has no
  argument", "unused argument", "cannot honor", "Cannot interpret ...
  argument"), unless brms's own pattern is about such a refusal;
- `expect_equal()`, `expect_identical()` or `expect_equivalent()` held
  on two NULLs where brms's expected value is not the literal `NULL`
  and one side is a `$` read that partial-matches or reads through an
  absent name (two correctly NULL `names()` still pass);
- `expect_null()`, `expect_true(is.null())` or `expect_length()` held
  on such a `$` read (a final absent name, which brms may mean, still
  passes);
- the assertion reads an object whose last assignment failed, including
  a failed complex assignment `x$a <-`, `x[[i]] <-`, `names(x) <-`, and
  a failed `assign("x", ...)`;
- its text reaches `brms::` or `brms:::` (never run).

Round 1's rule against `is.*` type checks on an all-`NA` value is GONE:
it also rejected an `NA` that is the right answer (recheck V3), and no
rule can tell the two apart. `expect_type()` and compound type checks on
an `NA` are not rejected either (U3, U4). A hollow one is marked by
hand, and the guards assert that these forms are not rejected.

27 runs that testthat alone counts as passes are not counted; one more
is marked hollow by hand (`brmsfit-methods:394`, `Evid.Ratio` is `NA` on
a fit and holds only as `is.numeric(NA_real_)`). `brmsfit-methods:1035`,
hand-marked in round 1, is now rejected by the partial-`$` rule.

**The own-words verdict (rule 1).** `brms_port_own()` evaluates brms's
assertion as written, which must NOT hold (else the row is a plain pass
and the verdict is stale), then the same call with the regexp replaced
by frmtmb's own pattern for that case (and `fixed` dropped), which must
hold under every rule above. The call is otherwise brms's, so the
refusal fires on brms's own case, and a different error cannot satisfy
it: `brm:81`'s "object 'sei' not found" is exactly the case the guard
constructs. It also rejects a row that reads a stale object, as
`brms_port()` does, and a pattern that matches the empty string or one
of seven unrelated control messages (`brms_own_controls`). The guards
check each of the 27 patterns against the message frmtmb gave for the
other 26 rows: none matches another row's message.

**What the tier asserts.** A pass must hold; every other verdict must
NOT hold, so the tier is green while the ledger says what is wrong, and
it goes red the day a verdict goes stale. `dev/brmsport-guards.R`,
`dev/brmsport-log/guards.txt`, 60 of 60: 53 harness cases each firing
or not firing as labeled (the reviewer's hollow constructions H1 to H7,
the recheck's U1 to U8 and V1, V2, their genuine counterparts, both
reviewer controls, and the own-words cases O2, O3 and five more); the
helper copy identical to core's while a copy with one body changed is
not; the helper-copy TEST passing on the real tree (33 expectations, 0
skips) and failing on the mutated copy (1 failure); the cross-row
specificity check with an unspecific control pattern seen to match; and
the gate (0 expectations with `FRMTMB_BRMS_FIT_TESTS` unset, 84 with it
set).

## 3. The fixtures

brms's `brmsfit_example1..6` are not in the tarball. The stand-ins
use brms's own DATA, read from its namespace, so every count the suite
asserts (40 rows, 10 patients, 4 visits, 8 subjects) is brms's.
`dev/brmsport-log/fixtures.txt`:

| fit | brms | frmtmb stand-in | changed | fit time, convergence |
|---|---|---|---|---|
| 1 | `count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) + (1 + Trt \| visit) + arma(visit, patient)`, `sigma ~ Trt`, student | same with `arma(visit, patient, cov = TRUE)` | `cov = TRUE` (D8) | 8.5 s; **does not converge**: nlminb code 1, false convergence (8), NaN standard errors, the `mo()` simplex flat in 3 directions |
| 2 | nonlinear gamma(identity) with `weights(AgeSD)` and `(1 \| ID1 \| patient)` | same | `Trt` enters the body as its numeric codes, which is brms's Stan data `C_1`, because frm() dies on the factor (L3) | 0.2 s, code 0 |
| 3 | `count ~ Trt * me(Age, AgeSD) + (1 + mmc(Age, volume) \| mm(patient, visit))` | `me(Age, AgeSD)` becomes `Age` | D8 | 0.1 s, code 0 |
| 4 | `rating ~ x1 + cs(x2) + (cs(x2) \|\| subject)`, `disc ~ 1`, sratio | `rating ~ x1 + cs(x2) + (1 + x2 \|\| subject)` | no `disc`; no `cs()` inside a group term | 0.2 s, code 0 |
| 5 | `count ~ Age + (1 \| gr(patient, by = gender))`, `mu2 ~ Age`, mixture(gaussian, exponential) | `(1 \| patient)` | D8 | 0.1 s, code 0 |
| 6 | `volume ~ Trt + gp(Age, by = Trt, gr = TRUE)` gaussian and `count ~ Trt + Age` poisson | `gp(Age)` | D8 | 0.2 s, code 0 |

**Fixture 1 does not converge, and no pass reads a value from it.** I
read every pass that touches `fit1` (and the draws sampled from it):
they assert dimensions, row counts, class, printed structure, that an
exp-link `sigma` is positive (`:326`), and that `scale = "linear"`
differs from the response scale (`:329`). None reads an estimate. A
future value assertion on fixture 1 would need a converged stand-in.

The draws stand-ins for 1, 2, 3 and 5 are
`frm_sample(chains = 1, iter = 75, warmup = 50, seed = 20260917)`,
brms's own design (`ndraws` 25, `niterations` 25): 13.0, 1.1, 0.7 and
1.1 s (`dev/brmsport-log/draws.txt`). The fits are built with
`do.call()` so the stored call carries formula and data, which
`update()` re-evaluates.

**The standata shim writes brms's zeros into `rcens`.** frmtmb keeps an
upper bound on every row; brms's `rcens` holds it only on interval rows
and 0 elsewhere, and the shim writes those zeros. So in `standata:163`
the zero half is the shim's, not frmtmb's; the reviewer's mutant with
999 on non-interval rows still passed, while a mutant shifting the
interval rows failed. The interval half is genuine. Recorded here
rather than removed, because frmtmb has no quantity for the zero half.

## 4. The ledger

`dev/brmsport-ledger.tsv`, one row per assertion: file, block, line,
tier, outcome, class, reason, id, held per package, brms's assertion
text, and frmtmb's message. Pasted verbatim from
`dev/brmsport-log/ledger-summary.md`:

<!-- BEGIN GENERATED: dev/brmsport-ledger.R -->

### Per file

| file | assertions | pass | defect | divergence | pending 2.6d | cannot transfer |
|---|---|---|---|---|---|---|
| `tests.brm.R` | 23 | 16 | 5 | 0 | 0 | 2 |
| `tests.brmsfit-helpers.R` | 3 | 0 | 0 | 0 | 0 | 3 |
| `tests.brmsfit-methods.R` | 223 | 61 | 92 | 29 | 12 | 29 |
| `tests.brmsformula.R` | 16 | 7 | 0 | 0 | 0 | 9 |
| `tests.brmsterms.R` | 5 | 0 | 0 | 0 | 0 | 5 |
| `tests.data-helpers.R` | 6 | 0 | 3 | 0 | 0 | 3 |
| `tests.emmeans.R` | 11 | 4 | 2 | 0 | 0 | 5 |
| `tests.families.R` | 84 | 45 | 1 | 1 | 0 | 37 |
| `tests.priors.R` | 36 | 17 | 4 | 4 | 0 | 11 |
| `tests.standata.R` | 87 | 42 | 5 | 0 | 0 | 40 |
| **total** | **494** | **192** | **112** | **34** | **12** | **144** |

### Outcome by class

| outcome | class | assertions |
|---|---|---|
| pass | - | 165 |
| pass | own-words | 27 |
| defect | accepts-refused | 3 |
| defect | argument | 23 |
| defect | different-error | 1 |
| defect | filed | 1 |
| defect | fit-data | 13 |
| defect | internal-error | 5 |
| defect | misparse | 1 |
| defect | naming | 2 |
| defect | output | 5 |
| defect | refuses-accepted | 28 |
| defect | shape | 22 |
| defect | silent | 4 |
| defect | spelling | 4 |
| divergence | class | 18 |
| divergence | hollow | 1 |
| divergence | no-draws | 8 |
| divergence | policy | 7 |
| pending 2.6d | shape | 12 |
| cannot transfer | absent | 106 |
| cannot transfer | brms-internal | 10 |
| cannot transfer | fixture | 1 |
| cannot transfer | mcmc | 4 |
| cannot transfer | no-draws | 1 |
| cannot transfer | stan | 22 |

Bin 1 passes: 192 of 494 (38.9%).
Against bins 1 and 2: 192 of 823 (23.3%); bin 2 was not ported.
Runs testthat alone would count as a pass and the harness does not (missing function or object, stale object, argument-name refusal, a NULL read through a partial $ match), over both packages: 27; hollow passes marked by hand: 1.

### The frmtmb.sample half

| file (tier) | assertions | cannot transfer | defect | divergence | pass |
|---|---|---|---|---|---|
| `tests.brmsfit-methods.R (sample)` | 61 | 13 | 12 | 0 | 36 |
| `tests.brmsformula.R (both)` | 16 | 9 | 0 | 0 | 7 |
| `tests.brmsterms.R (both)` | 4 | 4 | 0 | 0 | 0 |
| `tests.families.R (both)` | 84 | 37 | 1 | 1 | 45 |
| **total** | **165** | **63** | **13** | **1** | **88** |

frmtmb.sample passes 88 of its 165 runs (53.3%): 36 of the 61 sample-tier and 52 of the 104 both-tier assertions.
<!-- END GENERATED -->

**Classes.** Under `pass`: `own-words` is rule 1. Under `defect`:
`argument` is a brms argument refused by name; `refuses-accepted` and
`accepts-refused` say what they are; `shape` and `output` are rule 3 or
printed text; `fit-data` reads `fit$data`; `spelling` is a feature
frmtmb has under another name; `different-error` is a refusal of a
different thing; `silent` a wrong answer with nothing said. Under
`cannot transfer`: `absent` is a function, term or argument frmtmb does
not have (it could transfer if built), while `stan`, `brms-internal`,
`mcmc` and `fixture` cannot.

### Against the review's corrected totals

The review's generated totals (`dev/brmsport-log/rev-totals.txt`)
already include items 4, 6 and 8 of the punch list: 192 pass, 112
defect, 34 divergence, 12 pending 2.6d, 144 cannot transfer. This
ledger now has the same five figures. Punch round 1 had 117 defects and
139 cannot transfer, because it moved the five emmeans refusals of a
nonlinear or multivariate fit (`emmeans:27`, `:35`, `:38`, `:42`, `:50`)
to defect. Punch round 2 moved them back to cannot transfer, `absent`,
on the recheck's reading (R5): frmtmb's compatibility table declares
both cases refused for now (`R/compat.R:1356` and `:1612`, enforced by
`emm_mu_linpred()`), so they are a recorded limitation, not a refusal
of a supported model. `emmeans:21` and `:24` (`mo()`) are declared
nowhere and stay defects. That emmeans hides frmtmb's reason is filed
as its own defect, P2 in section 5. Within the unchanged totals, two
classes differ from round 1: `brmsfit-methods:1035` is now an ordinary
`fit-data` defect (13 rows) that a rule rejects, and
`brmsfit-methods:394` is a hand-marked `hollow` divergence (section 2).
The frmtmb.sample half agrees exactly (88, 13, 1, 63), since emmeans is
core only.

## 5. Defects, ranked

Silent wrong answers first, in the review's ranking; the numbers S1 to
S8 are the sections of `dev/brmsport-defects.R` and
`dev/brmsport-log/defects.txt`. frmtmb first and brms 2.23.0 on the same
call second; data seed 20260917.

### SILENT

- **S2. `ar()` and `ma()` accept an expression as the time term and
  fit.** With `cov = TRUE`: `ar(x + t, g)` logLik -30.33421 and
  `ar(t - 10 * x, g)` -31.00618, against -30.88513 for `ar(t, g)`;
  `ma(x + t, g)` -29.84347 (with an unrelated "singular convergence"
  warning) against -30.49819 for `ma(t, g)`. brms refuses all three,
  "Cannot coerce ... to a single variable name". Ledger `brm:106`.
- **S3. `ar(t, gr = g1/g2, cov = TRUE)` on numeric grouping codes
  groups by the QUOTIENT.** Series (1, 1) and (2, 2) both have ratio 1
  and merge; logLik -30.87554 against -30.88513, nothing said. brms
  refuses the term. `gr = g1:g2` on numeric codes is R's sequence
  operator, warns, then refuses on repeated times. Ledger `brm:108`.
- **S6. `variables()` of an ordinal fit omits the thresholds, on every
  ordinal family.** `cumulative`, `sratio`, `cratio` and `acat` on the
  same data each list 1 name for 3 parameters; fixture 4 lists 3 of 9
  (no thresholds, no `cs()` coefficients). The categorical control lists
  4 of 4. brms lists `b_Intercept[k]` and `bcs_x2[k]`. Not a bin-1
  assertion.
- **S5. A hypothesis with no relation is answered.**
  `hypothesis(fit3, "Trt1 + Age")` returns exactly the row of
  `"Trt1 + Age = 0"` (-2.224308, SE 2.336819); brms refuses, "Every
  hypothesis must be of the form 'left (= OR < OR >) right'". The other
  half is LOUD: `"b_Age x 0"` errors as the unknown parameter
  `b_b_Agex0`. Ledger `brmsfit-methods:417`.
- **S8. `fit$data` is `fit$data2`, by `$` partial matching.** A fit has
  no `data` element (`fit[["data", exact = TRUE]]` is `NULL`) and no
  `$.frmtmb_fit`, so `fit$data` returns `data2`: an empty list on
  fixture 1 (`identical(fit1$data, fit1$data2)` is TRUE), and
  `list(A = <6 x 6 matrix>)` on a fit with `data2 = list(A = A)`.
  brms keeps the 40 x 7 frame there. This is the hazard
  `frm_hazard_reads()` exists for, on a public object. 13 `fit-data`
  rows, `:1035` among them.
- **S1. `point_estimate` and `ndraws_point_estimate` are ignored on
  draws.** `posterior_epred(ds, point_estimate = "median",
  ndraws_point_estimate = 2)` returns 25 x 40 where brms returns 2 x 40.
  Ignoring an UNKNOWN argument is brms parity: brms also returns 25 x 40
  for `not_an_argument = 2`, so the first record's claim about the
  dots was wrong and is withdrawn. Ledger `:713`, `:714`.
- **S4. `fitted(ds)` and `residuals(ds)` on draws return `NULL`**;
  brms returns 40 x 4. Not a bin-1 assertion.
- **S7. Predictive errors on an ordinal fit.** `residuals(fit4)` is the
  response code minus the expected score, `y - sum(k P(Y = k))`: an
  IDENTITY, checked at full precision (max relative difference 0).
  `pp_check(fit4, "error_binned")` plots it. brms refuses both. A
  defined number, the least harmful of the eight. Ledger `:837`, `:703`.

### LOUD

| # | defect | ledger rows | construction |
|---|---|---|---|
| L1 | `categorical()` refuses integer responses with the FALSE claim "fewer than two categories"; `factor()` of the same codes is accepted; brms reports ncat 3 | 5 | defects.txt L1 |
| L2 | `frm(y ~ ., d)` dies with R's "'.' in formula and no 'data' argument" | 1 | L2 |
| L3 | a factor covariate in a nonlinear body dies with an RTMB advector error | 0 (fixture 2) | L3, `dev/brmsport-probe-nl.R` |
| L4 | newdata with a factor's numeric codes dies with "contrasts apply only to factors" | 3 | L4 |
| L5 | **rule 3, item 2.6f**: `fitted()`, `residuals()` return vectors (brms 40 x 4), `fixef()` a list per dpar (brms 9 x 4), `ngrps()` an integer vector (brms a list), `vcov()` 10 x 10 (brms 9 x 9), `summary()` has no `$fixed` or `$random`; four of these sit behind a refusal (`:314`, `:317`, `:825`, `:887`) | 22 shape, 5 output | L5 |
| L6 | 23 brms arguments refused by name, including `families:102`'s `order` | 23 | ledger class `argument` |
| L7 | `update()`: `data =` refits where brms demands `newdata`; three `bf()` update forms refused | 4 | defects.txt L6 |
| L8 | frmtmb's refusals of what brms answers: emmeans on a `mo()` term (2); `fitted()`/`residuals()` on a multivariate fit (4); the default Wald band on a nonlinear predictor (2); and more. emmeans on a nonlinear or multivariate fit is a declared limitation, cannot transfer (section 11) | 28 | ledger class `refuses-accepted` |
| L9 | **PRE-EXISTING, not a regression**: `default_prior()` and `get_prior()` validate the response and refuse brms's `rnorm` response under `Beta()`; brms returns 7 rows. The review measured the refusal at every commit back to base 0.58.0. priorform's P13 showed a table only because `dev/priorform-ledger.R` replaced brms's `y` with `runif(10)`, which neither its findings nor its review mention | `priors:74` | defects.txt L8, `dev/brmsport-rev-regress.R` |
| L10 | frmtmb.sample refuses `nsamples()` and `posterior_samples()` (rule 3) | 6 | ledger `:593` to `:635` |
| L11 | features under another spelling: `parnames()` (defined and refused by frmtmb.sample) and `conditional_smooths()` (documented as covered by `conditional_effects()`) | 4 | ledger class `spelling` |
| L12 | internal R errors: `set_rescor()` on the right-hand side, an unknown `pp_check` type, a missing `group` for `violin_grouped` (`:698`), a formula update with a variable the data lack, and 1 more | 5 | ledger class `internal-error` |
| L13 | `brm:81`: frmtmb evaluates the data first and stops on a missing variable, where brms refuses `se()` for weibull from the formula | 1 | `different-error` |
| L14 | prior-table naming `s(x).fx1` against brms's `sx_1`, where `variables()` already writes brms's name; no document decides it | `priors:55`, `:59` | ledger class `naming` |
| L15 | no "treated as continuous" warning on the ordinal fit4 | `:217` | ledger class `output` |
| P1 | `ma(x)` with brms's default `cov = FALSE` is refused on EVERY family, where brms fits the gaussian and student cases (Stan code of 1525 and 2411 characters) and refuses only the others ("Please set cov = TRUE when modeling MA structures for this family"). frmtmb documents the refusal (`R/autocor.R:102`, `dev/feature-gaps.md:228`: the residual-regression form is not implemented), so this is a known gap that the refusal names correctly; it matters here because `brm:110`'s own-words pass rests on the same message | 0 (`brm:110` is the poisson case, a correct refusal) | `dev/brmsport-log/punch2-defects.txt` P1 |
| P2 | emmeans hides frmtmb's refusal reason. `emmeans(fit2, "Age")` on the nonlinear fixture and `emmeans(fit6, "Age")` on the multivariate one both fail with "Perhaps a 'data' or 'params' argument is needed", which points the user at the wrong fix. frmtmb's `recover_data()` raises "emmeans support needs a linear mu predictor" and "emmeans support is univariate-only for now"; emmeans catches it with `try()`, prints it to stderr as `Error : ...`, and signals only its own message. Control: `y ~ x` reaches a grid of 1 row | 0 (`emmeans:27`, `:35`, `:38`, `:42`, `:50` are cannot transfer) | `dev/brmsport-log/punch2-defects.txt` P2 |

**Filed here, FIXED by lane `wt-shapes` in punch round 1.**
`R/confint.R` set `class(out) <- c("frmtmb_hypothesis",
"brmshypothesis")`, so `hypothesis()` output carried a brms class name
against rule 2. No bin-1 pass depended on it. The class is
`"frmtmb_hypothesis"` alone now; the stated reason for keeping it
(that `print()` and `plot()` would dispatch differently) was measured
false, since frmtmb exports both generics for its own class and that
class comes first.

## 6. Divergences, with their citations

34 assertions:

| class | n | citation |
|---|---|---|
| class | 18 | **user decision, 2026-09-17, rule 2**; includes `:373` (`is.brmsformula`), `:270` (`brms_conditional_effects`) and `priors:27` (`brmsprior`) |
| no-draws | 8 | `dev/brmsnames-findings.md` section B and "The hypothesis() object" |
| hollow | 1 | `:394`, `Evid.Ratio` NA on a fit (D4); marked hollow by hand, asserted to keep holding |
| policy | 7 | acat logit only (`?frmtmb-links`); priorform P1, P15, P26; `ranef(condVar =)` (`dev/brms-api-diff.md` (b)); `variables()` of a fit (D5, D10) |

`:399` to `:401` (`r_visit[4,Intercept]`) cite D4, which is about
`scope`; brmsnames "Found and NOT fixed" 16 is closer. Kept as
divergence, as the review did.

## 7. Tier runs

`sh dev/brmsport-tier.sh`, one R process per file, verdicts asserted,
`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`, log
`dev/brmsport-log/tier-final.txt`:

<!-- BEGIN GENERATED: tier -->
```
RESULT frmtmb test-brms-suite-brm.R pass=23 fail=0 err=0 skip=0 blocks=1
RESULT frmtmb test-brms-suite-brmsfit-helpers.R pass=3 fail=0 err=0 skip=0 blocks=2
RESULT frmtmb test-brms-suite-brmsformula.R pass=16 fail=0 err=0 skip=0 blocks=6
RESULT frmtmb test-brms-suite-brmsterms.R pass=5 fail=0 err=0 skip=0 blocks=2
RESULT frmtmb test-brms-suite-data-helpers.R pass=6 fail=0 err=0 skip=0 blocks=2
RESULT frmtmb test-brms-suite-emmeans.R pass=11 fail=0 err=0 skip=0 blocks=3
RESULT frmtmb test-brms-suite-families.R pass=84 fail=0 err=0 skip=0 blocks=7
RESULT frmtmb test-brms-suite-methods.R pass=162 fail=0 err=0 skip=0 blocks=23
RESULT frmtmb test-brms-suite-priors.R pass=36 fail=0 err=0 skip=0 blocks=12
RESULT frmtmb test-brms-suite-standata.R pass=87 fail=0 err=0 skip=0 blocks=19
RESULT frmtmb.sample test-brms-suite-brmsformula.R pass=16 fail=0 err=0 skip=0 blocks=6
RESULT frmtmb.sample test-brms-suite-brmsterms.R pass=4 fail=0 err=0 skip=0 blocks=1
RESULT frmtmb.sample test-brms-suite-families.R pass=84 fail=0 err=0 skip=0 blocks=7
RESULT frmtmb.sample test-brms-suite-helper-copy.R pass=33 fail=0 err=0 skip=0 blocks=1
RESULT frmtmb.sample test-brms-suite-methods.R pass=61 fail=0 err=0 skip=0 blocks=18
TIER ran 15 of 15 files
```

631 expectations over 15 files; 0 files with a failure, error or skip. Generated lines past 80 columns: 12.
<!-- END GENERATED -->

The sampler half RAN: frmtmb.sample's methods file 61 of 61 with 0
skips, on draws it sampled. The fifteenth file is the helper-copy test,
33 expectations.

`dev/release/run-gated.ps1`'s glob throws on an empty match:
`dev/brmsport-globcheck.ps1` lifts the block and runs it on the tree (15
jobs, the helper-copy test included) and on an empty directory (throws),
log
`dev/brmsport-log/globcheck.txt`.

**`dev/release/run-tests.R` reads `rellib-r3`**, which predates
`4e179c0`, while these verdicts are measured on 0.59.0 from
`brmsport-lib`. The release runner was not changed; the consolidating
session installs the round's reference library first.

The full package suites and `R CMD check` were NOT run: no package code
changed, the helper only defines functions, and the generated files
skip whole unless the gate variable is set.

## 8. What was not done

- `tests.log_lik.R`, bin 2 and bin 3.
- The 23 refusal assertions inside `tests.stancode.R`.
- No defect was fixed.
- A string literal longer than a line (a brms message or block label)
  cannot be broken without changing it; section 7 gives the generated
  count of lines past 80 columns.

## 9. Time

Round 0, from file times: reading and setup 15 min; the no-fit half 15;
fixtures and the divergence list 10; the fit half 20; guards, defects
and ledger 15; the record 15. Punch round 1: harness rules and guards
20; own-words verdict 15; verdict moves and ledger 15; constructions 15;
this record 15. Punch round 2 was interrupted by a machine shutdown and
finished by a second worker: harness rules, verdict moves and the copy
test before the shutdown; after it, establishing state 15, constructions
10, guards 15, record, ledger and tier 25, this record 15.

## 10. Punch round 1

| finding | status |
|---|---|
| MAJOR 1, the regression claim | FIXED. `priors:74` is a pre-existing loud defect (L9) naming priorform's `runif(10)` substitution; "regression" is gone from the record |
| MAJOR 2, S1 | FIXED. Narrowed to `point_estimate` and `ndraws_point_estimate`; the dots claim withdrawn with brms's own 25 x 40 as the control |
| MAJOR 3, S8 | FIXED. The mechanism is the `$` partial match of `data2`, with the `data2 = list(A = A)` construction; the 13 row reasons say so |
| MAJOR 4, hollow `:394`, `:698` | FIXED. Both are now rejected by harness rules (all-NA type check; missing-argument error), not by hand: `:394` divergence no-draws, `:698` defect internal-error. A hand-marked hollow verdict for `:394`, as `:1035` has, would assert that it keeps holding; the rule is stronger, since it can never be counted a pass. For `:394`, REVERSED in punch round 2: the all-NA rule is removed (section 11, R3 and R4), and `:394` is hand-marked hollow |
| MAJOR 5, seven hollow constructions | FIXED. H1 to H7 each fail, each has a genuine counterpart that passes, both reviewer controls still fail as they must (guards.txt) |
| MAJOR 6, citations | FIXED. D9 cites rule 2; D6 withdrawn and `:107` absent; `:213`, `:215` refuses-accepted, `:217` output, `priors:55`, `:59` naming |
| MAJOR 7, the three rules | FIXED. 27 own-words passes (25 wording plus `families:44`, `brm:110`), asserted with frmtmb's own pattern on brms's call; `brm:81` different-error; `families:102` argument; `:373`, `:270`, `priors:27` class; `:314`, `:317`, `:825`, `:887` shape; `:358`, `:359`, `:840`, `:841` refuses-accepted |
| MAJOR 8, spellings | FIXED. `:995`, `:269`, `:275`, `:277` defect spelling |
| MINOR 1, S2 wider | FIXED, with `ma(x + t, g)` and `ar(t - 10 * x, g)` on this lane's design |
| MINOR 2, S6 every ordinal family | FIXED, four families measured plus a categorical control |
| MINOR 3, S5 half silent | FIXED |
| MINOR 4, S7 computation | FIXED, stated as an identity with its residual |
| MINOR 5, glob fails open | FIXED, throws on an empty match; seen firing and not firing |
| MINOR 6, fixture 1 | RECORDED, and every pass reading `fit1` checked (section 3) |
| MINOR 7, `rcens` zeros | RECORDED (section 3) |
| MINOR 8, default-link passes | not changed: brms's test's weakness, and frmtmb honors non-default links |
| MINOR 9, helper copy drift | FIXED, function-by-function identity guard with a mutated copy |
| MINOR 10, emmeans | FIXED, all seven defect refuses-accepted. REVERSED for five rows in punch round 2 (section 11, R5) |
| MINOR 11, defects.txt numbering | FIXED, S1 to S8 and L1 to L8 match this record |
| MINOR 12, `:399` to `:401` citation | noted in section 6 |
| NIT, run-gated header and width | FIXED |
| NIT, `run-tests.R` reads `rellib-r3` | RECORDED (section 7), runner unchanged |
| NIT, `:203` and `:719` messages | now own-words passes under rule 1; the misleading "pass the original data via data =" for a misspelled effect remains a wording nit |
| NEW, `R/confint.R:2673` `brmshypothesis` | FILED (section 5); FIXED by `wt-shapes` punch round 1 |

## 11. Punch round 2

The recheck of round 1 (`dev/reviews/20260917-brmsport.md`, "Recheck,
round 1"). A machine shutdown interrupted the round. The second worker
compared file times against the logs, kept what held, and reran every
step after the last harness edit: the generator, the record (15 files,
one process each, `dev/brmsport-log/record-punch2.txt`), the ledger,
the generator again (no file changed), the guards, the tier and the
glob check.

| finding | status |
|---|---|
| R5, emmeans | FIXED. `emmeans:27`, `:35`, `:38`, `:50` (`R/compat.R:1612`) and `:42` (`:1356`) are cannot transfer, `absent`; `:21` and `:24` stay defects. Totals 192, 112, 34, 12, 144, as the recheck gave. NEW defect P2 filed: emmeans signals "Perhaps a 'data' or 'params' argument is needed" in place of frmtmb's reason (section 5) |
| R1, caveat on four rows | FIXED. `brm:110`, `standata:112`, `standata:179` and `brmsfit-methods:203` carry the caveat in `dev/brmsport-verdicts-own.tsv`, so it reaches the ledger reason |
| R1, `ma()` with `cov = FALSE` | FILED as P1 (section 5), with brms's Stan code for gaussian and student and its refusal for poisson. frmtmb documents the refusal as a gap, which the entry says |
| R1, specificity | FIXED. `brms_port_own()` rejects a pattern matching the empty string or one of seven control messages; guards O3 with `.`, `x*` and `is not supported` each fire. The guards also check every one of the 27 patterns against the messages frmtmb gave for the other 26 own-words rows: 0 matches, and the control pattern `response` is seen to match |
| R2, stale objects in own words | FIXED. `brms_port_own()` runs `brms_stale_reads()`; guard O2 fires on a stale object and passes on a fresh one |
| R6, helper copy | FIXED. `test-brms-suite-helper-copy.R` is a gated tier test in frmtmb.sample, 33 expectations. The guards run it on a temporary tree with the real copy (0 failures, 0 skips) and with one body changed (1 failure) |
| R3, U1, U2, U5 | FIXED. A new rule rejects `expect_null()`, `expect_true(is.null())` and `expect_length()` on a partial or intermediate `$` read. Control: a final absent name, `list(a = 1)$b`, and a real element still pass |
| R3, U3, U4 | RECORDED, not built. A type check on an `NA` cannot be told from an `NA` that is the right answer without knowing the right answer, which is the assertion's job. For the same reason round 1's all-`NA` `is.*` rule is removed. The guards assert that U3 and U4 are NOT rejected, so a later rule that starts rejecting them is seen. `brmsfit-methods:394`, the one live row of this shape, is marked hollow by hand, a divergence |
| R3, U6 | FIXED. "Cannot interpret ... argument" is an argument refusal |
| R3, U7 | FIXED. `assign("x", ...)` marks `x` as `<-` does; the control, a successful `assign()`, is fresh |
| R4, V1, V3 | FIXED. Two correctly NULL `names()` pass (V1); an `NA` that is the right answer passes (V3, by removing the rule) |
| R4, V2 | held since round 1: an assertion about an argument refusal passes |
| NIT, `run-gated.ps1` header | FIXED. The header reads as whole sentences |

Guards: 60 of 60 (`dev/brmsport-log/guards.txt`). Tier: 15 of 15 files,
631 expectations, 0 failures, errors or skips (section 7). Rows whose
verdict or class changed this round: the five emmeans rows,
`brmsfit-methods:394` (divergence, `no-draws` to `hollow`) and
`brmsfit-methods:1035` (defect, `hollow` to `fit-data`).
