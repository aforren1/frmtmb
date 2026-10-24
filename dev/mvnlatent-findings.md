# Lane mvnlatent: `mixture_mvn()` from core into frmtmb.latent

Round of 2026-10-23, base 19b162b9 (frmtmb 0.73.1), branch
`wt-mvnlatent`. Applies the fourth proposed move of
`dev/core-boundary.md`. Status: DONE, uncommitted. Core, frmtmb.latent
and frmtmb.sample are roxygenised and installed in
`C:/Users/adf44/source/r/wt-mvnlatent-lib`, fix round 1 included
(section 9). R CMD check --as-cran after fix round 1: frmtmb 1 NOTE,
frmtmb.latent 1 NOTE (both the V8 math-rendering NOTE of the HTML
manual); frmtmb.sample OK in the first pass, unchanged since (section
8).

Logs are in `dev/mvnlatent-out/` (gitignored); every number below is
pasted from the log named beside it.

## 1. Measurements before the move (guardrail 3)

### Lines of code in R/ that move

Script `dev/mvnlatent-measure.R`, run on `R/families.R` and
`R/compat.R` of a `git archive 19b162b9` snapshot; log
`dev/mvnlatent-out/measure.txt`. A block's count includes the comment
block right above it.

```
mvn_cov_models   lines 7229-7238: code 2, comment 8, blank 0
mvn_cov_spec     lines 7240-7336: code 83, comment 14, blank 0
mixture_mvn      lines 7338-7551: code 123, comment 88, blank 3
mvn_sim_rows     lines 7553-7586: code 27, comment 7, blank 0
TOTAL families.R: code 235, comment 117, blank 3, all 355
compat.R lines naming mixture_mvn: 13
compat.R rules naming mixture_mvn: 9
families.R total lines: 10236
R/ total lines: 86204 in 70 files
share of R/ lines: 0.41%
```

So 355 lines of `R/families.R` (235 of code) and 9 compatibility rules
(13 lines of `R/compat.R`, plus one entry of the structure vocabulary)
move. That is 0.41% of the lines in `R/`. (The reviewer counts 358:
the script leaves out the 3 blank lines between the four blocks.)

What core LOSES is less, because the seam and the new prose add lines
back. `git diff --numstat 19b162b9 -- R/`: +101 / -417, net -316,
after the first pass (the reviewer's count); +182 / -417, net -235,
after fix round 1, which adds `us_theta()`, the help text of the
component interface and the corrected `?us_sigma`.

### Test time, and its share of the core suite

Every suite of the base, one file per process, 8 at a time,
`NOT_CRAN=true`, on rellib-r11 (frmtmb 0.73.0) with the test trees of
19b162b9 (`bash dev/mvnlatent-par.sh before <snapshot> 8 ...`, runner
`dev/mvnlatent-run1.R`, which also prints each `test_that()` block's
`real` time). Log `dev/mvnlatent-out/before.log` and
`dev/mvnlatent-out/before/*.txt`.

```
RESULT test-mvn-mixture.R pass=277 fail=0 err=0 skip=0 warn=0 tests=9 elapsed=6.3
TIME  test-mvn-mixture.R  mixture_mvn matches direct ML                               0.40   8
TIME  test-mvn-mixture.R  mixture_mvn recovers the faithful clusters                  0.04   5
TIME  test-mvn-mixture.R  mixture_mvn class means can depend on covariates            0.08   5
TIME  test-mvn-mixture.R  mixture_mvn dpar overrides and covariate gating work        0.12   4
TIME  test-mvn-mixture.R  mixture_mvn validation and guards                           0.16  12
TIME  test-mvn-mixture.R  mixture_mvn covariance models size the template correctly   0.64  74
TIME  test-mvn-mixture.R  mixture_mvn covariance models agree at a shared parameter point 0.20 17
TIME  test-mvn-mixture.R  mixture_mvn class means take covariates under every model   0.59  40
TIME  test-mvn-mixture.R  mixture_mvn matches mclust for every covariance model       2.63 112
TIME  test-quadrature-defects.R  mixture_mvn() refuses REML and profile too           0.08   3
TIME  test-structure.R  a rowwise mixture-type family declares capabilities only      0.08  15
TIME  test-unpinned-seams.R  method = 'predict' refuses a family that draws the response whole 0.15 6
TIME  test-prior-location-dpar.R  a family brms does not have is refused by frmtmb's rule 0.02 3
TIME  test-simulators.R (frmtmb.sample)  mixture_mvn simulates its class covariances  0.12   8
```

(The last column is the number of expectations of the block.) The core
suite's per-file elapsed times sum to 2046 s over 262 files (same log;
the sum is of single-process times, measured with 8 processes running,
so it is a load-inflated figure for both numerator and denominator).
`test-mvn-mixture.R` takes 6.3 s of it, 0.31%; its nine blocks take
4.86 s of test time, and the four blocks in other core files 0.33 s
(of which the `test-structure.R` block is only half mixture_mvn).
So the move takes about 0.3% of the core suite's time out of it.

That share is ONE measurement under load, good as an order of
magnitude only. The same 262 core files summed to 1678.9 s in the first
after-arm run (`after.log` of the first pass), against 2045.6 s before:
a swing of 367 s between two runs of the same machine, 58 times the
6.3 s the move removes.

### Churn in the last five rounds

Range `4f5ea39f..19b162b9` (0.68.1 to 0.73.1: the 0.69, 0.70, 0.71,
0.72 and 0.73 rounds), 37 commits.

- `git log -L '/^mvn_cov_models <- /,/^fam_von_mises <- function/:R/families.R' 4f5ea39f..19b162b9`:
  **0 commits** touched the code that moves. Over the whole history
  the same range has 11 commits.
- `git log 4f5ea39f..19b162b9 -- tests/testthat/test-mvn-mixture.R`:
  **0 commits**.
- `git log -G 'mixture_mvn|mvn_cov_|mvn_sim_rows' 4f5ea39f..19b162b9 -- R/`:
  3 commits (7b2bf010, afb1cb27, 83708f40), and each one changed only
  a comment or a roxygen sentence that names the family elsewhere in
  core, not the family.

### Dependencies only it needs

`mclust`, in core's Suggests, used by `test-mvn-mixture.R` alone (the
other mention was the srr sentence in `tests/testthat/helper-reference.R`).
No Imports dependency. The code uses RTMB and stats, which core and
frmtmb.latent both import already.

### Ported brms rows that depend on it

**0.** `grep -c mixture_mvn dev/brmsport-verdicts.tsv` is 0, and no
file of `dev/brmsport-*` or of the generated `test-brms-suite-*.R`
names it. brms 2.23.0 has no multivariate normal mixture
(`dev/mvnlatent-brms-probe.R`, log `dev/mvnlatent-out/brms-probe.txt`):

```
exports matching mvn|mvnorm|multi_normal|mixture: dmulti_normal, mixture, pp_mixture, rmulti_normal
brmsfamily('mixture_mvn') : ERROR: mixture_mvn is not a supported family. ...
brmsfamily('mvnormal') : ERROR: mvnormal is not a supported family. ...
mvbind + mixture, rescor = FALSE : OK
mvbind + mixture, rescor = TRUE : ERROR: Currently, estimating 'rescor' is only possible in multivariate gaussian or student models.
```

brms's nearest spelling, `bf(mvbind(y1, y2) ~ 1)` with
`family = mixture(gaussian, gaussian)`, gives each response its own
class variable: a different model, which frmtmb fits in core and which
does not reach `mixture_mvn()` (`dev/mvnlatent-mvbind-probe.R`: its
log-likelihood, -1219.858, equals the sum of the two univariate
mixtures, -1219.858). With residual correlation brms refuses it.

### Reading of the numbers

The candidate is small (0.41% of `R/`), stable (no change to its code
in five rounds), needs no Imports dependency, carries no brms row, and
costs about 0.3% of the core suite. By guardrail 3's own sentence
("a feature that is small, stable, and needs no extra dependency is a
weak candidate, whatever its domain") this is a weak candidate on the
numbers. What argues for it is the rule's question 3 (a model class
for one domain) and that it already used the extension protocol, so
the move cost one export. The user decides at consolidation.

### No code branch in core named the family

The 38 mentions of `mixture_mvn` in core's `R/` were the family itself
(`R/families.R` 7229-7586), 13 lines of compatibility rules and one
vocabulary entry (`R/compat.R`), and roxygen or comment prose in nine
files. No branch of core code tested for the family: it already reached
core only through generic slots (`fam$mix`, `extra_pars`, `sim_ctx`,
`post$mean_fn`, `post$dpar_response`, `structure`), which `lca()` and
`mixture()` also use. So no new family slot was needed. The one
internal the family read was `us_sigma()`.

## 2. What changed

### Core (frmtmb)

- `R/families.R`: `mvn_cov_models`, `mvn_cov_spec()`, `mixture_mvn()`
  and `mvn_sim_rows()` removed (moved verbatim). `?mixture` gains a
  paragraph that names `mixture_mvn()` in frmtmb.latent and how to load
  it. The prose of `?frmtmb_family` (Structured simulators),
  `?mixture_probs` and `mixture_posterior()` names the family with its
  package.
- `R/compat.R`: the vocabulary entry `mixture_mvn` (kind `structure`)
  and its 9 rules removed; frmtmb.latent registers them.
- `R/covstruct.R`: **new export `us_sigma(theta, d, floored = FALSE)`**,
  the seam (below).
- Prose only: `R/brms-postfit.R` (`?pp_mixture`),
  `R/conditional-effects.R` (comment), `R/predict.R`
  (`?simulate.frmtmb_fit`, a `@noRd` block), `R/priors.R` (`@noRd`),
  `R/simulate-new.R` (`?frm_simulate`), `R/simulate-newdata.R`
  (`@noRd`), `R/structure.R` (`?frmtmb-extension-api`). Every remaining
  mention of `mixture_mvn` in core `R/` is roxygen or comment prose that
  names frmtmb.latent; no code line names it, which
  `test-core-boundary.R` now enforces.
- `DESCRIPTION`: mclust removed from Suggests. `NAMESPACE`:
  `export(mixture_mvn)` out, `export(us_sigma)` in. `man/`: regenerated
  (`mixture_mvn.Rd` deleted, `us_sigma.Rd` new). `_pkgdown.yml`:
  `mixture_mvn` out, `us_sigma` in after `frmtmb-robust-dpars`.
  `README.md`: the family is described under frmtmb.latent. `NEWS.md`:
  a development-version section, BREAKING.
- Tests: `test-mvn-mixture.R` deleted (moved). `test-quadrature-defects.R`
  and `test-structure.R` lose their mixture_mvn blocks (moved).
  `test-unpinned-seams.R` and `test-prior-location-dpar.R` keep their
  core seam pinned on a family built in the test (a `sim_ctx`-only
  family with a matrix mean; a two-location family brms does not have),
  and their mixture_mvn versions moved. `test-core-boundary.R` polices
  `mixture_mvn`, `mvn_cov_spec` and `mvn_sim_rows` in core code, its
  positive control sees them in frmtmb.latent, and it holds the core arm
  of guardrail 1. `test-us-sigma.R` is new. `helper-reference.R` drops
  mclust from its srr sentence.

### The seam: `us_sigma()`

(Fix round 1 corrected this paragraph's first claim and added the
inverse `us_theta()`; see section 9.) The covariance a `us()` parameter
vector encodes (log SDs, then the scaled-Cholesky correlation entries
in column-major order of `L`'s lower triangle, which is glmmTMB's order
only up to d = 3),
AD-safe. It is useful on its own: any family with a free covariance
matrix as an extra parameter needs exactly this, and a copy would fork
the parameterization that `us()` blocks and `VarCorr()` read (the new
test checks `us_sigma(theta, 2, floored = TRUE)` against `VarCorr()` of
a `(x | g)` fit). Documented on its own help page. Two changes to the
function that was internal:

- It refuses a `theta` of the wrong length, and a bad `d`, by name.
  Before, a wrong length was recycled with a warning and a matrix came
  back (`dev/mvnlatent-us-sigma-base.R`, log
  `dev/mvnlatent-out/us-sigma-base-length.txt`: lengths 5 at d = 3,
  4 at d = 2 and 7 at d = 3 each gave "number of items to replace is
  not a multiple of replacement length").
- `floored = TRUE` now floors a 1 x 1 block too; it returned
  `exp(2 theta)` whatever `floored` said. No caller reached it:
  `gr_cov` and `gr_prec` handle d = 1 before they call, and
  `mixture_mvn()` needs D >= 2.

Every caller passes the exact length (`gr_cov`, `gr_prec`, and the EEE
and VVV models of `mixture_mvn()`), so neither change moves a fit; see
section 4.

`test-us-sigma.R` was seen to fail on rellib-r11 behaviorally (the test
environment reaches the internal `us_sigma()` there), log
`dev/mvnlatent-out/seen-fail-base.txt`:

```
DETAIL FAIL [us_sigma() refuses a vector of the wrong length by name] Expected `us_sigma(numeric(5), 3)` to throw a error.
DETAIL FAIL [us_sigma() refuses a vector of the wrong length by name] Expected `us_sigma(numeric(4), 2)` to throw a error.
DETAIL FAIL [us_sigma() refuses a vector of the wrong length by name] Expected `us_sigma(0, 0)` to throw a error.
DETAIL FAIL [us_sigma() refuses a vector of the wrong length by name] Expected `us_sigma(numeric(3), 2.5)` to throw a error.
DETAIL ERROR [us_sigma() refuses a vector of the wrong length by name] Error in `if (d == 1L) return(RTMB::matrix(exp(2 * theta[1]), 1, 1))`: missing value where TRUE/FALSE needed
DETAIL FAIL [us_sigma(floored = TRUE) keeps a vanishing variance positive] Expected `as.matrix(us_sigma(-400, 1, floored = TRUE))[1, 1]` > 0. | Actual comparison: 0.0 <= 0.0 | Difference: 0.0 <= 0
RESULT test-zz-lane-us-sigma.R pass=12 fail=5 err=1 skip=0 warn=2 tests=5 elapsed=0.8
```

The new boundary assertions were seen to fail on the base tree (same
log):

```
DETAIL FAIL [no core file names a structured family or an ODE symbol] ... 'compat.R|mixture_mvn', 'families.R|mixture_mvn', 'families.R|mvn_cov_spec', 'families.R|mvn_sim_rows' ...
DETAIL FAIL [the boundary scanner still sees a token where one exists] Expected `any(grepl("\\|mvn_sim_rows$", names(control)))` to be TRUE.
```

A first version of the length check was wrong: `d * (d - 1L) %/% 2L`
binds as `d * ((d - 1) %/% 2)`, so a 2 x 2 block "took 2". The string
probe caught it on the first install ("us_sigma(): a 2 x 2 matrix
takes 2 parameters, not 3"), and the test pins the counts at d = 2
and 3.

### frmtmb.latent

- `R/mixture-mvn.R` (new): the four objects, code byte-identical. The
  roxygen of `mixture_mvn()` is reworded: no spaced hyphens; the stale
  "`simulate()` is not supported" replaced by what simulate does; links
  qualified with `frmtmb::`; examples lead with
  `family = mixture_mvn(...)`; a paragraph on brms's nearest spelling.
  `mixture_mvn_compat_rules()` holds the 9 rules, text unchanged.
- `R/zzz.R`: registers `c(mixture_mvn = "structure")` with its rules
  FIRST, because the `hmm()` and `lca()` rules name `mixture_mvn`, and a
  rule side must resolve when it is registered
  (`dev/mvnlatent-register-order.R`, log
  `dev/mvnlatent-out/register-order-probe.txt`: "the rule `REML x
  mixture_mvn` names 'mixture_mvn', and the registry has no such
  feature").
- `R/frmtmb.latent-package.R`: imports `us_sigma`. `NAMESPACE`:
  `export(mixture_mvn)`, `importFrom(frmtmb, us_sigma)`. `man/`:
  `mixture_mvn.Rd` new; `lca.Rd` (its link now resolves locally) and the
  package page regenerated.
- `DESCRIPTION`: Description names three families; mclust in Suggests.
  The Depends floor is NOT changed (no version numbers); see section 6.
- Tests: `test-mixture-mvn.R` (new, 14 blocks) is frmtmb's
  `test-mvn-mixture.R` verbatim except two lines, plus the moved blocks
  and two new ones (the family is this package's and not core's; its
  compatibility rows arrive with the package). `diff
  --strip-trailing-cr` against
  `git show 19b162b9:tests/testthat/test-mvn-mixture.R` (log
  `dev/mvnlatent-out/moved-test-diff.txt`), added lines not shown:

  ```
  2a3,43
  250c291
  <     v <- unlist(make_start(fit$frame, st), use.names = FALSE)
  305c346
  <     c(log(sqrt(diag(S))), us_theta_cor(stats::cov2cor(S)))
  381a423,540
  ```

  Those two lines read core internals. They call local helpers now:
  `mvn_start_vector()`, which refuses a start that leaves a component
  out and checks its length against `fit$obj$par`, and
  `mvn_theta_cor()`. `test-prior-dpar.R` gains the moved prior-rule
  block, on that file's existing helpers.
- `vignettes/latent.Rmd`: the title and the introduction name three
  families; a new section clusters Old Faithful and cross-checks the
  log-likelihood with mclust. It renders and tangles
  (`dev/mvnlatent-render.R`): frmtmb -1130.264, mclust -1130.264.
- `dev/out-of-tree-inventory.md`: a section for the third family.
- `NEWS.md`: a development-version section.
- `.github/workflows/check-frmtmb-latent.yaml`: comments only. It
  installs Suggests from DESCRIPTION, so mclust arrives by itself.

### frmtmb.sample

It does simulate and sample mixture_mvn fits: `test-simulators.R`
"mixture_mvn simulates its class covariances" runs `simulate()`,
`posterior_predict(ml_draws())` and `frm_simulate()`. frmtmb.latent was
already in its Suggests. The block now calls
`skip_if_not_installed("frmtmb.latent")` and
`frmtmb.latent::mixture_mvn()`, as its `hmm()` and `lca()` blocks do.
Its help links `[frmtmb::mixture_mvn()]` (three in
`R/methods-draws.R`) point to `frmtmb.latent::`; `man/` regenerated;
a `NEWS.md` bullet. No code change.

## 3. Guardrails

### Guardrail 1, no silent gaps

No string spelling reaches the family, before or after. Probe
`dev/mvnlatent-probe-strings.R`, logs
`dev/mvnlatent-out/probe-strings-base.txt` and `-lane.txt`:

| arm | `family = "mixture_mvn"` | `brmsfamily("mixture_mvn")` | `mixture_mvn(K = 2, D = 2)` |
|---|---|---|---|
| base, frmtmb only | not a supported family | not a supported family | OK |
| base, with latent | not a supported family | not a supported family | OK |
| lane, frmtmb only | not a supported family | not a supported family | could not find function |
| lane, with latent | not a supported family | not a supported family | OK |

So `extension_family_homes` gets no entry. brms has no such name. An
entry would not mend a gap; it would ADD spellings that never existed:
`family_from_name()` passes named arguments through, so
`brmsfamily("mixture_mvn", K = 2, D = 2)` would start to work, and
`family = "mixture_mvn"` would give the constructor's K/D refusal.
(The first pass of this record named only the second; the reviewer
pointed out the first.) Neither is brms vocabulary, so neither is
owed. Both arms are tested: core's `test-core-boundary.R` and latent's
`test-mixture-mvn.R`. Since fix round 1 the core block skips only its
last line, `exists()`, and only when frmtmb.latent is ATTACHED; the
two string refusals run in every session, the preload one included.
The user may want the core-boundary.md row for this move to read
"none; no string spelling reached it" instead of naming an
`extension_family_homes` entry (finding 3, the coordinator's).

**Decision: no exported stub in core.** The reasons:

1. No brms user types it. brms 2.23.0 has no multivariate normal
   mixture. Its nearest spelling (`mvbind()` with a `mixture()` family)
   is a different model, which stays in core and works
   (`dev/mvnlatent-mvbind-probe.R`). The wiener precedent is for a brms
   STRING that brms users do type; here there is none.
2. A stub would be masked. frmtmb.latent Depends on frmtmb, so
   attaching it puts its `mixture_mvn` ahead of core's and prints "The
   following object is masked from 'package:frmtmb'" at every
   `library(frmtmb.latent)`, for every user, to help a user who did not
   load it. `frmtmb::mixture_mvn` would then be a function that only
   refuses.
3. The precedent inside this project: `hmm()` and `lca()` left core with
   no stub.
4. Users still find it: `?mixture` names the family and its package,
   the README lists it under frmtmb.latent, and NEWS says which package
   to load.

### Guardrail 2, interface before move

`grep -rn 'frmtmb:::\|getFromNamespace\|asNamespace'` over
frmtmb.latent's `R/`, `NAMESPACE`, `DESCRIPTION` and `vignettes/`: no
output (log `dev/mvnlatent-out/grep-internals-latent-R.txt` is empty).
Over its `tests/` (log `grep-internals-latent-tests.txt`):

```
tests/testthat/test-hmm-starts.R:389:  f <- getFromNamespace("hmm_starts_scale", "frmtmb.latent")
tests/testthat/test-mixture-mvn.R:514:  expect_false(exists("mixture_mvn", envir = asNamespace("frmtmb"),
tests/testthat/test-prior-dpar.R:19:  des <- frmtmb:::prior_design(formula, data, family, list())
tests/testthat/test-prior-dpar.R:20:  r <- frmtmb:::resolve_priorlist(des, pl)
tests/testthat/test-prior-dpar.R:23:    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
tests/testthat/test-structure-latent.R:45:  ns <- asNamespace("frmtmb.latent")
tests/testthat/test-structure-latent.R:50:  core <- ls(asNamespace("frmtmb"), all.names = TRUE)
```

All were there before except line 514, which asserts that the family
is ABSENT from core. The moved block of `test-prior-dpar.R` uses that
file's pre-existing `pd_reached()`, which reads frmtmb internals.
`test-structure-latent.R`'s codetools check ("both families reach
frmtmb only through its exported seams") covers `R/mixture-mvn.R` and
passes. New export in core: `us_sigma()` only.

### Rule 4, nothing is lost

`frm_compat()` and `frm_compat_features()` dumped on both builds
(`dev/mvnlatent-compat-dump.R`), compared as sets of rows with each
pair oriented the same way (`dev/mvnlatent-compat-compare.R`, log
`dev/mvnlatent-out/compat-compare.txt`). The vocabulary order moves (a
contributed feature comes last), so the 19 pairs of mixture_mvn with a
feature after it in core's list change orientation. Nothing else moves;
statuses and notes are compared:

```
=== all (core and all seven extensions loaded)
table rows before 11350 after 11350
rows only before: 0
rows only after: 0
vocabulary before 164 after 164 only before 0 only after 0
SETS EQUAL: TRUE
=== latent (core and frmtmb.latent)
table rows before 7735 after 7735
rows only before: 0
rows only after: 0
vocabulary before 134 after 134 only before 0 only after 0
SETS EQUAL: TRUE
=== core (before minus mixture_mvn rows)
before rows naming mixture_mvn, dropped: 131
table rows before 7339 after 7339
rows only before: 0
rows only after: 0
vocabulary before 131 after 131 only before 0 only after 0
SETS EQUAL: TRUE
```

The export (`mixture_mvn`), its help topic and its tests have a new
home. The methods it uses (`fitted`, `predict`, `simulate`,
`mixture_probs`, `pp_mixture`, `conditional_effects`) stay in core and
work on its fits (section 4). Core had no vignette section on it; the
latent vignette gains one.

## 4. Fits are unchanged

`dev/mvnlatent-fitcompare.R` fits every mixture_mvn model of the moved
tests, on the data those tests build with the same seeds: 59 fits
(direct ML, faithful, covariates, overrides and gating, guards, 8
template sizes, 8 shared-point starts, 16 "covariates under every
model" with `restarts = 3`, 16 mclust-started, and the quadrature,
unpinned-seams and frmtmb.sample fits). For each it saves logLik, df,
`unlist(estimates)`, `opt$par`, `sqrt(diag(vcov()))`, `fixef()`,
`predict()`, `fitted()`, `mixture_probs()`,
`simulate(nsim = 2, seed = 11)` and the convergence code. It ran twice
on rellib-r11 (frmtmb 0.73.0, the family from frmtmb) and twice on the
lane build (the family from frmtmb.latent), each in its own process;
compared with `dev/mvnlatent-fitcompare-diff.R`, log
`dev/mvnlatent-out/fits-compare.txt`:

```
=== base1 vs lane1
fits 59 ; quantities compared 649 ; identical 649
character results equal in every fit: TRUE
logLik range: -2730.0258619  -249.0267112
=== base2 vs lane2
fits 59 ; quantities compared 649 ; identical 649
=== base1 vs base2
fits 59 ; quantities compared 649 ; identical 649
=== lane1 vs lane2
fits 59 ; quantities compared 649 ; identical 649
```

`identical()` holds on every quantity: bitwise equal. None of the
compared quantities was a refusal (`dev/mvnlatent-fitcompare-kinds.R`,
log `fits-kinds.txt`): se, fixef, predict, fitted, probs and simulate
are values in 59 of 59 fits, the 321 SEs are all finite, and the
convergence code is 0 in 59 of 59. The code is byte-identical and the
results are bitwise equal on reference BLAS, so no OpenBLAS pass was
run: the move adds no numerical path and no branch that rounding
decides.

## 5. Tests: before and after, every one accounted for

Every suite one file per process, 8 at a time, `NOT_CRAN=true`.
Before: rellib-r11 with the 19b162b9 test trees (log
`dev/mvnlatent-out/before.log`). After: the lane library ahead of
rellib-r11 with the worktree's test trees, rerun in full after fix
round 1 by the unedited driver (`bash dev/mvnlatent-par.sh after ...`,
log `dev/mvnlatent-out/after.log`, 423 of 423 RESULT lines). All 423
`lib:` lines name frmtmb 0.73.1 from `wt-mvnlatent-lib`. (The first
pass's `after.log` had one line, `test-core-boundary.R`, spliced in
from a separate rerun of that file; this run replaces that log
entirely, so no line of it is edited.) Generated by
`dev/mvnlatent-count.R` (log `dev/mvnlatent-out/count.txt`):

```
BEFORE
             pkg files  pass fail err skip warn tests elapsed
          frmtmb   262 18846    1   0  202    0  2742  2045.6
 frmtmb.coupling    11   551    0   0    5    0    77    25.0
      frmtmb.eam    32  1838    0   0    5    0   337   464.7
   frmtmb.latent    11   368    0   0    2    0    77    55.6
    frmtmb.learn    15   430    0   0   13    0   112    49.9
      frmtmb.ode    11   549    0   0    1    0   139    75.2
   frmtmb.sample    64  3022    0   0    9    0   454   393.4
   frmtmb.spline    16   594    0   0    1    0   115    51.8
AFTER
             pkg files  pass fail err skip warn tests elapsed
          frmtmb   262 18617    0   0  202    0  2741  1471.4
 frmtmb.coupling    11   551    0   0    5    0    77    25.0
      frmtmb.eam    32  1838    0   0    5    0   337   485.5
   frmtmb.latent    12   679    0   0    2    0    92    72.6
    frmtmb.learn    15   430    0   0   13    0   112    47.4
      frmtmb.ode    11   549    0   0    1    0   139    73.1
   frmtmb.sample    64  3022    0   0    9    0   454   411.1
   frmtmb.spline    16   594    0   0    1    0   115    51.8

files only before: frmtmb test-mvn-mixture.R
files only after: frmtmb test-us-sigma.R, frmtmb.latent test-mixture-mvn.R

files whose counts changed:
           pkg                         file tests.b tests.a pass.b pass.a
        frmtmb                test-compat.R      36      36    772    771
        frmtmb         test-core-boundary.R       2       3      4     10
        frmtmb test-dirichlet-multinomial.R      13      13     66     67
        frmtmb           test-mvn-mixture.R       9      NA    277     NA
        frmtmb    test-quadrature-defects.R       7       6     43     40
        frmtmb             test-structure.R      28      28    135    128
        frmtmb              test-us-sigma.R      NA       8     NA     52
 frmtmb.latent           test-mixture-mvn.R      NA      14     NA    308
 frmtmb.latent            test-prior-dpar.R       4       5     10     13
(fail, err, skip and warn unchanged in all of them, except
test-dirichlet-multinomial.R, fail 1 -> 0)

failing, erroring or warning files after: 0 rows
```

The core elapsed sum is 1471.4 s in this run, 1678.9 s in the first
pass and 2045.6 s before: load, not the change (section 1).

The one failure before is not this lane's: rellib-r11 is frmtmb 0.73.0
and the test tree is 0.73.1's, whose `fitted()` count clamp
(`R/predict.R`) that test pins ("Expected `ft[3, "Q97.5", 2]` to be
identical to 50 ... actual 49.9999999999999929").

Expectations, core: 18847 before (18846 pass and 1 fail), 18617 after,
-230 =

- -277, `test-mvn-mixture.R` moved to frmtmb.latent;
- -3, `test-quadrature-defects.R` "mixture_mvn() refuses REML and
  profile too", moved;
- -7, the mixture_mvn half of `test-structure.R` "a rowwise
  mixture-type family declares capabilities only", moved;
- -1, `test-compat.R` "rule patterns name real kinds, groups, and
  features" asserts once per distinct side of a CORE rule, and no core
  rule names `mixture_mvn` now (156 -> 155 in that block, from the
  TIME lines);
- +6, `test-core-boundary.R`: 2 in the positive control (the scanner
  sees `mixture_mvn` and `mvn_sim_rows` in frmtmb.latent) and 4 in the
  new guardrail-1 block;
- +52, `test-us-sigma.R`, new (18 in the first pass; fix round 1 adds
  the d = 4 order, the Jacobian, the `d` message, the floors and
  `us_theta()`).

`test-unpinned-seams.R` (30) and `test-prior-location-dpar.R` (55) keep
their counts: each block was rewritten on a family built in the test,
and the mixture_mvn version moved.

Expectations, frmtmb.latent: 368 -> 679, +311 = 308 in
`test-mixture-mvn.R` (277 + 3 + 7 + 6 moved from core files, plus 15
in its two new blocks) and +3 in `test-prior-dpar.R` (moved).
frmtmb.sample: 3022 -> 3022; `test-simulators.R` runs its mixture_mvn
block (it does not skip, since frmtmb.latent is installed). test_that()
blocks, core 2742 -> 2741 (-9 -1 +1 +8), latent 77 -> 92 (+14 +1).

So 296 expectations moved (277 + 3 + 7 + 6 + 3), 9 core expectations
were rewritten on built families (6 + 3), none was dropped, and 1 core
expectation disappeared with the rule it iterated over.

### One session per suite, three orders

`FRMTMB_RELLIB=<lane lib> NOT_CRAN=true bash dev/release/run-sessions.sh
mvnlatent-sess 8` (all eight packages), log
`dev/release/mvnlatent-sess.log`:

```
SESSION frmtmb forward files=252/262 pass=18583 fail=0 err=0 skip=202 warn=0
SESSION frmtmb reverse files=252/262 pass=18685 fail=0 err=0 skip=202 warn=0
SESSION frmtmb preload files=252/262 pass=18782 fail=0 err=0 skip=204 warn=0
SESSION frmtmb.coupling forward|reverse|preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward|reverse|preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward|reverse|preload files=12/12 pass=679 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward|reverse|preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward|reverse|preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward|reverse|preload files=59/64 pass=3022 fail=0 err=0 skip=9 warn=0
SESSION frmtmb.spline forward|reverse|preload files=16/16 pass=594 fail=0 err=0 skip=1 warn=0
MVNLATENT-SESS ran 24 of 24
```

The three lines of each extension are identical and are folded here.
fail 0 and err 0 in all 24. Every LEAK line carries `[with load of
...]`, so none is a test's own. One is new: frmtmb.sample's
`test-simulators.R` in reverse order, "compat features: 1 -> 4 entries
[with load of frmtmb.latent]", the latent load that block does, which
now registers three features where it registered two.

### The gated tier

The 44 files that skip in the plain tier (any package), rerun with
`FRMTMB_BRMS_FIT_TESTS=true NOT_CRAN=true` on the lane build, one file
per process (`dev/mvnlatent-gated.sh`, log
`dev/mvnlatent-out/gated-after.log`): 44 of 44 RESULT lines, pass 3386,
fail 0, err 0, warn 0, skip 29. The 29 skips are other tiers' gates,
read from the DETAIL lines: `test-drmtmb-agreement.R` 13
(`FRMTMB_DRMTMB_FIT_TESTS`), `test-fuzz.R` 1 (`FRMTMB_FUZZ`), and the
`test-scale.R` of every package 15 (`FRMTMB_SCALE_TESTS`). No before
arm was run for this tier, since nothing failed.

## 6. Floors, bumps, overlaps

- **Floor of frmtmb.latent on core: the core release that exports
  `us_sigma()`, which is this change.** Measured: this tree's
  frmtmb.latent installed alone into a scratch library over rellib-r11
  (frmtmb 0.73.0) fails at lazy loading, log
  `dev/mvnlatent-out/floor-install.txt`: "Error: object 'us_sigma' is
  not exported by 'namespace:frmtmb' ... ERROR: lazy loading failed for
  package 'frmtmb.latent'". Its DESCRIPTION still says
  `frmtmb (>= 0.72.0)`, because the lane writes no version numbers; the
  consolidation must raise it to the new core version.
  `dev/rel073-floor.sh` measures frmtmb.sample only; the next floor
  script needs a frmtmb.latent arm.
- **The converse breaks harder: frmtmb.latent 0.6.2 does not LOAD over
  this core.** Its `.onLoad()` registers `hmm()` rules that name
  `mixture_mvn`, which core no longer supplies (`dev/mvnlatent-old-latent.R`,
  log `dev/mvnlatent-out/old-latent.txt`): ".onLoad failed in
  loadNamespace() for 'frmtmb.latent' ... the rule `hmm x mixture_mvn`
  names 'mixture_mvn', and the registry has no such feature". So the new
  core and the new frmtmb.latent must be released and installed
  together; R has no way for core to declare that floor on its
  extension. A user who upgrades core alone loses `hmm()` and `lca()`
  as well as `mixture_mvn()` until frmtmb.latent is upgraded. The
  alternative, keeping the `mixture_mvn` vocabulary entry in core, was
  not taken: core would then describe a feature it does not have.
  NEEDS THE USER: accept the lockstep release, or have frmtmb.latent
  name `mixture_mvn` in `expects =` (which would not help 0.6.2,
  already released).
- frmtmb.sample: no floor change (its code is unchanged); its
  `test-simulators.R` needs a frmtmb.latent that has `mixture_mvn()` or
  it fails at that block rather than skipping. A version-aware skip was
  not written; the release installs both together.
- Bumps: core **minor** (an export removed, BREAKING; one export
  added). frmtmb.latent **minor** (new exported family). frmtmb.sample
  **patch** (tests and help links only).
- Overlaps with the other lanes of this round:
  - `R/families.R`: a 359-line deletion around old line 7229 and prose
    at lines 265, 5810, 7158 and 7186. `extension_family_homes` is NOT
    touched (splitsurv adds there).
  - `R/compat.R`: the structure vocabulary line 865 and 9 rule
    deletions (REML, quadrature, profile, cens()/trunc(), the
    mixture_mvn block near 1720, group:autocor). spatial's guardrail-1
    test over the compat vocabulary will find `mixture_mvn` only when
    frmtmb.latent is loaded; it is not a brms name.
  - `R/covstruct.R`: `us_sigma()` exported and its roxygen block (lines
    58-115). spatial builds a covariance-structure seam in this file and
    may export `us()`-parameterization helpers too: one name, one
    implementation, should be agreed at consolidation.
  - `R/predict.R`, `R/simulate-new.R`, `R/simulate-newdata.R`,
    `R/conditional-effects.R`: prose and comments only.
  - `NAMESPACE`, `DESCRIPTION` (mclust out of Suggests), `NEWS.md` (new
    top section), `man/`, `_pkgdown.yml` (two lines), README.
  - frmtmb.sample: `R/methods-draws.R` (three doc links),
    `tests/testthat/test-simulators.R`, `NEWS.md`, two `man/` pages.
  - Not touched: `.github/workflows/test-coverage.yaml`,
    `dev/release/run-sessions.sh` (no new package).

## 7. Not done, and for the user

- The decision itself: on the numbers of section 1 the move is a weak
  candidate (small, stable, no Imports dependency, no brms row, 0.3% of
  the core suite). It is done in full so that the user can keep or drop
  it at consolidation.
- No stub `mixture_mvn()` in core (section 3, argued).
- No OpenBLAS pass (section 4, argued: bitwise-equal fits, no new
  numerical path).
- The moved tests keep their absolute tolerances (`1e-5`, `1e-4`,
  `1e-10`, `0.3`) unchanged, so that the move is visibly verbatim;
  rewriting them as ratios is a separate change.
- `R/methods-draws.R` of frmtmb.sample had lines over 80 columns at
  the three edited links before the edit; they are still over 80.

## 8. R CMD check --as-cran

`bash dev/mvnlatent-check.sh <pkg>` (build with vignettes, then check,
inside `dev/mvnlatent-check/`, since deleted; `R_LIBS` = lane lib, then
rellib-r11, then the user library; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`;
pandoc and TinyTeX on PATH). Logs kept as
`dev/mvnlatent-out/rcmdcheck-<pkg>.txt`:

After fix round 1 (the runs of record for core and frmtmb.latent; the
logs `rcmdcheck-frmtmb.txt` and `rcmdcheck-frmtmb.latent.txt` are from
this run):

```
frmtmb: Status: 1 NOTE          (checking HTML version of manual: V8 unavailable; tests 19m)
frmtmb.latent: Status: 1 NOTE   (the same V8 NOTE; tests 63s)
```

First pass (frmtmb.sample was not changed by fix round 1, so its
check is not repeated):

```
frmtmb: Status: 1 NOTE          (the V8 NOTE; tests 20m)
frmtmb.latent: Status: 1 NOTE   (the V8 NOTE; tests 60s)
frmtmb.sample: Status: OK       (tests 311s)
```

Fix round 1 also reran core and frmtmb.latent in one session in three
orders (`run-sessions.sh mvnlatent-sess2 8 frmtmb frmtmb.latent`, log
`dev/release/mvnlatent-sess2.log`):

```
SESSION frmtmb forward files=252/262 pass=18617 fail=0 err=0 skip=202 warn=0
SESSION frmtmb reverse files=252/262 pass=18719 fail=0 err=0 skip=202 warn=0
SESSION frmtmb preload files=252/262 pass=18819 fail=0 err=0 skip=203 warn=0
SESSION frmtmb.latent forward|reverse|preload files=12/12 pass=679 fail=0 err=0 skip=2 warn=0
MVNLATENT-SESS2 ran 6 of 6
```

Every LEAK line carries `[with load of ...]` (frmtmb.sample, frmtmb.eam,
the sf stack), so none is a test's own. The preload session now skips
one test more than forward, not two: `test-core-boundary.R` runs all
10 of its expectations there (F7), and the one extra skip is the
pre-existing wiener test of `test-families-extreme.R`
(`isNamespaceLoaded("frmtmb.eam") is TRUE`).

A correction to the record: after the fix round a report said that
the core check and the core sessions had stopped partway through.
They had not. They were still running; the report read logs that were
being written, an empty `check-frmtmb.txt` (written only at the end),
and a `check-latent.txt` left from the first pass, which the latent
check overwrites only when it starts. The process list and free
memory (44109 of 64907 MB) at the time are in
`dev/mvnlatent-out/stop-evidence.txt`. No run was stopped or
restarted, and nothing points to another lane.

## 9. Fix round 1 (review `dev/reviews/2026-10-23-mvnlatent.md`)

- **F1, the order of `us_sigma()`.** `?us_sigma` no longer claims
  glmmTMB's parameterization. It states the order: the correlation
  entries fill `L`'s strict lower triangle by column (`L[2, 1]`,
  `L[3, 1]`, ..., `L[d, 1]`, `L[3, 2]`, ...), which equals glmmTMB's map
  up to d = 3 and differs from d = 4, where glmmTMB fills by row. New
  test "the correlation entries fill L's lower triangle by column":
  at d = 4 with six distinct nonzero entries, against a reference that
  assigns each cell by name (no `lower.tri()`), plus a guard that the
  by-row build differs, so the assertion can fail.
- **F2, the inverse.** New export `us_theta(S)` on the same help page
  (`@rdname us_sigma`): the `theta` whose `us_sigma()` is `S`, numeric
  only. It refuses a non-square, non-finite, asymmetric (judged
  against `1e-10` of the matrix's own scale) or non-positive-definite
  `S`, and drops dimnames. Test "us_theta() inverts us_sigma()" round
  trips both ways at d = 1 to 5 (seed 21), and a refusal test covers
  each case. latent's `test-mixture-mvn.R` drops its copy
  `mvn_theta_cor()` and starts the mclust fits at `us_theta(S)`; it
  passes (308). The hunk stays inside `us_sigma()`'s block of
  `R/covstruct.R`, so a per-file merge with the spatial lane's seam
  should resolve it. Its help page links to the extension pieces it
  serves; a link to the spatial lane's `?frmtmb_register_covstruct`
  is left for consolidation, since that topic is not in this tree.
- **F4, the record.** Section 1 gives the net change to core `R/`
  (+101 / -417, net -316, after the first pass; +182 / -417, net -235,
  now) and the time-share caveat (one loaded measurement; the same 262
  files summed to 2045.6 s, 1678.9 s and 1471.4 s in three runs).
  Section 5 is regenerated from a full rerun whose `after.log` no
  script or hand edited. The preload session's skip=204 against 202:
  from the first pass's session logs, the two extra skips were
  `test-core-boundary.R` "mixture_mvn() is not in core, and no string
  spelling reaches it" (`isNamespaceLoaded("frmtmb.latent") is TRUE`,
  this lane's, fixed by F7) and `test-families-extreme.R` "family =
  "wiener" names frmtmb.eam when it is not loaded"
  (`isNamespaceLoaded("frmtmb.eam") is TRUE`, pre-existing and correct:
  that refusal does depend on eam being absent). Logs
  `dev/mvnlatent-out/skip-forward-old.txt` and `skip-preload-old.txt`.
  Section 3 corrects the `extension_family_homes` argument (F3's half
  that is the lane's).
- **F6.** The `mclust` entry is removed from `codemeta.json` (12
  lines); the file still parses (`jsonlite::fromJSON()`), and no
  `mclust` is left in it.
- **F7.** The guardrail-1 block of `test-core-boundary.R` has no skip
  before its string refusals; only `exists("mixture_mvn")` is skipped,
  when `package:frmtmb.latent` is on the search path.
- **F8.** The tape test now compares the AD Jacobian of
  `as.vector(us_sigma(p, 3))` with central differences (h = 1e-6), as
  the ratio max |J - Jfd| / max |J| < 1e-6. `us_sigma(theta, 1,
  floored = TRUE)` now returns `sd2_floored(theta)`, `exp(2 theta) +
  1e-300`, which is what the 1 x 1 random-effect blocks use; it was
  `(exp(theta) + 1e-150)^2`. No caller in frmtmb passes d = 1 with
  `floored = TRUE` (`gr_cov`, `gr_prec` branch before), so no fit
  moves. The test pins the value with `expect_identical()` against
  `sd2_floored(-400)`, and the d = 2 diagonal against
  `sd_floored(-400)^2`. The help text says which floor applies when.
- **F9.** A numeric `d` that is not whole prints with 15 significant
  digits: `us_sigma(numeric(3), 2 + 1e-12)` says "not 2.000000000001";
  the test pins it.
- **F10.** `?frmtmb-extension-api` gains a section, "The mixture
  component interface": a `mix` element makes a family a mixture to
  frmtmb (what that switches on: `mixture_probs()`, `pp_mixture()`,
  the several-location prior rule, the `ar()` and `incl_thres`
  refusals), and the three fields `mixture_posterior()` reads: `K`,
  `log_pi(dpars)` and `comp_lpdf(y, dpars, aterms, k, extra)` with the
  fifth argument passed only when the family declares `extra_pars`.
  It says that `ord`, `ref`, `type`, `comp_dpars`, `comp_sim` and the
  bookkeeping are read only for `mixture()` itself. `?mixture_probs`
  points to it. The field list was read from core's code
  (`grep` of `[["mix"]]` reads in `R/`).
- **F11.** The 87-column test title is now "method = 'predict' refuses
  a draw of the whole response". The latent vignette drops "The family
  moved here from frmtmb, and its fits did not change in the move.",
  and its mclust assertion is a ratio: |logLik difference| / |mclust
  logLik| < 1e-7. Measured on the rendered fit
  (`dev/mvnlatent-vig-ratio.R`, log `dev/mvnlatent-out/vig-ratio.txt`):
  frmtmb -1130.2639601854, mclust -1130.2639601859, relative gap
  5.146e-13. The vignette renders and tangles. Both NEWS files say
  "on the same machine". Left as they are, as asked: the moved
  `test-prior-dpar.R` block reaches three core internals through that
  file's pre-existing `pd_reached()` (`frmtmb:::prior_design`,
  `resolve_priorlist`, `par_template_names`); and `predict(newdata =)`
  refuses this family with a message ("a hidden state sequence, a
  group-level latent class, or a correlated residual",
  `R/predict-brms.R`) none of whose three cases is a multivariate
  normal mixture, which was so before the move.
- Not done, as asked: F5 (the floor raise, consolidation's) and F3's
  document change (the coordinator's).
- Fits after the fixes: `dev/mvnlatent-fitcompare.R` on the lane build
  again, against the base save: 59 fits, 649 of 649 quantities
  identical (log `dev/mvnlatent-out/fits-compare-fix1.txt`).
