# Three recorded debts, settled

Items 2, 3 and 9 of "Follow-ups carried out of the 0.54.0 round" in
`dev/feature-gaps.md`. They are unrelated; each section stands alone.

Every number below was measured in this worktree against a private
library, one test file per R process, `NOT_CRAN=true`.

---

## Item 2. The capped importance-correction warning

### What was wrong

`fit$importance$capped` is TRUE for two fits that want OPPOSITE advice,
and until now both got the same sentence.

* A correction that is SHORT OF ROUNDS shrinks its move by a factor of
  three to ten each round. One or two more rounds land it, and "raise
  `frmtmb_control(importance_rounds =)`" is right.
* A STALLED correction takes the SAME step every round. Its total shift
  is that step times the round count, so a larger cap buys a
  proportionally larger number rather than a better one, and the advice
  is wrong.

Before, on both (frmtmb 0.54.0, private library):

```
== stalled (12 groups of 3 Bernoulli rows, sd = 0, 50 draws)
moves 0.90907826 0.90907292 0.90907839 0.90907821 0.90907291
WARNING: The importance correction used all 5 of its rounds and the
estimates were still moving by 0.909 at the last one. Raise
frmtmb_control(importance_rounds =), or raise the draw count so each
round lands in the same place

== short of rounds (scalar Bernoulli, cap 2, 200 draws)
moves 0.109326438 0.010590534
WARNING: The importance correction used all 2 of its rounds and the
estimates were still moving by 0.0106 at the last one. Raise
frmtmb_control(importance_rounds =), or raise the draw count so each
round lands in the same place
```

After:

```
== stalled
WARNING: The importance correction moved by the same amount, 0.909, in
every one of its 5 rounds, so it has not converged and its total shift
is that step times the round count rather than an estimate. The step is
a property of the draws and not of the data, so raising
frmtmb_control(importance_rounds =) would only move the estimates
proportionally further. A variance component the Laplace fit has
already collapsed does this, because nothing is left to reweight: check
VarCorr() before reading the corrected estimates

== short of rounds
WARNING: (unchanged)
```

### The step is a property of the draws, not of the model

The recorded debt calls the step a "cap". There is no step cap in the
code: `optimize_obj()` bounds nothing, and `imp_move_tol` is a stopping
tolerance, not a limit. The constant step is the fixed-point map
becoming a pure translation. When the Laplace fit has collapsed a
variance component the likelihood carries no information about the
random effect, so the frozen-proposal objective at anchor `s` is
maximized at `s` times a constant determined by the DRAWS, and every
round multiplies the standard deviation by that same constant.

Measured, and this is what the new sentence asserts:

| what varies | mean step over the 5 rounds |
|---|---|
| `importance_seed = 1`, 8 subjects, 50 draws | 0.9496142 |
| `importance_seed = 7`, same data | 0.3179309 |
| `importance_seed = 99`, same data | 0.2609419 |
| a DIFFERENT dataset, seed 1, same size | 0.9496088 |
| core Bernoulli, 8 groups of 3, seed 1, 50 draws | 0.9496129 |

The last two rows are the point. Two datasets at one seed differ by
5.4e-06 (5.7e-06 relative), and a core Bernoulli fit and a
`frmtmb.learn` `bandit2arm_delta()` fit, different packages, different
families, different data, differ by 1.3e-06 when they share eight
groups, fifty draws and one seed. Both are smaller than the spread
BETWEEN ROUNDS of a single fit (1.8e-05 relative here), so the step
varies less across datasets and packages than it does across the
rounds of one run. A dataset whose component is identified does not
stall at all: the same 8-subject design at 60 trials instead of 40
moves by 0.0155 with a relative spread of 3.16.

### The threshold, and both regimes measured

`imp_stalled()` compares the spread of the moves to their own size,
`(max - min) / mean`, so no parameter scale enters it.

Capped, and genuinely still moving (5 fits):

| design | rounds | moves | spread |
|---|---|---|---|
| probe 2x2 Bernoulli, cap 2 | 2 | 0.0958, 0.0257 | 1.15 |
| gaussian mu-and-sigma, cap 2 | 2 | 0.1087, 0.0212 | 1.35 |
| scalar Bernoulli, cap 2 | 2 | 0.1093, 0.0106 | 1.65 |
| probe 2x2 Bernoulli, cap 3 | 3 | 0.0958, 0.0257, 0.0072 | 2.06 |
| 20 groups of 3, sd = 0, 50 draws | 5 | 0.364 ... 0.287 | 0.322 |

The last row is not converging either, it is wandering; the moves are
not equal, so it keeps the old message, which is the honest one for it.

Stalled (18 fits: 6 reversal fits from `frmtmb.learn` at 100 and 400
draws, 4 `ln_small` fits at 24 to 200 draws, and 8 core Bernoulli fits
simulated with sd = 0), relative spread from **5.4e-06 to 2.6e-04**.

The two regimes are more than three orders of magnitude apart.
`imp_stall_tol = 1e-2` sits 39 times above the widest stalled spread
and 32 times below the narrowest capped-but-moving one. It is loose on
purpose: a fit whose moves really agreed to one percent while still
shrinking would have a per-round factor of 0.9975, and from a move of
0.36 that needs about 2400 more rounds to reach `imp_move_tol`, where
"raise the round count" is not honest advice either.

One round is never called stalled: a single move is trivially equal to
itself, and a cap of one round is exactly where raising the cap is
right.

### Changed

* `R/importance.R`: `imp_stall_tol` and `imp_stalled()`, both `@noRd`,
  beside `imp_move_tol` and `imp_ess_floor` where the measured
  constants live.
* `R/fit.R`: the `capped` branch of `check_convergence()` chooses
  between two sentences; `@param importance_rounds` says the two cases
  differ.
* `tests/testthat/test-importance.R`: three tests, 23 assertions.
* `extensions/frmtmb.learn/tests/testthat/test-factorization.R`: the
  8-subject fit is a stalled one, so its `expect_warning()` moves to
  the new text, plus an assertion that its moves agree to a part in ten
  thousand.

### Shown failing first

Against the unmodified installed 0.54.0, `test-importance.R` gives
`PASS=203 FAIL=0 ERROR=TRUE`, with all three new tests erroring on
`object 'imp_stalled' not found`, and the two message assertions could
not be reached. The message half is shown directly by the before and
after above, on the same fixture.

### Not done

`sd(alpha)` collapse detection is not added to the refusal path. The
warning names `VarCorr()` and stops there, because a collapsed
component is a property of the Laplace fit that the user can read
themselves, and a refusal on it would fire on the one 20-trial dataset
in six whose component survived.

---

## Item 3. The spline span refusal keyed on the wrong thing

### What the two calls actually held

In `frm_curve_feature()` the block was

```r
parts <- sp_curve_parts(sp$fit, stk, ...)   # stk: 5 points per root
if (length(parts$span)) {
  sp_span_stop(sp_span_on_grid(sp$fit, nd, ...))   # nd: the whole grid
}
```

`parts$span` holds the `ps()` span messages raised by `frm_lp_basis()`
on the five-point STENCIL, which is row 1 replicated with `var` set to
`roots + c(-e2, -e1, 0, e1, e2)`. `sp_span_on_grid()` holds the
messages raised by one `predict()` on the WHOLE grid `nd`. Those are
different questions about different point sets, and the block used the
first as the gate and the second as the message.

The two can disagree in both directions, and both were constructed on
the two-`ps()` fixture in `test-span.R` (`t` and `z`, grid ending
5e-06 of its range inside `t`'s span):

| root | grid | gate `parts$span` | message `sp_span_on_grid()` | 0.54.0 |
|---|---|---|---|---|
| at the bracket end | clean | `t`, the stencil fringe | none | returns 1.29803460501, se 0.14249 |
| at the bracket end | `z` out of span in row 10 | `t`, the stencil fringe | `z` | REFUSES, quoting `z` |
| in the middle | clean | none | none | returns 0.840355589493, se 0.0131381 |
| in the middle | `z` out of span in row 10 | none | `z` | returns 0.840355589493, se 0.0131381, in SILENCE |

Rows 2 and 4 are the defect. The same out-of-span row is refused or
ignored depending on whether the located root happens to land within
`e2` (a ten-thousandth of the grid's range) of a knot in a DIFFERENT
term, and the estimate is bit-identical either way, because the search
holds `z` at row 1's value and never predicts at row 10 at all. When it
does refuse, the sentence "the search bracket leaves a ps() term's knot
span" is false as applied: the bracket is entirely inside `t`'s span,
and `z` is not the bracket.

`frm_curve()` warns on the same grid, so the two doors disagreed about
it as well.

### What it should read, and why

The search evaluates nothing but row 1 with `var` moved: the scan, the
Newton steps and the stencil are all `row1[rep(1L, n), ]` with one
column overwritten. So:

1. If every column other than `var` is PINNED to row 1, the scan has
   already predicted at every covariate combination the grid holds. A
   grid that is dirty in `var` implies a dirty scan (for a crossing the
   scan's points ARE the grid's `var` values; for a stationary point it
   evaluates `x - e1` and `x + e1`, and a value outside the span on one
   side stays outside when moved further that way), so the FIRST
   refusal already covers it.
2. The only excursion neither the scan nor the stencil can see is one
   in another column of another row. That is what the second check
   exists for, and it is a question about the grid, not about the
   stencil.

So the gate is now "does a column other than `var` vary down the
grid", and the message is what the grid itself reports. Trigger and
message are the same question, and the answer no longer depends on
where the root landed.

After, same four cases: rows 1 and 3 unchanged, rows 2 and 4 both
refuse with the message about `z`.

The stencil's own fringe is still not asked about, on purpose: it is
the same doctrine `frm_curve_deriv()` follows in counting the grid and
not the stencil it widens the grid into, and it is what
`test-span.R`'s "a grid ending exactly on a knot is inside the span"
pins.

### Cost, and the false-alarm rate

`sp_grid_pinned()` is `duplicated()` on each column other than `var`,
so a grid whose other columns are constant, which is every grid in this
package's own examples, tests and vignettes, costs one vector pass and
no prediction. The extra `predict()` runs only when a column varies.

The refusal cannot fire on a grid that is inside every span: it refuses
if and only if `predict()` on the grid raises a span message, which is
the same condition `frm_curve()` warns on. On a model with no `ps()`
term nothing is ever raised, which `test-span.R` already pins.

### Changed

* `extensions/frmtmb.spline/R/curve-cov.R`: `sp_grid_pinned()`,
  `@noRd`, beside `sp_span_on_grid()` whose doctrine it extends.
* `extensions/frmtmb.spline/R/curve-feature.R`: the gate, the comment
  that says what it keys on and what it deliberately does not, and the
  "Past a `ps()` knot span" section of the help. The check also moved
  above the stencil build, so the refusal happens before the covariance
  work rather than after it.
* `extensions/frmtmb.spline/tests/testthat/test-span.R`: the two-`ps()`
  section renamed and extended; the middle-root test is new.

### Shown failing first

Against the unmodified installed frmtmb.spline, the new test fails:

```
── 1. Failure ('test-span.R:351:3'): a second ps() term's span refuses
   at a root in the middle
Expected `frm_curve_feature(...)` to throw a error.
```

with `PASS=55 FAIL=1`. After: `PASS=58 FAIL=0`.

### Found, not fixed

For `type = "crossing"` the `+/- e2` stencil points do not enter any
reported number: the estimate uses `blk(3)`, the standard error uses
`C0` and `f1`, and `f2` reaches only the `curvature` attribute. For a
stationary point they DO enter, through `denom <- f2`. So a maximum
located within `e2` of a knot has its curvature, and therefore its
standard error, computed across the partition-of-unity cliff. The
control case above shows what that looks like from outside: a crossing
3e-05 inside the span reports se 0.14249 where the same curve reports
0.0131 in the middle. Refusing it would re-introduce the bug the
grid-not-stencil doctrine exists to prevent, so nothing here changes
it; a future span change should decide whether the STATIONARY case
wants its own answer.

---

## Item 9. The hazard-container guard now runs in each extension

### The decision

The guard MOVES, and the list does not. `R/hazard-containers.R` holds
the container list, moved byte for byte out of the test file (diffed:
identical), and exports one function, `frm_hazard_reads()`. Each
extension asserts the rule on itself in its own
`tests/testthat/test-bracket-access.R`, eight lines of assertion under
a comment and no list of its own, so it fails its OWN `R CMD check`.
Core keeps the source sweep of every extension as a backstop and gains
an assertion that every extension in the tree carries the file.

### Why not a copy per package, measured

The file's own header said a cloned copy "would trade a single
container list for five that drift, which is the worse bargain". There
are seven extensions now, and the repository already runs the
experiment: `test-message-uniqueness.R` is duplicated.

| package | lines | md5 | condition constructors handled |
|---|---|---|---|
| frmtmb | 101 | 18f24106 | yes |
| frmtmb.coupling | 114 | 15eeddaa | yes |
| frmtmb.eam | 80 | 707a37cd | no |
| frmtmb.latent | 80 | 65aee423 | no |
| frmtmb.learn | 86 | 6083708b | no |
| frmtmb.sample | 72 | 2264e52c | no |
| frmtmb.spline | 114 | ced97e15 | yes |
| frmtmb.ode | absent | | |

Seven copies, seven distinct contents, two different parsers, and one
extension with no copy at all. The split is not cosmetic: four copies
do not collect templates raised through `errorCondition()` or
`warningCondition()` and do not skip the `class =` and `call =`
arguments, so a duplicate template raised that way is unasserted there
and nothing says so. Today only `frmtmb.spline` raises any (2 calls),
and its copy is one of the three that handle them, so the hole is
latent rather than active. That is the drift the bracket guard would
have inherited, on a file that also carries 22 reserved names.

### What was possible, tested rather than assumed

* An extension's test CAN call core: every extension has `frmtmb` in
  `Depends` or `Imports`, so `frmtmb::frm_hazard_reads()` resolves.
* Core CAN ship a testing aid. It is one exported function, documented
  as such, with the container list kept internal.
* The source-tree scanner cannot be the extension's guard. It reads
  `../../R`, which does not exist under `R CMD check` on a built
  tarball, where the tests run beside an INSTALLED package; core's own
  guard skips there by design. `frm_hazard_reads()` walks the syntax
  tree of every function in the NAMESPACE instead, so it runs wherever
  the package does.

### The two scanners agree, and finding out fixed one of them

Measured on the pre-sweep `frmtmb.coupling` sources (git 7a01ed2), the
tree the release tally caught:

```
source scanner on the pre-fix frmtmb.coupling R/: 34 reads
namespace scanner on the same code: 34 reads
same set of container$slot pairs? TRUE
```

The scan costs 0.03 to 0.22 s per extension (20 to 170 namespace
objects) and 1.78 s on core's 982, so an extension pays for its own
guard in a tenth of a second.

34 is the number `dev/feature-gaps.md` records. The first version of
the namespace walker found 33: it missed `coupling.R:185`,
`fit$frame$data`, because it only looked at the outer expression's
left-hand SYMBOL, and there the outer `$`'s left side is a call. The
rule this project documents is that a chain is judged one link at a
time by its own left-hand token, so `fit$frame$data` is a hit on
`frame` and `fit[["frame"]]$data` is not. The walker now takes the
inner `$`'s slot name as the left-hand name, which reproduces the
source scanner exactly on this tree, and both the equivalence and the
chained idiom are pinned in `test-bracket-access.R`.

### An extension fails its own check, demonstrated

`sp_debts_probe <- function(dpars) dpars$mu` added to
`extensions/frmtmb.spline/R/curve-cov.R`, spline reinstalled:

```
── Failure ('test-bracket-access.R:20:3'): no frmtmb.spline function
   reads a hazard container
Expected `hits` to be identical to `character(0)`.
`actual`:   "sp_debts_probe: dpars$mu"
Use [["name"]], or rename the local: these names are reserved.
```

and in the package's OWN `R CMD check` of its OWN tarball:

```
[ FAIL 1 | WARN 0 | SKIP 1 | PASS 328 ]
Error: ! Test failures.
Status: 1 ERROR
```

The probe was then removed (`git diff` on the file shows only the
26-line helper) and spline reinstalled and rechecked clean.

### Changed

* `R/hazard-containers.R` (new): the list, moved verbatim, and
  `frm_hazard_reads()`.
* `NAMESPACE`, `man/frm_hazard_reads.Rd`, `_pkgdown.yml`: one export.
* `tests/testthat/test-bracket-access.R`: reads the list from
  `frmtmb:::hazard_containers`; three new tests, including the
  every-extension-carries-the-file assertion.
* `extensions/*/tests/testthat/test-bracket-access.R` (7 new files, 25
  lines each, identical but for the package name).
* `dev/bracket-sweep.md`: "The lint" says what the arrangement is now.

### User-visible? Mostly not

The seven extension files and core's test changes are internal. The one
user-visible part is the new export, which is in core's NEWS.

### Left for consolidation

The seven extensions' `frmtmb (>= 0.53.0)` floors must rise to the
version that ships `frm_hazard_reads()`. They are NOT bumped here: core
in this worktree is still `0.54.0`, and a floor above the core in the
same tree would refuse to install and take every check with it. The
guard fails closed rather than skipping if it meets an older core,
which is the right way round for a guard but does mean the floors
matter.

`frmtmb.learn` needs the same floor for a second reason, from item 2:
its `test-factorization.R` now asserts the reworded capped warning, so
its suite is against the newer core too.

---

## Verification

### Test files, one per R process, private library, `NOT_CRAN=true`

| package | file | pass | fail |
|---|---|---|---|
| frmtmb | test-importance.R | 226 | 0 |
| frmtmb | test-bracket-access.R | 33 | 0 |
| frmtmb | test-message-uniqueness.R | 6 | 0 |
| frmtmb.spline | test-span.R | 58 | 0 |
| frmtmb.spline | test-curve.R | 62 | 0 |
| frmtmb.spline | test-deriv.R | 38 | 0 |
| frmtmb.spline | test-gratia.R | 13 | 0 |
| frmtmb.spline | test-royston-parmar.R | 67 | 0 |
| frmtmb.spline | test-rp-floored.R | 42 | 0 |
| frmtmb.spline | test-surface.R | 48 | 0 |
| frmtmb.spline | test-message-uniqueness.R | 4 | 0 |
| frmtmb.spline | test-bracket-access.R | 1 | 0 |
| frmtmb.learn | test-factorization.R | 35 | 0 |
| frmtmb.learn | test-surface.R | 146 | 0 |
| frmtmb.learn | test-bracket-access.R | 1 | 0 |
| the other 5 extensions | test-bracket-access.R | 1 each | 0 |

Before the change, for the record: `test-importance.R` gave 203 pass
with 3 errors, and `test-span.R` 55 pass with 1 failure.

### `R CMD check --as-cran`, on built tarballs

`_R_CHECK_CRAN_INCOMING_=false`, `_R_CHECK_FORCE_SUGGESTS_=false`,
pandoc and TinyTeX on PATH, the private library first in `R_LIBS`, no
`--no-manual`.

| package | status | testthat |
|---|---|---|
| frmtmb | 1 NOTE | FAIL 0, WARN 1, SKIP 140, PASS 8238 |
| frmtmb.coupling | OK | FAIL 0, WARN 4, SKIP 1, PASS 381 |
| frmtmb.eam | 1 NOTE | FAIL 0, WARN 1, SKIP 2, PASS 1369 |
| frmtmb.latent | 1 NOTE | FAIL 0, WARN 0, SKIP 1, PASS 224 |
| frmtmb.learn | OK | FAIL 0, WARN 0, SKIP 10, PASS 262 |
| frmtmb.ode | OK | FAIL 0, WARN 0, SKIP 0, PASS 212 |
| frmtmb.sample | OK | FAIL 0, WARN 0, SKIP 3, PASS 977 |
| frmtmb.spline | 1 NOTE | FAIL 0, WARN 0, SKIP 1, PASS 329 |

Core was checked twice, the second time on the final tree after a
comment and a documentation sentence changed; both runs give
`Status: 1 NOTE` and the same `FAIL 0 ... PASS 8238`.

Every NOTE is the same line, `checking HTML version of manual ...
NOTE / Skipping checking math rendering: package 'V8' unavailable`,
which is environmental. Core's one testthat WARN is a pre-existing
`singular convergence (7)` from the optimizer at
`test-prior-compat.R:505`, raised by the branch of
`check_convergence()` this lane did not touch: it needs
`fit$importance` to be NULL to be reached at all.

`roxygen2::roxygenise()` on core and on `frmtmb.spline`: two
consecutive runs give byte-identical `NAMESPACE` and
`man/frm_hazard_reads.Rd` (md5 compared). The only files roxygen
rewrote are `NAMESPACE`, `man/frm_hazard_reads.Rd` (new),
`man/frmtmb_control.Rd` and `man/frm_curve_feature.Rd`.
