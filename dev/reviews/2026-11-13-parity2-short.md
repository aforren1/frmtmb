# Short check of lane parity2's fix round (round of 2026-11-13)

Base ec8178d3 (0.82.0, rellib-r21). Lane tree `frmtmb-wt-parity2`,
library `wt-parity2-lib` (frmtmb and frmtmb.sample built 19:01 and
19:02 UTC, after the last source edit at 19:00 UTC). This checks the
fixes of the review `2026-11-13-parity2.md` (F1 to F3, LOW items) and
the "Fix round" section of `dev/parity2-findings.md`; it does not
re-review the lane. Scripts `dev/reviews/parity2-short-*`, logs
`dev/reviews/parity2-short-out/` (gitignored). No `git diff` or
`git status`; the NEWS diff was read from `git show ec8178d3:<path>`
copies outside the tree. Free space on C: 134.4 GB at start, 135 GB at
end; the OpenBLAS copy (463 MB) was deleted.

## Verdict: FIX

One small item, F3a. F1, F2, the LOW items and NEWS hold up.

### F3a. A numeric grouping variable counts as a factor at `re_formula = NA`

`ce_numeric_first()` (`R/conditional-effects.R:458`) counts every
variable of `ce_group_vars(x)` as a factor, whatever `re_formula` is.
brms 2.23.0 takes its `group_vars` from `update_re_terms(x$formula,
re_formula)` (`conditional_effects.R:245`, then
`get_group_vars(bterms)`). At the default `re_formula = NA` there are
none, so a numeric variable that is also a grouping factor is numeric
and the pair keeps its order.

`y ~ kk * x + (1 | kk)`, `kk` integer 1 to 6, data seed 21, n 200:

| arm | re_formula = NA (default) | re_formula = NULL |
|---|---|---|
| brms 2.23.0, fixed_param fit (`brms-fit.txt`) | kk x kk:x | |
| brms 2.23.0, its own steps on `brm(empty = TRUE)` (`kk-brms.txt`) | kk x kk:x | kk x x:kk |
| lane (`kk-frmtmb.txt`) | kk x **x:kk** | kk x x:kk |
| 0.82.0 (`a-base.txt`) | kk x kk:x | |

0.82.0 matched brms on this shape at the default; the lane does not.
An edge case, but a new parity difference that the fix itself
introduced, and the roxygen of `ce_numeric_first()` ("a grouping
variable whatever its type") states brms's rule wrongly. Expected: pass
the grouping variables of the terms that `re_formula` keeps (none at
`NA`, all at `NULL`), and add `y ~ kk * x + (1 | kk)` at both values to
`test-ce-default-effects.R`.

## 1. F1: holds

- `ce_snl_fit()` defaults to `"gaussian"` (`ce_snl(dots, default =
  "gaussian")`); `ce_snl()` keeps `"uncertainty"` as its own default,
  and frmtmb.sample calls `ce_snl(dots)`
  (`conditional-effects-draws.R:178`), so draws keep brms's default.
- `a-lane.txt` (seed 100 data of the review): `"uncertainty"` and
  `"old_levels"` refused by name on boot ("a refit holds each level's
  conditional mode ..."); `sample_new_levels = "gaussian"` on boot 20,
  seed 2, is `identical()` to the default: TRUE.
- `dev/parity2-snl-width.R` rerun (`snl-width.txt`): default boot 200,
  seed 2, width 4.031 3.690 4.264; boot 400 4.036 4.344 4.560; Wald
  3.967 3.840 3.940. Equal to the findings to the digit.

### old_levels on a fit

- 0.82.0 refused the argument on every band: "conditional_effects() has
  no argument `sample_new_levels`" (`a-base.txt`, all 11 calls). So
  refusing `"old_levels"` is not a regression, and NEWS rightly does not
  mark it BREAKING.
- Boot: the refusal is justified. At an observed level
  (`conditions = data.frame(g = l)`), the boot band carries no
  conditional variance of the level (`obslevel.txt`, boot 400, seed 2):
  level 10, Wald 3.391 3.143 3.169 against boot 2.077 1.727 1.807;
  level 3, Wald 3.524 3.327 3.389 against boot 2.078 1.812 1.942. A
  lent level would read the same narrow band.
- Wald: `fitted(sample_new_levels = "old_levels")` picks level 10
  (seed 7) and gives that level's fitted value and Est.Error 0.8650
  0.8019 0.8084, equal on every digit to `fitted()` at g = 10 and to the
  Wald CE at `conditions = data.frame(g = "10")` (width 3.391 3.143
  3.169, `se__` 0.8650 0.8019 0.8084). So the Wald band could honor
  `"old_levels"` by `fitted()`'s route (pick one seen level per call,
  evaluate at it). The refusal is loud and not wrong, but it is a gap,
  not a law the band cannot carry. LOW, file: honor `"old_levels"` on
  the Wald band, or name `conditions = data.frame(<group> = <level>)` in
  the refusal.

### Pre-existing, not this lane

The boot band at an OBSERVED level is 55 to 61 percent of the Wald band
on the same fit, the same on 0.82.0 (`obslevel-base.txt`, identical
numbers). The bootstrap of the level's mode omits its conditional
variance, so `band = "boot"` with `re_formula = NULL` at a seen level
undercovers in silence. No backlog item found. MEDIUM, file.

### The width test

`test-ce-new-levels-snl.R`, third block: data `set.seed(100)`, boot
`seed = 2`; both fixed. `parity2-short-width.R`:

- 20 ulp perturbations `y * (1 + 2e-15 * k)`: min ratio 0.961 on 20 of
  20, reference BLAS and OpenBLAS 0.3.32 with LAPACK,
  `OPENBLAS_CORETYPE=Haswell` (`width-ref.txt`, `width-ob332hw.txt`),
  0 warnings.
- 30 boot seeds at the test's data: min ratio 0.820 to 1.189, mean
  0.985, sd 0.086, 0 of 30 at or below 0.8; equal on both BLAS.
- So the test is rounding-safe. LOW: its comment says 0.8 is "three"
  Monte Carlo sd below 1; for the minimum over three x values it is
  about 2.1 sd below the mean (0.985, sd 0.086). The fixed seed carries
  it.
- Under OpenBLAS 0.3.32 Haswell (`ob332hw/summary.txt`):
  test-ce-new-levels-snl pass 14, test-ce-default-effects pass 13, fail
  0, warn 0.

## 2. F2: holds

`a-lane.txt`: `"uncertainty"` and `"old_levels"` on `band = "wald"` and
`"profile"` refused, "this band draws no new level; ...".
`sample_new_levels = "gaussian"` on Wald is `identical()` to the
default: TRUE. `"gaussian"` on profile with `re_formula = NULL` stops on
the band's own rule (profile exists only at `re_formula = NA`), as
before. At `re_formula = NA`, `"gaussian"` answers and `"uncertainty"`
is refused; brms ignores the argument there; the review called refusal
there optional.

## 3. F3: holds on every requested shape; F3a above

Names of the default displays, data seed 21, n 200
(`a-lane.txt`, `a-base.txt`, `brms-fit.txt` fixed_param fits,
`brms-internal.txt`):

| shape | brms 2.23.0 | lane | 0.82.0 |
|---|---|---|---|
| `y ~ f * x` | f x x:f | f x x:f | f x f:x |
| `y ~ x * f` | x f x:f | x f x:f | x f x:f |
| `y ~ f * h` (two factors) | f h f:h | f h f:h | f h f:h |
| `y ~ x * z` (two numerics) | x z x:z | x z x:z | x z x:z |
| `y ~ f * x + (1 \| grp)` | f x x:f | f x x:f | f x f:x |
| `y ~ s(x, by = f)` | x f x:f | x f x:f | x |
| `y ~ ch * x` (character) | ch x x:ch | ch x x:ch | ch x ch:x |
| `y ~ lg * x` (logical) | lg x x:lg | lg x x:lg | lg x lg:x |
| `y ~ kk * x + (1 \| kk)` | kk x kk:x | kk x x:kk | kk x kk:x |

Effects the caller names keep their order: `"f:x"` gives `f:x`, `"x:f"`
gives `x:f` in both brms and the lane (same for `ch`, `lg`), and
`c("x:f", "x")`, `c("x", "f:x")` keep the given order in the lane, as
brms's `effects[sort(matches)]` keeps it. In the lane the `x:f` display
varies `x` on the axis (`effect1__` numeric, 5 values), `f` second, as
brms.

## 4. LOW items: hold

- `` s(`my x`) `` refused by name in `dpar_frame_rhs()`
  (`R/frame.R:169`); pinned in `test-brms-surface-mixed.R:27-29` (32
  pass).
- Newton-step allowance: `test-ce-default-effects.R:43-46` puts
  `allowed = "A Newton step from the stop"` (text, fixed match, as
  `allow_warnings()` matches) on `bf(y ~ x + (1 + z | g), sigma ~ w)`
  alone. `dev/parity2-ob-warn.R` under OpenBLAS 0.3.32 LAPACK, Haswell
  and native core (`ob-warn-hw.txt`, `ob-warn-native.txt`): that fit
  warns once ("false convergence (8). A Newton step from the stop ..."),
  the other eight 0. Not `require`d, correct since reference BLAS does
  not warn.
- The laplace refusal (`methods-draws.R:225-230`) adds, when the plan
  lends a level, "sample_new_levels = \"gaussian\" draws it from each
  draw's covariance and needs none"; asserted at
  `test-laplace-draws.R:116`.
- NEWS: "on eleven formulas" is gone.
- `dev/upstream-bugs.md` brms-52: reproducer, expected/actual (stanc,
  `ordered_logistic_glm_lpmf(Y | X, b, 0)`), "Status: not reported". No
  other worktree has a brms-52.

## 5. NEWS: holds

Both packages open with a development-version section; no lane names,
no `dev/` paths, no version numbers. Core: the default-effects bullet is
BREAKING and now carries the numeric-first order (F3); the
`sample_new_levels` bullet is not BREAKING, correct, since 0.82.0
refused the argument. frmtmb.sample: the draws default
`"uncertainty"` is BREAKING. If F3a is fixed, the F3 sentence "A pair
puts its numeric variable first, as brms does" stays true.

## 6. Tests on the lane library

One file per process, `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`,
`lib:` lines show `wt-parity2-lib` (`tests/summary.txt`):

| file | pass | fail | err | skip | warn |
|---|---|---|---|---|---|
| test-ce-new-levels-snl | 14 | 0 | 0 | 0 | 0 |
| test-ce-default-effects | 13 | 0 | 0 | 0 | 0 |
| test-brms-surface-mixed | 32 | 0 | 0 | 0 | 0 |
| test-effects | 54 | 0 | 0 | 0 | 0 |
| test-ce-levels | 101 | 0 | 0 | 0 | 0 |
| sample test-laplace-draws | 120 | 0 | 0 | 0 | 0 |
| sample test-draws-brms-surface | 21 | 0 | 0 | 0 | 0 |
| sample test-postfit-draws | 89 | 0 | 0 | 0 | 0 |
| test-latres-sparse | 67 | 0 | 0 | 0 | 0 |

test-latres-sparse no longer hits "No space left on device".

## To file at consolidation

- LOW: `"old_levels"` on the Wald band of a fit, by `fitted()`'s route,
  or a refusal that names `conditions =`.
- MEDIUM (pre-existing): `band = "boot"` at an observed level omits the
  level's conditional variance (55 to 61 percent of Wald).
