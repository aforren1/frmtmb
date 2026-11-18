# Review of lane ordstarts, build round (backlog 6136, C0 + D), 2026-11-18

Reviewer's scripts, logs and tables: `dev/reviews/ordstarts-build-review-*`
in this worktree. Lane build from `C:/Users/adf44/source/r/wt-ordstarts-lib`,
"before" from `rellib-r22` (0.83.0). Base files from `git show b255350e:`
into `dev/reviews/ordstarts-build-review-base/`, compared with
`diff --strip-trailing-cr` (`ordstarts-build-review-diff*.txt`). No package
code changed, no git operation beyond status/log/show, at most two R
processes at once.

## Verdict: FIX

The swap does what the lane says on the designs it measured: my own
generator reproduces the lane's logLik to the last digit on 40 of 40 fits,
the refit paths run no swap, the private stream leaves the session's
random state alone, and the cost medians hold. But the MEDIUM fix widened
"flat" in a way that switches off the very remedy the new warning names,
and the swap's keep rule depends on the scale of constant case weights.
Items 1 and 2 must be fixed; 3 to 7 are small and should be fixed or
filed in this round.

## Findings that need a fix

1. **MEDIUM. `frmtmb_control(starts = 10)`, which the degenerate warning now
   names on every such fit, runs no start on 40% of them.**
   `mixture_starts()` runs jittered starts only when the default is not
   flat (`run_jit <- n > 0L && t1 != 1L && ...`). The MEDIUM fix made a
   default with a threshold standard error above 100, or a lost one, flat
   (`mixture_ord_start_class()`), so on those fits the advice does
   nothing: the refit returns the same point, still warned, and the text
   changes to "Compare several starts (frm(start = ))".
   - `ordstarts-build-review-p1.R`, section 6 (`mk(seed, 300)`,
     `y ~ x + z`, seeds 1 to 40): of the default fits whose warning
     names `starts = 10`, the call ran 0 starts on 9 of 20
     (cum_cum) and 8 of 23 (cum_sratio), every one ranked flat.
   - `ordstarts-build-review-p4.R` / `-p4-sum.txt` (build's `starts = 10`
     against 0.83.0's, no random effects, 120 fits plus the compat cell):
     the build ran no start on 30 of 120 fits (0.83.0: 9). Where 0.83.0's
     starts went higher, the build stays: test-compat.R cell 2
     (`mixture(cumulative(threshold = "sum_to_zero"), sratio())`) ends at
     -539.1765 with the warning (threshold SE 21900), and 0.83.0's starts
     reached -537.4069 with no warning at the fit or at `summary()`;
     cum_cum n300 seed 3: -385.5437 warned against -385.4217 silent; seed
     16: -382.2427 against -379.8336.
   - The new test "a degenerate fit names the starts to run, unless they
     ran" (test-ordinal-mixture.R) has this case in its own "absent" arm:
     `omx_flat_data(16)` with `starts = 10` runs 0 starts (p1, cs n300
     seed 16), so it checks only that the text changed, not that starts
     ran.
   - The "a flat default stays" rule was justified by a flag that the
     optimizer's travel decides (a gap of 1e-10 against 0, cum_cum seed 10
     under autoscale). That reason does not apply to a threshold with a
     standard error of 21900.
   Fix: run the jittered starts (and let the pick replace the default
   with a better proper point) when the default is flat by the
   standard-error or lost test only; keep "stays" for collapsed and
   run-off thresholds. Independently, name `starts = 10` in the warning
   only where the starts would run. Then rewrite test-compat.R's comment
   ("it is the data's, not the cell's" is not so: 0.83.0's starts found a
   clean point 1.77 higher).

2. **MEDIUM. Constant case weights change which maximum the default
   returns.** The swap-alone rule keeps a swap only when it is higher by
   more than `grad_tol` (1e-3) on the WEIGHTED objective, so weights 1/n
   divide every gain by n. `ordstarts-build-review-p1.R` section 1
   (cum_sratio, no random effects, n 1000, seeds whose swap the lane
   kept): with `weights(w)`, w = 1/1000, the swap is not kept on 5 of 6
   seeds (1, 2, 4, 10, 13; unweighted gains 0.14 to 0.93), and the
   returned estimates differ from the unweighted fit's by 2.55 to 4.39;
   with w = 10 it is kept on 6 of 6, w = 1 on 6 of 6. No warning says so.
   The same absolute `tol` is in `mixture_starts_pick()` and in
   `mixture_swap_holds()`.
   A second weights artifact: on test-ordinal-mixture.R's sound fit
   (cum_cum seed 3, n 500), w = 1/500 makes the standard-error check
   lose a threshold, the ranking calls the point flat, and
   `starts = 10` runs nothing (record NULL; w = 1: proper, 10 starts
   run). The lane filed this as LOW, "not measured"; it is now measured.
   This is the construction asked for in question 1: the ranking refuses
   a proper point because of the weights' scale.
   Fix: compare gains per unit of mean case weight, as
   `mixture_ord_degeneracy()` already does for the sharpening, in the
   swap rule, the pick and the holds check; and do not let a lost list
   that moves with the weights' scale make a point flat (for example,
   read the lost list of the covariance per unit of mean weight).

## Smaller items (fix or file in this round)

3. **LOW. "One kind" is not exact label switching, in both directions**
   (`ordstarts-build-review-p2.R` section 2, `-p6.R`).
   - Called one kind, not the same model: components with different
     priors, fixed weights (`theta1 = 0.8`), an offset in one component,
     or a fixed `disc2 = 2`. The hand-made label switch moves the
     objective by 7.7 to 119 on all 18 such fits (0 on the 6 controls),
     and a refit from it reached a better point on 2 (seed 1: prior
     N(-0.8, 0.2) on `mu1` x, 0.084; `theta1 = 0.8`, 0.240). The build
     runs no swap on any of them.
   - Called distinct, the same model: `mixture(sratio(), cratio())` with a
     logit or probit link (P(Y = j | Y >= j) = F(tau_j - eta) in both).
     The swapped start sits at exactly the fit's objective on 10 of 10
     fits (cloglog control: 30 to 73 away), and the swap costs 1.4 to 3.7
     times the fit for nothing.
   Fix: decide by the objective, which is the definition: evaluate
   `obj$fn` at the swapped start and skip it when it equals the fit's to
   rounding. One evaluation per assignment, exact for priors, offsets,
   fixed parameters and equivalent families alike. At least, correct the
   help and NEWS ("Components of one family, link and threshold layout
   are exchangeable").

4. **LOW. A theta formula that moves the weights' reference leaves the
   weights unswapped.** `mixture_ord_swap_starts()` assumes the reference
   is component K (`theta1 .. theta(K-1)`). brms's usual spelling for
   K = 2, `theta2 ~ ...`, moves it to component 1 (`mixture_theta_reference()`).
   p2 section 2b, seed 13: `theta2_(Intercept)` is -0.3091 in the fit and
   -0.3091 in the swapped start (with `theta1 ~ 1` it is negated). At
   K = 3 with `theta2 ~ 1, theta3 ~ 1` the weight keys differ and the
   weights are not moved either (`-p7.R`). No harm measured (seed 13
   still kept the swap, gain 0.8353), but the start is not the one the
   help describes. Fix: read the reference from the family.

5. **LOW. The holds check rejected nothing and, when it fires, returns
   the worse answer.** `ordstarts-build-review-p3.R` (a wrapper that runs
   the same four runs and returns the same verdict): 23 kept swaps on
   120 fits (cum_sratio, no random effects n 1000 and n 300, `(1 | g)` n
   300), 0 rejected, all 92 jittered ends within 3.6e-8 of the swap's
   point and all proper. It costs a kept swap 7.91 against 2.84 times the
   fit (item 8). When a jittered end does climb, the code drops the swap
   and returns the default, which is lower than the swap's point, even
   when the climbing end is proper (the lane's one motivating case, a
   group sd running to zero, gained 0.0099 over a better point than the
   default; dropping it returns the default), and it counts a climb to a step function,
   the supremum the pick exists to avoid, as a rejection. Fix: keep the
   best proper converged point among the swap's end and its jittered
   ends, and ignore ends that are not proper.

6. **LOW. Test edits.**
   - test-ordinal-mixture.R "a threshold with no standard error is said
     so": on the new default the fit ends at -399.3975 (code 1); `fixef()`
     has NaN standard errors for `mu2_Intercept[1]` and `[2]`, and no
     warning names them (only the optimizer note and the degenerate
     warning on component 1; p2 section 7a, same with `starts = 10`). The
     edited test now asserts only the first `expect_false()`; the block
     at -404.9916 is unreachable. The defect is filed, so the edit is not
     hiding an unfiled one, but the F3 contract is no longer tested. Fix:
     reach the -404.9916 point with `start =` so the "has no standard
     error" path keeps a test.
   - test-compat.R: the `allow_warnings()` covers all four cells, but
     only cell 2 warns (p2 section 7b: cells 1, 3, 4 have no degenerate
     warning under the default, `starts = 10` or `starts = 0`). Allow it
     on cell 2 only. See item 1 for the comment.
   - The removed "random effects run no starts" test was the only unit
     test that fitted an ordinal mixture with `(1 | g)`. Nothing now
     tests that the swapped start moves the group sd entries
     (`thpos[to] <- thpos[from]`). Add a cheap check on a `(1 | g)` fit
     with `starts = 0` and `mixture_ord_swap_starts()`.

7. **LOW. NEWS does not say what C0 costs.** The first BREAKING bullet
   gives the speed-up but not that one-family mixtures without random
   effects are now short more often (4 of 40 silent shortfalls at n 1000
   where 0.83.0 had none; degenerate warning 14 to 20 of 40 at n 300).
   `?frmtmb_control` (`starts`) and `?mixture` give the moves and name
   `starts = 10`, so users of the help are told; add one sentence to
   NEWS. (And item 1 decides whether `starts = 10` actually helps those
   fits.)

## The eight questions

1. MEDIUM fix: the inconsistency (ranking reads a lost threshold as Inf,
   the fit-end check reads the number) is safe in one direction: the
   ranking never keeps a point the check then warns about (no fit moved
   down from A0 in the lane's 320 + 160 fits or my 20). It is not
   safe in the other: a point the check passes, or warns about only as
   SE-flat, is ranked flat, and that switches off the starts the user
   asks for (item 1) and, under weights 1/500, starts on a sound fit
   (item 2). Constructed, see items 1 and 2.
2. Label switching: not exact. Item 3 (priors, fixed weights, offsets,
   fixed dpars called identical; sratio/cratio with a symmetric link
   called distinct) and item 4 (the weights' reference).
3. RNG: `.Random.seed` and `RNGkind()` are identical before and after a
   fit that runs the swap and its holds check, under Mersenne-Twister
   and L'Ecuyer-CMRG, with the same estimates (p1 section 3). Every
   `frm()` fit creates `.Random.seed` when none existed, in 0.83.0 too
   and for a gaussian model (`-p5.R`): not this lane's, and no user
   stream is changed. The holds check: item 5.
4. Refit paths (`-p2.R` section 4, a counter on
   `mixture_swap_candidates()`): 0 calls under `refit()`,
   `frm_bootstrap()`, `influence()`, `frm_allfit()` (both
   `start_from_mle` values), `start =` of one's own, and `starts = 0`; 1
   under `update()` and `update(newdata = )`, and 1 under
   `autoscale = TRUE` (verdict not "choose"). `update()` is a new fit
   from scratch, so running the swap there is consistent with `frm()`;
   the help could say so. Under an autoscale "choose" verdict both fits
   are cold, so the swap would run twice (not triggered; as with the
   starts in 0.83.0).
5. Reproduction (`-rep.R`, `-rep-sum.R`, `-rep-sum.txt`; own generator,
   10 seeds each): lane A0 equals 0.83.0's default to 0 on 20 of 20, and
   every A0 and D logLik equals the lane's record. cum_sratio int n1000:
   silent shortfalls before 2 (seeds 1, 3; +0.815, +0.916), after 0;
   swap kept 3, D never below A0, D above 0.83.0's 10 starts on 1 (seed
   5), below on 1 (seed 9, where B lost a threshold SE). cum_cum int
   n300: before 0, after 0; D = A0 on 10 of 10. The lane's per-seed
   record agrees (`sum3.txt`: seeds 1 and 3 among its 15).
6. C0 without random effects: documented in `?frmtmb_control` and
   `?mixture`, not in NEWS (item 7). The new hint is wrong on 17 of the
   43 warned fits measured (item 1); it is not noisy beyond the
   warning itself, which fires anyway.
7. Test edits: item 6. None hides an unfiled defect, but the F3 subject
   is no longer exercised, the compat allowance is four times too wide,
   and the hint test's "absent" arm runs no start.
8. Cost, cum_sratio int n1000, seeds 1 to 10, one process per library:
   D/A0 median 3.22; where the swap was kept 7.91 (n 3), where it ran and
   was not kept 2.84 (n 7); 0.83.0's 10 starts 18.1 times the default.
   The lane's own rows for the same seeds: 8.03 and 2.89. Its pooled
   medians over all random-effect designs (7.1 and 2.36) are consistent.
   cum_cum int n300: D/A0 median 0.99 (no swap).

## Also seen (no action)

- On cum_sratio int n1000 seed 2 the build's default warns (degenerate,
  threshold SE 513) at the same logLik as A0, which is silent: the
  swap's evaluations make the deferred check build the covariance. The
  warning is right; this is the filed deferred-check item, not new.
- Column scale without autoscale changes the default fit itself (x times
  0.001 moves the fit by up to 31 in logLik), so the swap's verdict there says
  nothing about the swap (`-p1.R` section 3b).
