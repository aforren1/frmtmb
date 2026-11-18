# Short check of lane ordstarts, build fix round 1, 2026-11-18

Checks the lane's answer ("Build fix round 1" in `dev/ordstarts-findings.md`)
to `dev/reviews/2026-11-18-ordstarts-build.md`. Lane build from
`C:/Users/adf44/source/r/wt-ordstarts-lib` (reinstalled); "before" from
`rellib-r22`, reusing the base arms of the build review. Scripts and logs:
`dev/reviews/ordstarts-short-*` (`-run.sh` drives `-rep`, `-p1`, `-p2`, `-p4`,
`-p6`, `-p7`, `-p8`; then `-p9`, `-p10`, `-p11`). No package code changed,
no git operation beyond status/log/show, at most two R processes at once,
nothing denied.

## Verdict: SHIP

All seven items are fixed as claimed, on my own probes. Three LOW items
remain; none changes a returned fit for the worse, and two are older than
this lane. File them (below); they do not block.

## The seven items

1. **Fixed.** `-p4-sum.txt` (build's `starts = 10` against 0.83.0's, no
   random effects): the build runs no start on 6 + 2 + 1 = 9 of 120 fits,
   0.83.0 on the same 9 (it was 30); test-compat cell 2 reaches -537.4069
   with no degenerate warning (0.83.0's point; `-p2.txt` 7b: 0 and NULL
   still -539.1765 warned); the build is below 0.83.0's starts by more
   than 1e-3 on 4 fits, none silent at 0.83.0. cum_cum n300 seeds 3 and
   16 run 10 starts (`-p1.txt` section 6: seed 3 moves +0.1220 to a clean
   point). The hint: of the 14 (cum_cum) and 21 (cum_sratio) n300
   defaults whose warning names `starts = 10`, the call runs 10 starts on
   all 35 (it ran none on 17 of 43 before). A default flat by travel
   (run off, collapsed) no longer names the call: the count of hint fits
   fell from 20 to 14 and from 23 to 21, the travel-flat ones. The hint
   rule (`step_any || (!travel && converged)`) matches `run_jit`
   (`!travel && (converged || t1 == 2)`) branch for branch.

2. **Fixed for the keep rule; one optimizer effect remains (LOW, not this
   lane's).** `-p1.txt` section 1: on the 6 seeds the swap is kept under
   w = 1/1000, 1 and 10 alike (it was dropped on 5 of 6 at 1/1000).
   `-p9.txt` (seeds 1 to 20, w = 1/1000 and 10 against 1, default D and
   `starts = 0` A0): 39 of 40 D pairs end at the same likelihood per unit
   weight within 1.3e-7. Parameters differ by up to 30 on 11 pairs, but
   A0 differs the same way on 9: flat directions, the optimizer's, at
   equal likelihood. The exception is seed 2 at w = 1/1000: D ends 0.226
   per unit below w = 1 (gain 0.0715 against 0.2978). `-p10.txt` shows
   why: plain `nlminb()` from the same swapped start stops at -1292.6624
   per unit at w = 1/1000 and at -1292.4361 at w = 1, both "relative
   convergence (4)". The keep rule is scale-free now; nlminb's stop is
   not. The lane's "w = 1/1000, 1, 10 give the same maximum" should say
   "on seed 13" (its test's seed), not in general.
   The lost-threshold rule: dividing H by the mean weight before
   `cov_from_hessian()` is the same as multiplying the check's absolute
   curvature cuts by the mean weight; relative cuts are unchanged. Could
   the ranking (unit H) and the fit-end check (raw H) disagree visibly? Only
   where the two lost lists differ at one point; on 18 fits (`-p11.txt`,
   cum_cum seeds 3, 5, 8, w = 1, 1/500, 1/5000, starts 0 and 10) they
   agreed every time. The weights effect there is again the optimizer's
   stop: seed 3 at w = 1/500 with `starts = 0` stops 0.0019 per unit
   below w = 1, with a lost `tau_raw1_3` (SE warning); with `starts = 10`
   the starts reach the w = 1 point and the warning goes. At w = 1/5000
   the 10 starts run and the fit stays at the lower point. Seed 8 at
   w = 1/500 stops at collapsed thresholds (degenerate warning) at the w =
   1 likelihood, where w = 1 has none. Both happen with `starts = 0`, which
   is 0.83.0's fit to the last digit, so they predate this lane. So the
   unit-H device is harmless but did not change any fit I ran; the
   residual weights dependence is nlminb's.

3. **Fixed.** `mixture_swap_distinct()` (`-p2.txt` section 2, `-p6.txt`):
   cum_cum controls run 0 swaps (6 of 6); a per-component prior (12),
   `theta1 = 0.8` (2), an offset in one component (2) and `disc2 = 2` (2)
   now run 1 each (seed 1 with `theta1 = 0.8` keeps it, +0.2396).
   sratio with cratio under logit and probit: 10 of 10 swapped starts at
   exactly the fit's objective, 0 run; cloglog: 5 of 5 run. The inner
   start (`-p8.txt`, 15 fits with `(1 | g)`: cum_sratio n1000 and cum_cum
   n300 and n1000, seeds 1 to 5, built against `mixture_swap_inner()` made
   a no-op in the same process): returned logLik differ by at most
   1.4e-12 and parameters by at most 6.7e-8; it changes speed only.

4. **Fixed.** `-p7.txt`: for every swapped start, component k gets
   component sigma(k)'s mixing weight to 1.1e-16, at K = 2 with the
   reference on component 2 or 1, and at K = 3 with it on 3, 1 or 2 (5
   starts each). `-p2.txt` 2b: `theta2_(Intercept)` -0.3091 becomes
   +0.3091.

5. **Fixed.** `mixture_swap_settle()` replaces the swap's point only with
   a converged end that is higher by more than the gain unit and proper,
   so it never returns below the swap's point, which is above the default
   by more than `tol`. No default fell below the `starts = 0` fit: 0 of 15
   (`-p8.txt`), 0 of 20 (`-rep-sum.txt`), and the lane's 320 + 160. A step
   function end is not counted (`proper()` gates it; the unit test passes
   `function(r) FALSE` and keeps the point).

6. **Fixed.** The F3 test reaches -404.9916 through `start =` (a start
   turns the swap off) and asserts "has no standard error", no quoted SE,
   and the `bcs` names; test-compat.R allows the degenerate warning on
   cell 2 only, with a corrected comment; "a swap moves the group sd
   entries with their component" checks the two sds and the modes block
   for block. The hint test's absent arm now asserts that 10 starts ran.
   Nit: "constant weights do not change the maximum returned" asserts
   `max(abs(par diff)) < 1e-3`, an absolute tolerance (lane-env forbids
   them) on a quantity that differs by up to 30 on flat directions at
   other seeds (`-p9.txt`); on seed 13 it holds on 21 of 21 roundings
   (lane's `ulps-weights.txt`). Prefer the per-unit logLik ratio alone.

7. **Fixed.** NEWS's first bullet gives C0's cost (4 of 40 silent at n
   1000, the warning 14 to 20 of 40 at n 300).

## Before and after (own reproduction, 10 seeds each)

`-rep-sum.txt` (lane arms rerun on the fixed build, base arms from the
build review):

- cum_sratio `(1 | g)` n1000: lane A0 = 0.83.0 default (max |dll| 0);
  silent shortfalls 2 before (seeds 1, 3), 0 after; swap kept 3; D below
  A0 0; D/A0 median 3.13, 6.91 where kept, 2.54 where it ran and was not
  kept (lane, pooled over designs: 6.89 and 2.23).
- cum_cum `(1 | g)` n300: A0 = 0.83.0 (0); 0 before, 0 after; D = A0 on 10
  of 10; D/A0 median 1.00.

Consistent with the lane's 19 to 3 of 320 and 0 moved down.

## To file (LOW)

- **nlminb's stop depends on the scale of constant weights** (older than
  this lane; `starts = 0` shows it). cum_sratio n1000 seed 2: the refit
  from the swapped start stops 0.226 per unit below at w = 1/1000
  (`-p10.txt`); cum_cum n500 seed 3 at w = 1/500 stops 0.0019 below with
  an SE warning, seed 8 at a collapsed threshold with the degenerate
  warning (`-p11.txt`). Same likelihood per unit at w = 1.
- **The objective test passes some label switches at unconverged
  points.** On one-family `(1 | g)` fits the swapped start of a code-1
  default can differ from the fit's objective by more than
  `sqrt(eps) * |f|`, so a futile swap runs: 14 of 160 one-family fits in
  the lane's own after-run (cum_cum int n300 5, slope n1000 5, slope n300
  4), cost up to 1.72 times, swap never kept. Whether it runs flips with
  the inner start (`-p8.txt`: cum_cum n300 seed 4 runs with it, seed 5
  without). Cost only.
- The absolute tolerance in the weights test (item 6).
