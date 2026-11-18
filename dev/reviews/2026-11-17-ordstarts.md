# Review of lane ordstarts (backlog 6136), 2026-11-17

Reviewer's scripts and logs: `dev/reviews/ordstarts-review-*` in this
worktree. Base build `rellib-r22` (frmtmb 0.83.0) for every run. No
package code read from anything but `R/` at b255350e; no package code
changed. At most two R processes at once.

## Verdict: FIX

The measurement is sound in design, and the main tables regenerate
exactly from the logs. But the headline count of "silent shortfalls" is
overstated, the MEDIUM is wider than filed (it reaches the default that
ships today), and the question to the user leaves out the user's own
precedent and two options that the numbers make relevant. Items 1 to 7
below must be fixed in `dev/ordstarts-findings.md`, in its question, and
in the backlog items before the user decides.

## What regenerates (my summarizer, not the lane's)

`dev/reviews/ordstarts-review-sum.R`, log `ordstarts-review-sum.txt`:

- 320 fits per arm, 7 errors per arm on the same rows. Moved 46, up 34,
  down 12, none with 1e-6 < |d| <= 1e-3. Silent up moves 33. Ran 188,
  reached 1358 of 1880. code1A 161 (code 1), 7 NA. Per-design table
  (err, code1A, class, ran, reached, moved, up, down, silent, degA, degB,
  hintA, stepB) identical to `sum.txt` row by row. Class and code
  transition tables identical. Up d 0.01276 to 2.477, median 0.4451; down
  -10.21 to -0.7548. Random intercept 42 moved, slope 4. No proper
  default moved down; no code 0 to 1.
- Cost: per-design sums and ratios identical; all fits 1408.14 s against
  10127.95 s (7.19); per-fit ratio where starts ran median 17.96 (5.22 to
  95.80); evaluations median 11.98.
- Deferred check: 148 of 313 default fits, 9 of 313 starts fits.
- Section 3: 24 reproduced (`ordstarts-review-sum2.R`), gain 0.0396 to
  1.69; the 3 degenerate-warned B points (cum_sratio int n1000 seed 24,
  n300 seeds 24 and 33).
- Section 4 (`ordstarts-review-paths-sum.R`, `ordstarts-review-paths-sum.txt`):
  boot 40 pairs 12 above 4 below; B vs B_starts0 39 finite of 40, 16
  above 4 below; influence 47/15, 50/24, 7/18; elapsed 1676.3/94.2,
  1586.6/94.1, 3185.0/52.3, 42.0; refit starts 0/10/0, 325.9 s against
  1.8 s; allfit no starts, never reaches B. All as the lane states.
- `moves.tsv` is exactly my moved set (max |d| difference 5.3e-15).

## Independent re-runs

`dev/reviews/ordstarts-review-rerun.R` (and `-rerun-b.R`, the same with a
guard for a non-finite Hessian), logs `ordstarts-review-rerun1.txt`,
`-rerun2.txt`, `-rerun3.txt`. My own data generator from the findings'
description; y, x and g identical to the lane's `mkg()` on every seed.
Seeds: cum_sratio int n1000 13, 1 (moved up, silent); cum_cum int n1000
9 (log-gap point), 1 (not moved); cum_cum int n300 27 (moved up);
cum_sratio int n300 24 (B warns, the MEDIUM), 9 (down move off a step
function); cum_sratio int n1000 2 and 20 (not moved, hidden threshold
SE). Every A and B logLik matches the lane's (for example seed 13:
-1294.423301 and -1292.730604; n300 seed 9: -363.372625 and -372.985093).
Eight refits from each point with jitter N(0, 0.5), five times the
lane's, found nothing above B on any of these seeds, and `starts = 30`
equaled `starts = 10` on seeds 13, 1 and 27. Arms, seeds and the moved
definition are sound.

## Findings

1. **The 24 "clean interior maxima from a silent default" are 19.**
   The lane's "clean" asks for no DEGENERATE warning on B, not no
   warning; its text says "no warning". On 5 of the 24 the B fit warns
   "Standard errors are not available ... tau_raw<k>_<j>: the likelihood
   is flat along it" (cum_cum slope n1000 seed 17; cum_sratio int n300
   seeds 4, 6, 14, 31): flat threshold points, not interior maxima.
   And "silent" reads warnings at the fit only, while the standard-error
   check of 148 of 313 defaults waits for the first use: `summary()` of
   the default warns on 13 of 158 `(1 | g)` fits that were silent at the
   fit (all n = 300), 4 of them up-move seeds (cum_sratio int n300 6, 12,
   14, 24) (`ordstarts-review-swap.R`, `ordstarts-review-sum2.txt`).
   With both corrections: 19 of 320, gain 0.0396 to 1.69 (cum_sratio
   int n1000 15, cum_sratio int n300 2, cum_cum int n300 2). 8 of the 19
   B points have a group sd at zero ("Boundary (singular) fit"): valid
   boundary maxima, but not "interior"; the near-zero Hessian
   eigenvalues (1e-7 to 1e-11) on those rows are that sd, not a flat
   mixture direction. Replace "24 of 320 / about 1 in 13" with "19 of
   320, 19 of 160 with `(1 | g)`" everywhere, and define silent as "no
   warning at the fit or at the first summary()".

2. **The MEDIUM reaches the default that ships now.** The lane wrote
   "code reading, not measured" for fits without random effects. I
   measured it (`ordstarts-review-nore.R`, `ordstarts-review-nore-sum.txt`;
   `mk(seed, n, mix = TRUE)`, `y ~ x + z`, seeds 1 to 40): with
   `starts = 0` a proper fit with no degenerate warning; with the
   default 10 starts the fit moves up to a point with the degenerate
   warning on 3 of 120 (cum_sratio n1000 seeds 13 and 34, threshold SE
   2.2 to 2705 and 12.0 to 2853; cum_sratio n300 seed 28, 1.37 to 5190).
   The pick rule's promise ("no fit the fit-end check passes gains its
   warning") is broken on the shipped default, independent of the
   decision on 6136. The item must say so and should be fixed first.
   Also, of the lane's 3 random-effect cases only cum_sratio int n1000
   seed 24 is a clean instance: at n300 seed 24 the default has a lost
   threshold (tau_raw1_3; `summary()` warns), and at n300 seed 33 the
   default has a threshold SE of 4014 that the deferred check hid. The
   proposed fix (read the threshold SE in `mixture_ord_start_class()`)
   would not catch the 5 points of item 1: `mixture_ord_thres_se()`
   reads a lost threshold's covariance entry as a number (1.8 at n300
   seed 24), so the fix must also count a threshold the standard-error
   check lost as flat.

3. **"6 of 148 hide a threshold SE above 100" is 4.** 6 deferred
   defaults have a forced threshold SE above 100, but 2 of them warn
   anyway as step functions (cum_sratio int n1000 seed 23, n300 seed
   22). The hidden ones: cum_sratio int n1000 seeds 2 and 20, n300 seeds
   33 and 38; `summary()` adds no warning on 2 and 20 (checked), so they
   stay silent at threshold SEs of 513 and 2520. On 2, 20 and 38 the
   starts arm gains the degenerate warning at the SAME logLik: the
   starts did not move the point; their evaluations made
   `se_check_at_fit()` build the covariance. Fix the LOW item's count,
   and say in the findings that part of degB's difference from degA is
   this side effect, not the starts' point.

4. **Time per fit.** "about 2.2 s and 40 s a fit (n = 1000 and n = 300
   alike)" is wrong for n = 300: medians 1.41 s and 20.4 s at n = 300,
   2.1 s and 41.7 s at n = 1000 (`ordstarts-review-sum.txt`). The ratio
   claims stand.

5. **B is not the best point; the realness check is local only.** On
   cum_cum int n300 seed 27, 2 of 8 refits jittered about the DEFAULT
   fit (sd 0.5, `set.seed(9100 + k)`) reach a proper point 1.52 above A
   and 0.49 above B: code 0, threshold SEs 0.37 and 3.87, sharpenings
   -23.4 and -0.88, no degenerate warning
   (`ordstarts-review-probe27.txt`). `starts = 10` and `starts = 30`
   both stop at B. The lane's jitter (sd 0.1 about B) shows B is a local
   maximum, not that it is the answer. Say that the shortfall counted is
   a lower bound, and that 10 starts about the cold start leave such
   points (the help already says "the best of the points these starts
   reach").

6. **Missed option: one swapped start.** The gain concentrates where the
   two components are different families (15 of the 19 on cum_sratio
   n1000; cum_cum int n1000 has 0 of 40). On seed 13, A puts the slope
   -0.42 on the cumulative component and 1.51 on the sratio one; B has
   1.63 and -0.51: the two assignments of the latent classes to the two
   families, which are different models. One refit from A with the
   components' parameters swapped (mu1 and mu2, tau_raw1 and tau_raw2,
   the group sds, the weight's logit negated; `ordstarts-review-swap.R`,
   `ordstarts-review-swap-sum.txt`) reaches B on 15 of 18 up moves of
   cum_sratio int n1000 and 4 of 10 at n300, and is above B (proper) on
   2 seeds where the starts did not move (cum_sratio n1000 seed 5,
   +0.18; n300 seed 16, +0.12). It costs a median 0.46 s against 1.8 s
   for the fit and 35 s for 10 starts. On cum_cum it finds nothing (0 of
   5), as expected for exchangeable components. This is a uniform rule
   (every mixture, every design), about 1.25 times the fit, and catches
   most of what the 18-times rule catches. Not a recommendation; it
   needs its own measurement (K > 2 means K! assignments), but the user
   should see it before paying 18 times.

7. **The question omits the user's precedent.** The user rejected
   `quadrature = TRUE` by default for huber() with group terms
   (2026-11-13) although Laplace left about half the fits short, because
   "a default that depends on the model's form surprises users and costs
   time they did not choose (7 to 49 times slower)". Recommendation C
   puts 5.2 to 95.8 times (median 18.0) on every ordinal mixture with
   random effects. The lane's reading, that the `!random` exclusion is
   itself a form-dependent default, is fair; but the user's stated
   reason is cost, which argues against C as much as for it. A fair
   question lists both uniform rules: C (10 starts everywhere), and 0
   starts everywhere by default with `starts = 10` documented (which
   reverses 0.82.0's choice without random effects, where the starts
   cost 5.1 to 6.9 times on my no-RE runs and moved up 3 to 7 of 40),
   plus the cheaper uniform rule of item 6 as a candidate to measure.

## Premises of C that hold (checked in code)

- `mixture_starts_count(control, cold, random)` returns 0 under
  `starts = NULL` whenever `cold` is FALSE. `refit()` passes
  `template = object$estimates`; `influence()` passes the subset frame's
  template, always non-NULL; `frm_bootstrap()` goes through `refit()`;
  `anova(refit = TRUE)` passes `template = fit$estimates`;
  `frm_allfit()` sets `starts = 0`. So removing the `random` condition
  does NOT make a default count inherited: only an explicit
  `starts = n` in the fit's control passes into refits (the refit
  record: B's refits ran 10, A's 0). Confirmed by the paths log.
- Not in the findings, and part of C's cost: `update()` is a new cold
  fit and would run the starts; under `autoscale` with the verdict
  "choose", `autoscale_choose()` runs the plain and the seeded fit, both
  cold, so the starts run twice (`R/fit.R` `starts_cold = cold_start`;
  already true today without random effects). Not triggered on the
  lane's N(0, 1) covariates.
- C must also change the degenerate warning's hint ("A fit with random
  effects runs no jittered starts by default", `R/families.R`
  `mixture_ord_degenerate_check()`) and the `starts` help.

## The other filed items

- LOW refit paths: numbers regenerate. Correct.
- LOW deferred check: fix the count (item 3).
- LOW 7 errors: the 7 rows reproduce; the message is "NA/NaN gradient
  evaluation". Correct as filed; not reproduced in a new process.
- LOW `(1 + x | g)` code 1 on 138 of 160: reproduces. Note the slope
  design therefore says little about the starts (42 of 160 ran).

## What the user must know before choosing A, B or C

- Fix the MEDIUM first, whatever the choice: it affects today's default
  without random effects (3 of 120), and its fix must count a lost
  threshold as well as a threshold SE above 100.
- The silent, clean shortfall with random effects is 19 of 320 (19 of
  160 with `(1 | g)`), 15 of them on one design (cum_sratio, n = 1000);
  for cum_cum with `(1 | g)` it is 2 of 80. The benefit is concentrated
  in mixtures of two different families, while C's cost (median 18
  times, 20 s to 42 s a fit here) falls on every such fit.
- The user's huber precedent weighed exactly this kind of cost.
- One swapped start may capture most of the gain at about 1.25 times
  the fit; worth a measurement before an 18-times default.
- Under every option the refit paths keep running no starts unless the
  fit set `starts = n`.

## Process note

One of my background launches ran a command from the main tree's working
directory: it created `dev/reviews/ordstarts-review-rerun2.log` (53
bytes, "command not found") in `C:\Users\adf44\source\r\frmtmb`. I deleted
that one file at once. The main tree's `git status` was clean before and
after (the file was gitignored); nothing else there was touched.
