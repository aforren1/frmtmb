# Short check of lane nlband, fix round 1, 2026-11-17

Reads the fix round of `dev/reviews/2026-11-17-nlband.md` (findings
section 10): `nl_hole_located()`, `nl_cold_par()` and the
`frmtmb_boundary_fit` filter in `nl_fit_probed()`. After arm
`wt-nlband-lib` (the installed `nl_hole_located`, `nl_cold_par` and
`nl_pole_escape` are identical to the worktree's source: `SRC ... TRUE`
in `after.txt`); before arm `rellib-r22` (0.83.0). Scripts
`dev/reviews/nlband-short-*.R` and `-*.sh`; logs in
`dev/reviews/nlband-short-out/` (gitignored). Every number below is
pasted from those logs. Code compared by `git show b255350e:<file>` and
`git merge-file -p` on copies in a scratch folder; no `git diff`. At most
two R processes ran at once.

## Verdict: FIX

Items 1 and 2 of the review are fixed as claimed, the band seeds still
reach the exact maximum, the located rule costs nothing measurable, and
the test file passes 122 of 122. Three small items remain, none a silent
wrong value that 0.83.0 did not also give:

1. **(MEDIUM) The kept hole fit loses its own boundary message.** The new
   filter drops the first fit's `frmtmb_boundary_fit` message, and
   `nl_refit_held()` runs the refit under `suppressMessages()`, so a
   kept fit at a limit says nothing about it. student() on the gaussian
   band data, seed 2 (`nlband-short-studg.R`, `studg-after.txt`):

   ```
   STUDG2 default     | par 0.99656 -2.0112 -1.6103 19.877 | fit msg frmtmb_nl_pole_start
   STUDG2 from (1,-2) | par 0.99656 -2.0112 -1.6103 20.397 | fit msg frmtmb_limit_fit
      MSG Boundary fit: nu_(Intercept) is at the end of its range (20.4 on its link scale, ...
   ```

   Same maximum (logLik 66.606612 both, `extra-after.txt`), nu at its
   limit in both, but only the direct fit explains the NaN standard
   error of nu. Hold the refit's `frmtmb_boundary_fit` messages in
   `nl_refit_held()` and say them with `pe$held` (as
   `autoscale_choose()` replays the kept fit's messages). Pin with this
   seed: the default fit must carry `frmtmb_boundary_fit` exactly when
   the fit from (1, -2) does.

2. **(MEDIUM) The located threshold refuses a better local maximum, in
   silence.** "SE below half the hole's width" sits at the noise floor;
   the separated fits it is meant to refuse have SE 1e131 (review item
   1). Sweep of 150 noisy fits (n 60, narrowest flagged hole, amp 0.15
   to 0.4, sd 1 and 2, seeds 1 to 15; `nlband-short-extra.R sweep`):

   ```
   SWEEP fits 150 | located calls 101 | refused 5 | refused that were the maximum 0 | max se/half 1.404
   REJ amp 0.15 sd 2 seed 8  | se 0.5441 half 0.408  | refused refit ll -122.565609 | returned ll -122.621985 code 0 ... warn - msg - | GLOB -117.674259 at -0.3572 in hole FALSE
   REJ amp 0.15 sd 2 seed 13 | se 0.3827 half 0.2725 | refused refit ll -127.942636 | returned ll -128.928989 ... warn - msg - | GLOB -122.683752 at -2.4170
   REJ amp 0.20 sd 2 seed 8  | se 0.4494 half 0.408  | refused refit ll -122.567193 | returned ll -122.821752 ... warn - msg - | GLOB -117.943800 at -0.3569
   REJ amp 0.20 sd 2 seed 13 | se 0.2948 half 0.2725 | refused refit ll -127.943073 | returned ll -129.405402 ... warn - msg - | GLOB -123.163770 at -2.4171
   REJ amp 0.20 sd 2 seed 14 | se 0.3427 half 0.3339 | refused refit ll -134.124191 | returned ll -134.287448 ... warn - msg - | GLOB -130.595307 at -3.4166
   ```

   No refused run was the global maximum (it lies at a pole between
   rows in all 5, the dense-pole item already filed HIGH), so the
   returned fit is 0.83.0's and no worse than before. But each refused
   refit is a code-0 maximum 0.056 to 1.462 above the fit frm() returns,
   with no message. The rule cannot reject a correct winner on these
   designs only because, at this noise, the winner is never in the hole.
   Either move the threshold away from the noise floor (any bound
   between 1.4 and 1e130 half widths separates both groups here; one tied
   to the data, such as the walls' span, is better than a constant), or
   test the supremum itself; and say a refusal (a message) rather than
   drop a higher fit silently. Also: the gate runs on the best run only,
   and a refusal returns `NULL` without trying the next run.

3. **(LOW) The student test uses an absolute tolerance.**
   `test-fit-restart.R` line 483 compares two log-likelihoods with
   `ref$control$grad_tol` (a gradient tolerance) as the bound. The
   round's rule asks for a ratio to something the run measures.

For the consolidation (no lane action): the filter conflicts textually
with lane gpselimit (below).

## 1. Separated bernoulli seeds 1 to 10, student seeds 1 to 3

`nlband-short-check.R repro`, `before.txt` and `after.txt`. The 10 SEP
lines are identical between arms (`diff` empty), for example:

```
SEP seed 2 | ll -178.76368 code 1 beta -1.20482 -0.0284362 | warn frmtmb_conv_code+frmtmb_conv_gradient msg frmtmb_nl_multistart
SEP seed 7 | ll -193.24878 code 1 beta -0.9999 -0.0336437 | warn frmtmb_conv_code+frmtmb_conv_gradient msg -
```

10 of 10 code 1 with both convergence warnings, as claimed. Student
(maximum = the fit from (1, -2) with `nl_starts = FALSE`):

```
before: STUD seed 1 | default ll -549.888859 ... msg frmtmb_limit_fit | from (1,-2) ll -42.621233 | short 507.267627
        STUD seed 2 | default ll -524.537292 ... msg frmtmb_limit_fit | short 481.455429
        STUD seed 3 | default ll -599.663089 ... msg frmtmb_limit_fit | short 521.071435
after:  STUD seed 1 | default ll -42.621233 code 0 par 0.991931 -2.00659 -1.58685 0.945136 | warn - msg frmtmb_nl_pole_start | short 0.000000
        STUD seed 2 | default ll -43.081863 ... | warn - msg frmtmb_nl_pole_start | short -0.000000
        STUD seed 3 | default ll -78.591654 ... | warn - msg frmtmb_nl_pole_start | short 0.000000
```

At the maximum on 3 of 3; the stale nu message is gone. As claimed.

## 2. The located rule as a gate

Can it refuse a correct winner? `nlband-short-check.R gate`: n 60 on
U(0, 4), the narrowest hole about the pole at 2 that `nl_pole_holes()`
flags (h 0.20 to 0.26, half widths 0.27 to 0.34), y = amp / (x - 2) +
N(0, sd), amp 1, 0.5, 0.3, 0.15, sd 1 and 0.2, seeds 1 to 4; global
maximum by the exact profile over b (`glob1()`, as the review's):

- 32 fits; 28 `nl_hole_located()` calls, 28 TRUE; largest SE / half
  width 0.639 (amp 0.3, sd 1, seed 2: se 0.2149, half 0.3361).
- 26 of 32 have the global maximum in the hole; all 26 return it
  (short between -0.0000 and 0.0000).
- The other 6 have it at a pole between rows (dense-pole item); on 2 of
  them the kept hole fit is 1.4858 and 2.3239 short with the pole
  message, on 4 the first fit is returned silently, 0.1374 to 3.2450
  short. Not caused by this round.

The wider sweep (item 2 of the verdict, sd up to 2): 5 refusals of 101
calls, none the maximum; the fallback is 0.83.0's fit and is silent
(`warn - msg -`).

Can it keep a wrong run? A supremum in a coefficient other than the
hole's: bernoulli `a / (x + b) + c * z`, y = 1 wherever z = 1
(`nlband-short-extra.R supz`, both arms). After keeps the hole fit (SE
of b 0.07254 to 0.1095, half widths 0.3003 to 0.3584) with c 22.2 to
22.8, SE of c 5413 to 6736, code 0, no warning. That is the separation
of c, and 0.83.0 is just as silent on it, both by default and from
inside the hole:

```
before: SUPZ seed 1 default                   | ll -180.473655 code 0 beta -0.35717 0.47783 22.366 se 0.2805 0.579 7710 | warn - msg -
before: SUPZ seed 1 hole start, nl_starts off | ll -146.224041 code 0 beta 0.95102 -2.0991 22.76 se 0.1334 0.09056 8203 | warn - msg -
after:  SUPZ seed 1 default                   | ll -146.224041 code 0 beta 0.95102 -2.0991 22.366 se 0.1334 0.09056 6736 | warn - msg frmtmb_nl_pole_start
```

So the gate keeps a run whose own coefficient is located while another
runs off; the silence is the SE check's (the filed LOW "SE check at
1e133" is the same gap), not this lane's. A supremum in the hole's own
coefficient was not found to pass: on separation its SE is 1e131.

## 3. Band seeds 1 to 3 and 201 to 205

`nlband-short-check.R band`, after:

```
BAND seed 1   | ll 47.731985 code 0 b -2.00375 | warn - msg frmtmb_nl_pole_start | GLOB 47.731985 at -2.00375 | short -4.10e-12
BAND seed 2   | ll 66.606612 ... | short -1.52e-12
BAND seed 3   | ll 57.327072 ... | short -2.13e-14
BAND seed 201 | ll 52.785978 ... | short -1.30e-12
BAND seed 202 | ll 56.630805 ... | short -3.25e-12
BAND seed 203 | ll 92.121313 ... | short -1.71e-13
BAND seed 204 | ll 98.938010 ... | short -1.42e-14
BAND seed 205 | ll 54.734956 ... | short -4.19e-12
```

8 of 8 at the maximum, the same values as the review's first-build run.

## 4. Cost

`nlband-short-check.R cost`: one process, the fixed build against an
emulation of the first build (`nl_hole_located()` TRUE, `nl_cold_par()`
NULL, as `dev/nlband-mockrun.R`); `ev` counts optimize_obj()
evaluations, `he` counts obj_he_cached() calls:

```
COST hole wins (band data) row 1  | fixed ev 414 he 4 ... | first ev 346 he 4 ... | ev ratio 1.197
COST hole wins (band data) row 2  | fixed ev 450 he 4 ... | first ev 385 he 4 ... | ev ratio 1.169
COST hole wins (band data) row 3  | fixed ev 438 he 4 ... | first ev 370 he 4 ... | ev ratio 1.184
COST hole loses (truth b 0.5) row 1 | fixed ev 344 he 2 | first ev 293 he 2 | ev ratio 1.174
COST hole loses (truth b 0.5) row 2 | fixed ev 384 he 2 | first ev 334 he 2 | ev ratio 1.150
COST hole loses (truth b 0.5) row 3 | fixed ev 369 he 2 | first ev 322 he 2 | ev ratio 1.146
COST hole wins, a ~ 1 + (1 | g) row 1 | fixed ev 539 he 0 | first ev 471 he 0 | ev ratio 1.144
COST hole wins, a ~ 1 + (1 | g) row 2 | fixed ev 9273 he 0 | first ev 9266 he 0 | ev ratio 1.001
```

`se_fit_analysis()` in `nl_hole_located()` is free: on all 38 calls in
`after.txt` the refit's own fit-end check had already cached the
analysis (`cached TRUE`, 0 `cached FALSE`), Hessian calls are equal
between arms, and the call took 0.000 to 0.020 s. The 1.15 to 1.20
extra evaluations come from the second (cold) hole run, which NEWS
already prices against 0.83.0. Not measured: a model where
`se_check_at_fit()` defers the Hessian (more than 10 outer parameters
and few evaluations); there `force = TRUE` builds the finite-difference
Hessian (at most 4 gradients per outer parameter), once per fit, cached
for later use.

## 5. Tests

`nlband-short-tests.sh`, `t-fit-restart.txt`:

```
lib: frmtmb 0.83.0 C:/Users/adf44/source/r/wt-nlband-lib; frmtmb 0.83.0 C:/Users/adf44/source/r/wt-nlband-lib
RESULT test-fit-restart.R pass=122 fail=0 err=0 skip=0 warn=0
```

## 6. Merge with lane gpselimit (`R/par-template.R`)

Both trees change the same line of `nl_fit_probed()`'s pole-path filter
(base line 1619): nlband adds `"frmtmb_boundary_fit"`, gpselimit adds
`"frmtmb_gp_basis"`. `git merge-file -p` on copies (CRLF stripped):
exit 1, one conflict block over that hunk. gpselimit's other two hunks
(`hold_gp` in the saddle restart, base line 866; the sign-probe filter,
base line 1629) merge cleanly.

They compose: gpselimit's message is
`class = c("frmtmb_gp_basis", "frmtmb_boundary_fit")`
(`frmtmb-wt-gpselimit/R/gp-start.R` line 1967), so nlband's entry
already drops it on the pole path; its warning form is dropped as a
`"warning"`. Resolve by keeping both entries and both comments:

```r
      !inherits(cnd, c("warning", "frmtmb_nl_local", "frmtmb_nl_restart",
                       "frmtmb_boundary_fit", "frmtmb_gp_basis"))
```

Note: on the sign-probe path only gpselimit filters (`frmtmb_gp_basis`);
a first fit's other boundary messages are still said there after the
probe replaces the fit (0.83.0 behavior, not this round's).
