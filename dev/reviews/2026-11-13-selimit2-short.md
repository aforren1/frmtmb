# Short review of lane selimit2's fix round (2026-11-13)

Build: wt-selimit2-lib, installed 17:36:23, after the last edit of R/
(se-check.R 17:36:22). Base: rellib-r21. Scripts
`dev/reviews/selimit2-short-*` (and the first review's
`selimit2-review-*`, rerun); logs `dev/reviews/selimit2-short-out/`. At
most 2 R processes; the largest peak measured by gc() was 598 MB.

## Verdict: FIX (one item, F4b; small)

F1, F2, F3 and F4 are fixed as claimed. One defect remains on the 3370
item's own model class: above `joint_root_max` (4000 joint parameters)
joint_cov_white_inner() returns NULL and the old route gives a silent
standard error of exactly 0 (F4b). It is pre-existing (base gives the same
numbers), but NEWS states the fix without a size bound. Minimum: file it
HIGH and bound the NEWS bullet; better: a path that does not need the dense
factor. Everything else can ship.

## F1: fixed

`selimit2-review-grad2.R lane` (`grad2-lane.txt`): start (-25, 24.5),
iter.max 1: code 1, gradient 94.2 94.2, group b's mean 85 percent off, the
two base warnings ("would still gain 21.6", "Large maximum absolute
gradient (94.2"); iter 2 (gradient 2.57, gain 0.0304): the same two; iter
3 (gradient 0.00614, gain 1.7e-7): the code-1 note and the limit words.
The other starts behave the same way.
Attacks, by reading se_boundary_stop_rest() (`R/se-check.R:4572`): a
rising single end is left out of the projected span, so its gradient
stays in the test (a stop on the way up still warns). A residual gradient
across the end near grad_tol decides the warning as the plain gradient test
does on base; the projection adds no new rounding branch beyond a QR of at
most a few columns. No defect found.

## F2: fixed for the families it certifies

`selimit2-review-fp.R lane pois zi` (`fp-lane.txt`): delta 1e-3 no word,
SEs 5131.9 342.1; delta 1e-4 and 1e-5 "flat" (as base); one count at the
top still the joint limit; `zi ~ f`, `zi ~ 0 + f`, `hu ~ f` with level a
all zero still limits (sep1); 1 or 2 nonzero rows of 40: no sep1, finite.
`selimit2-short-cert.R` (`cert-lane.txt`, `cert-base.txt`), every
finite-MLE case identical to base:

    bernoulli, reversed pair 1e-4 apart (finite)     SEs 816.37 81.64, no word
    poisson y ~ g, level a one count of 200          SEs 1.000 1.003, no word
    poisson y ~ x, a count mid-range and at the top  SEs 2.369 0.304, no word
    zip zi ~ f, level a 3 nonzero of 40              zi_fb end only (quiet), zi_(Intercept) SE 0.603
    poisson y ~ g, level a all zero                  lane joint limit; base "flat"
    bernoulli quasi (a row at the boundary)          separation (both); lane also records joint TRUE

Reading se_limit_joint_cert1(): the per-family rules are exact recession
conditions for the log/logit links (count mean: no row up, rows down all 0;
bernoulli; hu; zi; zoi; coi among y in {0, 1}). The QR hold is not a
rounding-decided false positive: the projected direction leaves the held
rows exactly in place and is then certified by the same rules, so a
certified end is a true recession direction whatever side of 1e-3 a row
fell; the threshold only decides a false negative ("flat"). The zi far
check is sound: a zero row moved down is held when it moves by 1e-3 of the
largest or less, else put 40 units past its predictor, and the objective
level there to 1e-6 means its zi was already far below the count density's
zero, so the rest of the ray loses less than that.

Open, not this round's regression:
- The mixing-weight far check does not hold or saturate rows moved by less
  than 1e-3 of the largest, so a row near the mixing boundary can lose
  only beyond the far point. I built no joint case. The single-coefficient
  rule (0.82.0, unchanged) has the false positive:
  `selimit2-short-theta.R <arm> 1e-8` (`theta-lane-b.txt`, `theta-base.txt`),
  two gaussians 16 sd apart with a reversed pair 1e-8 apart in z:
  "theta1_z (at -821.2, toward -Inf): the data separate the rows ... the
  maximum likelihood estimate does not exist" on both builds, while theta
  times 1e4 raises the objective by 0.082 and times 1e6 by 43 (finite MLE).
  At delta 1e-6 and 1e-4 the fit stops at its finite maximum, no word.
  File MEDIUM: certify a mixing weight's end (single or joint) by holding
  the rows it moves by less than `se_sep_move` and checking every row's
  side, as for the mean families.
- Test-ordinal-mixture.R #28/#29 hold: traced run on the lane build
  (`files/omx-lane.trace`): 7 LIM lines (base 5), the two joint ends with
  the recorded estimates (15423.5, 14087.3) and (4015.4, 3744.6); 99/0.

## F3: fixed

`selimit2-review-span.R lane` (`span-lane.txt`): rows (1, 1), (1, 0.8),
(1.5, 0.9) NaN with "3 of 5 predictions move along a direction the fit does
not determine". Rows (0, 0.5) and (3, 0) keep 0.3270 and 0.1434, and they
do not move: over zi_x from -30 to -6 (objective within 5.4e-7) their means
read 1.98834 and 4.06386 at every point printed. eam test-defects.R 63/0
(below). Cost of moving every predictor (`selimit2-short-cost.R`,
`cost-*.txt`): poisson `y ~ f + x`, 60 levels of 1000 rows, 45 all zero (45
ends): fitted() on 60000 rows 3.12 s, 354 MB on the lane against 0.58 s,
346 MB on base (5.4 times, bounded by the 50-direction cap); a ZIP of 50000
rows with no limit 0.09 against 0.08 s.

## F4: fixed up to joint_root_max; F4b above it

`selimit2-review-3370.R <arm> cost 3000` (`3370-cost-*-3000.txt`), run
side by side: fitted() 0.71 s, max used 598 MB (base 0.39 s, 459 MB; the
first review's build 70.5 s, 701 MB); SEs 0.030084767 0.008550877
0.031072024, the reference's to 1e-9 (base 0.035269882 0.008372948
0.035603539).

F4b. Above `joint_root_max` joint_cov_white_inner() returns NULL
(`R/predict.R`, `if (nrow(Qg) > joint_root_max) return(NULL)`) and
joint_cov_repair() inverts `Qrr` in natural coordinates, the 3370 defect.
`selimit2-short-3370ref.R <arm> 4500 <seed>` (4510 joint parameters,
degree-7 raw polynomial, sd lost "boundary"), warnings printed, none
raised by the fit or by fitted():

    seed 2  conditional SE 0.026547773 0.006987372 0.024168042
            fitted() SE    0.000000000 0.006780297 0.000000000   (lane = base)
    seed 3  conditional SE 0.026203606 0.006970559 0.025911831
            fitted() SE    0.028095623 0.007124937 0.027467583   (7 percent high; lane = base)

A standard error of exactly 0 with no warning on the item's own model
class. Seed 1 (sd not lost, the whole route, not the lane's): fitted()
0.023272941 against a conditional lower bound of 0.024359899 (4.5 percent
low), lane = base; also pre-existing, file with it. Expected: the
lane's digits at any size, or NaN with a warning. Fix: take `a' Qrr^-1 a`
for the predicted rows from the sparse Cholesky (solve with `M' a`, no
dense factor), or at least refuse to give a number above the limit; and
bound the NEWS bullet until then.

## Other checks

- test-id-kron.R record 10, five runs per build through the lane's traced
  runner (`files/idkron-*.trace`): all ten read
  `lost=theta_3:flat se=...,0.9613283437,NaN`, identical; 61/0 each. The
  lane's final run also reads 0.9613; the NaN for theta_2 is in its older
  base and lane runs (`selimit2-suite-base`, `selimit2-suite-lane`). Not
  a regression; a flake of that record in both builds.
- 5423 block: `skip_if_not()` replaced by `expect_true(any(lost ==
  "nocurv"))` and the words chosen by whether every lost reason is
  "nocurv"; it runs on every BLAS. LOW: the fit is still wrapped in
  `suppressWarnings()`, which the round's rules replace with
  `allow_warnings()`.
- NEWS: development section, no lane name, no dev/ path, no version; the
  4824 bullet is BREAKING and says the gradient across the end is still
  tested; the F3 bullet is rewritten (every row tested; the limit taken
  first). The 3370 bullet needs the size bound (F4b).
- Merge: `git merge-file` of R/se-check.R and R/fit.R (selimit2, then
  quadseam2, then slopescale, over ec8178d3, LF, in a scratch directory):
  0 conflicts each; both merged files parse; se_limit_joint_cert() and
  se_sd_gain_up0() are in the merged se-check.R, quad_reanchored_fn() in
  fit.R. The quadseam2 and slopescale files changed since the first
  review only in R/fit.R.
- Reference BLAS, lane build (`files/lane--*.txt`, every `lib:` line the
  lane's): test-se-limit.R 378/0, test-se-tier.R 112/0,
  test-mi-families.R 194/0, test-se-reparam-reml.R 60/0,
  test-ordinal-mixture.R 99/0, frmtmb.eam test-defects.R 63/0, warn 0.
