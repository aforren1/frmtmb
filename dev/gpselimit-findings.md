# Lane gpselimit: an approximate field at its basis's white limit

Round of 2026-11-16, base b255350e (frmtmb 0.83.0). Item 5571 of
`dev/test-backlog.md` (MEDIUM, gpfield review). Brief:
`round-briefs/r084-gpselimit.md`.

Scripts (all in `dev/`): `gpselimit-lib.R` (helpers; arms `lane` =
wt-gpselimit-lib over rellib-r22, `base` = rellib-r22, `r21` =
rellib-r21, `gpfield` = wt-gpfield-lib over rellib-r21),
`gpselimit-sweep.R` (the review's 30 fits, same data, optional ulp
jitter), `gpselimit-cmp.R` (bitwise base against lane),
`gpselimit-show.R` (what a user reads), `gpselimit-notflat.R` and
`gpselimit-notflat2.R` (fits where the basis fires and the ridge is not
level), `gpselimit-api.R`, `gpselimit-hyp.R`, `gpselimit-hyp2.R` (the
determined product through `hypothesis()`), `gpselimit-exact.R`,
`gpselimit-exact2.R` (exact `gp()`), `gpselimit-install.sh`,
`gpselimit-par.sh` (test driver). Logs are in `dev/gpselimit-out/`
(local, gitignored); the numbers below are pasted from them.

## 0. Fix round 1 (review dev/reviews/2026-11-17-gpselimit.md, FIX)

The review found three defects in the first hand-back and one error in
this record; section 1 below is corrected, and section 9 holds the fix
round's changes and numbers. In short:

- The certificate did not look outward: a stop inside the white plateau
  below a resolved maximum got the boundary message (5 fits, 3.6 to 14.3
  short). Fixed by the way-back rule.
- Section 1 said the item was misdiagnosed. That holds for single fields
  only. For `gp(by = f)` white sub-GPs, 0.83.0 kept finite SEs, as the
  item said (review finding 3; 42 of 42 draws, row 21 sdgp 197.7 with
  an SE of 4916 on the log scale), and this lane's `gp_white_dirs()`,
  which section 2 called a safety net, is what fixed them.
- A fit with a convergence warning got the message. Now it keeps the
  warning.
- Held: the message in nl_fit_off_saddle() and a replaced fit's message
  in nl_fit_probed().

## 1. The item is right for gp(by = ) only: single fields had NaN SEs

Correction (fix round 1): what follows holds for a single `gp(x, k)`
field. For the white sub-GPs of `gp(x, by = f, k)` the item was right:
0.83.0 kept finite SEs on 4 of 4 white hyperparameters on sweep row 21
and on 2 of 4 on row 3, on 21 of 21 roundings each, with no word
(dev/reviews/gpselimit-review-out/byjit-base.txt). The first hand-back
did not test `gp(by = )` and called the whole item misdiagnosed; that
was wrong.

The item says that on the 13 fits where `frmtmb_gp_basis` fires, sdgp
and lscale "keep finite standard errors (`se_lost 0` on all 13)". For
these single fields they do not, on any build:

- rellib-r22 (0.83.0), all 30 rows (`sweep-base.txt`): the 13 white rows
  have `se(theta) NaN,NaN` and `se_lost` `theta_1=boundary,
  theta_2=boundary`; the 17 resolved rows have finite SEs and lose
  nothing.
- rellib-r21 (0.82.0) and wt-gpfield-lib over it (the build the review
  measured), rows 1, 3, 9, 21, 25 (white) and 6, 7 (resolved)
  (`sweep-r21.txt`, `sweep-gpfield.txt`): the same, NaN and "boundary"
  on every white row.
- Under 20 ulp perturbations of the response (`y * (1 + 2e-15 * j)`,
  j = 1..20) on rellib-r22, all 30 rows (`sweep-base-jit.txt`): 260 of
  260 white fits NaN with "boundary", 340 of 340 resolved fits finite.
  The white fits' length scales are 2.75e-5 to 1.59e-4 of the spacing.

The review's column `se_lost` was
`any(grepl("Standard errors are not available", q$warn))`
(dev/reviews/gpfield-review-basis.R): it read the SE warning's text,
and the SE check is silent about this verdict by design (a smooth, gp()
or hsgp() term at its own limit loses its SEs without a message,
`?frmtmb_control`, `check_se`). The verdict came from
`se_boundary_names()`: tier 3 removed the ridge direction (eigenvalue
3.65e-9 of the largest on row 3, `show-base-3.txt`), and
`se_block_ridge_edge()` (lane setier) called the block's hyperparameters
"limit" kind of "boundary". The summary on row 3 already read
`sdgp(gpx) 154.65 NaN NaN NaN`, with the generic boundary footnote ("a
variance at zero, or a correlation, mixing or range parameter at its
limit").

What was left of the item, and what this lane did:

1. The verdict rested on `se_block_ridge_edge()`: one step of 2 along
   the block's flattest Hessian direction, level within `grad_tol`
   (1e-3) on one side. That is not a proof that the ridge is level out
   to the limit, and it passed a stop short of it (section 3).
2. The summary's footnote did not describe the ridge (sdgp at 155 is
   not "a variance at zero").
3. The warning stayed a warning on fits the SE check had already
   handled as a boundary.

## 2. The change

`R/gp-start.R`:

- `gp_white_ridge(fit, bk)`: the certificate. For an hsgp block under
  `gp_basis_white()`'s gate, the objective (`se_probe_fn()`) along each
  ridge of `gp_white_ridge_dirs()` (isotropic: log sdgp up D/2, log
  length scale down 1; anisotropic: one per length scale, log sdgp up
  1/2, that length scale down 1) at `se_edge_step` (2) and
  `se_edge_step + se_sd_far` (12) must stay within `se_edge_tol` (1e-6)
  of the fit's. Why two points suffice: every kernel's spectral density
  departs from its white limit by a factor whose log is proportional to
  the squared length scale near 0, so the gap falls by exp(-2) per unit
  of the move; at the gate it is at most 5 percent of a basis
  function's variance, at t = 12 4e-11 of that (the limit to rounding).
  Cached per block in `fit$cache$gp_white` with the point; the
  objective's state is put back. Cost: four objective evaluations per
  white block (more for anisotropic), none on a fit under no gate.
- `gp_white_dirs(fit, np)`: the certified ridges as joint directions for
  tier 3. On single fields tier 3 removed them by itself (260 of 260
  white fits); on `gp(by = f)` it did not, and this is the fix there
  (corrected in fix round 1; it is not a mere safety net).
- `gp_basis_message()` returns up to two conditions: a message of class
  `c("frmtmb_gp_basis", "frmtmb_boundary_fit")` for certified terms
  ("Boundary fit: ... is at the white limit of its basis ...", with the
  `hypothesis()` string of the product, built by `gp_white_product()`
  from `gp_brms_terms()`), and the old warning for the rest.
  `gp_end_warn()` routes the message through `se_boundary_act()`, so
  `check_se = "ignore"` silences it and `"stop"` stops, as for the SE
  check's own boundary messages.

`R/se-check.R` (three small hunks):

- `cov_from_hessian()`: `gp_white_dirs()` appended to `joint`.
- `se_boundary_names()`: for a block under the gate, kind "gpwhite" only
  when `gp_white_ridge()` holds; otherwise no kind (the parameters stay
  "flat"). Other gp/hsgp blocks keep `se_block_ridge_edge()`.
- `se_lost_clauses()`: a "gpwhite" clause for the summary footnote and
  the SE warning ("an approximate gp() term at the white limit of its
  basis: the likelihood stays level from the estimate to a length scale
  of 0 ... reads only sdgp^2 times the product of the length scales,
  which hypothesis() reports: a boundary fit").

The reason stays "boundary" (decision 1 below). Docs: `?frmtmb-gp`
(R/gp-by.R) and `?frmtmb_control` (`check_se`, R/fit.R). NEWS: two
bullets under `# frmtmb (development version)`, one BREAKING (warning to
message).

## 3. Measurements

### The 30 fits, before (rellib-r22) and after (lane)

Generated from `sweep-base.txt` and `sweep-lane.txt` (separate
processes):

<!-- generated from dev/gpselimit-out/sweep-table.md -->
| row | n | ls0 | k | ls/spacing, sdgp | base: se(theta); lost; warning; message | lane: se(theta); lost; warning; message |
|---|---|---|---|---|---|---|
| 01 | 60 | 0.7 | 5 | 0.000159 sd 33.67 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 02 | 60 | 0.7 | 8 | 4.19e-05 sd 104.5 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 03 | 60 | 0.7 | 10 | 3e-05 sd 154.7 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 04 | 60 | 0.7 | 15 | 2.75e-05 sd 94.52 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 05 | 60 | 0.7 | 20 | 5.47e-05 sd 59.98 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 06 | 60 | 0.7 | 30 | 1.06 sd 0.9725 | 0.2467,0.1582; none; none; none | 0.2467,0.1582; none; none; none |
| 07 | 100 | 1.5 | 5 | 0.606 sd 0.9665 | 0.4782,0.7679; none; none; none | 0.4782,0.7679; none; none; none |
| 08 | 100 | 1.5 | 8 | 0.696 sd 9.856 | 0.4597,0.4391; none; none; none | 0.4597,0.4391; none; none; none |
| 09 | 100 | 1.5 | 10 | 3.59e-05 sd 96.25 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 10 | 100 | 1.5 | 15 | 4.56e-05 sd 41.79 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 11 | 100 | 1.5 | 20 | 0.869 sd 0.5102 | 0.3039,0.3197; none; none; none | 0.3039,0.3197; none; none; none |
| 12 | 100 | 1.5 | 30 | 1.63 sd 0.7968 | 0.3337,0.1833; none; none; none | 0.3337,0.1833; none; none; none |
| 13 | 100 | 3.0 | 5 | 0.447 sd 1.657 | 0.3791,1.0211; none; none; none | 0.3791,1.0211; none; none; none |
| 14 | 100 | 3.0 | 8 | 3.76e-05 sd 16.37 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 15 | 100 | 3.0 | 10 | 0.893 sd 0.4943 | 0.4568,0.3954; none; none; none | 0.4568,0.3954; none; none; none |
| 16 | 100 | 3.0 | 15 | 1.13 sd 0.3033 | 0.3636,0.2879; none; none; none | 0.3636,0.2879; none; none; none |
| 17 | 100 | 3.0 | 20 | 2.14 sd 0.7077 | 0.4138,0.2421; none; none; none | 0.4138,0.2421; none; none; none |
| 18 | 100 | 3.0 | 30 | 3.38 sd 0.8417 | 0.393,0.218; none; none; none | 0.393,0.218; none; none; none |
| 19 | 200 | 1.0 | 5 | 0.823 sd 3.59 | 0.5699,0.4791; none; none; none | 0.5699,0.4791; none; none; none |
| 20 | 200 | 1.0 | 8 | 0.747 sd 86.57 | 0.3330,0.3407; none; none; none | 0.3330,0.3407; none; none; none |
| 21 | 200 | 1.0 | 10 | 0.000106 sd 182.5 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 22 | 200 | 1.0 | 15 | 0.000113 sd 75.77 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 23 | 200 | 1.0 | 20 | 0.822 sd 0.7789 | 0.2815,0.2755; none; none; none | 0.2815,0.2755; none; none; none |
| 24 | 200 | 1.0 | 30 | 1.64 sd 1.073 | 0.2848,0.1422; none; none; none | 0.2848,0.1422; none; none; none |
| 25 | 50 | 2.0 | 5 | 4.65e-05 sd 52.78 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 26 | 50 | 2.0 | 8 | 4.75e-05 sd 23.68 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 27 | 50 | 2.0 | 10 | 3.25e-05 sd 26.31 | NaN,NaN; boundary x2; frmtmb_gp_basis; none | NaN,NaN; boundary x2; none; frmtmb_gp_basis |
| 28 | 50 | 2.0 | 15 | 0.981 sd 1.05 | 0.3604,0.2680; none; none; none | 0.3604,0.2680; none; none; none |
| 29 | 50 | 2.0 | 20 | 2.16 sd 1.229 | 0.4888,0.2510; none; none; none | 0.4888,0.2510; none; none; none |
| 30 | 50 | 2.0 | 30 | 1.56 sd 0.3709 | 0.3961,0.3255; none; none; none | 0.3961,0.3255; none; none; none |

Bitwise (`gpselimit-cmp.R`, `cmp-base-lane.txt`): on all 30 rows the
estimates (`fit$estimates`), the outer point, the log-likelihood, every
SE (NaN as NaN) and the lost reasons are `identical()` between base and
lane; only the conditions differ, on the 13 white rows.

Under 20 ulp perturbations (`sweep-lane-jit.txt` for j = 1..14 and
`sweep-lane-jit2.txt` for j = 15..20; `cmp-jit.txt`): 260 of 260 white
fits get the message (no warning), 340 of 340 resolved fits no
condition, and 600 of 600 rows are `identical()` to base in estimates,
point, log-likelihood, SEs and lost reasons. (A first run of j = 15..20
read rellib-r22, because this lane reinstalled its library while it ran;
its `lib:` lines show it, and it was rerun. That run's lines are not
counted.)

### A fit where the basis fires and the ridge is not level

`gpselimit-notflat2.R` on sweep row 3 (`notflat2-base.txt`,
`notflat2-lane.txt`). Ridge columns: the objective at t = -2, 2, 12
along the ridge minus the fit's.

| case | ls/spacing | ridge t=-2, 2, 12 | base | lane |
|---|---|---|---|---|
| no prior | 3.0e-5 | 2.25e-7, -4.12e-9, -4.2e-9 | NaN, boundary, warning | NaN, boundary, message |
| prior lognormal(log(0.01 sp), 0.05) | 0.01 | 800, 800, 2.88e4 | finite (0.341, 0.050), warning | same, warning |
| prior at 0.03 sp | 0.03 | 800, 800, 2.88e4 | finite, warning | same, warning |
| prior at 0.05 sp | 0.05 | 801, 800, 2.88e4 | finite, warning | same, warning |
| prior at 0.09 sp | 0.0900 | 807, 800, 2.88e4 | finite, warning | same, warning |
| iter.max 15 (stop short, code 1) | 0.0097 | 0.024, -4.31e-4, -4.39e-4 | NaN, boundary, basis and code warnings | NaN, flat, basis and code warnings |

So the prior fits keep the warning (their SEs read the prior: lscale's
0.050 is the prior's sd), and the stop 4.3e-4 short of the limit is no
longer called a boundary. The first attempt at a data-made interior
maximum under the gate (`gpselimit-notflat.R`: basis coefficients set to
their prior sd at 0.02 to 0.08 of the spacing, with and without an
intercept, K 10 and 30) ended at the white limit every time
(`notflat-base.txt`): at that gate a basis function's variance varies by
at most 5 percent across the basis, and k coefficients cannot hold that
against the intercept and the noise. I did not find such a fit without a
prior.

### Higher dimensions and Matern

`gpselimit-hyp.R` (`hyp-lane.txt`): two-dimensional fields drawn from
the white limit itself (coefficients iid N(0, 1) on the k = 4 basis,
noise 0.1, seeds 5161 to 5163), exp_quad and matern32, isotropic and
anisotropic. Seed 5161 ends white in all four (length scales 1.8e-5 to
3.8e-5 of the spacing): each gets the message, "gpwhite" for every
hyperparameter, and its `hypothesis()` string evaluates to -2.6368 with
SE 0.39439. Seeds 5162 and 5163 end resolved (0.25 to 1.99 of the
spacing) with no condition. One anisotropic fit (exp_quad, seed 5163)
has one length scale at 0.036 of the spacing and the other at 0.36:
`gp_basis_white()` asks all dimensions, so no condition, finite SEs, as
on base.

That exp_quad and matern32 give the same product is the limit, not a
bug (`gpselimit-hyp2.R`, `hyp2-lane.txt`): both fits put every basis
coefficient's prior sd at 0.670665, and the two kernels' spectral
constants at D = 2 agree (log 2 pi = 1.8379 against matern32's
1.8378, arithmetic in the section's notes), so `2 log sdgp + 2 log
lscale` is -2.636847 in both.

### The determined product through the API

`gpselimit-api.R` (`api-base.txt`), base: `hypothesis(fit, "2 *
log(sdgp_gpx) + log(lscale_gpx) = 0", class = "")` on row 3 gives
-1.72 with SE 0.68 (z -2.53), while `hypothesis(fit, "sdgp_gpx = 0")`
gives NaN with the `frmtmb_se_lost_prediction` warning. So the API
could already report the product; what was missing was saying so. The
message now names that string; `hyp-lane.txt` evaluates it from the
message text: row 3 -1.7194 (SE 0.68044). `fitted()` and
`frm_linpred(se.fit = TRUE)` keep their SEs on the white fit (0.373,
0.207, 0.200 on rows 1 to 3 of the data), as they did.

## 4. Exact gp() (brief item 5)

The same ridge cannot occur. An exact field's covariance at length
scale l is `sdgp^2 K(l)` with `K(l) -> I` as l falls under the smallest
gap: the limit is an iid effect per position with variance `sdgp^2`, so
sdgp stays determined (with replicates per position; with one row per
position it trades against sigma instead) and only the length scale is
at the end of its range. `gpselimit-exact.R` (`exact-base.txt`): an iid
effect per position plus noise, 40 positions, 1 or 3 rows each, seeds
1 to 4:

- rep 3 seed 4: length scale 0.040 of the smallest gap, the objective
  level toward 0 (0 at -2 and -12): `theta_2` "boundary" (kind
  "limit"), sdgp SE 0.114, silent. Right.
- rep 1 seed 1: sdgp ran to 0.0066 (sigma took the variance): "zero" and
  "undefined". Right.
- rep 3 seed 2: length scale 0.18 of the smallest gap, the objective
  level toward 0 (-6.1e-9 at -2 and -12): in this process, after the
  script's own probes of the objective, `sdr_of()` lost `theta_2` as
  "concave" with no condition; in a fresh process (`gpselimit-exact2.R`,
  `exact2-base.txt`) the fit-time analysis stood at tier 1 and
  `summary()` gave lscale an SE (0.12 shown, CI 0 to Inf), with no
  condition. Filed (below); not this lane's to fix.

## 5. Decisions and questions for the user

1. **Reason "boundary", not "limit".** The brief asks for reason "limit"
   "as 0.82.0 does for a coefficient at the end of its range". In the SE
   check's vocabulary "limit" is the verdict of `se_edge_limit()` on
   linear-predictor coefficients, and its consumers
   (`se_limit_quiet()`, `se_limit_draw_held()`, `se_null_drop_held()`,
   `se_limit_ranges()`, `jc_limit_extrap()`) are built for coefficients
   run off; a covariance parameter at the edge of its space ("a ...
   range parameter at its limit") is "boundary", which is what the 13
   fits already had. Moving the gp hyperparameters to "limit" would
   route them through those consumers for no change a user can read,
   and change the predict() draw law on these fits (drawdir's area). I
   kept "boundary" with a kind of its own ("gpwhite") for the text.
   Recommendation: keep it. If the user wants the label "limit", it is
   one line in `se_relabel()` plus the consumers' review.
2. **Message, not warning, where proved.** The fit is the white limit's
   model (k coefficients with one common prior variance), at the
   likelihood's maximum: the review profiled the 13 (within 0.012 at
   0.05 of the spacing, falling monotonically above). That is a boundary
   fit, which this package tells with a message. The two cases where the
   proof fails (prior, stop short) keep the warning. This is BREAKING,
   and marked so.

## 6. Tests

New file `tests/testthat/test-gp-white-limit.R` (three tests, core only,
about 3 s):

- a white fit (sweep row 3): one `frmtmb_gp_basis` condition, a message
  of class `frmtmb_boundary_fit`; sdgp and lscale SE NaN, the rest
  finite; `se_lost` "boundary" twice; the summary's footnote says "at
  the white limit of its basis"; the message's `hypothesis()` string
  evaluates, its estimate equal to `2 log sdgp + log lscale` from
  `gp_brms_values()`.
- guard absent: a prior holding lscale at 0.05 of the spacing keeps the
  warning and finite SEs; the white fit moved back along its own ridge
  to 0.01 of the spacing fails `gp_white_ridge()` while the fit itself
  passes.
- a resolved field (row 6, k 30) gives no condition and finite SEs;
  `check_se = "ignore"` silences the message and `"stop"` makes it an
  error.

Seen to fail on rellib-r22 (`new-base.log`): `pass=13 fail=7 err=2
skip=0 warn=4` (the condition is a warning, the footnote generic, no
hypothesis string, `gp_white_ridge` not found, `check_se` does not
govern it). On the lane build: `pass=24 fail=0 err=0 skip=0 warn=0`
(`new-lane.log`).

Changed test sites. The certified fits now give a message, which printed
into the test output where a test allowed the basis warning; each such
site keeps `allow_warnings()` (a rounding may fail the proof and warn)
and gains `allow_boundary()`: test-brms-lik-gp-exact.R (2 sites),
test-brms-likelihood.R (2), test-gp-by.R (4), test-gp-multidim.R (1, two
calls), test-marginaleffects-re.R (1), test-disc-no-fixed.R (1, a
`suppressWarnings()` fit). test-gp-field.R's "a field under its basis's
spacing warns" required the warning text, and failed on the lane build
(`t1-lane.log`, `pass=21 fail=1`); it now counts one `frmtmb_gp_basis`
condition, message or warning.

## 7. Suites, OpenBLAS, check

Targeted core run, gated (`FRMTMB_BRMS_FIT_TESTS=true`), every file that
writes `gp(` plus every test-se-*.R, 45 files, one per process
(`t1-lane.log`): 44 clean, test-gp-field.R `fail=1` (above). After the
test edits, the 8 changed files again (`t2-lane.log`):

    test-gp-white-limit.R pass=24 fail=0 err=0 skip=0 warn=0
    test-gp-field.R pass=22 fail=0 err=0 skip=0 warn=0
    test-brms-lik-gp-exact.R pass=74 fail=0 err=0 skip=0 warn=0
    test-brms-likelihood.R pass=556 fail=0 err=0 skip=0 warn=0
    test-disc-no-fixed.R pass=12 fail=0 err=0 skip=0 warn=0
    test-gp-by.R pass=86 fail=0 err=0 skip=0 warn=0
    test-gp-multidim.R pass=67 fail=0 err=0 skip=0 warn=0
    test-marginaleffects-re.R pass=59 fail=0 err=0 skip=0 warn=0

and no message text in any of their logs.

OpenBLAS 0.3.26, LAPACK mode (`dev/ciharden-openblas.sh 0.3.26
lapack`, `OPENBLAS_NUM_THREADS=4`; `gpselimit-openblas-run.sh`,
`ob0326-summary.txt`), ungated:

    test-gp-by.R pass=71 fail=0 err=0 skip=2 warn=0
    test-gp-field.R pass=21 fail=0 err=0 skip=0 warn=0
    test-gp-kernels.R pass=39 fail=0 err=0 skip=0 warn=0
    test-gp-multidim.R pass=67 fail=0 err=0 skip=0 warn=0
    test-gp-start.R pass=251 fail=0 err=0 skip=0 warn=0
    test-gp-white-limit.R pass=24 fail=0 err=0 skip=0 warn=0
    test-marginaleffects-re.R pass=59 fail=0 err=0 skip=0 warn=0
    test-se-check.R pass=134 fail=0 err=0 skip=0 warn=0
    test-se-limit.R pass=379 fail=0 err=0 skip=0 warn=0
    test-se-tier.R pass=114 fail=0 err=0 skip=0 warn=0

(The two skips are the brms gate.) The 30-fit sweep under the same copy,
jitter 0 to 5 (`sweep-lane-ob0326.txt`): 78 of 78 white fits get the
message, 102 of 102 resolved fits nothing. Not run under 0.3.32: the
certificate's gate (1e-6) sits two to three orders above every value
measured on a white fit (`f(t) - f0` between -8.7e-8 and -4.1e-9 at
t = 2 and 12 in `notflat-base.txt` and `notflat2-*.txt`), so no fit
seen is near the rounding-decided side. `dev/ciharden-out/` deleted.

Whole suites, ungated (`NOT_CRAN=true`), one file per process, two at a
time, every log's `lib:` line reading wt-gpselimit-lib for frmtmb
(`all-lane.log`, 550 of 550 files ran; totals summed from the RESULT
lines). R CMD check --as-cran of core, built with vignettes, inside
`dev/gpselimit-check/` (`dev/gpselimit-check.sh`): `Status: 1 NOTE`, the
expected V8 math-rendering note on the HTML manual; tests `[ FAIL 0 |
WARN 0 | SKIP 279 | PASS 25407 ]` in 55 min, vignettes re-built OK.

One session per suite (`dev/gpselimit-sessions.sh` over
`dev/release/run-session.R`, `sessions.log`), lane library:

    frmtmb forward files=337/350 pass=25670 fail=0 err=0 skip=243 warn=0
    frmtmb.sample forward files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.sample reverse files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.sample preload files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0

27 LEAK lines, every one tagged `[with load of ...]` or `[extension
session state]` (none from a file this lane touched). `files=` counts
files with a result row; the pass totals equal the one-file-per-process
run's. Core reverse and preload were not run (each is a 350-file session
of about an hour on this load); frmtmb.sample is the only extension
whose tests write `gp(`.

The one-file-per-process totals (`all-lane.log`):

    frmtmb files=350 pass=25670 fail=0 err=0 skip=243 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1851 fail=0 err=0 skip=5 warn=0
    frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.spatial files=8 pass=551 fail=0 err=0 skip=6 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

The skips are the gates of an ungated run (brms fits, Stan). One file
printed the new message: frmtmb.sample's test-brms-draws-order.R, a
`suppressWarnings(frm(y ~ gp(x, by = f, k = 5)))` whose two sub-GPs end
white; it now muffles `frmtmb_boundary_fit` there, and reran `pass=39
fail=0 err=0 skip=0 warn=0` with no message (`t3-lane.log`). No
extension's R code changed, so none was reinstalled; each loaded the
lane's core.

R CMD check --as-cran of frmtmb.sample (one test file changed;
`dev/gpselimit-check-sample.sh`): `Status: 1 NOTE` (the V8 note); tests
`[ FAIL 0 | WARN 0 | SKIP 26 | PASS 4272 ]`. Both 00check.log files are
kept in `dev/gpselimit-out/`; `dev/gpselimit-check/` is deleted.

Bump and floors: minor (a BREAKING behavior change, warning to message).
No extension needs a new frmtmb floor: none calls a changed function, and
frmtmb.sample's change is to a test only.

Tool denials: a `git diff --no-index` on copies of the base files, and
an `awk` line-length scan right after it, were denied by the permission
classifier; neither was retried with git. The 80-column check was done
with the Grep tool instead (no line over 80 in this lane's hunks).

## 8. Not done

- No change for exact `gp()` (section 4): filed.
- `gp_basis_white()` asks every dimension; an anisotropic field with one
  white dimension is not covered (one such fit in `hyp-lane.txt`, with
  finite SEs on base and lane alike). Not filed: no evidence its SEs
  mislead.

## 9. Fix round 1

Review: dev/reviews/2026-11-17-gpselimit.md (verdict FIX). The review's
case, by-jitter and nl2 scripts were copied to `dev/gpselimit-rv-*.R`
(output to `dev/gpselimit-out/rv/`, so the review's logs stay) and run
on the fixed build.

### Changes

- FIX 1, `gp_white_ridge()` (R/gp-start.R): the way back. Along each
  ridge outward (toward larger length scales), every `se_edge_step` until
  the length scale reaches the gate (0.1 of the spacing; the smallest
  spacing for an isotropic block), the objective must not fall by more
  than `se_edge_tol`; otherwise FALSE with attribute `short`. Asked even
  where the way to 0 already failed, so the warning can say the fit
  stopped short. Cost: up to five more evaluations per ridge on a white
  fit at 3e-5 of the spacing; none on a fit under no gate.
- The warning for a short stop says "the likelihood still rises toward
  larger length scales, so the fit stopped short of its maximum: refit
  from a larger length scale (start =)", in place of "The fitted field
  and its predictions hold".
- FIX 3: `gp_basis_message()` gives the message only where the fit has
  no `frmtmb_conv_*` verdict. `check_convergence()` (R/fit.R) keeps its
  verdict classes in `fit$cache$conv_classes`, and `gp_end_warn()` moved
  from before the convergence check to after `latres_laplace_report()`
  in `fit_assembled()`, so the verdict exists when it runs. The gp()
  conditions now follow the convergence warnings instead of preceding
  them. (alapse also edits `check_convergence()`'s return line and the
  `se_explained` test; my two hunks are the lines just before the
  warning loop and after `latres_laplace_report()`.) The short check of
  the fix round (dev/reviews/gpselimit-short.md) measured the order on
  plain fits: the gp() conditions (nugget, Laplace, basis) now come after
  sdreport(), the convergence check, the fit-end checks and
  latres_laplace_report(), and before the SE check; before, they came
  ahead of sdreport(). No test depends on the order.
- Code-read gap, filed LOW (short check): with separation at a nonzero
  code, check_convergence() returns before it sets `conv_classes`, so a
  certified white field there could get the boundary message beside the
  separation warning. Not reproduced; no code change in this round.
- FIX 4 (R/par-template.R; nlband edits R/fit-restart.R only, and no
  lane edits par-template.R): `nl_fit_off_saddle()` holds a
  `frmtmb_gp_basis` message as it holds `frmtmb_limit_fit`, and
  `nl_fit_probed()` drops a replaced first fit's `frmtmb_gp_basis` when
  a sign probe or pole escape returns another fit; the returned fit says
  its own through `gp_end_warn()` in `frm()`. Not reproduced by a test:
  the review's 52 nonlinear fits never printed a stale message, and
  neither did the 24 of this round's nl2 run; I could not construct a
  replaced white first fit in the time I gave it. The hunk is a
  code-read fix.
- FIX 2: `gp_white_dirs()`'s comment now says it is load-bearing for
  `gp(by = )`; NEWS has a BREAKING bullet for it; section 1 and the
  backlog CLOSED line are corrected (one line).
- Tests narrowed: the wrappers this lane added are `allow_gp_basis()`, a
  new helper in tests/testthat/helper-warnings.R that muffles only the
  `frmtmb_gp_basis` message, in place of `allow_boundary()`, which would
  also hide a singular fit's message. frmtmb.sample's inline handler
  muffles `frmtmb_gp_basis` only.
- Filed: a LOW item for `gp(x1, k = 10) + gp(x1, k = 5)` sharing brms
  names.

### New tests (tests/testthat/test-gp-white-limit.R)

- "a stop in the white plateau short of a maximum keeps the warning":
  row 6 (resolved) started at 3e-5 and at 1e-3 of the spacing; where it
  stays white, one `frmtmb_gp_basis` warning saying it stopped short, a
  log-likelihood under the default fit's, and no "boundary" on theta.
  This also tests the "flat" of NEWS bullet 3.
- "the white sub-GPs of gp(by = f) have no standard error": row 3 with
  `f` alternating, k 10: every white hyperparameter NaN and "boundary",
  every other SE finite, the condition a message.
- "a white field at a doubted stop keeps the warning": the review's nl2
  seed 2 (code 1 beside a white field); where a `frmtmb_conv_*` warning
  occurs, every `frmtmb_gp_basis` condition is a warning.

Seen to fail: on the first hand-back's build (`fx1-oldlane.log`):
`pass=37 fail=5`, the plateau test (message, no "stopped short", theta
"boundary") and the doubted-stop test (message beside the convergence
warnings). On rellib-r22 (`fx1-base.log`): `pass=24 fail=14 err=2
warn=5`, including the gp(by = f) test (finite SEs, no "boundary", a
warning). On the fixed build (`fx1-lane.log`): `pass=42 fail=0 err=0
skip=0 warn=0`.

### Measurements on the fixed build

- The 30-fit sweep under 21 ulp draws (j = 0..20, `sweep-lane-fx.txt`,
  `cmp-fx.txt`): 273 of 273 white fits get the message, 357 of 357
  resolved fits nothing, and 630 of 630 are `identical()` to base in
  estimates, point, log-likelihood, every SE and the lost reasons.
- The review's cases (`rv/cases-lane.txt`, 81 lines, against the
  review's `cases-lane.txt`): only the five plateau stops (rows 6, 12,
  18, 24, 30 started at 3e-5) changed, from the message and "boundary"
  to the warning and "flat" (with the SE warning). The 36 local refits
  of group B that return to the limit keep the message; every other
  case (priors, 2-D, `by = f`, `sigma ~`, two terms, REML, poisson) is
  unchanged in conditions, SEs and lost reasons.
- gp(by = f) under 21 ulp draws (`rv/byjit-lane.txt`): rows 3 and 21,
  42 of 42 fits with all four white hyperparameters NaN and the message.
- nl2 seeds 1 to 6 (`rv/nl2-lane.txt`): seeds 2 and 4 of b0 -0.8, k 20
  (code 1, `frmtmb_conv_code`, seed 2 also `frmtmb_conv_gradient`) now
  get the `frmtmb_gp_basis` warning; the 8 other white fits get one
  message each beside the restart message.
- OpenBLAS 0.3.26, LAPACK (`ob0326-fx-summary.txt`), the 10 files of
  section 7, one process at a time: all `fail=0 err=0 warn=0`
  (test-gp-white-limit.R `pass=42`, test-gp-by.R `pass=71 skip=2`, the
  rest as in section 7). The sweep under the copy, j = 0..5
  (`sweep-lane-obfx.txt`): 78 of 78 white fits the message, 102 of 102
  resolved nothing. gp(by = f) under the copy (`rv/byjit-lane-ob0326.txt`):
  42 of 42 NaN with the message. `dev/ciharden-out/` deleted.

### Suites, sessions and checks on the fixed build

`dev/gpselimit-fx-chain.sh`, two R processes at most, every `lib:` line
reading wt-gpselimit-lib for frmtmb.

All 550 files ungated, one per process (`fx-all.log`, 550 of 550 ran):

    frmtmb files=350 pass=25688 fail=0 err=0 skip=243 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1851 fail=0 err=0 skip=5 warn=0
    frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.spatial files=8 pass=551 fail=0 err=0 skip=6 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

(core +18 over the first hand-back: the three new tests.) No log holds
the white-limit message text. The 45 gp() and SE files gated
(`fx-t1.log`): 45 of 45, every one `fail=0 err=0 skip=0 warn=0`, no
message text.

One session each (`sessions.log`):

    frmtmb forward files=337/350 pass=25688 fail=0 err=0 skip=243 warn=0
    frmtmb.sample forward files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.sample reverse files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
    frmtmb.sample preload files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0

27 LEAK lines, all tagged `[with load of ...]` or `[extension session
state]`, 0 untagged.

R CMD check --as-cran: frmtmb `Status: 1 NOTE` (V8, HTML manual), tests
`[ FAIL 0 | WARN 0 | SKIP 279 | PASS 25425 ]`; frmtmb.sample `Status: 1
NOTE` (V8), tests `[ FAIL 0 | WARN 0 | SKIP 26 | PASS 4272 ]`
(`00check-frmtmb-fx.log`, `00check-frmtmb.sample-fx.log`).
`dev/gpselimit-check/` deleted.

Bump unchanged: minor. Floors: none.
