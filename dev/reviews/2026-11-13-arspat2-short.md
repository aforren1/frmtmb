# Short check of lane arspat2's fix round (round of 2026-11-13, base ec8178d3)

Scope: the fix round's items only (F1, F2, F3, F4 of
`2026-11-13-arspat2.md`). Lane library `wt-arspat2-lib` (core and
frmtmb.spatial; every log's `lib:` line says so), base `rellib-r21`.
Scripts `dev/reviews/arspat2-short-*.R|sh`, logs
`dev/reviews/arspat2-short-out/*.txt` (rds files deleted). Base compared
by `git rev-parse ec8178d3:<path>` against `git hash-object`; no `git
diff`, no git writes. At most 2 R processes; largest probe 2.88 GB (base
fcor n = 400, m = 100). Free space on C: 134.8 GiB at start, 136.0 GiB at
end.

## Verdict: FIX (NEWS text only)

The code is right on every check: exact to 6.7e-12 relative or better,
read values stay latent and are right, no-mi() fits are bit-identical,
the interval keeps 0.82.0's coverage, and the six test files pass.
Three NEWS sentences are false or overstate. The fixes are text only (or
one small code change, optional, at N2):

- N1. The core NEWS says the margin is "the fit over the observed rows
  with the term's matrix cut to them". This is true for `fcor()` only. For
  `sar()` it is false (see item 7).
- N2. "a slot that disagrees with `whiten()` is refused by the term's
  name". This overstates: the check is a spot check, and item 3 shows a
  wrong slot that the check accepts.
- N3. "0.3 to 3.1 nats low" and "1.1 to 2.0 nats low" are the lane's
  cases only. The measured range runs to 4.66 nats (reviewer, fcor
  m = n/2) and 45.17 nats (sar n = 400, m = 300, item 4). For gaussian
  with RE in sigma, it runs down to 0.21 nats.

## 1. Exactness: PASS

`arspat2-short-exact.R`, log `exact-lane.txt`. The objective of the `mi()` fit is
compared with the observed-rows `fcor(M_oo)` fit (fcor), or with a dense
reference written in the script (sar). The reference is the t margin with
`Sigma_oo = sigma^2 [(A'A)^-1]_oo` and mean `(B (X beta + Z b))_o`, with b
integrated by RTMB's Laplace. Three outer points; nu set to 3, 6, 21 (or
2, 5, 20). AD gradient against `numDeriv::grad` at the optimum + 0.1
(scaled).

| case | m | unread / latent | relative diff (3 points) | AD vs FD |
|---|---|---|---|---|
| fcor (1\|h), student | 9 | 9 / 0 | 1.5e-16, 0, 1.5e-16 | 1.2e-9 |
| fcor (1\|h), profile | 9 | 9 / 0 | 1.5e-16, 0, 1.5e-16 | 2.0e-9 |
| fcor, nu 1.02 / 1.1 / 1.5 (data df 1.2) | 7 | 7 / 0 | 1.2e-16, 2.4e-16, 0 | 7.0e-10 |
| fcor m = n/2 = 30, (1\|h) | 30 | 30 / 0 | 2.5e-16, 2.5e-16, 2.6e-16 | 2.3e-10 |
| sar lag, binary rook 7x7, (1\|h) | 15 | 15 / 0 | 5.2e-12, 6.7e-12, 4.3e-16 | 4.2e-9 |
| sar error, row-std rook 8x8, (1\|h) | 16 | 16 / 0 | 1.9e-12, 1.7e-13, 1.3e-14 | 6.4e-10 |

Every fit has convergence 0. R/fit.R is the same blob as at ec8178d3
(b68a5bdb). `autocor_seam_t_missing`, `autocor_seam_t_laplace`,
`autocor_seam_scale_integrated` and `mi_inner_beta` occur nowhere in
core or extension R/.

## 2. A value another density reads: PASS

`arspat2-short-read.R`, logs `read-lane.txt` and `read-base.txt`. The base arm
evaluates its own objective at the lane's points. Design: `y | mi() ~ x +
fcor(M)` and `z | mi() ~ mi(y) + x` (gaussian), rescor off, n = 40. y is
missing at 3, 11, 19 and 27 (z observed there, so read) and at 6, 14 and 33 (z
missing there too). z is also missing at 9 and 22.

- Routes, lane: y latent rows 3, 11, 19, 27; y unread 6, 14, 33; z
  unread 6, 9, 14, 22, 33; `length(miss)` 4. Base: all 7 values of y
  latent.
- B1 gaussian: the lane equals two references to 4.1e-16. The first is the
  Laplace over the read values of [y density on the rows other than 6, 14,
  33, M cut to them] x z's density. The second is the dense joint normal
  margin of (y_o, z_o). Lane - base objective: +1.4e-14, 0, -1.4e-14 (base
  holds all values latent, which is exact on gaussian). logLik is the
  same to 10 digits.
- B2 student y: the lane equals the Laplace reference (exact margin over
  the unread rows, Laplace over the read ones) to 7.5e-13, 6.0e-15 and
  1.1e-12. Base is 2.1e-3 to 3.2e-3 relative from that reference (lane -
  base -0.21, -0.19, -0.14). The 0.82.0 route is not the reference there.
  The two logLiks agree at the optimum only because nu runs toward the
  gaussian limit there.
- B3: a registered reader (`z | mi() ~ mi(y) + x + fcor(M)`) keeps every
  y row latent (3, 6, 11, 14, 19, 27, 33; none unread). z's own missing
  values are unread. Lane - base: within 1.4e-14.
- B4: `ar(t, cov = FALSE)` reads y's lag. All 7 values stay latent, and
  lane - base is 0 at all 3 points.
- B5: `ar()` and `fcor()` on one response are refused on both builds
  ("carries 2 residual correlation terms"). A lagged ar reader and a
  registered term cannot be on one response.

## 3. The margin slot: PASS. The check is a SPOT CHECK, not a proof

`arspat2-short-slot.R`, log `slot-lane.txt`.

- (a) A slot that forgets sigma in G is refused by the term's name: "The
  'myfix' term's margin() slot is wrong for myfix(V): ...". The same slot
  without a missing value fits, because the slot is not read.
- (b) The fallback is the term without a slot (`myfix`, one `whiten()` per
  missing value; the slot is NULL on the entry). It is compared with the
  observed-rows `fcor(M_oo)` fit, n = 60, 9 missing, (1|h), on gaussian and
  student, with one sigma and with `sigma ~ x`. The largest relative
  difference over the 12 points is 3.2e-16, and logLik is equal to 8
  digits.
- (c) The build-time check compares G only at entries {1, m} x {1, m} and
  v only at {1, m}, at one point. A sar slot that writes
  `W_mm + W_mm'` as `2 W_mm` is right only for a symmetric W. The test
  uses it on a row-standardized 7x7 rook W, missing rows 1, 2, 8, 9, 17,
  24, 25, 49 (first and last not neighbors; interior asymmetry 0.167). It
  is ACCEPTED, and the objective at the right fit's optimum is off by
  0.0110 (logLik -59.0374 against -59.0265). The roxygen ("checked at two
  columns and a generic point") states this scope. The NEWS sentence
  (N2) does not. The two in-tree slots are exact (item 1). Optional
  code fix: check every column on numbers once per objective. This
  costs m `whiten()` calls off the tape, and it makes the check exact
  at that point.

## 4. Cost: PASS

`arspat2-short-cost.sh` / `-cost.R` (the lane's design, one fit per
process, sequential, lane and base interleaved), log `cost.txt`. The fcor
base run at n = 400, m = 200 is not run (5.3 GB, above this check's 4 GB
limit), and n = 1000 is not run.

| n | m | term | base s | lane s | ratio | base peak GB | lane peak GB | base logLik | lane logLik |
|---|---|---|---|---|---|---|---|---|---|
| 400 | 100 | sar | 0.53 | 0.61 | 1.15 | 0.37 | 0.38 | -539.2286 | -531.8415 |
| 400 | 200 | sar | 1.09 | 1.82 | 1.67 | 0.43 | 0.50 | -365.1112 | -338.4426 |
| 400 | 300 | sar | 2.01, 2.10, 2.13 | 4.07, 4.14, 4.04 | 1.90 to 2.06 | 0.49 | 0.64 | -230.5469 | -185.3720 |
| 400 | 100 | fcor | 9.30 | 2.32 | 0.25 | 2.88 | 0.64 | -539.2497 | -531.8628 |
| 400 | 200 | fcor | (lane's record 26.20, 5.34 GB) | 3.89 | | | 0.80 | | -339.4116 |
| 400 | 300 | fcor | (lane's record 48.97, 7.06 GB) | 7.20 | | | 1.12 | | -185.9245 |

Every lane fit has convergence 0 and a finite gradient. The machine
carried less load than in the lane's run: absolute times are about 0.6x
the lane's, and the ratios agree. The borderline cell, sar m = 300, is
confirmed and is slightly worse here than the lane's figures: 1.90 to
2.06x time (lane 1.73 to 1.95x) and +31% peak, +0.15 GB (lane +22 to 27%).
This is acceptable. The cost is about 2 s and 0.15 GB on a fit whose 0.82.0
logLik is 45.17 nats wrong. At m = 100 to 200, the overhead is 1.15 to
1.67x. fcor is 4x faster and 4.5x smaller.
The NEWS memory figures (1.1 GB, 0.6 against 0.5 GB) match. The NEWS
seconds are from a busier machine (12 s against 7.2 s here for fcor
m = 300). The ratios hold.

## 5. F1, rho's interval: PASS

`arspat2-short-cover.R`, logs `cover-lane.txt` and `cover-base.txt`. The
lattice is a row-standardized rook 5x5, lag, true rho 0.9, seeds 1 to 100.

- Type `"clamped"` on the lane (`"raw"` on base). `summary()$cor_pars`
  equals `confint_varcorr()` on both builds.
- Coverage is 0.960 on the lane and 0.960 on base. The upper end is past
  hi in 68 of 100 base fits, and in those 68 the lane's upper end is
  exactly hi. The lane interval equals `pmin(pmax(base raw, lo), hi)`
  with max |diff| 0 (lwr and upr). Estimates and Est.Error are
  identical.
- `frm_multiple()` with 4 imputations (seed 3, y + N(0, 0.1)): type
  clamped, pooled estimate = mean of the rho estimates (to 1e-10). The
  Rubin interval on rho's own scale is (0.629439, 1.052236), and the
  reported interval is (0.629439, 1.000000) = clamp(raw), diff 0.
- test-sar.R "rho's interval is the raw Wald interval held inside its
  range" pins type `"clamped"`, `upr == hi` where the raw end passes it
  (seed 3), and the raw interval inside (seed 1). summary() and the SE
  are checked against the delta method. No test covers `frm_multiple()`'s
  clamp. This is not required.

## 6. Fits without mi() on a registered term: PASS (bit-identical)

`arspat2-short-ident.R`, both arms pinned to logical CPU 0 (`start
/affinity 1`), logs `ident-lane.txt` and `ident-base.txt`. The fits are fcor
gaussian; fcor student + (1|h); sar lag gaussian; sar error student + (1|h);
gaussian REML `y | mi() ~ x + (1 | h)`; and mv `y ~ x + fcor(M)` + `z | mi() ~ x`.
For all 6, the following are `identical()`: `opt$par`, logLik,
`cov.fixed`, the objective and gradient at 3 points, `fixef()`, and the
messages. `confint_varcorr()` is identical except the sar rows' type
(clamped / raw). Their numbers are identical.

## 7. NEWS: FIX (N1 to N3), otherwise PASS

- BREAKING bullets are present: `fit$estimates$miss` no longer holds the
  values that leave the fit, and the interval type changes from
  `"raw"` to `"clamped"`. No lane names, `dev/` paths or version numbers
  (core and spatial).
- N1: `arspat2-short-news.R`, log `news-lane.txt`. sar, 7x7 row-std, 12
  missing, gaussian. The mi() fit equals the dense margin at its optimum
  (lag -53.637158, error -54.820154). `sar(W_oo)` on the observed rows
  gives -55.203983 (lag, rho 0.10 against 0.50) and -55.083296 (error,
  rho 0.21 against 0.39). Reword: "the fit over the observed rows; for
  `fcor()`, the fit with `M` cut to them". The fcor compat note and
  `?frmtmb-autocor` say this for `M` only and are right.
- N2, N3: see the verdict.
- Cost numbers: true within machine load (item 4).
- LOW: the spatial NEWS bullets change `estimates$miss` and the
  `lagsar`/`errorsar` interval type, but they are not marked BREAKING
  (core marks both).
- Spatial floor: the F4 section of the findings says consolidation MUST raise
  `Depends: frmtmb (>= ...)` (DESCRIPTION still says `>= 0.82.0`), and
  the spatial NEWS says it needs the frmtmb that takes `margin =`. This is
  right for consolidation.

## 8. Tests on the lane library: PASS

`arspat2-short-test.sh` (`NOT_CRAN=true`, one file per process, two at
a time). Logs are `test-test-*.txt`. Each `lib:` line names
wt-arspat2-lib for both packages.

| file | pass | fail | error | skip | warn |
|---|---|---|---|---|---|
| test-autocor-se-mi.R | 85 | 0 | 0 | 0 | 0 |
| test-autocor-seam.R | 90 | 0 | 0 | 0 | 0 |
| test-autocor.R | 198 | 0 | 0 | 0 | 0 |
| test-mi-families.R | 194 | 0 | 0 | 0 | 0 |
| frmtmb.spatial test-sar.R | 177 | 0 | 0 | 4 | 0 |
| frmtmb.spatial test-car.R | 132 | 0 | 0 | 0 | 0 |

These match the lane's fix-round table.
