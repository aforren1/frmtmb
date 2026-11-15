# Short check of lane gpfield's fix round (2026-11-13)

Build: wt-gpfield-lib install14 (frmtmb.rdb 17:32:58, after the last
R/ edit at 17:31:35). Scripts `dev/reviews/gpfield-short-*`, plus the
review's `gpfield-review-predict.R` and `gpfield-review-cost.R`; logs
`dev/reviews/gpfield-short-out/`. At most 2 R processes at a time, no git
writes, no lane file edited. The OpenBLAS copy of R was built in the
session scratch directory from `dev/ciharden-openblas.sh` and deleted
afterwards. No run failed on disk space or allocation.

## Verdict: SHIP

F1, F2 and F3 are fixed. Three numbers in the record and NEWS are a
little optimistic (n1 to n3 below). Correct them at consolidation; none
blocks the merge.

## F1: the gp() warnings on the fit frm() returns

`gp_end_warn(fit)` runs in the `announce_start` block (R/fit.R:1984) and
around `nl_fit_probed()` in `frm()` (R/fit.R:1143). The flag is
`fit$cache$gp_end_said`. `gpfield-short-f1.R` (`f1.txt`):

    F1 sweep:mm_V:bernoulli:40:1    conv 0 | gp warnings: frmtmb_gp_laplace x1   (the review's reproducer; probe refit)
    F1 sweep:exp_a:bernoulli:120:2  conv 0 | gp warnings: frmtmb_gp_laplace x1   (nl, first fit kept)
    F1 gate:approx_a                conv 0 | gp warnings: frmtmb_gp_basis x1     (nl, approximate)
    F1 ng:student:step:80:1         conv 0 | gp warnings: frmtmb_gp_laplace x1   (no body: fit-end block only)

Every warning is said exactly once on every route.

- The non-nl route returns `run(start, TRUE)` before the wrapper.
- On the nl route the wrapper sees the flag the fit-end block set.
- The autoscale choice runs both arms announced, inside its capture, and
  replays only the kept arm. The flag lives in the kept fit's cache
  environment, which list copies share.

saveRDS() keeps the flag (`read fit flag: TRUE`), and that is correct:
the saved fit already warned. It does not suppress a legitimate warning
on a refit made from the saved fit:

    F1 update(readRDS fit) | gp warnings: frmtmb_gp_laplace x1      (update() builds a new fit and cache)
    F1 update(readRDS st80) | gp warnings: frmtmb_gp_laplace x1
    F1 gp_end_warn twice on read fit: 0 warnings

Paths that pass neither call:

- **Quadrature and importance fits.** `gp_end_warn()` returns early. This
  is by design: 0.82.0's nugget warning had the same gate (`is.null(imp)`).
- **refit().** `refit()` (R/sugar.R:261) calls `fit_assembled()`
  unannounced, so it gives no gp() warning:

      F1 refit(readRDS fit, same y) | gp warnings: none | messages on it now: 1

  This is the settled design for refit(), influence and bootstrap (as
  for `nl_flat_message()` and the nugget warning), not a regression. I
  note it so nobody reads refit() output as checked.

## F2: predict()'s cost on an exact field

**Kriging law, rerun on install14** (`law-eq-4000-s7.txt`,
`law-m32-2000.txt`): identical, line for line, to the review's. exp_quad,
x 10.5 to 60: z mean -0.010 to +0.006, sd 0.992 to 0.998, KS p 0.649 to
0.989, pair correlation 0.997 against 0.997. matern32: sd 0.966 to
0.998, pair 0.973 against 0.975.

**Cost** (base and lane in 2 processes at once, same seed; per-draw time
from the second call; `cost-*-g500.txt`, `cost-*-p300.txt`):

| fit | call | base s/draw | lane s/draw | ratio |
|---|---|---|---|---|
| gaussian n 500 | newdata 6 rows | 0.0625 | 0.0553 | 0.88 |
| gaussian n 500 | no newdata | 0.0232 | 0.0306 | 1.32 |
| poisson n 300 | newdata 6 rows | 0.0112 | 0.0176 | 1.57 |
| poisson n 300 | no newdata | 0.0054 | 0.0116 | 2.15 |

The predictions' sds equal the review's tape-route values in every
column.

**Newton steps against the tape's inner mode** (`gpfield-short-f2.R`,
`f2-*.txt`). For 40 replicates per fit: the steps the Newton route needs,
whether it falls back, and |z_newton - z_tape| in posterior sds of z (the
inverse inner Hessian at the replicate):

| fit | fell back | steps 1/2/3/4/5/6/7 | max err | median err |
|---|---|---|---|---|
| gaussian n 200 | 0 | 40/0/0/0/0/0/0 | 5.3e-8 | 2.1e-8 |
| poisson n 300 | 0 | 0/4/33/3/0/0/0 | 2.0e-3 | 2.0e-5 |
| poisson low counts n 200 | 0 | 0/3/35/2/0/0/0 | 5.2e-4 | 3.8e-6 |
| bernoulli n 200 | 0 | 0/3/25/11/1/0/0 | 2.8e-4 | 7.3e-6 |
| bernoulli n 60 | 0 | 0/2/35/3/0/0/0 | 1.5e-4 | 1.4e-6 |
| binomial(20) n 80 | 8 | 0/1/7/6/15/1/2 | 1.2e-3 | 1.9e-5 |

**The stopping rule.** The binomial field (trials 20, `3 sin(x)`) does
not settle in 8 steps on 8 of 40 replicates (15 of 40 on another seed,
`f2b.txt`), and those fall back to the tape. The fallback is the tape
route exactly (`gpfield-short-f2b.R`, `f2b.txt`): with the steps capped
at 1, every replicate falls back, and the draws are identical to those of
a run with the Newton plan removed:

    F2B bern200    steps<=1: 40 calls, 40 fell back | steps<=1 == tape TRUE
    F2B binom80    default: 40 calls, 15 fell back | steps<=1 == tape TRUE
    F2B student80  default: 40 calls, 40 fell back | steps<=1 == tape TRUE   (curvature not positive at outliers)

**The factor cache.** It does not leak. After predict() no block of the
fit carries `gp_factor_cache`, and `object.size(fit)` is unchanged on all
5 fits. predict() with the same seed is identical before and after a
predict() with other ndraws and newdata (`F2REPRO ... TRUE TRUE`, 5
fits). The cache key is the exact theta, on the replicate's copy of the
block only.

Notes for the record (not blocking):

- **n1.** The accuracy claim "at most 5.5e-4 posterior sds" is from one
  poisson fit and 20 draws. Here the maximum is 2.0e-3 (poisson n 300)
  and 1.2e-3 (binomial). Still under 0.2 percent of a posterior sd.
  Restate it as about 2e-3.
- **n2.** NEWS says "0.9 to 1.8 times as long as before". The poisson
  no-newdata row read 2.15 here (0.58 s against 0.27 s for 50 draws,
  small and noisy). Say "0.9 to about 2".
- **n3.** On binomial fields with many trials the Newton route falls back
  to the tape on 20 to 38 percent of replicates. Those replicates cost
  what install7 cost. Record it beside the cost.

## F3 and the minor items

- **F3.** Both warning bullets in NEWS are now BREAKING.
- **Basis wording.** The warning text, `?frmtmb-gp` and NEWS all say
  "every basis function gets the same prior weight". NEWS gives the rate
  (13 of 30). The only "white noise over the basis" left is the internal
  `@noRd` comment of `gp_basis_white()` (R/gp-start.R:1722), which is
  harmless.
- **m1.** The tape route of `gp_field_draw()` and `gp_end_messages()` use
  `obj_state_save()` / `obj_state_restore()`.
- **m3.** The inner-Hessian branch now runs only where every random
  effect sits in a predictor that holds an exact field
  (`all(ak %in% gk)`, R/gp-start.R:1899). `?frmtmb-gp` states each
  branch's scope.
- **m4.** test-gp-by.R `fm` allows only "spacing its basis resolves".
- **m6.** The 5319 item (dev/test-backlog.md:4824) points at parity2's
  `draws_latres_cols()` and silent2's `draws_brms_order()`, after both
  merge.
- **The new OpenBLAS allowance in test-gp-field.R** ("Optimizer did not
  report", F1 test). It is allowed by text only; `require` is
  "first-order correction". It is true under OpenBLAS 0.3.32 in LAPACK
  mode with `OPENBLAS_CORETYPE=Haswell` (`gpfield-short-obf1.R`):

      obf1-0332.txt: OBF1 conv 1 | W frmtmb_conv_code: Optimizer did not report convergence: false convergence (8). A Newton step from the stop gains 1e-09 ... | W frmtmb_gp_laplace: ... 60.2 ...
      obf1-ref.txt:  OBF1 conv 0 | W frmtmb_gp_laplace: ... 60.2 ...

  The code-8 note is the accurate at-the-maximum note. The Laplace
  warning is said once on both BLAS. test-gp-field.R under that copy
  (`ob0332/test-gp-field.txt`): `pass=21 fail=0 err=0 skip=0 warn=0`.

## Tests on reference BLAS (gated tier, one file per process, 2 at a time)

`tests/*.txt`, every `lib:` line wt-gpfield-lib:

    RESULT test-gp-field.R pass=22 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-kernels.R pass=39 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-predict-draw.R pass=7 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-by.R pass=86 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-start.R pass=251 fail=0 err=0 skip=0 warn=0
