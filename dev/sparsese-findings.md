# Lane sparsese: prediction SEs above joint_root_max by sparse solves

Round opened 2026-11-16, base b255350e (frmtmb 0.83.0, rellib-r22).
Item 4952 (HIGH) of `dev/test-backlog.md`.

## Status

Done. Item 4952 fixed; brief items 1 to 5
done, 5e4 groups not run (over 4 GB, section 5). One user decision filed
(section 7).

## 1. Reproduction on rellib-r22

`Rscript dev/sparsese-3370ref.R base 4500 <seed>` (REML, degree-7 raw
polynomial, `(1 | g)` with 4500 groups of 2 rows, 4510 joint parameters;
logs `dev/sparsese-out/ref-base-4500-s<seed>.txt`):

Standard errors at x = 1, 1.5, 2 (the full lines are in section 4):
- Seed 1, nothing lost. Conditional 0.02435989850 0.00688185325
  0.02522057232; fitted() 0.02327294104 0.00689557714 0.02427105025
  (-4.46, +0.20, -3.76 percent); a' V a from vcov() the same as fitted().
- Seed 2, theta_1 lost ("boundary"). Conditional 0.026547772782
  0.006987372076 0.024168041760; fitted() NaN NaN NaN with the warning
  "above 4000 joint parameters"; a' V a from vcov() 0 0.006617954497 0,
  and every diagonal entry of vcov() negative, with no warning.
- Seed 3, theta_1 lost. Conditional 0.026203605934 0.006970559395
  0.025911831185; fitted() NaN NaN NaN with the same warning; a' V a
  from vcov() 0.028089953209 0.006961031483 0.027369759315 (+7.2,
  -0.14, +5.6 percent), with no warning.

A fitted SE adds the variance parameters' propagation to the conditional
one, so it must not fall below it; seed 1 does by 4.46 percent at x = 1,
the item's number. vcov() of the same fits is wrong too: it reads
joint_cov_whole(), whose factor sdr_unwhiten() does not form above the
constant, and inverts the natural-coordinate precision.

## 2. The change

Route: a REML fit on a whitened design (`fit$beta_whiten`) reads its
joint covariance by sparse solves in the whitened coordinates whenever it
has more than joint_root_limit() joint parameters
(`options(frmtmb.joint_root_max)`, default 4000), or whenever its stored
factor carries some rows only (jcs_white_wanted()).

- R/joint-sparse.R. jcs_joint_cov() takes `white` (default
  jcs_white_wanted()) and `whole`. With `white` it factors the inner block
  `Qg[r, r]` of the precision sdr_unwhiten() was handed
  (`sdr$joint_precision_white`; Matrix::Cholesky, supernodal, permuted),
  forms `Jg = -Qg_rr^-1 Qg_rf` by solves, and keeps `M`, the whitener at
  joint positions (natural = `M` gamma; joint_whitener()). `Vf` is the
  check's repaired covariance where it kept one (and its factor `Rf`
  where se_reparam_cov() kept one), else the inverse of the Schur
  complement, which is the whole inverse's outer block. `whole` skips the
  repair, for joint_cov_whole(). Every reader carries its loads to the
  stored coordinates first, `A M[pos, ]` (jcs_white_load()), so
  `a' V a = (M' a)' V_gamma (M' a)`: this is the item's `a' Qrr^-1 a`
  through `M' a`, plus `J Vf J'` through `Jg`. `[.frm_joint_sparse` forms
  a natural block as `M[i, ] V_gamma[ci, cj] M[j, ]'` (jcs_block() is the
  old body). The size refusal (jcs_check_size()) words "this REML fit"
  instead of "a latent residual" for such an object.
- R/predict.R. get_joint_cov() takes the sparse route for such a fit
  before the dense ones. lp_basis_out() keeps the sparse object for a
  whitened fit (lb_rowvar() then reads it through the loads), where it
  formed the natural block for fewer than 200 inner positions.
  joint_cov_white_inner() returns `too_big` (the 0.83.0 NaN rows) only
  where jcs_white_wanted() and the sparse route returned NULL, so the
  dense attempt is not repeated on a block that did not factor.
- R/se-reparam.R. sdr_unwhiten() compares with joint_root_limit() (the
  option, read when the fit forms its sdreport, at its first standard
  error) instead of the constant. Above it, it keeps the rows of the
  factor at every position but the group effects (`beta` and the outer
  parameters), from the sparse Cholesky factor of the whole `Qg`
  (joint_root_rows(): rows of `P' L^-T` by solves, made square by
  square_root_rows(), in a matrix with a row per joint parameter, NaN at
  the rows not carried, attribute `rows`). Those rows are what
  hypothesis(), the draws (fit_draw_space()) and interop's
  coef_vcov_root() read as the factor; `diag.cov.random`'s beta entries
  and the "lossy" attribute come from them as below the limit.
  joint_cov_whole() returns the sparse whole-inverse object as `V`
  beside those rows, and vcov() reads that `V` (corrected in fix round
  1: this paragraph first said vcov() read the rows).
  joint_whitener() is the sparse whitener sdr_unwhiten() built inline
  (now also with `T` for the readers).
- R/se-check.R. se_pred_warn()'s "size" warning no longer advises
  raising the option (that reaches the same failed factorization) and
  says why: the inner block does not factor in the design's coordinates.

## 3. Agreement with the dense route (3510 joint parameters)

`Rscript dev/sparsese-agree.R 3500 <seed> <group sd>`: one data set fitted
twice in one process, under `frmtmb.joint_root_max = 1e6` (the dense
factor and joint_cov_white_inner()) and `= 100` (the sparse route); the
estimates are identical() in every run. Largest relative difference per
reader, then the time of each reader's first call, dense / sparse (logs
`dev/sparsese-out/agree-3500-s<seed>-g<0|3>.txt`):

```
agree-3500-s1-g0: lost []
  fitted_newdata_noRE n=11 rel 3.577e-10 dense 8.50 s sparse 0.02 s
  fitted_newdata_RE n=11 rel 4.943e-11 dense 0.41 s sparse 0.00 s
  fitted_insample n=7000 rel 1.962e-10 dense 157.52 s sparse 0.06 s
  linpred_se n=11 rel 4.943e-11 dense 0.16 s sparse 0.02 s
  vcov n=64 rel 7.852e-15 dense 8.21 s sparse 0.00 s
  vcov_full n=4 rel 1.994e-14 dense 7.83 s sparse 0.00 s
  draw_space_V n=100 rel 1.994e-14 dense 15.58 s sparse 0.00 s
  draw_space_rootrows n=100 rel 2.229e-13 dense 16.09 s sparse 0.02 s
  hypothesis n=3 rel 7.260e-14 dense 8.07 s sparse 0.14 s
agree-3500-s1-g3: lost []
  fitted_newdata_noRE n=11 rel 2.728e-10 dense 8.05 s sparse 0.01 s
  fitted_newdata_RE n=11 rel 7.822e-12 dense 0.36 s sparse 0.03 s
  fitted_insample n=7000 rel 1.402e-11 dense 172.82 s sparse 0.06 s
  linpred_se n=11 rel 7.822e-12 dense 0.20 s sparse 0.01 s
  vcov n=64 rel 8.912e-15 dense 8.47 s sparse 0.00 s
  vcov_full n=4 rel 4.569e-15 dense 8.53 s sparse 0.02 s
  draw_space_V n=100 rel 8.912e-15 dense 16.94 s sparse 0.01 s
  draw_space_rootrows n=100 rel 1.162e-14 dense 18.37 s sparse 0.06 s
  hypothesis n=3 rel 2.254e-14 dense 10.77 s sparse 0.16 s
agree-3500-s2-g0: lost []
  fitted_newdata_noRE n=11 rel 2.315e-10 dense 8.47 s sparse 0.02 s
  fitted_newdata_RE n=11 rel 9.333e-11 dense 0.36 s sparse 0.02 s
  fitted_insample n=7000 rel 1.970e-10 dense 157.37 s sparse 0.05 s
  linpred_se n=11 rel 9.333e-11 dense 0.19 s sparse 0.00 s
  vcov n=64 rel 1.067e-14 dense 8.27 s sparse 0.00 s
  vcov_full n=4 rel 3.358e-13 dense 7.78 s sparse 0.01 s
  draw_space_V n=100 rel 3.490e-13 dense 15.54 s sparse 0.01 s
  draw_space_rootrows n=100 rel 7.858e-13 dense 16.00 s sparse 0.02 s
  hypothesis n=3 rel 2.271e-14 dense 8.06 s sparse 0.19 s
agree-3500-s2-g3: lost []
  fitted_newdata_noRE n=11 rel 2.627e-10 dense 10.72 s sparse 0.03 s
  fitted_newdata_RE n=11 rel 5.243e-12 dense 0.50 s sparse 0.01 s
  fitted_insample n=7000 rel 2.042e-11 dense 216.19 s sparse 0.08 s
  linpred_se n=11 rel 5.243e-12 dense 0.22 s sparse 0.01 s
  vcov n=64 rel 5.753e-15 dense 10.67 s sparse 0.01 s
  vcov_full n=4 rel 6.199e-15 dense 10.79 s sparse 0.00 s
  draw_space_V n=100 rel 6.366e-15 dense 20.42 s sparse 0.02 s
  draw_space_rootrows n=100 rel 1.369e-14 dense 19.89 s sparse 0.02 s
  hypothesis n=3 rel 7.139e-14 dense 10.00 s sparse 0.16 s
agree-3500-s3-g0: lost []
  fitted_newdata_noRE n=11 rel 5.522e-10 dense 7.80 s sparse 0.01 s
  fitted_newdata_RE n=11 rel 3.237e-11 dense 0.37 s sparse 0.03 s
  fitted_insample n=7000 rel 1.367e-10 dense 172.70 s sparse 0.06 s
  linpred_se n=11 rel 3.237e-11 dense 0.22 s sparse 0.02 s
  vcov n=64 rel 1.361e-14 dense 8.62 s sparse 0.02 s
  vcov_full n=4 rel 2.871e-13 dense 8.71 s sparse 0.00 s
  draw_space_V n=100 rel 2.964e-13 dense 16.91 s sparse 0.02 s
  draw_space_rootrows n=100 rel 9.163e-14 dense 18.61 s sparse 0.11 s
  hypothesis n=3 rel 3.982e-14 dense 11.05 s sparse 0.10 s
agree-3500-s3-g3: lost []
  fitted_newdata_noRE n=11 rel 2.107e-10 dense 10.63 s sparse 0.03 s
  fitted_newdata_RE n=11 rel 1.257e-12 dense 0.48 s sparse 0.02 s
  fitted_insample n=7000 rel 7.966e-12 dense 216.78 s sparse 0.06 s
  linpred_se n=11 rel 1.257e-12 dense 0.24 s sparse 0.02 s
  vcov n=64 rel 4.780e-15 dense 11.24 s sparse 0.02 s
  vcov_full n=4 rel 3.808e-15 dense 10.76 s sparse 0.02 s
  draw_space_V n=100 rel 4.780e-15 dense 20.06 s sparse 0.02 s
  draw_space_rootrows n=100 rel 1.134e-14 dense 19.53 s sparse 0.03 s
  hypothesis n=3 rel 6.558e-14 dense 9.70 s sparse 0.12 s
```

In all six runs nothing was lost (the whole-inverse route). Largest
relative difference over the four standard-error readers (fitted() on new
rows with and without group effects, in-sample fitted(), frm_linpred()):
5.522e-10; over vcov(), vcov(full = TRUE), fit_draw_space() and
hypothesis(): 7.858e-13 (both generated from the logs above with
`awk '{print $6}' | sort -g | tail -1`). The brief's bound is 1e-8.

The repaired route (sd lost) at 110 joint parameters
(`Rscript dev/sparsese-agree.R 100 <1|2> 0`, logs
`dev/sparsese-out/agree-100-s<1|2>-g0.txt`):
fitted() 2.3e-10 and 6.1e-10, in-sample 7.2e-10 and 1.1e-9, vcov() 1.1e-15
and 2.0e-15. At 4510 (section 4) the sparse route is within 3.2e-9 of the
exact conditional standard errors, where the dense route is not formed.

Cost of the dense route itself, measured on the way: at 3510 joint
parameters it took 7.8 to 11.2 s for every vcov() call (joint_cov_whole()
forms `root root'`, 3510^2, on each call), 157 to 217 s for an in-sample
fitted(), and 8.5 s for the first fitted() on new rows; the sparse route
took 0.00 to 0.19 s for each. See section 7 for what this suggests.

## 4. The item's three seeds, before and after (4510 joint parameters)

`Rscript dev/sparsese-3370ref.R <base|lane> 4500 <seed>`, default limit.
Lines copied from the logs `dev/sparsese-out/ref-<arm>-4500-s<seed>.txt`:

```
seed 1 base:
  lost:  
  conditional SE: 0.02435989850 0.00688185325 0.02522057232 
  fitted(): 1.14 s; R heap max used 458 MB; working set before 570 MB, peak 
    before 822 MB, peak after 873 MB
  fitted SE: 0.02327294104 0.00689557714 0.02427105025 
  SE / conditional - 1: -0.044621  0.001994 -0.037649 
  vcov SE (a' V a): 0.02327294104 0.00689557714 0.02427105025 
seed 1 lane:
  lost:  
  conditional SE: 0.02435989850 0.00688185325 0.02522057232 
  fitted(): 0.04 s; R heap max used 206 MB; working set before 517 MB, peak 
    before 821 MB, peak after 821 MB
  fitted SE: 0.024365044545 0.006882459922 0.025220634049 
  SE / conditional - 1: 2.113e-04 8.816e-05 2.448e-06 
  vcov SE (a' V a): 0.024365101009 0.006895274644 0.025253973546 
seed 2 base:
  lost: theta_1 boundary 
  conditional SE: 0.026547772782 0.006987372076 0.024168041760 
  fitted W: 3 of 3 predictions read coefficients whose covariance this REML 
    fit cannot give accurately above 4000 joint parameters, so their standard er
  fitted(): 1.59 s; R heap max used 653 MB; working set before 545 MB, peak 
    before 822 MB, peak after 951 MB
  fitted SE: NaN NaN NaN 
  SE / conditional - 1: NaN NaN NaN 
  vcov SE (a' V a): 0.000000000000 0.006617954497 0.000000000000 
seed 2 lane:
  lost: theta_1 boundary 
  conditional SE: 0.026547772782 0.006987372076 0.024168041760 
  fitted(): 0.03 s; R heap max used 207 MB; working set before 571 MB, peak 
    before 822 MB, peak after 822 MB
  fitted SE: 0.026547772746 0.006987372078 0.024168041838 
  SE / conditional - 1: -1.330e-09  3.119e-10  3.214e-09 
  vcov SE (a' V a): 0.026547605879 0.006989158954 0.024159400268 
seed 3 base:
  lost: theta_1 boundary 
  conditional SE: 0.026203605934 0.006970559395 0.025911831185 
  fitted W: 3 of 3 predictions read coefficients whose covariance this REML 
    fit cannot give accurately above 4000 joint parameters, so their standard er
  fitted(): 1.64 s; R heap max used 653 MB; working set before 504 MB, peak 
    before 822 MB, peak after 951 MB
  fitted SE: NaN NaN NaN 
  SE / conditional - 1: NaN NaN NaN 
  vcov SE (a' V a): 0.028089953209 0.006961031483 0.027369759315 
seed 3 lane:
  lost: theta_1 boundary 
  conditional SE: 0.026203605934 0.006970559395 0.025911831185 
  fitted(): 0.03 s; R heap max used 207 MB; working set before 571 MB, peak 
    before 821 MB, peak after 821 MB
  fitted SE: 0.026203605850 0.006970559408 0.025911831161 
  SE / conditional - 1: -3.224e-09  1.850e-09 -9.390e-10 
  vcov SE (a' V a): 0.026203834376 0.006975243774 0.025935794274 
```

Reading: seeds 2 and 3 (sd lost, so the conditional standard error is
the exact one) move from NaN with the size warning to within 3.214e-9 of
the exact values. Seed 1 (nothing lost) moves from 4.46 percent below the
conditional lower bound to 2.1e-4, 8.8e-5 and 2.4e-6 above it, which is
the variance parameters' propagation; the dense route on the same fit with
the limit raised before the fit (`cost-dense-4500-s1.txt`) reads
0.024365044547 0.006882459922 0.025220634044 against the sparse
0.024365044545 0.006882459922 0.025220634049. vcov() moves from negative
variances (seed 2) and +7.2 percent (seed 3) to the factor's values; its
"vcov SE (a' V a)" line is a' V a formed from the returned natural
matrix, which loses digits on this design by construction, and vcov() now
says so with its "a product such as a' V a" warning, as it did below the
limit.

## 5. Time and memory

`Rscript dev/sparsese-3370ref.R <arm> <groups> <seed> [limit]`, one process
per line, at most two at once. "dense" is the lane build with the limit
raised to 1e6 before the fit (sdr_unwhiten()'s dense factor and
joint_cov_white_inner()), "sparse" the lane build with the default limit
(or 100 at 3500 groups), "base" rellib-r22. fitted() is the first
fitted(re_formula = NA) on 3 new rows, which forms the joint covariance;
"heap" is gc()'s max used after it, "peak" the process's peak working set
(ps::ps_memory_info()) after it, which includes the fit's own peak.
Generated from the logs:

```
ref-base-4500-s1       p=4510   lost=[]
  fit 3.84 s, peak working set 822 MB
  fitted 1.14 s heap 458 MB peak 873 MB; vcov 0.43 s, peak 1007 MB
  SE/cond-1 -0.044621  0.001994 -0.037649 
ref-base-4500-s2       p=4510   lost=[theta_1 boundary]
  fit 4.22 s, peak working set 822 MB
  fitted 1.59 s heap 653 MB peak 951 MB; vcov 0.37 s, peak 991 MB
  SE/cond-1 NaN NaN NaN 
cost-base-10000-s1     p=10010  lost=[]
  fit 10.16 s, peak working set 1318 MB
  fitted 4.24 s heap 1536 MB peak 2430 MB; vcov 3.05 s, peak 3190 MB
  SE/cond-1 -0.0820360 -0.0002583 -0.0703402 
cost-dense-3500-s1     p=3510   lost=[]
  fit 3.21 s, peak working set 734 MB
  fitted 9.61 s heap 374 MB peak 763 MB; vcov 9.20 s, peak 763 MB
  SE/cond-1 0.0003948 0.0003718 0.0004138 
cost-dense-3500-s2     p=3510   lost=[]
  fit 3.23 s, peak working set 734 MB
  fitted 9.65 s heap 374 MB peak 778 MB; vcov 9.33 s, peak 778 MB
  SE/cond-1 1.278e-06 1.193e-04 3.113e-06 
cost-dense-4500-s1     p=4510   lost=[]
  fit 4.29 s, peak working set 789 MB
  fitted 20.20 s heap 499 MB peak 966 MB; vcov 21.70 s, peak 966 MB
  SE/cond-1 2.113e-04 8.816e-05 2.447e-06 
cost-dense-4500-s2     p=4510   lost=[theta_1 boundary]
  fit 5.03 s, peak working set 789 MB
  fitted 1.65 s heap 1120 MB peak 1416 MB; vcov 20.08 s, peak 1416 MB
  SE/cond-1 -1.320e-09  3.060e-10  3.128e-09 
cost-sparse-10000-s1   p=10010  lost=[]
  fit 9.77 s, peak working set 1318 MB
  fitted 0.09 s heap 226 MB peak 1318 MB; vcov 0.03 s, peak 1318 MB
  SE/cond-1 9.261e-06 2.174e-05 2.613e-05 
cost-sparse-10000-s2   p=10010  lost=[]
  fit 10.01 s, peak working set 1318 MB
  fitted 0.06 s heap 226 MB peak 1318 MB; vcov 0.03 s, peak 1318 MB
  SE/cond-1 7.646e-05 1.988e-05 8.850e-05 
cost-sparse-30000-s1   p=30010  lost=[theta_1 boundary]
  fit 36.04 s, peak working set 3335 MB
  fitted 0.16 s heap 301 MB peak 3335 MB; vcov 0.05 s, peak 3335 MB
  SE/cond-1 8.872e-09 -3.379e-09 -5.678e-10 
cost-sparse-3500-s1    p=3510   lost=[]
  fit 3.22 s, peak working set 733 MB
  fitted 0.03 s heap 202 MB peak 733 MB; vcov 0.01 s, peak 733 MB
  SE/cond-1 0.0003948 0.0003718 0.0004138 
cost-sparse-3500-s2    p=3510   lost=[]
  fit 3.28 s, peak working set 734 MB
  fitted 0.06 s heap 202 MB peak 734 MB; vcov 0.00 s, peak 734 MB
  SE/cond-1 1.278e-06 1.193e-04 3.113e-06 
ref-lane-4500-s1       p=4510   lost=[]
  fit 3.33 s, peak working set 821 MB
  fitted 0.04 s heap 206 MB peak 821 MB; vcov 0.00 s, peak 821 MB
  SE/cond-1 2.113e-04 8.816e-05 2.448e-06 
ref-lane-4500-s2       p=4510   lost=[theta_1 boundary]
  fit 3.88 s, peak working set 822 MB
  fitted 0.03 s heap 207 MB peak 822 MB; vcov 0.01 s, peak 822 MB
  SE/cond-1 -1.330e-09  3.119e-10  3.214e-09 
```

Reading:
- The sparse route's first fitted() takes 0.03 to 0.16 s from 3510 to
  30010 joint parameters, vcov() 0.00 to 0.05 s, and neither raises the
  process's peak above the fit's own (733 MB at 3510, 1318 MB at 10010,
  3335 MB at 30010). The dense route at 3510 took 9.6 s and 9.2 s, at
  4510 1.65 to 21.7 s with a peak up to 1416 MB.
- rellib-r22 at 10010 (seed 1, nothing lost) took 4.24 s and a peak of
  2430 MB for fitted(), 3190 MB for vcov(), and read 8.2 percent below the
  conditional lower bound (the item's 4.5 percent grows with size).
- 5e4 groups (50010 joint parameters) was NOT run: the fit itself peaked
  at 3335 MB at 30010 and grows by about 100 MB per 1000 groups
  (1318 MB at 10010), so a 5e4 fit would pass 4 GB before any standard
  error. That memory is TMB's fit, not this route. 30010 was run in its
  place. Neither base nor dense ran at 30010: each forms joint-size
  dense matrices, 7.2 GB each at 30010.

## 6. Tests

tests/testthat/test-se-reparam-reml.R: the 0.83.0 test "a whitened REML
fit above the size limit refuses its SEs" is replaced by four tests:
- "a whitened REML fit above the size limit reads sparse solves" (limit
  20 before the fit, 30 joint parameters, sd lost; 20 responses moved by
  2e-15 k; fitted() and vcov() SEs against the exact GLS values, 1e-6),
- "the sparse route agrees with the dense one on every reader" (one data
  set, sd not lost, fitted() with and without group effects, in-sample,
  frm_linpred(), vcov(), vcov(full = TRUE), frm_joint_cov()),
- "the item's whole-route fit above 4000 parameters is exact" (the item's
  seed 1 at the default limit, skip_on_cran; against the GLS value plus
  the variance parameters' propagation, or the GLS value where a rounding
  loses the sd),
- "a whitened REML fit whose sparse route fails keeps NaN SEs" (the kept
  0.83.0 path, jcs_joint_cov() mocked to NULL).

Before and after, one process each (`dev/sparsese-par.sh reml-<arm>`,
logs `dev/sparsese-suite-reml-<arm>/`):
```
base: frmtmb RESULT test-se-reparam-reml.R pass=115 fail=64 err=1 skip=0 warn=21
lane: frmtmb RESULT test-se-reparam-reml.R pass=181 fail=0 err=0 skip=0 warn=0
base, by test:
      1 ERROR the item's whole-route fit above 4000 parameters is exact
     60 FAIL a whitened REML fit above the size limit reads sparse solves
      1 FAIL a whitened REML fit whose sparse route fails keeps NaN SEs
      2 FAIL the item's whole-route fit above 4000 parameters is exact
      1 FAIL the sparse route agrees with the dense one on every reader
     20 WARN a whitened REML fit above the size limit reads sparse solves
      1 WARN a whitened REML fit whose sparse route fails keeps NaN SEs
DETAIL FAIL [the item's whole-route fit above 4000 parameters is exact] 
    Expected `max(abs(unname(fp[, "Est.Error"])/sqrt(vref) - 1))` < 1e-04. | 
    Actual comparison: 0.0448 >= 0.0001 | Difference: 0.044
DETAIL FAIL [the item's whole-route fit above 4000 parameters is exact] 
    Expected `is.null(attr(root, "rows"))` to be FALSE. | Differences: | 
    `actual`:   TRUE  | `expected`: FALSE | 
```

On rellib-r22 the item's seed 1 reads 4.48 percent off and every one of
the 20 small fits has NaN SEs with the size warning; the agreement test
fails there only on the `rows` attribute (rellib-r22 ignores the option
when it forms its factor, so both arms are dense there), and the mocked
test on the warning's new words. It is a guard of the new route, not a
pin of the defect.

## 7. Brief items 3 and 4, and what needs the user

- Item 3 (sdr_unwhiten() above the constant): done, not filed. Above the
  limit it keeps the factor's rows at `beta` and the outer parameters
  (joint_root_rows()). hypothesis(), fit_draw_space() and interop's
  coef_vcov_root() read those rows as the factor; vcov() and
  vcov(full = TRUE) read the sparse whole inverse that joint_cov_whole()
  returns beside them (jcs_joint_cov(whole = TRUE)), not the rows (fix
  round 1 corrects this sentence). All agree with the dense factor to
  7.9e-13 at 3510 (section 3). It now reads the option, not the constant, so a caller can lower
  it; the option is read when the fit forms its sdreport (its first
  standard error, usually inside frm()).
- Item 4 (where the NaN and warning stay): only where jcs_white_wanted()
  holds and jcs_joint_cov() returns NULL on a fit whose standard-error
  check repaired the covariance (a lost parameter): the inner block of
  the whitened precision does not factor, or `Jg` is not finite. A fit
  at a maximum has a positive definite inner block, so no case was
  constructed; the test mocks jcs_joint_cov(). The warning no longer
  advises raising the option, which reaches the same failed
  factorization. The whole route has no such refusal (filed LOW).
- Behavior beside the item: vcov() of such a fit now gives the
  "a product such as a' V a" warning, as below the limit, because the
  stored factor carries its "lossy" attribute; rellib-r22 gave no
  warning there (and wrong values).
- Needs the user (filed MEDIUM, updated in fix round 1): the dense
  route below the limit is 2 to 4 orders slower than the sparse one at
  3510 joint parameters, one order at 1510, and faster at 510 (review,
  section 6); the two agree to 5.5e-10. Recommendation: lower the
  default of `frmtmb.joint_root_max` to about 1000, where the two cost
  the same, now that F1 and F2 are fixed, or keep it and memoize the
  dense `V` in joint_cov_whole(). Not done here: it moves every
  whitened fit's route, and the brief scopes the change to above the
  limit.

### Suites, one file per process, two at a time (reference BLAS)

`FRMTMB_RELLIB=<lane lib> bash dev/sparsese-par.sh <tier> <jobs> 2`; every
log's `lib:` line names wt-sparsese-lib for frmtmb.
```
core (dev/sparsese-suite-core-lane/RESULTS.log):
  SUITE ran 349 of 349
  TOTAL pass=25762 fail=0 err=0 skip=243 warn=0
```
```
extensions (dev/sparsese-suite-ext-lane/RESULTS.log):
  SUITE ran 200 of 200
  TOTAL pass=9537 fail=0 err=0 skip=50 warn=0
  frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
  frmtmb.eam files=33 pass=1851 fail=0 err=0 skip=5 warn=0
  frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
  frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
  frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
  frmtmb.sample files=92 pass=4286 fail=0 err=0 skip=18 warn=0
  frmtmb.spatial files=8 pass=551 fail=0 err=0 skip=6 warn=0
  frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0
```

### OpenBLAS with LAPACK (dev/ciharden-openblas.sh <ver> lapack, dev/sparsese-ob.sh)

```
  0.3.26 Haswell Core: Haswell RESULT test-se-reparam-reml.R pass=181 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-se-reparam-reml.R pass=181 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-latres-sparse.R pass=67 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-latres-sparse.R pass=67 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-latres-newrow.R pass=39 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-latres-newrow.R pass=39 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-marginaleffects-re.R pass=59 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-marginaleffects-re.R pass=59 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-latres.R pass=305 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-latres.R pass=305 fail=0 err=0 skip=0 warn=0
  0.3.32 Haswell Core: Haswell RESULT test-se-reparam-reml.R pass=181 fail=0 err=0 skip=0 warn=0
  0.3.32 Sandybridge Core: Sandybridge RESULT test-se-reparam-reml.R pass=181 fail=0 err=0 skip=0 warn=0
  0.3.32 Haswell Core: Haswell RESULT test-latres-sparse.R pass=67 fail=0 err=0 skip=0 warn=0
  0.3.32 Sandybridge Core: Sandybridge RESULT test-latres-sparse.R pass=67 fail=0 err=0 skip=0 warn=0
```

Every log's lib line names wt-sparsese-lib. dev/ciharden-out/ deleted.

### One session per suite, three orders (dev/release/run-sessions.sh sparsese-sess 2)

```
  SESSION frmtmb forward files=336/349 pass=25762 fail=0 err=0 skip=243 warn=0
  SESSION frmtmb reverse files=336/349 pass=25927 fail=0 err=0 skip=243 warn=0
  SESSION frmtmb preload files=336/349 pass=26031 fail=0 err=0 skip=245 warn=0
  SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
  SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
  SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
  SESSION frmtmb.eam forward files=33/33 pass=1851 fail=0 err=0 skip=5 warn=0
  SESSION frmtmb.eam reverse files=33/33 pass=1851 fail=0 err=0 skip=5 warn=0
  SESSION frmtmb.eam preload files=33/33 pass=1851 fail=0 err=0 skip=5 warn=0
  SESSION frmtmb.latent forward files=13/13 pass=700 fail=0 err=0 skip=2 warn=0
  SESSION frmtmb.latent reverse files=13/13 pass=700 fail=0 err=0 skip=2 warn=0
  SESSION frmtmb.latent preload files=13/13 pass=700 fail=0 err=0 skip=2 warn=0
  SESSION frmtmb.learn forward files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
  SESSION frmtmb.learn reverse files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
  SESSION frmtmb.learn preload files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
  SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
  SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
  SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
  SESSION frmtmb.sample forward files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
  SESSION frmtmb.sample reverse files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
  SESSION frmtmb.sample preload files=87/92 pass=4286 fail=0 err=0 skip=18 warn=0
  SESSION frmtmb.survival forward files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
  SESSION frmtmb.survival reverse files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
  SESSION frmtmb.survival preload files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
  SESSION frmtmb.spatial forward files=8/8 pass=551 fail=0 err=0 skip=6 warn=0
  SESSION frmtmb.spatial reverse files=8/8 pass=551 fail=0 err=0 skip=6 warn=0
  SESSION frmtmb.spatial preload files=8/8 pass=551 fail=0 err=0 skip=6 warn=0
  SPARSESE-SESS ran 27 of 27
  DETAIL lines: 0; LEAK lines: 70
```

files=336/349 (core) and 87/92 (frmtmb.sample): the files with no
expectation outside the gated tier; the one-file run has 13 such core
files at pass=0 fail=0 err=0 skip=0.

## 8. R CMD check --as-cran, core

`bash dev/sparsese-check.sh` (built and checked in `dev/sparsese-check/`,
deleted after; logs `dev/sparsese-out/check-core.txt`,
`dev/sparsese-out/check-core-testthat.Rout`), R_LIBS = lane library,
rellib-r22, user library; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; pandoc
and TinyTeX on PATH; NOT_CRAN unset, so skip_on_cran() tests skip:

```
  Status: 1 NOTE   (the V8 math-rendering NOTE of the HTML manual)
  Running 'testthat.R' [29m] OK
  [ FAIL 0 | WARN 0 | SKIP 564 | PASS 20848 ]
```

No extension changed, so no extension was checked; their suites ran
against this core (section 6).

## 9. Files touched, bump

- R/joint-sparse.R, R/predict.R (get_joint_cov(), joint_cov_repair()'s
  comment, joint_root_limit()'s doc, joint_cov_white_inner(),
  lp_basis_out()), R/se-reparam.R (sdr_unwhiten(), joint_whitener(),
  joint_root_rows(), joint_cov_whole()), R/se-check.R (se_pred_warn()'s
  "size" warning only), NEWS.md, tests/testthat/test-se-reparam-reml.R,
  dev/test-backlog.md (CLOSED line under 4952; section "Filed in the
  round of 2026-11-16", subsection "From lane sparsese").
- New: dev/sparsese-*.R, dev/sparsese-*.sh, dev/sparsese-*.jobs, this
  file.
- No roxygen output changed (every new function is @noRd; man/ and
  NAMESPACE are as before).
- Bump: minor (a BREAKING bullet: NaN and its warning become values,
  vcov() values and its warning change above the limit). No floor
  changes; Matrix and TMB as already declared.
- No ledger row, compat row or equivalence pair is touched.
- The working tree of R/*.R is CRLF, as checked out; NEWS.md, the test
  file and dev/test-backlog.md came out LF after a sed edit and were put
  back to CRLF, so a file compare shows only the edits.

## 10. Fix round 1 (review `dev/reviews/2026-11-17-sparsese.md`)

### F1: frm_lp_basis()'s V_root and frm_curve() on the sparse route

lb_dense() (the exported frm_lp_basis()) now gives a whitened sparse fit
`V_root = M[pos, ci] L`, with `V_gamma[ci, ci] = L L'` over the whitened
positions the rows load (jcs_root_at(); jcs_psd_factor() takes a Cholesky
factor in the scale of the block's standard deviations, a zero row for a
position with zero variance, and the eigen-directions above 1e-12 of the
largest where the rest is singular), and `V_rest` zero. No dense factor
of the whole joint is formed: the block is over the rows' positions only.

### F2: fit_fd_se() (cumulative and categorical fitted()) on the sparse route

fd_joint_root() now returns, for a whitened sparse fit, jcs_fd_root():
`V[pos, pos] = Vrest + K K'` with a column of `K` per direction of the
drawn parameters (their whitened block's factor `Ld`, carried by `M`, and
`V_bd Ld'^+` on the group effects) and `Vrest` the group effects'
covariance given the drawn parameters. Why this and not a factor over
every position: fit_fd_se() spends two evaluations of the output per
column of `K`, so a column per group effect costs two per level (the
dense route's root has one column per joint parameter, which is its 62 s
at 4508); with columns for the drawn parameters only, the group effects
keep fit_fd_se()'s chain and batch routes, as on the dense route, and
the blocks come from the whitened coordinates, where they keep their
digits. The Schur complement `Vrest` is the same in either coordinates
of the drawn parameters, so it is formed in the whitened ones.

### Numbers

`Rscript dev/sparsese-fix1.R <lib> <seed> [group multiplier]`: the
review's base design (gaussian, `y ~ f + x1..x7 + (1 | g)`, 41 joint
parameters) and a cumulative fit (`yc ~ x1..x5 + (1 | g)`, 38), each
fitted under limit 1e6 and 20 in one process. Before is the round's first
lane build (`fix1-before-s1.txt`), after the fix (`fix1-after-*.txt`):

```
fix1-before-s1:
  V_root recipe: max |recipe / se.fit - 1| 0 
  frm_curve ERROR: frm_curve(): the assembled covariance of this grid disagr
  V_root recipe: max |recipe / se.fit - 1| NA 
  new rows Est.Error: max |sparse / dense - 1| 1.035e+02; max sparse 4.1280 
  in-sample Est.Error: max |sparse / dense - 1| 1.238e+02; max sparse 7.3926
fix1-after-s1:
  V_root recipe: max |recipe / se.fit - 1| 0 
  V_root recipe: max |recipe / se.fit - 1| 5.15e-11 
  frm_curve .se sparse / dense - 1: 6.05e-11 
  new rows Est.Error: max |sparse / dense - 1| 9.275e-09; max sparse 0.0984 
  in-sample Est.Error: max |sparse / dense - 1| 1.568e-08; max sparse 0.1949
fix1-after-s2:
  V_root recipe: max |recipe / se.fit - 1| 0 
  V_root recipe: max |recipe / se.fit - 1| 7.35e-11 
  frm_curve .se sparse / dense - 1: 7.34e-11 
  sparse : joint 38  whitened TRUE  lost [ nonfinite ]
  new rows Est.Error: max |sparse / dense - 1| 7.542e-08; max sparse 5.2859 
  in-sample Est.Error: max |sparse / dense - 1| 1.857e-07; max sparse 5.5931
fix1-after-s3:
  V_root recipe: max |recipe / se.fit - 1| 0 
  V_root recipe: max |recipe / se.fit - 1| 3.3e-11 
  frm_curve .se sparse / dense - 1: 6.88e-11 
  new rows Est.Error: max |sparse / dense - 1| 3.176e-08; max sparse 0.1409 
  in-sample Est.Error: max |sparse / dense - 1| 3.423e-08; max sparse 0.2032
fix1-after-s1-g0:
  V_root recipe: max |recipe / se.fit - 1| 0 
  V_root recipe: max |recipe / se.fit - 1| 3.1e-10 
  frm_curve .se sparse / dense - 1: 3.11e-10 
  new rows Est.Error: max |sparse / dense - 1| 9.850e-09; max sparse 0.1037 
  in-sample Est.Error: max |sparse / dense - 1| 1.117e-08; max sparse 0.1847
fix1-after-s2-g0:
  V_root recipe: max |recipe / se.fit - 1| 0 
  V_root recipe: max |recipe / se.fit - 1| 6.6e-10 
  frm_curve .se sparse / dense - 1: 4.99e-10 
  sparse : joint 38  whitened TRUE  lost [ nonfinite ]
  new rows Est.Error: max |sparse / dense - 1| 5.765e-08; max sparse 6.9911 
  in-sample Est.Error: max |sparse / dense - 1| 3.953e-07; max sparse 7.1884
```

In each block the first "V_root recipe" line is the dense arm, the second
the sparse arm. Reading: frm_curve() refused on the sparse route before
(0.00132 relative) and now agrees with the dense route to 3.1e-10 (6.6e-10
at a lost sd); the V_root recipe matches se.fit to 6.6e-10. Cumulative
fitted() was 104 to 124 times off before and now agrees to 4e-7 or
better. Seeds 2 (group multiplier 1 and 0) lose theta_1 as "nonfinite";
both routes then give Est.Error up to 5.6 to 7.2, the same pre-existing
defect on rellib-r22 (`Rscript dev/sparsese-cumlost.R <base|lane> 2 4000`,
`dev/sparsese-out/cumlost-*-s2.txt`: the fit warns, summary() has NaN and
NA, fitted() gives 289 of 900 entries above 0.5 with no warning of its
own); filed MEDIUM.

At 4508 joint parameters (the review's
`dev/reviews/sparsese-review-cum4500.R lane default`,
`dev/sparsese-out/cum4500-lane-default-fix1.txt`): fitted() 0.13 s, its
Est.Error equal to the review's dense run (limit 1e6, 62.17 s,
`dev/reviews/sparsese-review-out/cum4500-lane-1e6.txt`) in every printed
digit (0.050191 0.013987 0.036861 ... 0.052222 0.013501 0.039437), where
the first lane build matched rellib-r22 and was up to 1.2e-3 off.

### Tests

New in tests/testthat/test-se-reparam-reml.R: "frm_curve() and
frm_lp_basis()'s V_root read the sparse route" and "a cumulative fit's
fitted() reads the sparse route" (data helper sparse_route_data()). Seen
failing: with the two fixes switched off in place (fd_joint_root()'s
`white` set to FALSE and lb_dense()'s V_root branch to `if (FALSE)`),
built into the lane library, then restored and rebuilt
(`dev/sparsese-suite-reml-fix1-disabled/`): pass=185 fail=3 err=1, the
frm_curve() refusal (0.00132) and the cumulative ratios 103.5 and 123.8
and a largest Est.Error of 7.4. With the fixes
(`dev/sparsese-suite-reml-fix1/`): pass=193 fail=0 err=0 skip=0 warn=0.

### F5 wording

Corrected: the item-seed test's comment, the sdr_unwhiten() doc, section
2 and section 7 (vcov() reads the sparse whole inverse from
joint_cov_whole(), not the root rows), and section 3 (the agree-100 logs
are in dev/sparsese-out/). NEWS names what reads the sparse route,
including frm_lp_basis(), frm_curve() and cumulative fitted(), now that
F1 and F2 are fixed.

### Backlog

My section gains a MEDIUM (cumulative fitted() above 1 with a
"nonfinite" sd, pre-existing) and three LOW (the review's P1, P2, P3),
and the MEDIUM on the default limit now carries the review's timings
(the dense route is faster at 510, one order slower at 1510), its point
that a lower default was unsafe until F1 and F2 were fixed, and the
alternative of memoizing the dense `V`. The default is not changed.

### Fix round 1 runs

- tests/testthat/test-se-reparam-reml.R, reference BLAS: 193/0/0/0/0
  (`dev/sparsese-suite-reml-fix1/`).
- OpenBLAS with LAPACK (`dev/sparsese-ob.sh`, every lib line the lane's;
  `dev/ciharden-out/` deleted after):
```
  0.3.26 Haswell Core: Haswell RESULT test-se-reparam-reml.R pass=193 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-se-reparam-reml.R pass=193 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-curve-illcond.R pass=4 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-curve-illcond.R pass=4 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-ordinal-fitted.R pass=129 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-ordinal-fitted.R pass=129 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-predict-lp-basis.R pass=30 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-predict-lp-basis.R pass=30 fail=0 err=0 skip=0 warn=0
  0.3.26 Haswell Core: Haswell RESULT test-latres-sparse.R pass=67 fail=0 err=0 skip=0 warn=0
  0.3.26 Sandybridge Core: Sandybridge RESULT test-latres-sparse.R pass=67 fail=0 err=0 skip=0 warn=0
  0.3.32 Haswell Core: Haswell RESULT test-se-reparam-reml.R pass=193 fail=0 err=0 skip=0 warn=0
  0.3.32 Sandybridge Core: Sandybridge RESULT test-se-reparam-reml.R pass=193 fail=0 err=0 skip=0 warn=0
```
- Core and frmtmb.sample, one file per process, two at a time
  (`dev/sparsese-suite-fix1-lane/RESULTS.log`, every lib line the lane's):
```
  SUITE ran 441 of 441
  TOTAL pass=30060 fail=0 err=0 skip=261 warn=0
  frmtmb files=349 pass=25774 fail=0 err=0 skip=243 warn=0
  frmtmb.sample files=92 pass=4286 fail=0 err=0 skip=18 warn=0
```
- One session, three orders, core (`dev/release/run-sessions.sh sparsese-fix1-sess 2 frmtmb`):
```
  SESSION frmtmb forward files=336/349 pass=25774 fail=0 err=0 skip=243 warn=0
  SESSION frmtmb reverse files=336/349 pass=25939 fail=0 err=0 skip=243 warn=0
  SESSION frmtmb preload files=336/349 pass=26043 fail=0 err=0 skip=245 warn=0
  SPARSESE-FIX1-SESS ran 3 of 3
  DETAIL lines: 0; LEAK lines: 20, untagged: 0
```
  Every LEAK line is tagged "[with load of frmtmb.<ext>]" (test-arg-refusal.R,
  test-cox-home.R, test-families-extreme.R loading an extension).
- R CMD check --as-cran, core (`bash dev/sparsese-check.sh`, deleted after;
  logs `dev/sparsese-out/check-core-fix1.txt`,
  `dev/sparsese-out/check-core-fix1-testthat.Rout`):
```
  Status: 1 NOTE   (the V8 math-rendering NOTE of the HTML manual)
  Running 'testthat.R' [27m] OK
  [ FAIL 0 | WARN 0 | SKIP 564 | PASS 20860 ]
```

Files touched in fix round 1: R/joint-sparse.R (jcs_psd_factor(),
jcs_root_at(), jcs_fd_root()), R/brms-shapes.R (fd_joint_root(), one
hunk; new to this lane), R/predict.R (lb_dense()), R/se-reparam.R
(sdr_unwhiten()'s doc), NEWS.md, tests/testthat/test-se-reparam-reml.R,
dev/test-backlog.md (my section), this file; new dev/sparsese-fix1.R,
dev/sparsese-cumlost.R, dev/sparsese-fix1.jobs,
dev/sparsese-fix1.suite.jobs. No man/ or NAMESPACE change. Bump and
floors unchanged (minor; none).

## 11. Short check of the fix round (`dev/reviews/sparsese-short.md`): SHIP

Record edits only; no code change and no reruns.
- The short check's correction to the review: the review's independent
  reference (`dev/reviews/sparsese-review-indep.R`, line 100) held
  `names(of) %in% names(lost)`, but se_lost names the parameter
  `theta_1` where opt$par names it `theta`, so the reference never held
  the lost sd (it kept a theta variance near 1.16e8). The b-row
  differences the review gave at 4500 groups, seeds 2 and 3 (7.45e-7,
  3.38e-7 fitted with group effects; 2.87e-6, 1.93e-6 in-sample), came
  from that mismatch and not from "the choice of the held direction".
  With the index stripped those rows agree to 8.0e-9 and 8.7e-9, and the
  b variances agree exactly. A held-direction choice is real only for a
  lost correlation (the short check's slope design, seed 2).
- Backlog (my section): the cumulative lost-sd item is raised from
  MEDIUM to HIGH with the short check's reference (0.1492 in-sample and
  0.0709 for new rows, so the package is 141 and 113 times too high, a
  silent wrong value on rellib-r22 and both routes); new LOW items for
  jcs_fd_root()'s `qr.solve()` tolerance (stops at sd ratios of 1e8 or
  more; `tol = 0` exact to 1.9e-16; not reached on a real fit) and for
  vcov(full = TRUE) of a REML fit printing finite values for a lost
  parameter (1.157e8 for a boundary sd, -5.85e5 for a boundary
  correlation); the P2 LOW extended to cumulative fits (7.8e-4 and 5.0e-4
  scaled, effect on Est.Error up to 1.4e-4 relative, both routes).
