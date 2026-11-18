# Lane drawdir: a per-direction draw check (round of 2026-11-16)

Base b255350e (frmtmb 0.83.0), branch wt-r84-drawdir, library
`C:/Users/adf44/source/r/wt-drawdir-lib`, base build rellib-r22. Brief
`round-briefs/r084-drawdir.md`: backlog items 5584 (MEDIUM, the ratio
rule misses a weak direction diluted by well-identified parameters) and
5644 (LOW, the check reads the median). Logs named here are local
(gitignored); every number this record rests on is pasted here.

STATUS: fix round 1 done (section 12), with option (b). Before it: built,
measured, all suites (one file per process, gated
tier, one session three ways) and `R CMD check --as-cran` of core.

## 1. Reproduction on rellib-r22

Script `dev/drawdir-repro.R` (a copy of `dev/reviews/predraw2-short-a.R`
with its library from `DRAWDIR_LIB`), log `dev/drawdir-repro-base.txt`.
Reproduced exactly; the review's numbers hold on 0.83.0.

Part S (each of the 18 mixture fits that fire alone, its own e and m at
the check's test draws, plus a separable block of K coordinates). The
first K at which the check stops firing, quadratic block (awk over the
log):

    mixture_gaussian seed 20 first unfired: K8
    mixture_gaussian seed 23 first unfired: none to K256
    mixture_gaussian seed 25 first unfired: none to K256
    mixture_gaussian seed 29 first unfired: K16
    mixture_gaussian seed 36 first unfired: K32
    mixture_gaussian seed 40 first unfired: K64
    mixture_poisson seed 11 first unfired: K8
    mixture_poisson seed 13 first unfired: K128
    mixture_poisson seed 18 first unfired: K64
    mixture_poisson seed 20 first unfired: K16
    mixture_poisson seed 24 first unfired: K32
    mixture_poisson seed 25 first unfired: K32
    mixture_poisson seed 26 first unfired: K32
    mixture_poisson seed 31 first unfired: K16
    mixture_poisson seed 33 first unfired: K32
    mixture_poisson seed 36 first unfired: K128
    mixture_poisson seed 39 first unfired: K16
    mixture_poisson seed 5 first unfired: K64

K 8 on 2, 16 on 4, 32 on 5, 64 on 3, 128 on 2, none on gaussian seeds 23
and 25: the item's table. With the 200-level poisson block (alone:
excess 2.23, ratio 0.0237), 16 of 18 do not fire.

Part R (the same data with K = 5, 10, 20 noise covariates in both
components; 72 fits with K = 0): 29 fire (18 of them K = 0), 36 unfired
with a summed excess of 1 or less, and the 7 the ratio decides:

    mixture_poisson  seed 13 K 20 p 45 ---- ex   1.071 rel  0.049 n 20 | mae lane 0.4665 z  29.04 | forced mae 0.4062 z  12.55 said frmtmb_draw_pareto
    mixture_poisson  seed 20 K 10 p 25 ---- ex   1.101 rel  0.090 n 40 | mae lane 0.2621 z  18.31 | forced mae 0.2829 z   8.74 said frmtmb_draw_pareto
    mixture_poisson  seed 25 K 20 p 45 ---- ex   1.147 rel  0.052 n 20 | mae lane 0.5662 z  26.50 | forced mae 0.5352 z  -7.15 said frmtmb_draw_pareto
    mixture_poisson  seed 26 K 20 p 45 ---- ex   1.624 rel  0.074 n 20 | mae lane 0.5508 z  27.80 | forced mae 0.5957 z  30.33 said frmtmb_draw_pareto
    mixture_poisson  seed 36 K 10 p 25 ---- ex   1.353 rel  0.111 n 40 | mae lane 0.2201 z  19.71 | forced mae 0.2667 z   5.30 said frmtmb_draw_pareto
    mixture_poisson  seed 36 K 20 p 45 ---- ex   1.555 rel  0.070 n 20 | mae lane 0.4375 z  29.95 | forced mae 0.4269 z   3.71 said frmtmb_draw_pareto
    mixture_poisson  seed 39 K 20 p 45 ---- ex   1.136 rel  0.051 n 20 | mae lane 0.4686 z  32.56 | forced mae 0.4998 z  15.40 said frmtmb_draw_pareto

(`z`: predict()'s mean over 1000 draws against fitted()'s, in Monte
Carlo standard errors; `mae`: rows' predictive means against the
truth; `forced`: the check forced on.) Note that the unfired poisson
fits with an excess of 1 or less read z of 6.7 to 27.9 too (for
example seed 11 K 20, z 27.85): a z far from 0 is not by itself the
mark of a wrong law here, since exp() of a draw with a weakly identified
second component has a mean above fitted() under any law with that
spread.

## 2. Measured before building: which directions, which distance

Probe scripts read rellib-r22 (0.83.0) through its own draw law and
`draw_logratio()` (the summed check's per-draw excess `e = drop - m / 2
- log Jacobian`), at +-c along a set of directions of the deviates:

- `dev/drawdir-probe.R` (logs `dev/drawdir-probe-*.txt`,
  `dev/drawdir-probe1-*.txt`; c = 1.5, 2, 3 and then 1, 1.5, 2): three
  sets, `eig` (eigen-directions of the law's covariance), `cor`
  (eigen-directions of its correlation matrix) and `par` (one per
  parameter: the parameter moved c of its standard deviations and the
  others by their regression on it, column j of the law's factor).
- `dev/drawdir-probe2.R` (logs `dev/drawdir-probe2-*.txt`): per
  parameter, `eA` along its `par` direction and `eB` along its own
  slice (the others at the estimates, c of its conditional standard
  deviations, no side law: what `draw_edge_law()` reads).
- Fits: `dev/drawdir-cases.R`. mix40 (the three mixture designs of
  predraw2 section 3, seeds 1 to 40), partR (the short check's part R,
  72 fits), kd (the K designs, 162 fits with draws), gglm (128 gaussian
  GLMs: n 20, 50, 200; 1, 5, 15 covariates; and `sigma ~ x`), gglmm (144
  gaussian GLMMs: `(1 | g)` with 5 to 30 groups of 4 or 10, sd 0.3 or 1;
  `(1 + x | g)` with 8 to 40 groups of 8, slope sd 0.2 or 0.5), and
  predraw2's corpus (wb 310, zip 30, sweep 30, target 6).

### 2a. The plain per-direction excess fires on gaussian GLMMs

The largest excess over directions, as `e / (c^2 / 2)` (the drop over
the law's drop, less 1), unfired fits, `par` set
(`dev/drawdir-sum.R`):

    par c 2.0                          fits | unfired max
    gri  (1 | g)                         96 |     7.61
    grs  (1 + x | g)                     48 |     10.5
    G6, R0.1, R0.3                  40/50/50 | 3.36 / 2.68 / 3.6
    S0.1, S0.2                        25/25 |      Inf / 8.58
    bglmm, zip                        30/30 | 4.84 / 7.8
    gglm (108), K designs (all)             | 0.249, 1.19
    mix40 fired (44)                        | min of max 6.67

A group log sd near 0 is steep on its upper side (`+theta_1`, excess
up to 21 at c = 2 on `grs`, 15.2 on `gri` seed 1) and flat on its lower
side. The side law reads that slice, and leaves it to the normal law
where the slice falls by 2 later than 1 conditional standard deviation
out. The plain reading puts these GLMMs inside the mixtures' range at
every c: at c = 1 the GLMM maximum was an excess of 1.72 and the two
weakest of the 18 fired mixtures read 1.29 (poisson seed 11) and 1.41
(seed 20). The `eig` and `cor` sets overlap the same way (`cor` c = 2:
GLMMs up to 5.0, fired mixtures from 4.14).

### 2b. The direction beyond its own slice separates them

S2 = max over parameters and sides of `eA - max(eB, 0)` (the
direction's excess beyond what the parameter's own slice shows), S1 =
max of `eA` (`dev/drawdir-an2.R` over all probe2 files; excess units):

    == c 1.5: max over unfired fits | min over fired fits
    gri                                  n  96 | S1     6.62 S2     2.07
    grs                                  n  48 | S1     7.31 S2     2.79
    mix40:mixture_cumulative:interior    n  29 | S1     29.3 S2     29.1 | fired 26 S1     4.25 S2     3.88
    mix40:mixture_gaussian:boundary      n  39 | S1       17 S2       17 | fired  6 S1     16.3 S2     16.3
    mix40:mixture_poisson:boundary       n  30 | S1     20.9 S2     20.8 | fired 12 S1     5.41 S2     5.41
    S0.1                                 n  25 | S1     35.5 S2     31.2
    S0.2                                 n  25 | S1     8.65 S2     4.01
    == c 2.0: max over unfired fits | min over fired fits
    gri                                  n  96 | S1     15.2 S2     4.38
    grs                                  n  48 | S1     21.1 S2     5.61
    mix40:mixture_cumulative:interior    n  29 | S1      116 S2      115 | fired 26 S1       20 S2     19.4
    mix40:mixture_gaussian:boundary      n  39 | S1     42.7 S2     42.6 | fired  6 S1     50.3 S2     50.2
    mix40:mixture_poisson:boundary       n  30 | S1     63.1 S2     62.4 | fired 12 S1     13.3 S2     13.3
    beta                                 n   5 | S1    0.146 S2  0.00767
    bglmm                                n  30 | S1     9.68 S2     0.34
    cumul                                n  10 | S1    0.365 S2    0.151
    G6                                   n  40 | S1     6.73 S2    0.333
    gamma                                n   5 | S1    0.143 S2  0.00399
    glm_gaussian                         n  10 | S1    0.203 S2 1.71e-13
    glm_logit                            n  10 | S1    0.445 S2    0.101
    glm_poisson                          n  10 | S1    0.147 S2   0.0488
    negbin                               n   5 | S1    0.219 S2 0.000134
    pglmm                                n  30 | S1    0.235 S2    0.216
    R0.1                                 n  50 | S1     5.35 S2    0.273
    R0.3                                 n  50 | S1     7.19 S2    0.538
    S0.2                                 n  25 | S1     17.2 S2     6.99
    student                              n   5 | S1    0.202 S2    0.202
    zip                                  n  30 | S1     15.6 S2     3.08
    gglms:30 / gglms:100                 n  10 | S1     1.17 S2    0.952 / 0.213 0.147
    gglm (9 designs, 108 fits)                 | S1 <= 0.498, S2 <= 4e-12
    kd (19 designs, 162 fits)                  | S1 <= 2.39, S2 = 0
    sweep (non-mixture, 26 fits)               | S2 <= 1.01

(The last four lines condense several printed rows each. S0.1 at
c = 2 has a direction the objective cannot evaluate, so its row prints
NaN there; per fit below.) For the mixtures S2 equals S1: their excess
is all joint, no single parameter's slice shows it, which is why the
side law does not narrow them. For a GLMM's group sd, most of the
excess is the slice's.

The largest S2 among 522 GLMs and GLMMs: 462 gaussian (gglm 128,
gglmm 144, wb's G6 40, R0.1 and R0.3 100, S0.1 and S0.2 50) and wb's 30
poisson and 30 bernoulli GLMMs (`dev/drawdir-top.R` with the regex
`^(gri|grs|S0|R0|G6|bglmm|pglmm|gglm)`; wb's 10 glm_gaussian fits are
in the table above, S2 at most 1.7e-13):

    fits 522 | S2 > 6: 4 | > 10.5: 3 | NA: 0
                 id  fire np     S2     par side     eA sratio
            S0.1:14 FALSE  6 46.080 theta_3    - 46.140  1.470
            S0.1:16 FALSE  6 29.030 theta_3    + 29.810  1.739
            S0.1:24 FALSE  6 12.680 theta_2    + 15.570  1.746
             S0.2:8 FALSE  6  6.986 theta_2    + 17.160  1.206
       grs:40:0.2:4 FALSE  6  5.609 theta_3    -  7.417  1.182
             S0.2:6 FALSE  6  4.669 theta_2    + 14.160  1.135
             S0.1:5 FALSE  6  4.390 theta_3    -  5.364  1.256
     gri:30:4:0.3:1 FALSE  4  4.377 theta_1    + 15.210  1.114

(`sratio`: the parameter's standard deviation over its conditional
one.) The four above 6 are random-slope fits with slope sd 0.1 or 0.2
whose slope sd and correlation are weakly identified together. They
are not false alarms (2c).

### 2c. On those random-slope fits the normal law is wrong

`dev/drawdir-force.R` (log `dev/drawdir-force-S.txt`): predict() at a
new level (`x = 2`), 4000 draws, seed 1, on rellib-r22 as it comes and
with the check forced on; the exact probability that the true
predictive law puts inside the interval, and its width (true 90% width
4.26 and 98% width 6.03 at slope sd 0.1; 4.41 and 6.24 at 0.2):

    S0.1:14    as is  cov90 0.953 w90 5.238 cov98 1.000 w98 10.724 | forced cov90 0.916 w90 4.670 cov98 0.990 w98 6.897 [k=0.778]
    S0.1:16    as is  cov90 0.967 w90 5.548 cov98 1.000 w98 12.276 | forced cov90 0.927 w90 4.653 cov98 0.989 w98 6.653 [k=0.941]
    S0.1:24    as is  cov90 0.930 w90 4.960 cov98 0.997 w98 8.201  | forced cov90 0.884 w90 4.190 cov98 0.977 w98 6.030 [k=0.738]
    S0.1:5     as is  cov90 0.937 w90 5.201 cov98 0.994 w98 8.059  | forced cov90 0.914 w90 4.944 cov98 0.984 w98 6.949 [k=0.773]
    S0.2:8     as is  cov90 0.977 w90 6.110 cov98 1.000 w98 11.218 | forced cov90 0.945 w90 5.243 cov98 0.994 w98 7.501 [k=0.829]
    S0.2:6     as is  cov90 0.940 w90 5.055 cov98 0.997 w98 8.189  | forced cov90 0.917 w90 4.667 cov98 0.983 w98 6.465 [k=0.701]
    S0.1:1     as is  cov90 0.899 w90 4.247 cov98 0.981 w98 6.097  | forced cov90 0.894 w90 4.200 cov98 0.979 w98 6.016 [k=0.823]
    S0.1:2     as is  cov90 0.858 w90 3.991 cov98 0.971 w98 5.820  | forced cov90 0.877 w90 4.287 cov98 0.969 w98 5.977 [k=0.826]

(k: the forced arm's Pareto k, each with the warning
`frmtmb_draw_pareto`; the side-law messages are cut.) Where S2 is
large (seeds 14, 16, 24 and S0.2:8) the normal law's 98% interval was
1.36 to 2.04 times the true width, and the reweighted one 1.00 to 1.20
times. ("True" here and below is the predictive law with the
parameters KNOWN, which a posterior predictive should over-cover. The
review, `dev/reviews/2026-11-17-drawdir.md` section 2, set a brms-prior
posterior beside it on seeds 14 and 16: the normal law's 98% width is
1.43 and 1.67 times the posterior's, its 90% width 0.99 and 1.07
times; the reweighted 90% width is 0.89 times the posterior's on both.
The defect is in the far tail, and the reweighting overshoots it.) So
on these fits the reading finds a real defect of the law;
S0.1:14 is backlog 5644's fit. predraw2 called these designs nominal
from coverage averaged over the 25 seeds, which hides the few that
over-cover. Where S2 is small (S0.1:1, :2) forcing changes little.

### 2d. Choices

- Directions: one per parameter (`par`). With a separable block the
  law's factor is block-diagonal, so each direction stays inside its
  block and a block of well-identified parameters adds nothing to any
  other direction's reading (an identity, not a measurement). The
  `eig` and `cor` sets mix blocks whose eigenvalues are close, need an
  SVD, and separated the GLMMs from the mixtures no better (2a). The
  `par` set also names a parameter for the message.
- Distance: c = 2 (`draw_edge_k`, the side law's own reach). At c = 1
  the joint excess of the weakest fired mixture (0.98, poisson seed 11)
  is barely above the GLMM maximum (0.92); at c = 1.5 the gap to the GLMMs
  is 1.9 times (5.41 against 2.79), at c = 2 it is 2.4 times (13.3
  against 5.61); c = 3 leaves the objective's domain on S0.1 (Inf).
- Statistic: S2, the direction's excess less its own slice's positive
  excess. (Corrected in the fix round, after the review's section 1.
  This said "beyond its own slice", as if S2 were the part no single
  parameter shows. It is not: the slice is read at 2 CONDITIONAL sds,
  nearer than the direction's 2 marginal sds, so for a skew in the
  parameter alone that grows as the cube of the distance, S2 = eB
  (sratio^3 - 1): 0 for an uncorrelated parameter, most of the error
  for a strongly correlated one. The review's 2-D case: an error of 10.7
  at rho 0 does not fire, one of 6.7 at rho 0.9 does. S2 is a heuristic
  that keeps the GLMMs' group sds out.) The first build read a
  direction the objective cannot evaluate, with a finite slice, as a
  joint wall (S2 = Inf); the fix round reads it again nearer the
  estimate (section 12).
- Threshold: 6, borrowed from the side law's `draw_edge_c` (it narrows
  a side that falls by `k^2 / 2` within `draw_edge_c` of `k` standard
  deviations). For a quadratic misfit it is a law twice as wide as the
  likelihood along the direction (a drop of 8 at k = 2 where the law
  says 2); for any other shape a density error of exp(6) at 2 sds. A
  heuristic, not a calibrated false-alarm rate (the statistic is
  deterministic, and its value on GLMMs is not 0). It sits between the
  largest gaussian GLM or GLMM reading that is not a random-slope fit
  with a measured defect (5.61, `grs:40:0.2:4`) and the smallest of the
  18 fired mixtures (13.3): it was set on those fits.
- Cost: the slice is read only for a direction whose excess passes 6
  (S2 is at most eA), so most fits pay 2 evaluations per parameter.
- Flat directions: not read. A flat side is the side law's (held at
  the estimate); a direction whose drop is a small share of the law's
  is the lower side of every group log sd near 0 (`par` min -0.99 of
  the law's drop on grs at c = 2), which the side law leaves to the
  normal law on purpose.
- The reading runs only where the summed check does not fire, and is
  cached with the verdict.

## 3. What was built (`R/predict-brms.R`)

- `draw_check()`: where the summed verdict does not fire, calls
  `draw_check_dirs()` and fires where its `stat` is above
  `draw_dir_tol`. The cached verdict gains `dir`, `dir_par` and
  `dir_n` (the evaluations it took); `dir` is NA where the summed check
  fired first.
- `draw_check_dirs()`: the `par` directions at +-`draw_edge_k` through
  the side law, read by `draw_logratio()` (2 evaluations per
  parameter); for each direction and side whose excess passes the
  threshold, the parameter's own slice at `draw_edge_k` conditional
  standard deviations (`draw_dirs_csd()`, a QR of the factor's
  columns), one evaluation each, every inner solve started from the
  estimates' modes. `stat` is the largest joint excess; where a
  direction did not pass the threshold, its own excess stands in (an
  upper bound of its joint excess).
- `draw_dir_tol` = `draw_edge_k^2 / 2 * (1 / draw_edge_c^2 - 1)` = 6.
- `draw_check_said_what()`: the message and warning now say which
  reading fired; a per-direction one names the parameter, for example
  "along one direction of the normal law of vcov(), the one that moves
  `theta1_(Intercept)` and the parameters correlated with it, the
  log-likelihood falls 77.3 further than that law says at 2 standard
  deviations, beyond what that parameter alone shows".
- `predict_par_drawer()`: passes the parameter labels to the check (one
  line).
- When it fires, the reweighting and Pareto k path of 0.83.0 runs
  unchanged.
- `?predict.frmtmb_fit`: a third bullet for the reading, what it costs,
  that a random-slope model can fire it, and that a weakness spread
  over many directions can still go unflagged.

## 4. The mixture designs against the truth, before and after

`dev/drawdir-mixcov.R` (predraw2's `dev/predraw2-mixcov.R` on
rellib-r22, recording `dir`), seeds 1 to 40 per design, 1000 draws,
predict() seed 1, base (rellib-r22) and lane arms in separate
processes; summary `dev/drawdir-mc-sum.R`, log
`dev/drawdir-mc-sum.txt`:

    == mixture_gaussian: 40 seeds with both arms
      lane only   n  9 | cov90 base 0.9111 lane 0.9014 | cov98 base 0.9838 lane 0.9792 | w98 base 4.9518 lane 4.7330 | mae base 0.0712 lane 0.0616 plug 0.0549 | |z| base 4.9093 lane 3.7152
      neither     n 25 | cov90 base 0.8984 lane 0.8984 | cov98 base 0.9766 lane 0.9766 | w98 base 4.7230 lane 4.7230 | mae base 0.0739 lane 0.0739 plug 0.0745 | |z| base 3.6446 lane 3.6446
      fired both  n  6 | cov90 base 0.9055 lane 0.9055 | cov98 base 0.9800 lane 0.9800 | w98 base 4.7311 lane 4.7311 | mae base 0.0645 lane 0.0645 plug 0.0628 | |z| base 2.9182 lane 2.9182
      lane-only fits: seed dir k | mae base lane | cov90 base lane
         1    35.09 0.892 | 0.0991 0.0750 | 0.933 0.924
         8    13.03  1.12 | 0.0369 0.0253 | 0.921 0.901
        12     6.86 0.961 | 0.0877 0.0968 | 0.89 0.887
        13    42.60  0.89 | 0.0473 0.0474 | 0.903 0.894
        14     8.44  1.17 | 0.0510 0.0628 | 0.922 0.921
        16    30.15  1.01 | 0.1292 0.0999 | 0.91 0.891
        17    25.27 0.893 | 0.0418 0.0391 | 0.917 0.907
        28    18.81  1.29 | 0.0512 0.0421 | 0.901 0.893
        30    21.26 0.869 | 0.0969 0.0655 | 0.902 0.895
      neither: 25 of 25 with mae equal to 1e-9 base/lane

    == mixture_poisson: 40 seeds with both arms
      neither     n 17 | cov90 base 0.9586 lane 0.9586 | cov98 base 0.9945 lane 0.9945 | w98 base 8.3642 lane 8.3642 | mae base 21.2682 lane 21.2682 plug 0.1165 | |z| base 2.8450 lane 2.8450
      lane only   n 11 | cov90 base 0.9668 lane 0.9615 | cov98 base 0.9966 lane 0.9951 | w98 base 8.2679 lane 7.7320 | mae base 0.1665 lane 0.1395 plug 0.1349 | |z| base 5.5226 lane 5.5903
      fired both  n 12 | cov90 base 0.9609 lane 0.9609 | cov98 base 0.9944 lane 0.9944 | w98 base 7.6337 lane 7.6337 | mae base 0.1132 lane 0.1132 plug 0.1063 | |z| base 3.6629 lane 3.6629
      lane-only fits: seed dir k | mae base lane | cov90 base lane
         2    13.86  1.06 | 0.1703 0.1510 | 0.966 0.963
         3    45.85 0.866 | 0.1247 0.0928 | 0.963 0.956
         7    37.44 0.887 | 0.2536 0.1558 | 0.968 0.958
        10    36.44  1.39 | 0.1927 0.2424 | 0.963 0.957
        14    17.51 0.947 | 0.0731 0.0748 | 0.971 0.968
        15    14.26  0.96 | 0.1713 0.1481 | 0.969 0.966
        21    46.60 0.59. | 0.2588 0.1568 | 0.967 0.963
        23    19.39 0.717 | 0.1721 0.1339 | 0.965 0.962
        27    11.93 0.989 | 0.1030 0.1138 | 0.966 0.963
        29    62.36  1.38 | 0.1789 0.0818 | 0.972 0.963
        35    61.94 0.977 | 0.1331 0.1837 | 0.964 0.956
      neither: 17 of 17 with mae equal to 1e-9 base/lane

    == mixture_cumulative: 40 seeds with both arms
      neither     n 11 | ... mae base 0.0585 lane 0.0585 plug 0.0585
      fired both  n 26 | ... mae base 0.0321 lane 0.0321 plug 0.0318 | |z| base 6.0001 lane 6.0001
      lane only   n  3 | ... mae base 0.0304 lane 0.0309 plug 0.0322 | |z| base 8.8542 lane 5.4340
      lane-only fits: seed dir k | mae base lane
        12    18.60   1.1 | 0.0330 0.0284
        22    18.53 0.719 | 0.0290 0.0340
        30   115.02 0.978 | 0.0294 0.0303
      neither: 11 of 11 with mae equal to 1e-9 base/lane

(`k` "0.59." is the message's own full stop, k = 0.59.) The reading
adds 9, 11 and 3 fired fits to the summed check's 6, 12 and 26. On
them the predictive means' error against the truth falls on the
gaussian (0.0712 to 0.0616, 6 of 9 better) and poisson designs (0.1665
to 0.1395, 7 of 11 better) and is unchanged on the ordinal one (0.0304
to 0.0309). The 90% coverage falls on all 20 interval fits: the
gaussian mean from 0.911 to 0.901 (nominal 0.90; seeds 13, 28 and 30
move from just above 0.90 to 0.893 to 0.895, seed 12 from 0.890 to
0.887), the poisson
mean from 0.967 to 0.962 (the plug-in itself covers about 0.96, from the
counts' discreteness). The 98% width falls. 22 of the 23 warn k > 0.7 (poisson seed 21 gets the
message, k 0.59). As in predraw2, reweighting a weakly identified
mixture mostly ends in the warning, and the draws are still better
than the normal law's on average. Every fit neither reading fires is
equal to the base build (53 of 53). The heavy-tailed poisson seeds 1
and 28 (backlog 5617, the "neither" row's mean error of 21.3) do not
fire: their joint excess is 2.70 and 3.32 at c = 2 (probe1, par2), and
the near-empty component's slope is flat, not steep. The reading does
not see 5617, as the item foresaw.

## 5. False alarms: the built check on the null designs

`dev/drawdir-verdict.R` on the lane library (predict(), 200 draws,
seed 5; logs `dev/drawdir-v-<group>.txt`), summarised without R by
`dev/drawdir-v-sum.sh` ("max dir unfired" is the cached `dir`, which
is the direction's own excess, an upper bound of the joint excess,
where no slice was read; the exact joint excesses are section 2b's):

    gglm                                 fits 108 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.498
    gglms                                fits  20 noverdict  0 | summed  0 | dir  0 | max dir unfired 1.17
    gri                                  fits  96 noverdict  0 | summed  0 | dir  0 | max dir unfired 5.2
    grs                                  fits  48 noverdict  0 | summed  0 | dir  0 | max dir unfired 5.61
    glm_gaussian                         fits  10 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.203
    glm_poisson                          fits  10 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.147
    glm_logit                            fits  10 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.445
    cumul                                fits  10 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.365
    negbin                               fits   5 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.219
    gamma                                fits   5 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.143
    beta                                 fits   5 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.146
    student                              fits   5 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.202
    G6                                   fits  40 noverdict  0 | summed  0 | dir  0 | max dir unfired 4.01
    R0.1                                 fits  50 noverdict  0 | summed  0 | dir  0 | max dir unfired 5.35
    R0.3                                 fits  50 noverdict  0 | summed  0 | dir  0 | max dir unfired 5.1
    S0.1                                 fits  25 noverdict  0 | summed  0 | dir  3 | max dir unfired 5.85 S0.1:14(46.1) S0.1:16(29) S0.1:24(12.7)
    S0.2                                 fits  25 noverdict  0 | summed  0 | dir  1 | max dir unfired 5.51 S0.2:8(6.99)
    pglmm                                fits  30 noverdict  0 | summed  0 | dir  0 | max dir unfired 0.235
    bglmm                                fits  30 noverdict  0 | summed  0 | dir  0 | max dir unfired 4.74
    kd (19 designs, 190 fits, 162 with a verdict)        | summed 0 | dir 0 | max dir unfired 2.39
    zip                                  fits  30 noverdict  0 | summed  0 | dir  0 | max dir unfired 3.99
    sweep (30 fits, 26 with a verdict)                   | summed 3 | dir 0 | max dir unfired 1.03
      (summed: the gaussian and poisson mixtures at the boundary, the ordinal one inside)
    target (zib3269, A1, A4)                             | summed 0 | dir 0 | max dir unfired 0.0348

(The kd and sweep rows condense the printed per-design rows; fits with
no verdict fall back to the plug-in, a 0.82.0 behavior.) Over 472
gaussian GLMs and GLMMs the reading fires on 4 (0.85%), the
random-slope fits of section 2c where the normal law's 98% interval
for a new level was 1.36 to 2.04 times the true width. It fires on none
of the 128 gaussian GLMs, the 144 gaussian GLMMs of `gglmm`, the 140
random-intercept fits of wb, the 60 poisson and bernoulli GLMMs, the
162 K-design fits (FIX 1's concern) or the zero-inflated and hurdle
fits. The verdicts agree with the probe's S2 (section 2b) fit for fit.

These 472 are the fits the threshold was set on (section 2d), so the
count is in-sample, not a false-alarm rate. The review's survey of new
designs (`dev/reviews/drawdir-review-survey.R`, logs
`dev/reviews/drawdir-review-sv-main.*` and `-sv-b.*`, 286 fits; lane
build before the fix round) is the out-of-sample one:

    bfew4    bernoulli (1|g), 4 groups of 20   fits 15 dir 0 | max dir unfired 3.583
    bfew6    bernoulli (1|g), 6 groups of 20   fits 15 dir 0 | max 0.817
    binfew   binomial(8) (1|g), 5 groups of 6  fits 10 dir 0 | max 0.233
    xg       gaussian (1|a)+(1|b), 6 x 8 x 2   fits 15 dir 0 | max 5.166
    xp       poisson (1|a)+(1|b), 5 x 12       fits 10 dir 0 | max 5.018
    xg3      gaussian crossed 3 x 4 x 4        fits 10 dir 0 | max 2.021
    ordg     cumulative (1|g), 8 groups        fits 15 dir 0 | max 3.039
    ordf     cumulative, 6 categories, n 60    fits 10 dir 0 | max 2.587
    nbglm    negbinomial, n 40, shape 1        fits 10 dir 0 | max 0.907
    nbglmm   negbinomial (1|g), 8 groups       fits 10 dir 0 | max 2.880
    nbpois   negbinomial on poisson data       fits 10 dir 0 | max 3.396
    dsig     gaussian (1|g), sigma ~ x         fits 10 dir 0 | max 5.935
    dphi     Beta, phi ~ x                     fits  5 dir 0 | max 0.249
    prs      gaussian (1+x|g) with priors      fits 15 dir 0 | max 2.540
    prb      bernoulli, prior normal(0,1.5) b  fits 10 dir 0 | max 0.888
    sigid    gaussian, identity-link sigma ~ x fits  5 dir 1 | max 3.276 sigid:4(9.70)
    prsp     poisson (1+x|g), slope sd 0.2     fits 15 dir 2 | max 4.147 prsp:1(12.1) prsp:10(57.5)
    brs      bernoulli (1+x|g), slope sd 0.3   fits 15 dir 0 | max 4.371
    ordrs    cumulative (1+x|g)                fits 10 dir 1 | max 3.301 ordrs:9(6.17)
    sleep    lme4::sleepstudy (Days|Subject)   fits  1 dir 0 | max 0.323
    slp      sleepstudy-like, 18 subjects      fits 30 dir 0 | max 5.783
    slpw     the same, slope sd 2 (not 5.9)    fits 20 dir 4 | max 5.741 slpw:12(Inf) :14(11.3) :15(27.1) :18(16.9)
    slp10    sleepstudy-like, 10 subjects      fits 20 dir 1 | max 3.306 slp10:10(11.6)

(copied from the review.) Without a random slope: 1 of 160 fires
(sigid:4, a correct fire: the normal law put an identity-link sigma at
or below 0 in 31 of 4000 draws, and predict()'s means were NaN). With
one: 8 of 126 (prsp 2 of 15, ordrs 1 of 10, slpw 4 of 20, slp10 1 of
20; sleepstudy itself and 30 sleepstudy-like fits 0). ordrs:9 is barely
over, at 6.17, along a correlation already at 1.0 to double precision.
slpw:12's Inf was the non-finite end of section 12.

## 6. Backlog 5644: a two-sided reading of the summed check

The probe recorded the quantiles of the excess over the summed check's
own 40 test draws. A two-sided reading (fire where the upper quartile
is above 1, or also the lower one below -0.5), counted with awk over
`dev/drawdir-probe-gglm.tsv` and `dev/drawdir-probe1-*.tsv`:

    gglm                               fits 108 summed   0 | q75>1  24 | q75>1 & q25<-0.5  24
    grs                                fits  48 summed   0 | q75>1   2 | q75>1 & q25<-0.5   1
    S0.1                               fits  25 summed   0 | q75>1   2 | q75>1 & q25<-0.5   0
    S0.2                               fits  25 summed   0 | q75>1   3 | q75>1 & q25<-0.5   0
    kd:gaus_K50_m3                     fits  10 summed   0 | q75>1  10 | q75>1 & q25<-0.5  10
    kd:gaus_K200_m3                    fits  10 summed   0 | q75>1  10 | q75>1 & q25<-0.5  10
    kd:gaus_K500_m3                    fits  10 summed   0 | q75>1  10 | q75>1 & q25<-0.5  10
    kd:gamma_K50_m3                    fits  10 summed   0 | q75>1  10 | q75>1 & q25<-0.5  10
    kd:pois_K100_m10_mu1.5             fits  10 summed   0 | q75>1  10 | q75>1 & q25<-0.5   9
    (every other K design: q75 > 1 on all its fits but the bernoulli and
    gaus_K10 ones; gglms, gri, G6, R, pglmm, bglmm, zip, glm_*: 0)
    S0.1:14 q25 -0.396 q75 4.34

Not adopted. The summed excess is a sum over directions, so its
spread grows with the number of parameters as its median does: the
upper quartile passes 1 on 24 of 108 gaussian GLMs (all 12 fits each
of n = 20 and n = 50 with 15 covariates) and on nearly every K design,
and both quartiles together still fire on those 24 gaussian GLMs and on
76 K-design
fits, while they miss the item's own fit at 40 draws (q25 -0.40).
The per-direction reading does see that fit: S0.1:14 fires by
direction (joint excess 46.1 along the slope sd and correlation), and
its reweighted new-level interval is near the true width (section 2c).
A law wrong on both sides fires there through its steep side; its
flat side (the slope sd toward 0) is not read, for the reason of 2d.

## 7. The item's K table and the 7 fits, after

Part S (a separable block added to each of the 18 fired mixtures) is
decided by an identity, not a run: with a block-separable objective the
law's covariance is block-diagonal, so is its Cholesky factor, and each
`par` direction moves only its own block. The slice of a coordinate
is inside its block too. So the reading of the joined fit is the larger
of the two blocks' readings. The mixture block's own joint excess is
13.3 to 1782 (section 2b, mix40 fired: gaussian 50.2 and up, poisson
13.3 and up); the quadratic block's is 0 at every direction; the
200-level poisson block's is 0 (`kd`: no direction's own excess passes
the threshold, at most 2.39). So all 18 fire at every K of the table
(0 to 256) and with the 200-level block: 18 of 18, where the summed
check alone fires on (counted from `dev/drawdir-repro-base.txt`)

    K4 18 of 18, K8 16, K16 12, K32 7, K64 4, K128 2, K256 2

and on 2 of 18 with the 200-level block.

Part R, `dev/drawdir-repro.R` on the lane library (log
`dev/drawdir-repro-lane.txt`, lines with K > 0 in
`dev/drawdir-repro-lane-R.txt`). Of the 54 fits with covariates, 20
fire (11 summed, as before, and 9 by direction: gaussian seed 36 K 5;
poisson 5 K 10, 13 K 5, 26 K 5, 33 K 10, and 4 of the item's 7). The
7:

    seed 13 K 20  ---- dir  3.81 (mu1_z5)        z 29.04 (unchanged)
    seed 20 K 10  ---- dir  2.21 (mu2_z6)        z 18.31 (unchanged)
    seed 25 K 20  FIRE dir 77.27 (theta1)        mae 0.5662 -> 0.5352, z 26.50 -> -7.15
    seed 26 K 20  FIRE dir 22.20 (mu1_z3)        mae 0.5508 -> 0.5957, z 27.80 -> 30.33
    seed 36 K 10  ---- dir  4.97 (mu2_z5)        z 19.71 (unchanged)
    seed 36 K 20  FIRE dir 10.77 (mu2_z14)       mae 0.4375 -> 0.4269, z 29.95 -> 3.71
    seed 39 K 20  FIRE dir  8.52 (mu1_(Intercept)) mae 0.4686 -> 0.4998, z 32.56 -> 15.40

(dir values from `dev/drawdir-v-partR.txt`; the lane's "mae lane" of a
fired fit equals the base run's "forced mae" on all four, as it must:
the same reweighting.) 4 of 7 now fire, all with the warning (k > 0.7).
The error against the truth improves on 2 and worsens on 2, as the
review found with the check forced on. The other 3 stay silent: their
largest direction excess is 2.2 to 5.0, under the threshold and inside
the gaussian GLMMs' range (5.6). With 10 or 20 covariates in each
component the second component's coefficients are all weakly
identified, and the weakness is spread over those directions, each
only a little off the law. No single-direction reading can flag them
without flagging GLMMs; the summed reading cannot without flagging the
K designs (the review's point). Recorded, not built further; the docs
say it.

Among the 5 other newly fired fits, the error improves on 3 (13 K 5
0.2176 to 0.1607, 26 K 5 0.2429 to 0.2093, gaussian 36 K 5 0.1450 to
0.1236) and worsens on 2 (5 K 10 0.5966 to 0.8494, 33 K 10 0.2151 to
0.2410).

## 8. Cost (`dev/drawdir-cost.R`, log `dev/drawdir-cost-lane.txt`)

predraw2's cost fits and two large-p GLMs, lane library, single runs on
a machine shared by eleven lanes (a block under 1 s is the least of
3); "dirs" is `draw_check_dirs()` alone, "side law+drawer" the first
`predict_par_drawer()` with the reweighting stubbed out:

    F1 gaussian y ~ f(200) + x           np  202  one fn 0.0000s  side law+drawer 0.72s check(20) 0.03s  dirs 0.33s (404 evals, stat 0.0134) | first predict 1.67s  second 0.61s  fire FALSE
    F2 poisson f(60) + x + (1 | g:2000)  np   62  one fn 0.0000s  side law+drawer 6.01s check(20) 1.16s  dirs 5.94s (124 evals, stat 0.0548) | first predict 10.12s  second 0.33s  fire FALSE
    F3 categorical 8 x 7                 np   49  one fn 0.0000s  side law+drawer 0.91s check(20) 0.08s  dirs 0.59s (98 evals, stat 0.0329) | first predict 2.36s  second 0.76s  fire FALSE
    F4 zinb (1 + x | g:300), zi ~ x      np    8  one fn 0.0000s  side law+drawer 2.05s check(20) 2.15s  dirs 1.71s (16 evals, stat 0.564) | first predict 5.70s  second 0.37s  fire FALSE
    F5 ordinal mixture (fires)           np    9  one fn 0.0000s  side law+drawer 0.04s check(20) 0.00s  dirs 0.00s (32 evals, stat 29.1) | first predict 1.72s  second 1.08s  fire TRUE
    F6 poisson y ~ 0 + f(500)            np  500  one fn 0.0000s  side law+drawer 0.20s check(20) 0.03s  dirs 0.55s (1000 evals, stat 2.39) | first predict 1.71s  second 1.00s  fire FALSE
    F7 gaussian y ~ 0 + f(1000)          np 1001  one fn 0.0000s  side law+drawer 1.17s check(20) 0.36s  dirs 4.78s (2002 evals, stat 0.0269) | first predict 8.67s  second 2.77s  fire FALSE

("one fn" reads 0 because TMB returns its cached value at the point it
last evaluated, as in predraw2.) The reading costs 2p evaluations once
per fit, as the item estimated: on the 2000-level poisson GLMM 5.9 s
(the item's estimate was 3.7 s from 0.03 s per evaluation; this run's
own 20-evaluation check took 1.16 s, 0.058 s each, so the machine was
about twice as loaded), on the zinb GLMM 1.7 s, on GLMs with 200 to
500 coefficients 0.3 to 0.6 s. At p = 1001 it is 4.8 s, more than 2p
times a GLM evaluation: the dominant part there is the QR of the
1001-column factor inside `draw_logratio()` and the 1001 x 2002 matrices
of deviates (about 16 MB each). F5 fires on the summed check, so the
reading did not run inside predict(); the 32 evaluations are the
direct call's.

Not capped. Where a cap would go: memory grows as 2p^2 doubles (p =
5000 is 400 MB for the deviates, on top of the p^2 factor predict()
already holds), and time as 2p inner solves on a fit with random
effects. A cap on p (skip the reading above, say, p = 2000) would lose
the reading on exactly the fits where a weak direction is most diluted;
a time budget would make the verdict depend on the machine. I recommend
no cap now, and filing it if a large fit is seen to pay for it (the
same question as backlog 4924 for the side law).

## 9. Rounding (`dev/drawdir-ulp.R`)

The three fits of the new tests, the data times 1 + 2e-15 k, k = 0 to
20 (x for the poisson mixture, y for the gaussian GLMMs), reference
BLAS, log `dev/drawdir-ulp-ref.txt` (summarised with awk):

    grp   dir 4.37664 to 4.37664   | rel -0.05987 to -0.05987 | distinct logLik 1 | evals 9
    mixk  dir 77.27069 to 77.27073 | rel 0.05197 to 0.05197   | distinct logLik 1 | evals 118
    slope dir 46.07685 to 46.07685 | rel -0.04943 to -0.04943 | distinct logLik 1 | evals 14

Every fit keeps one verdict on all 21 roundings (grp unfired, mixk and
slope fired by direction). The nearest assertion to a threshold is
grp's `dir < 6` at 4.38 (0.32 log units) with a slice read (9 > 8
evaluations), and the `3 * draw_dir_tol` = 18 assertions at 46 and 77.

OpenBLAS with LAPACK, core Haswell, 4 threads (`dev/drawdir-ob.sh`,
outputs `dev/drawdir-ob-out/`; 0.3.26: Rblas 74 symbols and Rlapack 616
to OpenBLAS; 0.3.32: 76 and 617), the same ulp script and both draw test
files:

    0.3.26  grp 21 roundings, fire FALSE dir 4.37664 to 4.37664
    0.3.26  mixk 21 roundings, fire TRUE dir 77.27069 to 77.27073
    0.3.26  slope 21 roundings, fire TRUE dir 46.07685 to 46.07685
    0.3.32  (the same three lines)
    0.3.26  test-predict-draw-check.R pass=50 fail=0 err=0 skip=0 warn=0
    0.3.26  test-predict-draw-edge.R  pass=45 fail=0 err=0 skip=0 warn=0
    0.3.32  test-predict-draw-check.R pass=50 fail=0 err=0 skip=0 warn=0
    0.3.32  test-predict-draw-edge.R  pass=45 fail=0 err=0 skip=0 warn=0

`dev/ciharden-out/` deleted after.

## 10. Tests

`tests/testthat/test-predict-draw-check.R`, five new tests:

- "a joint wall fires the direction check, a slice's own does not": a
  synthetic 2-D objective (no fit). A wall along a + b with a and b
  correlated 0.9 fires (joint excess 14.1, above 12; 6 evaluations: 4
  directions and 2 slices); the same wall along a alone with a and b
  independent reads 0 after its slice (the guard's condition shown
  present: 5 evaluations, so a's direction passed the threshold first).
- "the conditional standard deviations need no inverse":
  `draw_dirs_csd()` against `solve()`, ratio within 1e-10.
- "a weak direction among many coefficients fires the check": the
  part R fit poisson seed 25, K = 20 (summed ratio 0.052, under 0.15;
  joint excess 77.3); fires, says "along one direction", and its draws
  differ from the normal law's.
- "a random-slope law wrong on both sides fires by direction": the
  5644 fit S0.1 seed 14 (summed excess -0.14; joint excess 46.1); fires,
  and the new level's 98% interval over 2000 draws is under 0.85 times
  the normal law's (section 2c: 6.9 against 10.7 at 4000 draws).
- "a group sd steep along its own slice does not fire": gri seed 1
  (direction excess 15.2, joint 4.38); no fire, a slice was read, and
  the draws are the normal law's bit for bit. A guard of the absent
  case: it cannot fail behaviorally on 0.83.0, which has no reading.

Seen to fail on rellib-r22 (`DRAWDIR_BASE=1 dev/drawdir-run-tests.R`,
log `dev/drawdir-base-out/test-predict-draw-check.txt`): pass 34,
fail 2, err 5. The behavioral failures:

    DETAIL FAIL [a weak direction among many coefficients fires the check] Expected `chk$fire` to be TRUE. | `actual`: FALSE
    DETAIL FAIL [a random-slope law wrong on both sides fires by direction] Expected `chk$fire` to be TRUE. | `actual`: FALSE

The errors are the weak form (`draw_check_dirs`, `draw_dirs_csd`,
`draw_dir_tol` absent). `test-predict-draw-edge.R` on rellib-r22: pass
45. On the lane build, reference BLAS: test-predict-draw-check.R pass
50, fail 0, err 0, skip 0, warn 0 (31 before the lane).

## 11. Suites

Every suite, one file per process, one process at a time (the other
slot ran the measurements, then the sessions), runner
`dev/drawdir-run-tests.R` (logs a FIRE line, with `dir`, whenever the
check fires), driver `dev/drawdir-run-par.sh`, job list
`dev/drawdir-jobs-all.lst`, log `dev/drawdir-all.log`: ALL ran 549 of
549. Per package (summed from the RESULT lines):

    frmtmb           files 349 pass 25665 fail 0 err 0 skip 243 warn 0
    frmtmb.coupling  files  11 pass   551 fail 0 err 0 skip   5 warn 0
    frmtmb.eam       files  33 pass  1851 fail 0 err 0 skip   5 warn 0
    frmtmb.latent    files  13 pass   700 fail 0 err 0 skip   2 warn 0
    frmtmb.learn     files  16 pass   474 fail 0 err 0 skip  13 warn 0
    frmtmb.ode       files  11 pass   549 fail 0 err 0 skip   1 warn 0
    frmtmb.sample    files  92 pass  4286 fail 0 err 0 skip  18 warn 0
    frmtmb.spatial   files   8 pass   551 fail 0 err 0 skip   6 warn 0
    frmtmb.survival  files  16 pass   575 fail 0 err 0 skip   0 warn 0

Every `lib:` line names core from wt-drawdir-lib (the extensions from
rellib-r22: none calls a changed function, so none was reinstalled).
The check fired in three files: this lane's two draw test files and
`test-gr-by.R`, where `(1 + x | gr(g, by = f))` with 15 groups per
level fires by direction (11.8, along `theta_3`, a by-level's
correlation) inside the test "new data that cannot be routed is
refused by name": predict() builds the draw law before it refuses the
new level, so the reweighting message prints and the refusal follows.
The test passes; the message is the random-slope class of section 2c.
That predict() reads its draws before it refuses new data is 0.83.0's
order, not this lane's.

The skips are the gates (FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS,
the scale and fuzz tiers).

Gated tier: the 47 files whose plain run skipped on those two gates (36
core, 1 eam, 1 learn, 7 sample, 2 spatial; `dev/drawdir-jobs-gated.lst`),
both gates set, log `dev/drawdir-gated.log`: GATED ran 47 of 47, pass
4736, fail 0, err 0, skip 0, warn 0; no FIRE line; every `lib:` line
the lane library.

One session per suite, three ways (`dev/release/run-sessions.sh
drawdir 1`, FRMTMB_RELLIB = the lane library; log
`dev/release/drawdir.log`, sessions `dev/release/drawdir-sessions/`):
DRAWDIR ran 27 of 27, fail 0, err 0, warn 0 in every session:

    SESSION frmtmb forward files=336/349 pass=25665 fail=0 err=0 skip=243 warn=0
    SESSION frmtmb reverse files=336/349 pass=25830 fail=0 err=0 skip=243 warn=0
    SESSION frmtmb preload files=336/349 pass=25934 fail=0 err=0 skip=245 warn=0
    SESSION frmtmb.coupling (3 ways) files=11/11 pass=551 skip=5
    SESSION frmtmb.eam      (3 ways) files=33/33 pass=1851 skip=5
    SESSION frmtmb.latent   (3 ways) files=13/13 pass=700 skip=2
    SESSION frmtmb.learn    (3 ways) files=16/16 pass=474 skip=13
    SESSION frmtmb.ode      (3 ways) files=11/11 pass=549 skip=1
    SESSION frmtmb.sample   (3 ways) files=87/92 pass=4286 skip=18
    SESSION frmtmb.survival (3 ways) files=16/16 pass=575 skip=0
    SESSION frmtmb.spatial  (3 ways) files=8/8 pass=551 skip=6

(the extension rows condense three identical printed lines each, all
fail=0 err=0 warn=0). 70 LEAK lines, every one tagged `[with load of
...]` or `[extension session state]`; none untagged. The 13 core and 5
sample files with no result are gated whole.

`R CMD check --as-cran` of core (`dev/drawdir-check.sh`, built and
checked in `dev/drawdir-check/`, then deleted; R_LIBS = the lane
library, rellib-r22, the user library; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; key lines kept in
`dev/drawdir-check-status.txt`): Status: 1 NOTE, the V8 math-rendering
NOTE on the HTML manual. Tests inside the check: OK (31 min). The
extensions are unchanged and were not checked.

The rendered Rd (`tools::Rd2txt()` of `man/predict.frmtmb_fit.Rd`, into
`dev/drawdir-rd.txt`) shows the new bullet whole.

## 12. Fix round 1 (review `dev/reviews/2026-11-17-drawdir.md`, FIX)

### FIX 1: a non-finite end is not a wall

`draw_check_dirs()` read a direction whose end the objective could not
evaluate as an infinite excess, and with a finite slice it fired with
`dir` Inf ("falls Inf further"). On slpw:12 that direction was flatter
than the law (excess -0.64 at t = 0.75) and the objective failed in the
inner solve past a log slope sd of about -65 (review section 4).

Now such a direction is read again at `draw_dir_back` = 0.75 and then
0.5 of the distance (through the side law, `draw_logratio()`), its
excess divided by t^2 to put it on the full distance's scale (exact for
a quadratic misfit), and its slice read at the same t. A direction
with no finite reading, or whose slice cannot be evaluated, is left
unread. `stat` is never infinite; where no direction is readable the
reading returns NULL and `dir` stays NA.

Test "a direction the objective cannot evaluate at its end is no wall"
(a 2-D objective, half the law's drop, NaN past a point along a + b;
the same with a wall before the NaN). On the lane build before the fix
(`dev/drawdir-fix1-before.txt`): fail 4, the behavioral ones
`is.finite(flat$stat)` FALSE and `flat$stat` Inf, and the wall's
`is.finite()` FALSE. After: pass.

slpw:12 on the fixed build (`dev/drawdir-fires12.R`, log
`dev/drawdir-fires12.txt`): 12 directions read, 12 finite; it still
fires, `dir` 202.7 along `-theta_3` (the direction that moves the log
slope sd UP by tens of units, which hid behind the Inf before), and on
`+theta_1` 48.4, `+theta_2` 26.0, `-sigma_(Intercept)` 48.0.

### FIX 3 and 6: descriptions

Corrected in `?predict.frmtmb_fit`, the roxygen of `draw_check_dirs()`
and `draw_dir_tol`, the message (`draw_check_said_what()`: "an excess
of X after the excess of that parameter moved alone, nearer the
estimate, is taken off"), NEWS, and sections 2c, 2d, 5 and 8 here:
what S2 measures (eB (sratio^3 - 1) for a one-parameter cubic skew),
the threshold as a heuristic, "true width" as the known-parameter law,
the 472 as the fits the threshold was set on, the review's survey
beside them, and the cost's attribution to `crossprod`. The
`draw_dir_tol` comment now says 1.36 to 2.04.

### FIX 7: the warning's advice at a new level

`predict()` passes whether the call predicts a grouping level the fit
did not see (`new_lv`, from `predict_new_level_spec()`) through
`predict_par_drawer(new_levels = )` to `draw_reweight_say()`. There the
k > 0.7 warning adds ", but its posterior_predict() cannot yet predict
a grouping level the fit did not see" after naming `frm_sample()`. The
same function serves the 0.83.0 summed path, so both are fixed. Test
"the warning says frm_sample() cannot predict a new level yet" (with and
without `new_levels`); before: the argument did not exist (weak form).

### FIX 4 and 5: findings and backlog

Section 5 now carries the review's 286-fit survey. Backlog: the 3
silent fits' item is LOW with the forced arm's errors (worse on 2 of
3); two LOW items filed (the warning depends on `ndraws`; the deviates'
product formed twice).

### The facts for the user's decision on random-slope fires (review item 2)

Facts, on the FIX 1 build before option (b), 20000 draws, seed 1
(`dev/drawdir-fires12.R`, log `dev/drawdir-fires12.txt`; "lane" = the
reweighted draws); "mcmc" from the review's
`dev/reviews/drawdir-review-fires.txt` (`frm_sample()`, brms's default
priors). In-sample MAE of the predictive mean against the true
conditional mean; new level's 90% coverage and width against the
known-parameter law:

    id        firing directions (component)                 MAE lane  norm    mcmc  | new level w90 lane / norm (cov90)
    slpw:12   +theta_1 48.4, +theta_2 26.0, -theta_3 202.7   7.406   7.276   6.611 | 114.8 (0.839) / 19471 (1.000); mcmc 120.1 (0.858)
              (theta), -sigma_(Intercept) 48.0 (betad)
    slpw:14   +theta_2 9.0, -theta_3 11.3 (theta)             8.329   7.132   7.189 | 116.2 (0.808) / 133.2 (0.864); mcmc 125.9 (0.844)
    slpw:15   +theta_2 27.1 (theta)                           7.489   7.870   7.650 | 124.9 (0.866) / 116.8 (0.837); mcmc 135.1 (0.894)
    slpw:18   +theta_2 16.9, -theta_3 16.0 (theta)            9.090   7.301   6.864 | 145.0 (0.924) / 156.9 (0.945); mcmc 160.7 (0.951)
    slp10:10  +theta_1 8.5, -theta_3 11.6 (theta)            12.201  10.273  10.214 | 182.4 (0.833) / 162.1 (0.778); mcmc 204.7 (0.879)
    prsp:1    +theta_2 12.1 (theta)                           0.679   0.700   0.730
    prsp:10   +theta_2 51.3, +theta_3 57.5 (theta)            0.918   0.566   0.525
    ordrs:9   +theta_3 6.17 (theta)                           0.269   0.189   (per category; not comparable)
    S0.1:14   +theta_2 42.2, -theta_3 46.1 (theta)            0.353   0.318     -   |   4.62 (0.918) /  5.20 (0.948); oracle 4.26
    S0.1:16   +theta_2 25.7, +theta_3 29.0 (theta)            0.384   0.314     -   |   4.66 (0.928) /  5.54 (0.968); oracle 4.26
    S0.1:24   +theta_2 12.7, -theta_3 11.6 (theta)            0.328   0.302     -   |   4.13 (0.874) /  4.89 (0.928); oracle 4.26
    S0.2:8    +theta_2 7.0 (theta)                            0.342   0.248     -   |   5.06 (0.939) /  6.11 (0.977); oracle 4.41

(98% widths, lane / norm: S0.1:14 6.64 / 10.69, :16 6.72 / 12.17, :24
6.00 / 8.22, S0.2:8 7.33 / 12.04; oracle 6.03 and 6.24. The review's
brms posterior on :14 and :16: w90 5.27 and 5.20, w98 7.52 and 7.37.)

Reading it:

- In-sample, the reweighted draws are worse than the normal law's on 10
  of 12 (by 2% to 62%), better on 2 (slpw:15, prsp:1). Where the review
  has the posterior, the normal law's MAE is within 10% of it on 6 of 7.
- At a new level, the normal law's 98% tail is too wide (1.36 to 2.02
  times the oracle on the four S fits at 20000 draws, 1.43 and 1.67 times the posterior
  on two), its 90% interval is close to the posterior's; the reweighted
  90% interval is narrower than the posterior's (0.89 times on S0.1:14
  and :16). On slpw:12 alone the normal law is absurd at a new level
  (w90 19471) and the reweighted one is near the posterior's.
- Every firing direction on these 12 fits is a covariance parameter
  (component `theta`: a log sd or a correlation parameter), except
  slpw:12's `-sigma_(Intercept)` (a distributional coefficient,
  component `betad`).

What option (b) would look like in code. `fit_draw_space(fit)$map$comp`
names each draw position's template component. Variance components
are the covariance components: `theta` (group-level log sds and
correlation parameters, and the `gp()` and smooth scales, which live in
`theta` too: `gp_brms_natural()`), `thetaac` (autocorrelation),
`thetar` (residual correlation), `thetame` and `logsdme` (`me()`).
`beta` and `betad` (the coefficients of mu and of distributional
parameters, sigma among them) are not. `draw_check_dirs()` would return,
beside `stat`, the components of every direction above the threshold;
`draw_check()` would set a third verdict, "warn only", where every
firing direction is a covariance component; `draw_reweight()` would
then return the draws as they came with a warning (class
`frmtmb_draw_pareto` or a new one) saying the normal law's far tail is
too wide along that parameter and naming `frm_sample()` (with the
new-level caveat). On these 12 fits, (b) would warn without reweighting
on 11 and still reweight slpw:12 (its sigma direction fires), which is
the one fit where the normal law's new-level interval is absurd. It
composes with 0.83.0's summed path without touching it: the summed
check has no direction, so a summed fire reweights as now, and the
reading runs only where the summed check did not fire. The cost is one
`comp` lookup per firing direction. Not built, pending the decision.

### Decided and built: option (b) (user decision of 2026-11-17)

(c) filed MEDIUM in the backlog with the review's numbers.

Built (`R/predict-brms.R`):

- `draw_vc_comps` = `theta`, `thetaac`, `thetar`, `logsdme`, `thetame`:
  the template's covariance components. `theta` holds every group-level
  log sd and correlation parameter and the `gp()` and smooth scales
  (`R/frame.R`: "theta - all covariance / smoothing parameters";
  `gp_brms_natural()` maps gp sds and length scales into `theta`);
  `thetaac` the autocorrelation parameters; `thetar` the residual
  correlations; `logsdme`, `thetame` the `me()` latent sds and
  correlations. Why these: each is a scale or correlation of a
  covariance, with a side toward no variance (or no correlation) along
  which the likelihood flattens, where weights by the likelihood have
  no proper mass. Not in the set: `beta`, `betad` (the coefficients of
  every linear predictor, a distributional sigma's and a mixture
  weight's among them; on slpw:12 the sigma direction's reweighting
  mended a new-level interval the normal law put at 162 times the
  posterior's) and `meanme` (a location).
- `draw_check_dirs()` returns `vc_only`: every direction above the
  threshold belongs to a component in the set (`fit_draw_space()`'s
  `map$comp`, passed in by `predict_par_drawer()` as `dw$comps`).
- `draw_check()` caches `warn_only` = fired by the reading and
  `vc_only`.
- `draw_reweight()`: where `warn_only`, returns the draws as they came
  and calls `draw_tail_say()`, a warning of class
  `frmtmb_draw_wide_tail` on every call (as the Pareto warning), which
  names the parameter, says the normal law's far tail along it is too
  wide and that the draws are not reweighted, and names `frm_sample()`
  with the new-level caveat where the call predicts a new level.
  Every other fire, and every summed fire, reweights as before.

On the 12 fires (`dev/drawdir-fires12.R`, log `dev/drawdir-fires12b.txt`,
20000 draws): 11 warn (`frmtmb_draw_wide_tail`) and their in-sample MAE
and new-level coverage and widths are those of the normal law to every
printed digit (for example S0.1:14 MAE 0.3176, cov90 0.948, w98 10.69 on
both arms); slpw:12 still reweights (its `-sigma_(Intercept)` direction
is a `betad` coefficient): MAE 7.4062, w90 114.78, k 0.8, as before.

Tests (`tests/testthat/test-predict-draw-check.R`):

- "a random-slope law wrong on both sides warns, not reweights" (S0.1
  seed 14, new level, 500 draws): fires, `warn_only` TRUE, the warning
  "far tail" with "cannot yet predict a grouping level", and the draws
  identical to those with `draw_reweight()` mocked out. Seen failing on
  the fix-round build before (b) (`dev/drawdir-fixb-before.txt`): fail
  4, the behavioral ones: no "far tail" warning, no new-level caveat,
  `warn_only` FALSE, and `p1` not identical to `p0` (reweighted:
  0.5836 against 3.2831 in the first draw).
- "a weak direction among many coefficients fires the check" (the
  mixture, whose direction is the mixing weight's coefficient, `betad`):
  now also asserts `warn_only` FALSE; its draws still differ from the
  normal law's (reweighted).
- Unchanged where nothing fires, bit for bit: "a fit the normal law
  describes keeps its draws and stream" and "a group sd steep along its
  own slice does not fire" (both `expect_identical` against the mocked
  arm) pass.
- `tests/testthat/test-gr-by.R` (outside this lane's files, one hunk):
  its refusal test predicts on a `(1 + x | gr(g, by = f))` fit that
  fires the reading on a by-level's correlation (`theta_3`), before
  predict() refuses the new level; the new warning is allowed there
  (`allow_warnings(..., "far tail")`, not required).

After (b): `test-predict-draw-check.R` pass 62, `test-predict-draw-edge.R`
45, `test-gr-by.R` 71, fail 0, err 0, warn 0 (logs
`dev/drawdir-fixb-*.txt`).

OpenBLAS 0.3.26 with LAPACK (`dev/drawdir-ob2.sh`, outputs
`dev/drawdir-ob2-out/`; Rblas 74 symbols and Rlapack 616 to OpenBLAS),
final fix-round build: `test-predict-draw-check.R` pass 62,
`test-predict-draw-edge.R` pass 45, fail 0, err 0, warn 0; the ulp
script over 21 roundings: grp unfired at 4.37664, mixk fired at
77.27069 to 77.27073, slope fired at 46.07685, one verdict each.
`dev/ciharden-out/` deleted.

### Fix-round suites (final build)

The core test files that read predictions (218 files that name
`predict`, `posterior_epred`, `pp_check`, `conditional_effects` or
`loo_predict`) and all of frmtmb.sample, one file per process
(`dev/drawdir-jobs-fix1.lst`, log `dev/drawdir-fix1.log`): FIX1 ran 310
of 310.

    frmtmb         files 218 pass 18297 fail 0 err 0 skip 179 warn 0
    frmtmb.sample  files  92 pass  4286 fail 0 err 0 skip  18 warn 0

Every `lib:` line the lane library. FIRE lines only in the two draw
test files and `test-gr-by.R`.

Core in one session three ways (`dev/release/run-sessions.sh drawdirf 1
frmtmb`, FRMTMB_RELLIB = the lane library; log
`dev/release/drawdirf.log`): DRAWDIRF ran 3 of 3.

    SESSION frmtmb forward files=336/349 pass=25677 fail=0 err=0 skip=243 warn=0
    SESSION frmtmb reverse files=336/349 pass=25842 fail=0 err=0 skip=243 warn=0
    SESSION frmtmb preload files=336/349 pass=25946 fail=0 err=0 skip=245 warn=0

20 LEAK lines, every one tagged; none untagged.

`R CMD check --as-cran` of core (`dev/drawdir-check.sh`, in
`dev/drawdir-check/`, deleted after; key lines in
`dev/drawdir-check-status2.txt`): Status: 1 NOTE, the V8 math-rendering
NOTE. Tests inside the check OK (27 min). The rendered Rd
(`dev/drawdir-rd2.txt`) shows the corrected bullet and the
variance-component sentence whole.

Process: before option (b) was built, the fix round's first rerun and
core sessions were running on the FIX 1 build; they were stopped
(TaskStop, then Stop-Process on this lane's own processes, matched on
the worktree's path) so that the library was not reinstalled under a
reading process. Their partial logs were overwritten by the reruns
above. At most two of this lane's R processes ran at any time.

### After the short check (`dev/reviews/drawdir-short.md`, SHIP)

Record and doc edits only:

- Backlog, three LOW items filed in this lane's section: the message
  on a mixed fire names the largest direction (slpw:12: `theta_3`, a
  variance component) and not the one that made the draws reweight
  (`-sigma_(Intercept)`); a direction non-finite at 1, 0.75 and 0.5 is
  left unread, so a hard wall within half the distance goes unreported;
  predict() warns about the draw law before refusing an unrouted new
  level on a `gr(by = )` fit.
- `?predict.frmtmb_fit` states the set as the short check found it in
  `R/frame.R`: every covariance-structure parameter (sds and
  correlations, ranges, Matern smoothness, decay rates, ar
  correlations, `rr()` loadings, spatial structures, the gp and smooth
  scales, the latent-residual sd and its autocorrelation,
  autocorrelation, residual-correlation and `me()` parameters), not
  only variances; and that where another direction fires too the draws
  are reweighted. The roxygen of `draw_vc_comps` says the same. The
  paragraph after it was reflowed (one line was 96 columns).
- `tests/testthat/test-gr-by.R` allows "normal law's far tail" instead
  of "far tail" (an `me()` warning also says "far tail").
- Reruns, lane library: `test-predict-draw-check.R` pass 62,
  `test-gr-by.R` pass 71, fail 0, err 0, warn 0 (`dev/drawdir-fix2-*.txt`);
  the rendered Rd `dev/drawdir-rd3.txt` shows the paragraph whole.
  `R CMD check --as-cran` of core (`dev/drawdir-check.sh`, key lines in
  `dev/drawdir-check-status3.txt`, check directory deleted): Status: 1
  NOTE, the V8 math-rendering NOTE; tests inside the check OK (26 min).

## Not done, and why

- The 3 silent part R fits (a weakness spread over many directions):
  filed MEDIUM. A subspace reading was not tried.
- No cap on the reading's cost: filed LOW, with the recommendation
  not to cap now (section 8).
- Flat directions are not read (section 2d).
- 5617 (the near-empty component's heavy tail): the reading does not
  see it (section 4); nothing more built, as the brief said.
- The reading uses the same reweighting as 0.83.0. Where it fires,
  the Pareto k is above 0.7 on 22 of the 23 newly fired mixtures and on
  the random-slope fits, so the usual outcome is the warning naming
  `frm_sample()`.
- brms 2.23.0 was not run: the change is to frmtmb's own draw law, which
  brms does not have (brms samples), and the outcome measures here are
  against the simulated truth.

## Files

- `R/predict-brms.R`: `predict_par_drawer()` (one line: the labels to
  the check), `draw_reweight()` (the no-loo warning's text through
  `draw_check_said_what()`), `draw_check()` (calls the reading where the
  summed check does not fire; three new fields), new
  `draw_check_dirs()`, `draw_dirs_csd()`, `draw_check_said_what()`,
  `draw_reweight_say()` (its head through `draw_check_said_what()`), new
  constant `draw_dir_tol`; the roxygen of `draw_reweight()` and a third
  bullet in `?predict.frmtmb_fit`'s section "When the normal law does
  not describe the likelihood".
- `man/predict.frmtmb_fit.Rd` (roxygenised).
- `NEWS.md`: a new top section `# frmtmb (development version)`,
  "## Predictions", one BREAKING bullet.
- `tests/testthat/test-predict-draw-check.R`: five new tests; in fix
  round 1 two more (a non-finite end; the new-level caveat), the
  random-slope test rewritten for option (b), one assertion added to
  the mixture test.
- Fix round 1 in `R/predict-brms.R`: the call site in `predict()` (the
  new-level flag), `predict_par_drawer(new_levels = )` and `comps`,
  `draw_reweight()` (the warn-only branch), `draw_check()`
  (`warn_only`), `draw_check_dirs()` (reading back from a non-finite
  end; `vc_only`), `draw_reweight_say(new_levels = )`, new
  `draw_tail_say()`, `draw_dir_back`, `draw_vc_comps`; the roxygen and
  `?predict.frmtmb_fit` corrected.
- `tests/testthat/test-gr-by.R` (fix round 1, outside this lane's
  area): one test allows the new `frmtmb_draw_wide_tail` warning.
- `dev/test-backlog.md`: CLOSED lines under 5584 and 5644; section
  "Filed in the round of 2026-11-16", "### From lane drawdir", two items.
- Lane scripts `dev/drawdir-*.R`, `dev/drawdir-*.sh`, job lists
  `dev/drawdir-jobs-*.lst`, and this record.

No extension changed. Bump: minor (a BREAKING behavior change in
predict(): more fits are reweighted and warned on). Floor: none new;
loo stays in Suggests.

## Process note

Once, for about a minute, three of this lane's R processes ran at once:
a run of `test-predict-draw-check.R` started while the two mixture
coverage arms (section 4) were running. Every later run kept to two.
No action was denied by the permission system.
