# Lane loglik: log_lik() for residual covariance terms, and the 13 bin-2 sar() and fcor() rows

Base 70476dc3 (frmtmb 0.75.0), worktree `frmtmb-wt-loglik`, branch
`wt-loglik`, library `C:/Users/adf44/source/r/wt-loglik-lib` (core,
frmtmb.sample and frmtmb.spatial from this tree) in front of
`rellib-r13`. Nothing is committed. Logs are under `dev/loglik-*.txt`,
`dev/loglik-log/` and `dev/release/loglik-*` (all gitignored); every
number below is pasted from them, with its script and seed.

The brief named `C:\Users\adf44\source\r\frmtmb-wt-sarfcor\dev\sarfcor-brmsrows.R`.
That worktree was removed during the round; the tracked copy
`dev/sarfcor-brmsrows.R` in this tree is the one meant (the
coordinator's correction).

## Part A. brms's leave-one-out conditional density

### A1. What brms 2.23.0 does (run, not recalled)

brms's `R/log_lik.R` (`dev/brms-suite/brms/R/log_lik.R`, lines 246 to
412 and 1078 to 1100): for `ar()`, `ma()`, `arma()` (`cov = TRUE`),
`cosy()` and `unstr()` one call per time-series GROUP
(`log_lik_gaussian_time()`, `begin_tg:end_tg`); for `fcor()` and both
SAR types one call over all rows (`stopifnot(i == 1)`). Per draw, a
dense `solve()` of the group's covariance `C`, then
`g = C^-1 e`, `cbar = diag(C^-1)`, `yloo = y - g / cbar`,
`sdloo = 1 / sqrt(cbar)`; for student `dfloo = nu + k - 1` and
`sdloo^2 = (nu + e'g - g^2 / cbar) / (nu + k - 1) / cbar`
(`student_t_cov_factor()`, Shah et al. 2014 eq. 6), `k` the group size.
For the lag SAR, `e = y - (I - rho W)^-1 mu`. `get_cov_matrix_ac()`
reads lags by position within the group (only `unstr()` carries
`Jtime`), and `sigma ~ x` enters as `outer(sigma, sigma) * cor`.

brms was run on every case (`dev/loglik-brms-probe.R`, seed 1, log
`dev/loglik-brms-probe.txt`; grouped student cases again in
`dev/loglik-brms-probe2.R`, `dev/loglik-brms-probe2.txt`): a brms fit
whose draws are two chosen parameter vectors (Fixed_param, one chain
each; `dev/loglik-brms-lib.R`), and `log_lik()` against the conditional
computed directly from the joint covariance by the partition formula,
with no precision matrix. 24 rows in 5 groups of sizes 6, 5, 4, 6, 3,
rows shuffled:

    ar gaussian                draws 2  rows 24  max|brms - dense| / max|dense| = 3.981e-16
    ma gaussian                draws 2  rows 24  max|brms - dense| / max|dense| = 2.876e-16
    arma gaussian              draws 2  rows 24  max|brms - dense| / max|dense| = 3.597e-16
    cosy gaussian              draws 2  rows 24  max|brms - dense| / max|dense| = 2.774e-16
    fcor gaussian              draws 2  rows 24  max|brms - dense| / max|dense| = 6.253e-16
    fcor student               draws 2  rows 24  max|brms - dense| / max|dense| = 6.570e-16
    lag gaussian               draws 2  rows 24  max|brms - dense| / max|dense| = 4.111e-16
    lag student                draws 2  rows 24  max|brms - dense| / max|dense| = 4.396e-16
    err gaussian               draws 2  rows 24  max|brms - dense| / max|dense| = 3.378e-16
    err student                draws 2  rows 24  max|brms - dense| / max|dense| = 3.768e-16
    ar sigma~x gaussian        draws 2  rows 24  max|brms - dense| / max|dense| = 5.196e-16
    unstr gap gaussian         draws 2  rows 22  max|brms - dense| / max|dense| = 2.670e-16
    ar gap gaussian: brms vs by-position 2.955e-16, vs by-level 6.752e-02
    ar gaussian draw 1: sum of conditionals -25.055181, joint -26.097025

and, with the reference taken group by group (probe2):

    ar student                 draws 2  rows 24  max|brms - dense| / max|dense| = 5.532e-16
    ma student                 draws 2  rows 24  max|brms - dense| / max|dense| = 4.569e-16
    arma student               draws 2  rows 24  max|brms - dense| / max|dense| = 5.160e-16
    cosy student               draws 2  rows 24  max|brms - dense| / max|dense| = 4.352e-16
    ar sigma~x student         draws 2  rows 24  max|brms - dense| / max|dense| = 5.991e-16
    unstr gap student          draws 2  rows 22  max|brms - dense| / max|dense| = 5.879e-16

So brms's `log_lik()` is the exact conditional for every term and both
families; no brms defect here. The first probe's grouped student rows
differed by 0.10 to 0.17: its reference was one multivariate t over all
24 rows, a different model from brms's one t per group (each group
draws its own chi-square), so those numbers say nothing about brms and
probe2 replaces them. brms's `log_lik()` works on the error SAR, where
its `posterior_epred()` fails (brms-39, already filed). The ARMA(1,1)
covariance brms uses has gamma0 = (1 + 2 phi theta + theta^2) / (1 -
phi^2) at unit innovations, which the probe's reference assumed and
which agreed.

### A2. The route, and why it is exact and cheap

`R/autocor-loglik.R`, `autocor_row_loglik(fit, dpv, resp)`, exported
in the sampling API (`R/sampling-api.R`). Everything is read off the
precision `Q` of the response; neither `Q` nor the covariance of more
than one group is formed:

- the `ar()` family (`ar`, `ma`, `arma` with `cov = TRUE`, `cosy`,
  `unstr`): `Q = D^-1 R^-1 D^-1` per group, `D = diag(sigma)`. One
  `chol2inv(chol(R_p))` per PATTERN of present time levels (the
  patterns the objective already uses), shared by every group of the
  pattern, applied to the `k x G` matrix of standardized residuals in
  one product. Unequal group sizes, missing time points and `sigma ~ x`
  come from the patterns and the per-row sigma;
- a registered term (`fcor()`, `sar()`, any `frmtmb_register_autocor()`
  term): the whitening map is affine, `u = K (y - m)`, so `Q = K'K`,
  `g = K'u`, `cbar = colSums(K^2)`, `e'Qe = u'u`. `K` comes from a new
  documented slot, `frmtmb_register_autocor(jacobian =)`; `fcor()`
  gives `L^-1 D^-1`, frmtmb.spatial's `sar()` gives the sparse
  `(I - rho W) / sigma` (both types: they differ in the mean only).
  Without the slot, `K = B^-1` with `B` from `unwhiten()` of the
  identity: exact, dense;
- the density is written in the residual `g / cbar` itself (which is
  `y - yloo`), so no `y` is subtracted from a nearly equal `yloo`.

The review's estimate (`dev/reviews/2026-10-24-sarfcor.md`, finding 12)
was this route; brms's grouping and its per-group t were the parts to
verify, and A1 verified them.

Cost, measured (`dev/loglik-cost.R`, seed 5, log `dev/loglik-cost.txt`;
arms interleaved in one process, each block repeated past 1.2 s, min of
5 rounds, a control arm that repeats the first; run alone on the
machine's other load):

    sar n = 1600  nnz(W) = 6240
    agreement: sar vs brms 9.107e-16, sarfb vs sar 6.960e-15 (max abs / max)
    sar, per draw, min over 5 rounds:
      sar        1.924 ms per draw (reps 1024)  alloc    1809120 bytes
      sarfb   1730.000 ms per draw (reps    1)  alloc  207899968 bytes
      brms    4120.000 ms per draw (reps    1)  alloc  102560792 bytes
      ctrl       1.895 ms per draw (reps  512)  alloc    1809120 bytes
      brms / sar = 2141.6, sarfb / sar = 899.2, ctrl / sar = 0.98
    ar agreement with brms's algorithm: 6.193e-16
    ar(cov = TRUE), 400 groups x 10, per draw:
      ar         0.635 ms per draw (reps 2048)  alloc     729536 bytes
      brms      17.813 ms per draw (reps   64)  alloc   15536048 bytes
      ctrl       0.605 ms per draw (reps 2048)  alloc     729536 bytes
      brms / ar = 28.1, ctrl / ar = 0.95
    fcor agreement with brms's algorithm: 1.843e-15
    fcor(), 600 rows, per draw:
      fcor       2.637 ms per draw (reps  512)  alloc    5844960 bytes
      brms      97.500 ms per draw (reps   16)  alloc   14482512 bytes
      ctrl       2.695 ms per draw (reps  512)  alloc    5844960 bytes
      brms / fcor = 37.0, ctrl / fcor = 1.02

`sar` is the term with its jacobian slot, `sarfb` the same map
registered without it (the fallback), `brms` brms's algorithm
transcribed (`log_lik_gaussian_lagsar()` and `log_lik_gaussian_time()`,
dense). Allocation is `Rprofmem()` bytes per call. The control ratios
(0.95 to 1.02) bound the clock's noise; the ratios that matter are 28
to 2142. The fallback is 900 times the slot on the 1600-row lattice,
which is why `sar()` fills the slot; a term registered without one
still gets exact values. `fcor()` allocates two `n x n` matrices per
draw (`K` and `K * K`); a precomputed `colSums(L^-1 ^2)` would save one
of them, not done because the term is dense by nature and is 37 times
brms already.

### A3. Agreement with brms at the same draws

`dev/loglik-brms-agree.R` (seed 4, log `dev/loglik-brms-agree.txt`):
frmtmb.sample's `log_lik()` on hand-set draws (no sampler) against
brms's `log_lik()` on a brms fit whose three chains sit at the same
three parameter vectors. 43 rows in 8 groups of sizes 7, 5, 6, 4, 7,
3, 6, 5, rows shuffled; the `cosy()` and `unstr()` rows on the same
data with two time points removed (41 rows):

    ar gaussian            3 x 43  max|frmtmb - brms| / max|brms| = 2.033e-15
    ar student             3 x 43  max|frmtmb - brms| / max|brms| = 8.183e-16
    ma gaussian            3 x 43  max|frmtmb - brms| / max|brms| = 6.853e-16
    ma student             3 x 43  max|frmtmb - brms| / max|brms| = 7.070e-16
    arma gaussian          3 x 43  max|frmtmb - brms| / max|brms| = 1.185e-15
    arma student           3 x 43  max|frmtmb - brms| / max|brms| = 7.988e-16
    ar sigma~x gaussian    3 x 43  max|frmtmb - brms| / max|brms| = 6.163e-16
    ar sigma~x student     3 x 43  max|frmtmb - brms| / max|brms| = 4.288e-16
    cosy gaps gaussian     3 x 41  max|frmtmb - brms| / max|brms| = 6.024e-16
    cosy gaps student      3 x 41  max|frmtmb - brms| / max|brms| = 4.040e-16
    unstr gaps gaussian    3 x 41  max|frmtmb - brms| / max|brms| = 9.479e-16
    unstr gaps student     3 x 41  max|frmtmb - brms| / max|brms| = 5.220e-16
    fcor(V) gaussian       3 x 43  max|frmtmb - brms| / max|brms| = 1.210e-15
    fcor(V) student        3 x 43  max|frmtmb - brms| / max|brms| = 9.476e-16
    sar(W) gaussian        3 x 43  max|frmtmb - brms| / max|brms| = 6.347e-16
    sar(W) student         3 x 43  max|frmtmb - brms| / max|brms| = 5.366e-16
    sar(W, type = 'error') gaussian 3 x 43  max|frmtmb - brms| / max|brms| = 5.928e-16
    sar(W, type = 'error') student 3 x 43  max|frmtmb - brms| / max|brms| = 4.822e-16

The translation between the two parameterizations is the one the
log-density tier uses: brms's `ar()`, `ma()` and `arma()` sigma is the
innovation sd, frmtmb's the marginal one (`sigma_brms = sigma /
sqrt(gamma0)`, and for `sigma ~ x` the intercept shifts by
`log sqrt(1 - ar^2)`); brms's `Lcortime` is the Cholesky factor of
frmtmb's `unstr()` matrix; `lagsar`/`errorsar` are rho. The ar() rows
use data with no interior gap: with one, brms (lags by position) and
frmtmb (lags by time level, `?frmtmb-autocor`) fit different models,
and the probe's "ar gap" line shows brms's values are the by-position
conditional and differ from the by-level one by 6.8e-2. `cosy()` has
one correlation for every pair and brms's `unstr()` carries the time
index, so those agree with gaps.

The same comparisons are tests (A5). The ML-estimate check
(`dev/loglik-check.R`, seed 3, `dev/loglik-check.txt`) compares
`autocor_row_loglik()` at each fit's estimates with the dense reference
built from the fit's own covariance, for 10 formulas x 2 families on
ragged groups with an interior gap: 1.3e-16 to 8.7e-16, except
`unstr()`, whose estimates on 7 levels over 8 groups are near singular
(false convergence; the reference itself is NaN), so the tests set
`unstr()`'s correlations by hand. The same log shows the row sums are
not the joint log density, e.g. `fcor` gaussian row sum -64.1559,
joint -62.6940.

### A4. The seam slot, and core's log_lik() on an ML fit

`frmtmb_register_autocor(jacobian = function(sigma, theta, aux))`:
documented in a new section "The Jacobian" of its help page. When the
frame is built a given slot is checked once on numbers against the map
(`K (y - m) = whiten(y)$u`, relative 1e-6), and a wrong one is refused
by the term's name. frmtmb.spatial fills it through the exported
function alone: no `frmtmb:::` (grep of `extensions/frmtmb.spatial/R`
for `frmtmb:::` finds none).

Core's `log_lik()` on a `frmtmb_fit` keeps its refusal: it refuses for
EVERY model ("the pointwise likelihood at every posterior draw ... this
is a maximum-likelihood fit"), and making it answer for the covariance
terms alone would make the one family of models where a fit has a
pointwise density the exception. The computation at one parameter
vector is core's exported `autocor_row_loglik()`, which on an ML fit is
the value at the estimates (the A3 check runs it that way). The
coordinator decides whether the brief meant more.

### A5. Tests

New:

- `tests/testthat/test-autocor-loglik.R` (core): the dense conditional
  for ar, ar(p = 2), ma, arma, cosy, ar + (1 | g), ar with sigma ~ x
  (gaussian and student), unstr at hand-set correlations, fcor with
  and without sigma ~ x and student; the row sum against the joint
  density (it differs by more than 1 percent at ar = 0.8, and equals it
  at ar = 0, where each conditional is the row's own density; the joint
  is checked against `build_objective()` too); the jacobian slot and the
  fallback against fcor(); a wrong jacobian refused, with fcor()'s own
  slot as the absent case; the refusal on a response without a
  covariance term.
- `extensions/frmtmb.sample/tests/testthat/test-loglik-autocor.R`:
  `log_lik()` of hand-set draws against the dense conditional for ar,
  ma, arma (with the correlation from `stats::ARMAacf()`, not from
  frmtmb), cosy and unstr with gaps, fcor with sigma ~ x, gaussian and
  student, each column identical to core's function at the draw; a
  registered term without a jacobian slot; `loo()`, `waic()`,
  `loo_epred()`, `loo_predict()` on a covariance block; and, gated,
  brms's `log_lik()` at the same draws for every row of A3 but sar.
- `extensions/frmtmb.spatial/tests/testthat/test-sar.R`, two new
  blocks: the dense conditional for both types, both families, and the
  row sum against the joint; gated, brms's `log_lik()` at three
  parameter vectors per type and family.

Changed (the refusals removed, each pin now an agreement test):

- `extensions/frmtmb.sample/tests/testthat/test-loo.R`, "log_lik()
  refuses likelihoods with no per-observation column": the
  `ar(cov = TRUE)` arm asserted the refusal; it now checks the columns
  against the dense conditional of each subject at four draws and that
  `loo()` is finite. The mixture arm is unchanged.
- the same file, "loo() runs on a cov = FALSE fit and cov = TRUE still
  refuses", renamed "... and on a cov = TRUE one": the three refusal
  messages it pinned (ar, unstr, and that unstr's message offers no
  cov = FALSE form) are replaced by: each cov = TRUE column is core's
  function at the draw, it differs from the cov = FALSE column at
  ar = 0.5 and equals it at ar = 0, `loo()` is finite, and unstr gives
  finite columns.

Seen to fail on rellib-r13 (logs `dev/loglik-log/base-*.txt`):

    RESULT frmtmb test-autocor-loglik.R pass=0 fail=0 err=7 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=0 fail=0 err=6 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=132 fail=0 err=2 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=93 fail=0 err=2 skip=0 warn=0

frmtmb.sample's errors are the behavioral refusal ("log_lik() needs a
likelihood that factors into one term per observation, and the
response 'y' does not: it carries ar(t, g, cov = TRUE) ..."), for every
block, the registered term without a jacobian among them. Core's and
frmtmb.spatial's are the weak form, "could not find function
autocor_row_loglik"; the behavioral failure for `sar()` is the same
frmtmb.sample refusal, which the registered lag term of
test-loglik-autocor.R shows.

On the lane build (one file per process, gated):

    RESULT frmtmb test-autocor-loglik.R pass=47 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=93 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=146 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=109 fail=0 err=0 skip=0 warn=0

## Part B. The 13 bin-2 sar() and fcor() rows

Adapter (`dev/brmsport-bin2-adapters.R`, section "Residual correlation
terms"; spec `dev/brmsport-bin2-adapters-spec.R`): brms's
`log_lik_<family>_<term>(1, prep)` and
`posterior_predict_<family>_<term>(1, prep)` for `lagsar`, `errorsar`
and `fcor` become `brms_adapt_log_lik()` and
`brms_adapt_posterior_predict()`, which build a frmtmb model of the
prep's shape (response `Y`, one coefficient per row set to brms's `mu`
of the draw, the term over brms's matrix, brms's sigma, nu and rho per
draw) and return frmtmb's `autocor_row_loglik()` or `sim_draw()`.
ADAPTATION, recorded in the spec header: the SAR blocks'
`Msar = diag(10)` becomes `brms_adapt_sar_path(10)`, the
row-standardized path graph over the same 10 rows, because frmtmb
refuses `diag(10)` by design (`I - rho I` makes rho and the scale one
parameter). The substitution applies to the setup line, so brms's own
function in the value checks reads the same matrix. Blocks added to
`bin2_adapted`: "log_lik for SAR models runs without errors", "log_lik
for FCOR models runs without errors", "log_lik for FCOR student-t models
works as expected", "loglik for SAR models runs without errors" (brms's
label in tests.posterior_predict.R), "posterior_predict for FCOR models
runs without errors".

`Rscript dev/brmsport-bin2-gen.R`: "wrote 6 files; 238 assertion ids;
13 without a verdict; 149 adapted; 71 value checks". Only
`dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R` (+199 lines)
and `-posterior_predict.R` (+90) changed. The 13 without a verdict are
the 13 rows; the ledger and `dev/brmsport-bin2-verdicts.tsv` were not
touched, so the rows are "unclassified" until consolidation.

Recorded on the lane build (`FRMTMB_PORT_LIB` = the lane library,
`R_LIBS_SITE` = rellib-r13 for the extensions the other blocks load;
logs `dev/loglik-log/run-frmtmb-*.txt`, records
`dev/loglik-log/rec-frmtmb-*.tsv`):

    RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=97 fail=0 err=0 skip=0 blocks=15
    RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=61 fail=0 err=0 skip=0 blocks=14
    assert frmtmb log_lik:150 unclassified TRUE
    value  frmtmb log_lik:150  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 1.44e-16
    assert frmtmb log_lik:152 unclassified TRUE
    value  frmtmb log_lik:152  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 3.12e-16
    assert frmtmb log_lik:157 unclassified TRUE
    value  frmtmb log_lik:157  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 1.44e-16
    assert frmtmb log_lik:159 unclassified TRUE
    value  frmtmb log_lik:159  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 3.06e-16
    assert frmtmb log_lik:174 unclassified TRUE
    value  frmtmb log_lik:174  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 1.97e-16
    assert frmtmb log_lik:176 unclassified TRUE
    value  frmtmb log_lik:176  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 1.62e-16
    assert frmtmb log_lik:230 unclassified TRUE
    assert frmtmb posterior_predict:111 unclassified TRUE
    assert frmtmb posterior_predict:113 unclassified TRUE
    assert frmtmb posterior_predict:118 unclassified TRUE
    assert frmtmb posterior_predict:120 unclassified TRUE
    assert frmtmb posterior_predict:133 unclassified TRUE
    assert frmtmb posterior_predict:135 unclassified TRUE

(the fifth column is `held`). Every one of the 13 holds, so each would
be a pass at consolidation; the 6 value checks hold against brms's own
`log_lik_*()` on the same prep at 3.1e-16 relative or better, and
log_lik:230 is itself a value assertion (brms's naive multivariate-t
conditional, `expect_equal(ll, ll_student)`) and holds. The 4 SAR
posterior_predict rows and the 4 SAR log_lik rows pass through the
path-graph adaptation; the 2 fcor posterior_predict rows and the 3
fcor log_lik rows pass as written.

The base arm (rellib-r13, `dev/loglik-log/base-rec-frmtmb-*.tsv`): the
7 log_lik rows do not hold ("'autocor_row_loglik' is not an exported
object"), the 6 posterior_predict rows hold, as the sarfcor replay
found at 0.75.0. So the port moves 13 rows from gap to pass, 7 of them
because of part A.

Not ported, and now portable: the two "log_lik for ARMA models runs
without errors" rows (log_lik:130, 133; survey: "brms's row is the
leave-one-out conditional density ... which frmtmb does not compute")
read brms's prep of per-draw covariance matrices, which an adapter
would have to map onto a time-series fit; outside the brief's 13.

## Open items and defects

- No new brms defect. brms-39 (`posterior_epred()` on an error SAR)
  stands; brms's `log_lik()` on the error SAR works.
- Gaps under `ar()`, `ma()`, `arma()`: frmtmb's conditional is the
  conditional of its own (by-level) model, which differs from brms's
  by-position one exactly as the two likelihoods differ; documented in
  `?frmtmb-autocor`, "The pointwise log-density".
- `fcor()` allocates two `n x n` matrices per draw (A2).

## Suites

Every suite the change can reach: core and all eight extensions, the
lane's core, frmtmb.sample and frmtmb.spatial in front of rellib-r13
(the `lib:` lines show `wt-loglik-lib` for those three and rellib-r13
for the other six).

One file per process, gated (`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`,
`dev/release/run-par.sh loglik-gated`, 8 at once; per-file logs
`dev/release/loglik-gated-files/`): "LOGLIK-GATED ran 452 of 452",
totals from the per-file RESULT lines:

    frmtmb           files 280 pass 22243 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling  files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam       files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent    files 12 pass 681 fail 0 err 0 skip 2 warn 0
    frmtmb.learn     files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode       files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample    files 68 pass 3423 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial   files 8 pass 544 fail 0 err 0 skip 0 warn 0
    frmtmb.survival  files 15 pass 552 fail 0 err 0 skip 0 warn 0

The 29 skips are in other tiers' gates: test-drmtmb-agreement.R (13,
"set FRMTMB_DRMTMB_FIT_TESTS=true"), each package's test-scale.R,
test-fuzz.R and test-curve-scale.R. None is on a path this lane
changed.

One session per suite, three orders (`NOT_CRAN=true`,
`dev/release/run-sessions.sh loglik-sess 8`): "LOGLIK-SESS ran 27 of
27", fail 0 and err 0 in every session, 46 LEAK lines, every one tagged
"[with load of ...]", none untagged:

    SESSION frmtmb forward files=270/280 pass=19807 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=270/280 pass=19957 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=270/280 pass=20050 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3125 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3125 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3125 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=461 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=461 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=461 fail=0 err=0 skip=6 warn=0

frmtmb.sample's 5 files with no RESULT in a session are its five
`test-brms-suite-*.R` files, which skip as a whole outside the suite
tier; core's 10 are the same kind.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
"Rblas: 74 to OpenBLAS, 2 to Rblasref", "Rlapack: 616 to OpenBLAS, 10
to Rlapackref"; `OPENBLAS_NUM_THREADS=4`; gated; logs
`dev/loglik-log/ob-*.txt`):

    RESULT frmtmb test-autocor-loglik.R pass=47 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-fcor.R pass=57 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=93 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=146 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=109 fail=0 err=0 skip=0 warn=0

0.3.32 was not run: the new code has no branch that rounding decides
(no threshold on a computed value), and every assertion is a ratio to
the reference's own size with margins of 1e5 or more over the observed
differences. The OpenBLAS copy (`dev/ciharden-out/`) is deleted.

## R CMD check --as-cran

Once per changed package, built and checked in `dev/loglik-check/`
(`dev/loglik-rcheck.sh`; `R_LIBS` the lane library, rellib-r13, the
user library; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, `NOT_CRAN=true`; logs
`dev/loglik-log/check-*.txt`):

    frmtmb          Status: 1 NOTE
    frmtmb.sample   Status: 1 NOTE   (tests 344 s, OK)
    frmtmb.spatial  Status: 1 NOTE   (tests 29 s, OK)

Each NOTE is the expected "Skipping checking math rendering: package
'V8' unavailable". `dev/loglik-check/` and the bin-2 runner's staging
directory are deleted.

## Bumps and floors

- frmtmb: minor (a new export in the sampling API,
  `autocor_row_loglik()`, and a new argument of an exported seam,
  `frmtmb_register_autocor(jacobian =)`).
- frmtmb.sample: minor (`log_lik()` and the leave-one-out functions on
  residual covariance terms); floor the new core version, for
  `autocor_row_loglik()`.
- frmtmb.spatial: minor (`sar()` reaches `log_lik()`); floor the new
  core version, for `jacobian =` (an older core refuses the unknown
  argument when the package loads).

## Overlaps with the other lanes

`R/autocor-seam.R` (acprior touches `frmtmb_autocor_pars()` there; my
edits are the registration's `jacobian` argument, the map check and
`fcor_jacobian()`), `R/autocor.R` (doc sections only),
`R/sampling-api.R` (one alias and one export, each on its own line),
`NAMESPACE`, `NEWS.md`, `man/frmtmb-autocor.Rd`,
`man/frmtmb-sampling-api.Rd`, `man/frmtmb_register_autocor.Rd`,
frmtmb.sample (`R/loo.R`, `NEWS.md`, `man/sample-log_lik.Rd`,
`man/sample-loo.Rd`, `tests/testthat/test-loo.R`), frmtmb.spatial
(`R/sar.R`: the registration call gains `jacobian =` on its own line;
`NEWS.md`, `man/sar.Rd`, `tests/testthat/test-sar.R`), and
`dev/brmsport-bin2-adapters*.R` with the two regenerated bin-2 files.

## Fix round 1

Review `dev/reviews/2026-10-27-loglik.md`, verdict FIX. Numbered as the
review.

### 1. The Jacobian check at a generic point

Seen first on the build under review: the reviewer's
`dev/review-loglik-seam.R` (seed 81), run unedited on the lane library
before the fix (`dev/loglik-log/rv-seam-before.txt`):

    rvnosig sigma ~ x                              ACCEPTED, rel err vs dense = 1.047e-01
    rvnosig constant sigma 1.7                     ACCEPTED, rel err vs dense = 3.040e-01
    rvrows sigma ~ x                               ACCEPTED, rel err vs dense = 4.357e-01
    rvrows constant sigma 1.7                      ACCEPTED, rel err vs dense = 2.544e-16
    rvnorho rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 2.322e-01
    rvtrans rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 8.844e-02
    rvright rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 1.391e-15

The check is now `autocor_seam_check_jacobian()` (R/autocor-seam.R).
It requires `K (y - m) = whiten(y)$u` for two directions of `y - m`
(`sin(i)` and `cos(1.7 i)`), at two points. The first is the start
with sigma = 1. The second is a generic point: sigma 1.37, or
`exp(0.3 sin(i))` on row `i` under `sigma_varies = TRUE`, and every
parameter moved from the start by `0.4 k`, with k its position (0.04,
and then 0, when the map cannot be evaluated there). The same script
after the fix (`dev/loglik-log/rv-seam-after.txt`; each refusal reads
"The '<name>' term's map is wrong: jacobian() must give the n x n
matrix K with K (y - m) = whiten(y)$u, m the mean unwhiten(0), at every
sigma and parameter value; ..."):

    rvnosig sigma ~ x                              REFUSED
    rvnosig constant sigma 1.7                     REFUSED
    rvrows sigma ~ x                               REFUSED
    rvrows constant sigma 1.7                      REFUSED
    rvnorho rho = tanh(0.8)                        REFUSED
    rvtrans rho = tanh(0.8)                        REFUSED
    rvright rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 1.391e-15

"rvrows constant sigma 1.7" is refused too: the term is registered
with `sigma_varies = TRUE`, and the slot is wrong for that whatever
sigma the model has.

Why one generic point is enough: a slot passes only if its error
matrix `K_slot - K` annihilates both directions at that point. The
slips that matter are functions of sigma, of the parameters, or of the
orientation of `K`:

- sigma: caught by a sigma that is not 1 and differs by row, which
  also tells rows from columns;
- the parameters: caught by moving every parameter off the start;
- the orientation: a transpose differs on two generic directions
  unless `K` is symmetric, and then `K' = K` and the values are right.

A slot that is right at both points and wrong elsewhere (a branch on
the parameter value) is still possible. No finite set of points
excludes it, and such a slot is not a function of the map. The help
section "The Jacobian" states the points. The help example now says
that `Linv / sigma` holds for one sigma, and that a sigma per row
divides the columns.

New test block in `tests/testthat/test-autocor-loglik.R`, "a jacobian
slot wrong away from the start point is refused". It refuses the four
slips by name: no division by sigma; rows instead of columns under
`sigma_varies = TRUE`; rho ignored; the transpose on an asymmetric
row-standardized W. As the absent cases, it accepts the right slot of
each map and checks it against the dense conditional; fcor()'s and
sar()'s own slots pass in every fit.

### 2. The 8 SAR rows: class "adapted-input"

`dev/brmsport-bin2-adapters-spec.R` gains `bin2_adapted_input`: the 8
ids (log_lik:150, 152, 157, 159; posterior_predict:111, 113, 118, 120),
each with the reason "pass on an adapted input: brms's prep uses Msar =
diag(10), and the port replaces it by the row-standardized path graph
over the same 10 rows (brms_adapt_sar_path()). frmtmb refuses brms's
literal input by design: I - rho I makes rho and the scale one
parameter, so sar() needs a weights matrix with a negative and a
positive eigenvalue. ...".

`dev/brmsport-bin2-ledger.R` gives a held row on that list verdict
"pass", class "adapted-input" and that reason. The reason goes into
the verdicts it writes, so the generated tier file carries it. The
script stops on an entry whose row does not hold, and the summary
states how many passes are adapted. For a dry run, the script now
takes `FRMTMB_BIN2_REC_DIR` and `FRMTMB_BIN2_LEDGER_OUT`.

Dry run: all six bin-2 files were recorded on the lane build
(`dev/loglik-bin2-record.sh`, records in `dev/loglik-bin2-log/`), and
the ledger was written into `dev/loglik-ledger-out/` (log
`ledger-run.txt`). `git status` shows no tracked verdict, ledger or
summary file changed. The 13 rows:

    id	outcome	class	kind	value_check
    log_lik:150	pass	adapted-input	shape	holds
    log_lik:152	pass	adapted-input	shape	holds
    log_lik:157	pass	adapted-input	shape	holds
    log_lik:159	pass	adapted-input	shape	holds
    log_lik:174	pass		shape	holds
    log_lik:176	pass		shape	holds
    log_lik:230	pass		value
    posterior_predict:111	pass	adapted-input	shape
    posterior_predict:113	pass	adapted-input	shape
    posterior_predict:118	pass	adapted-input	shape
    posterior_predict:120	pass	adapted-input	shape
    posterior_predict:133	pass		shape
    posterior_predict:135	pass		shape

The dry-run summary:

    | 2P | pass | - | - | 162 |
    | 2P | pass | adapted-input | - | 8 |
    8 of the bin-2 passes hold on an adapted input (class adapted-input: the port replaced brms's input, which frmtmb refuses by design); 220 are literal.
    | pass | 392 | 44 | 184 | 620 |

The tracked summary has `| pass | 392 | 44 | 171 | 607 |`. The gain is
13 rows: 5 literal passes (the fcor rows log_lik:174, 176, 230 and
posterior_predict:133, 135) and 8 adapted. Every other row of the dry
run matches the tracked verdicts; the ledger stops on any difference.

### 3. The dense fallback warns above 1000 rows

Measured (`dev/loglik-fallback-cost.R`, seed 7,
`dev/loglik-fallback-cost.txt`; sar()'s lag map registered without its
slot, over a path graph; minimum of 5 calls, peak memory of one call):

    n =   500  call 0.030 s  peak Mb 21  4000 draws 0.03 h  silent
    n =  1000  call 0.500 s  peak Mb 48  4000 draws 0.56 h  silent
    n =  1500  call 0.720 s  peak Mb 105  4000 draws 0.80 h  warns
    n =  2000  call 1.670 s  peak Mb 155  4000 draws 1.86 h  warns
    n =  3000  call 4.950 s  peak Mb 345  4000 draws 5.50 h  warns

The reviewer's map measured 4.3 s and 306 Mb at 2000; its whiten and
unwhiten are dense too.

I chose a warning, not a refusal: the values are exact, and 345 Mb per
draw at 3000 rows is not a memory hazard on a machine that fits such a
model. Above `autocor_loo_dense_max` = 1000 rows, core raises
`frmtmb_autocor_dense_fallback`, which names the term and `jacobian =`.
frmtmb.sample's `log_lik()` passes on the first one per call and
muffles the rest. Tests:

- core, "a term without a jacobian slot warns above 1000 rows": 1001
  rows warn, 1000 do not;
- frmtmb.sample, "log_lik() says once per call that a term takes the
  dense route": 3 draws at 1001 rows give exactly 1 warning.

### 4. A test that sar() uses its slot

test-sar.R, "sar() reads its precision off its sparse jacobian slot,
not the dense fallback": `sar_jacobian()` returns a `sparseMatrix`; a
1024-row sar() raises no fallback warning; and the same map registered
without the slot warns. Seen to fail on a mutant
(`dev/loglik-mutant-sarslot.R`, `dev/loglik-mutant-sarslot.txt`): sar()
registered again without the slot, with test-sar.R run against it:

    slot after the mutant registration: TRUE (TRUE = no slot)
    RESULT pass=91 fail=1 err=0
    FAILING:
      sar() reads its precision off its sparse jacobian slot, not the
               dense fallback

### 5. The gap difference in `?frmtmb.sample::log_lik`

One paragraph added to "Residual covariance terms". Under `ar()`,
`ma()` and `arma()`, a group that skips a time point inside its range
differs from brms: frmtmb reads lags by time level and brms by
position. On one such group in the reviewer's case the difference is
17 percent (`dev/review-loglik-log/brms.txt`, "rel = 1.721e-01"),
against 7.8e-16 on the groups without a gap. `cosy()` and `unstr()`
agree.

### 6. The cross-reference

sample-loo's link now points to `[sample-log_lik]`.

### 7. The loo_*() functions against the loo package

test-loglik-autocor.R, "loo(), waic() and the loo_*() functions on a
covariance block are the loo package's on the same log_lik() matrix".
For `ar(cov = TRUE)` and `cosy()`, with 40 hand-set draws each:

- `loo()` equals `loo::loo(ll, r_eff = relative_eff(exp(ll), one
  chain))`, in the estimates and the pointwise elpd;
- `waic()` equals `loo::waic(ll)`;
- `loo_epred()` equals `loo::E_loo(posterior_epred(ds), psis(-ll))`,
  and also the weighted mean `colSums(w * posterior_epred(ds))` written
  out;
- `loo_predict()` equals `E_loo()` of `posterior_predict()` drawn from
  the same random state.

There is no brms reference for these. brms's values would need brms's
draws through brms's `loo_*()`, which read the same `log_lik()` (A3
agrees with brms to 2.0e-15) through the same loo package.

### 8. test-loo.R

The `ar(cov = TRUE)` arm is now its own block, "log_lik() of an ar(cov
= TRUE) block is the conditional of each row given the rest of its
subject", at the end of the section "log_lik() is the conditional
density, exactly". The refusal block keeps the mixture arm.

### 9. theta_map

`autocor_seam_loo_parts()` applies `autocor_seam_theta_on()` once, for
`whiten()` and for the slot. The fallback gets the raw parameters,
because `autocor_seam_B()` applies the map itself.

### Not done

- Under `ar()`, `ma()` or `arma()`, `frm()` gives no message when a
  group skips a time level inside its range. `autocor_warn_gaps()`
  warns only on non-consecutive GLOBAL levels, so a ported brms script
  fits the by-level model where brms fits the by-position one. The
  defect predates this lane; filed here, not fixed.
- brms accepts `se()` and `mi()` beside `ar(cov = TRUE)`, where frmtmb
  refuses them (review, "What I checked"). This also predates the lane.

### The reruns

The first full gated run of fix round 1 (`dev/release/loglik-fr1-gated`,
"LOGLIK-FR1-GATED ran 452 of 452") found 2 failures in core's
test-message-uniqueness.R: "duplicated error templates: The '' term's
map  for". The new `autocor_seam_check_jacobian()` defined a second copy
of the map check's refusal. It now receives the map check's own `bad()`.
Core was reinstalled, and its 280 files and its three sessions were run
again. The extension results below come from the first run; the fix
changes no code that they reach other than core's message wiring.

One file per process, gated (`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`;
core from `dev/release/loglik-fr1b-gated`, "LOGLIK-FR1B-GATED ran 280
of 280"; the extensions from `loglik-fr1-gated`):

    frmtmb           files 280 pass 22257 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling  files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam       files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent    files 12 pass 681 fail 0 err 0 skip 2 warn 0
    frmtmb.learn     files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode       files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample    files 68 pass 3434 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial   files 8 pass 547 fail 0 err 0 skip 0 warn 0
    frmtmb.survival  files 15 pass 552 fail 0 err 0 skip 0 warn 0

The skips are the other tiers' gates, as before.

One session per suite, three orders, `NOT_CRAN=true` (core:
`loglik-fr1b-sess`; sample and spatial: `loglik-fr1-sess`). No LEAK
line is untagged; 20 and 37 LEAK lines are all "[with load of ...]":

    SESSION frmtmb forward files=270/280 pass=19821 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=270/280 pass=19971 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=270/280 pass=20064 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3136 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3136 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3136 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=464 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=464 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=464 fail=0 err=0 skip=6 warn=0

OpenBLAS 0.3.26 with LAPACK ("Rblas: 74 to OpenBLAS, 2 to Rblasref",
"Rlapack: 616 to OpenBLAS, 10 to Rlapackref"; `OPENBLAS_NUM_THREADS=4`;
gated; logs `dev/loglik-log/ob-fr1-*.txt`; the copy is deleted again):

    RESULT frmtmb test-autocor-loglik.R pass=61 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-fcor.R pass=57 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=104 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=146 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=112 fail=0 err=0 skip=0 warn=0

The two bin-2 files on the final build (`dev/loglik-log/fr1-run-*.txt`,
records `fr1-rec-*.tsv`):

    RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=97 fail=0 err=0 skip=0 blocks=15
    RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=61 fail=0 err=0 skip=0 blocks=14

All 13 rows hold (held = TRUE), and the 6 value checks hold at 3.1e-16
relative or better, as in Part B.

R CMD check --as-cran (`dev/loglik-log/fr1-check-*.txt`): frmtmb,
frmtmb.sample and frmtmb.spatial each "Status: 1 NOTE". The NOTE is
"Skipping checking math rendering: package 'V8' unavailable".
`dev/loglik-check/` is deleted.

## Fix round 2

The reviewer's "Fix round 1 check" (`dev/reviews/2026-10-27-loglik.md`)
says SHIP, with three gaps. All three are fixed.

### 1. The generic point steps parameters both ways

The gap: the generic point moved the parameters upward only (0.4, then
0.04, then 0). The reviewer's `rvoneside` has rho = theta on (-1, 0],
starts at 0, its map stops for theta > 0, and its slot ignores rho. It
was accepted, and `log_lik()` was 3.465e-01 off at rho = -0.5
(`dev/review-loglik-log/fr1-seam2.txt`).

New test blocks in `tests/testthat/test-autocor-loglik.R`:

- "the generic point also steps a parameter down from its start": the
  reviewer's construction with the slot that ignores rho, which must be
  refused. The absent case is the same map with the right slot, which
  is accepted and equals the dense conditional at rho = -0.5.
- "a slot whose map cannot move off its start is checked there only,
  and says so": a map defined at its start alone, with a right slot.

Seen to fail on the build before this round
(`dev/loglik-log/fr2-before-autocor-loglik.txt`):

    RESULT frmtmb test-autocor-loglik.R pass=62 fail=2 err=1 skip=0 warn=0
    ── 1. Failure ('test-autocor-loglik.R:487:3'): the generic point also steps a pa
    Expected `ll_unfit(y ~ x + oneside(W), data = d, data2 = list(W = W))` to throw a error with class <frmtmb_error>.
    ── 2. Failure ('test-autocor-loglik.R:543:3'): a slot whose map cannot move off
    Expected `w` to inherit from "frmtmb_autocor_jacobian_unchecked".

The first failure shows the slip accepted; the absent case passed on
that build too.

The fix, in `autocor_seam_check_jacobian()`:

- For a term with parameters, the steps are 0.4, -0.4, 0.04 and -0.04
  times each parameter's position, in that order.
- If no moved point can be evaluated, the start is checked with the
  generic sigma, and a warning of class
  `frmtmb_autocor_jacobian_unchecked` names the term. It says that the
  slot was checked only at the start of its parameters, that a slot
  which ignores a parameter would pass, and that `log_lik()` reads its
  precision off that slot alone.
- It warns and does not refuse, because the slot may be right: refusing
  a correct model is worse than a check that says what it could not do.

The help section "The Jacobian" and the function's own doc state the
steps and the warning.

The reviewer's script `dev/review-loglik-fr1-seam2.R` (seed 101),
rerun unedited on the fixed build (`dev/loglik-log/fr2-rv-seam2.txt`):
every case got its expected outcome, and `rvoneside` is now refused.

    rvoneside (ignores rho) asymmetric W         expected REFUSED  got REFUSED  The 'rvoneside' term's map is wrong: jacobia
    rvok asymmetric W                            expected ACCEPTED got ACCEPTED rel err vs dense = 1.459e-16
    rvsym symmetric W                            expected ACCEPTED got ACCEPTED rel err vs dense = 2.927e-16
    rvtsym symmetric W                           expected ACCEPTED got ACCEPTED rel err vs dense = 2.927e-16

The other slips (`rvsign`, `rvcache`, `rvfirst`, `rvsym` and `rvtsym` on
an asymmetric W, `rvfrozen0`, `rvfrozen8` and `rvfrozen14`) stay
refused.

### 2. The gap paragraph

`?frmtmb.sample::log_lik`, "Residual covariance terms", now says what
was measured. The case has one such group, 22 rows and 3 draws. The
largest difference over the matrix is 17 percent of brms's largest
value, and the rows of the groups without a gap agree to 8e-16
(`dev/review-loglik-log/brms.txt`, "rel = 1.721e-01" and "rel on
groups with no interior gap = 7.844e-16").

### 3. The sar-uses-its-slot test reads the threshold

test-sar.R builds the lattice as `k = ceiling(sqrt(lim + 1))`, where
`lim` is `frmtmb:::autocor_loo_dense_max`. It asserts `n > lim` before
it asserts that no warning is raised, so a change to the threshold
cannot make the block pass on a case below it.

### Reruns

The changed files, one per process, in both tiers (logs
`dev/loglik-log/fr2-{plain,gated}-*.txt`). The plain tier's skips are
the brms-gated blocks.

    plain  RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    plain  RESULT frmtmb test-autocor-loglik.R pass=66 fail=0 err=0 skip=0 warn=0
    plain  RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
    plain  RESULT frmtmb.sample test-loglik-autocor.R pass=76 fail=0 err=0 skip=2 warn=0
    plain  RESULT frmtmb.spatial test-sar.R pass=93 fail=0 err=0 skip=4 warn=0
    gated  RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb test-autocor-loglik.R pass=66 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb.sample test-loglik-autocor.R pass=104 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb.spatial test-sar.R pass=113 fail=0 err=0 skip=0 warn=0

Core one-session runs (`NOT_CRAN=true`, `dev/release/loglik-fr2-sess`,
"LOGLIK-FR2-SESS ran 3 of 3"). There are 20 LEAK lines, all of them
tagged:

    SESSION frmtmb forward files=270/280 pass=19826 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=270/280 pass=19976 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=270/280 pass=20069 fail=0 err=0 skip=206 warn=0

OpenBLAS 0.3.26 with LAPACK ("Rblas: 74 to OpenBLAS", "Rlapack: 616 to
OpenBLAS"; gated; `dev/loglik-log/ob-fr2-*.txt`; the copy is deleted):

    RESULT frmtmb test-autocor-loglik.R pass=66 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=104 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=113 fail=0 err=0 skip=0 warn=0

R CMD check --as-cran (`dev/loglik-log/fr2-check-*.txt`): frmtmb.sample
and frmtmb.spatial "Status: 1 NOTE"; frmtmb "Status: 1 NOTE". Each NOTE
is the V8 math-rendering note. `dev/loglik-check/` is deleted.
