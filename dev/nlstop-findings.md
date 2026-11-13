# Lane nlstop: nonlinear and ridge fits that stop short in silence

Round opened 2026-11-11, base e8abed0a (frmtmb 0.81.0). Worktree
`frmtmb-wt-nlstop`, branch `wt-r82-nlstop`, library `wt-nlstop-lib`, base
library `rellib-r20`. Status: DONE after the fix round, nothing committed. Numbers are
pasted from the logs named beside them (`dev/nlstop-out/`, gitignored).
Before-arm runs set `NLSTOP_ARM=before` (rellib-r20 alone); after-arm runs
put `wt-nlstop-lib` first. Every log's first line names its library.

## Items, in the brief's order

| line | item | verdict |
|---|---|---|
| 4850 | `y ~ x` at 1e6 with `normal(0, 1000)` on the slope takes the lower mode | FIXED: a restart at the prior's location |
| 4960 | a sign-probe win returns a quieter lower maximum (bernoulli mm 55) | MISDIAGNOSED in part; its silence FIXED by a warning |
| 2584 | F11 seed 8 stops 0.1418 short at code 0 | NOT FIXED: needs a refit per sd; filed |
| 4705 | two weak-prior ridge stops lose the gradient warning | NOT FIXED: cause found to differ; filed |
| 2569 | huber() discontinuous objective | FILED with measurements |
| 4916 | `nl_local_message()` reads `stats::logLik()` | FIXED (NIT) |
| 4902 | sig4 bernoulli seed 55 refit | rerun: 2.8e-4 short, below grad_tol; the mm part is 4960 |
| 4842, 4912, 4793, 4856, 4931, 4897, 5091 | not silent | not reached; 4931 left because `frmtmb_control()`'s signature is shared |

## What changed

New file `R/fit-restart.R`; hunks in `R/fit.R` (one call site) and
`R/par-template.R` (nl_fit_probed(), nl_local_message()).

1. **prior_conflict_escape()** (4850). After the first descent of a
   Laplace fit (and gp_second_descent(), nl_basis_descent()), from the
   default start only: every fixed coefficient whose estimate lies more
   than `prior_conflict_z` = 4 scales of its normal, Student t or
   logistic prior from the prior's location (entries on the internal
   scale, no link, no offset, so a centered-intercept prior is not read;
   nonlinear coefficients excluded, since prior_nl_starts() already
   starts them there) is set to the location, the rest at the tape's own
   start, and the optimizer runs once more. The better optimum is kept
   when it gains more than nl_probe_gain(); on the user's call a message
   of class `frmtmb_prior_multistart` names the coefficients and the two
   penalized values. Skipped for a custom optimizer, a non-nlminb
   optimizer, an early-stopping optCtrl (iter.max = 0 included), a
   whitened beta, REML and profile (beta not outer). The restart runs
   only on conflict, so a fit without one pays nothing.
2. **nl_limit_message()** (4960, 4902). In nl_fit_probed(), on the fit
   frm() returns, when that fit raised no warning: the Hessian of the
   objective over the body's nonzero coefficients by central differences
   of the gradient (2 k gradients), taken to log magnitudes; along its
   flattest direction a move by `nl_limit_factor` = 1000 in the
   coefficient that moves most, both ways (2 objective evaluations). A
   side that makes some coefficient larger and lowers the log-likelihood
   by no more than grad_tol gives a warning of class `frmtmb_nl_limit`.
   The growth condition keeps a coefficient near 0 out when the move
   shrinks it. A second condition was added after the full suite: every
   coefficient the move carries must be at least one conditional
   standard error (1 / sqrt(H_ii)) from 0. The build before it warned on
   test-nl-start.R's "a prior's curvature enters a coefficient's unit"
   (`a * xs`, xs about 1e-20, `normal(0, 1e-3)` on a: the prior holds a
   near 0, a 1000-fold move changes nothing, and that is the prior's
   answer, not a ridge). On mm 55 Vm and K are each 1.8 conditional
   standard errors from 0 (Vm 12734.6 x sqrt(2.04e-8), K 177216 x
   sqrt(1.06e-10), the Hessian of `dev/nlstop-dbg3.R`).
3. **nl_local_message()** (4916) prints logLik_penalized() with the word
   "penalized" under a prior with a density, as nl_saddle_warning() does.

Tried and removed (not in the tree): `sd_down_escape()` for 2584 and a
prior gain along Ritz-flat directions in headroom_off_flat() for 4705;
the sections below give the measurements that removed them.

## 4850: a prior in conflict with the data

Reproduced on rellib-r20 (`dev/nlstop-prior.R`, the 63 rows of
`dev/review-nlstart-prior.R` with the penalized log-likelihood; log
`prior-before.txt`): the three `apbx 1e6 1000` rows read linear
-2932.62887486, -2938.69597870, -2923.87972866 against the body
spelling's -2462.95723096, -2466.21957741, -2457.55671429, each at code 0
with no warning. No other row differs by more than 1e-4.

Cause (`dev/nlstop-4850-probe.R`, log `4850-probe-before.txt`, seed
6201): two posterior modes. The fit's slope 4.9789688e+04 with log sigma
6.8786 (sigma 971), objective 2921.577; from either the tape's start or
the fit with the slope set to 0, the optimizer reaches objective
2463.231227 at slope 4.38e+03, log sigma 10.736 (sigma 4.6e4). A normal
likelihood with a free scale and a normal prior has a second mode only
when the estimate lies far beyond the prior's scale; a log-concave
likelihood under a normal prior has one.

After (`prior-after2.txt`): 63 rows, 0 with |gap| > 1e-4; 3 messages,
for the three rows above. The 20 ulp draws of seed 6201's response
(`dev/nlstop-4850-ulp.R`, `4850-ulp-after.txt`): every draw at
-2463.231227, the message on 20 of 20.

## 4960 and 4902: the "quieter lower maximum" is a linear limit

`dev/nlstop-row.R fresh bernoulli mm 55` and `dev/nlstop-mm55.R` on
rellib-r20 (logs `row-bern-mm-55-before.txt`, `mm55-before.txt`):

- The returned fit: -206.2711973, code 0, no warning, Vm = 12734.61, K =
  177215.71, SEs 2155044 and 29990076.
- It is the linear limit: `y ~ 0 + xp` reads -206.2711798 with slope
  0.0718582; Vm / K at the fit is 0.07185935. Scaling Vm and K together
  by s moves the log-likelihood by -1.746e-05 (s = 0.5), 8.730e-06 (2),
  1.571e-05 (10), 1.728e-05 (100), 1.746e-05 (1e4).
- The "higher maximum" of the reference (-204.429415, start 13) ends at
  K = -5 exactly, code 1 with two warnings: xp runs from 0.126 to 5, one
  row has xp = 5 and y = 1, and the body has a pole there. Not a maximum.

So the 1.84 shortfall is against a pole, and the defect is that the fit at
the limit says nothing. After (`dev/nlstop-limit.R`, `limit1.txt`):

- The row warns with `frmtmb_nl_limit`. Vm, K times 1000: +1.744e-05;
  divided by 1000: -1.767e-02.
- 20 ulp draws of xp (`xp * (1 + 2e-15 k)`): 17 of 20 give the limit
  warning at -206.271197; 3 (k = 1, 2, 4) stop at -206.277640 further
  along, where the SE check's own warning (frmtmb_warning) is given and
  the limit check therefore does not run. Every draw warns.
- Absent cases, no limit warning on any: `b1 * exp(b2 * x)` with true b2
  = 0 (seeds 1 to 5, b2 from -0.0106 to 0.00968); gaussian and poisson mm
  rows, seeds 11 to 55 (no warning of any class); the 450-row scan's sat
  data, seeds 1 to 5.

4902's second case, bernoulli `sig4` seed 55 of the fresh2 rows: -189.6964,
2.8e-4 below the best start, no warning, on both builds. That is below
grad_tol (1e-3), the package's convergence tolerance; left.

## The 450 rows, the fresh rows, the bigx grid

Driver `dev/nlstop-queue.sh` (job lists `jobs-before.txt`,
`jobs-after2.txt`, `jobs-after4.txt`). The after arm of the scan runs no
reference starts (`NLSTOP_REF=0`); the reference is the before arm's.
`after2` is the build before the conditional-SE gate, `after4` the final
build (`scansum-after4.txt`, `freshsum-after4.txt`, `prior-after4.txt`);
every count below reads the same on both, the limit warning on the same
two rows. The bigx grid ran on after2 only: it fits no nonlinear body,
so the gate cannot reach it.

450 rows (`dev/nlstop-scan.R`, `dev/nlstart-scansum.R`, log
`scansum-after2.txt`):

```
scan-before: rows 450, default NA 0, short > 1e-4: 0 (silent 0, warned 0), short > 1: 0 (silent 0)
scan-after2: rows 450, default NA 0, short > 1e-4: 0 (silent 0, warned 0), short > 1: 0 (silent 0)
scan-after2 minus scan-before default logLik: higher by > 1e-4: 0  lower by > 1e-4: 0  within 1e-4: 450
```

Warned rows: 99 before, 100 after. The new one is bernoulli seed 1 slope
0.6 `a + b1 * b2 * x`, -173.776165, `frmtmb_nl_limit`; the other 89 prod
rows warned already (only the product b1 b2 is determined).

Fresh rows (`dev/nlstop-fresh.R` and `-fresh2.R`, the review's 425 rows,
335 fit from the default start; the 90 rate-free sigmoid rows of the first
set fail at the default start on both builds, as in the record; summary
`dev/nlstop-freshsum.R`, log `freshsum-after2.txt`):

```
short > 1e-4: after 11 (silent 4), before 11 (silent 5)
short > 1: after 6 (silent 2), before 6 (silent 3)
after - before: higher > 1e-4 0, lower < -1e-4 0, within 335
multistart message on 87 rows; of them gain over probe-off < 1e-2: 0
```

The silent row that now warns is bernoulli mm 55. The 4 still silent:
binomial mm 11 (4.32) and 44 (0.92), whose fits stop with K between -0.2
and -0.1, a pole inside xp; cumulative pow 22 (1.03); bernoulli sig4 55
(2.8e-4). Filed. `frmtmb_nl_limit` appears on 1 of the 425 rows.

bigx grid (`dev/nlstop-bigx-grid.R`, `dev/bigx-summary.R`,
`dev/nlstop-bigxcmp.R`): 240 of 240 right on both builds, 0 warned;
cells whose dll, slope_z, code or warning count differ: 0.

No row of any set is worse.

## 2584: F11 seed 8, not fixed

`dev/nlstop-f11.R 8` (logs `f11-before.txt`, `f11-after1.txt`): default
-589.023726, code 0, no warning, theta (-7.99197, -0.79877); from theta_1
= -10 the same; from theta_2 = -10 (the sd of s(x2), in sigma's formula)
-588.881888, where theta_2 stays at -10. Gap 0.141838. The boundary is
higher than an interior local maximum.

Tried: an escape like sd_up_escape() downward, one objective evaluation
per interior sd at its log sd - 10 with the rest held. It does not see it
(`dev/nlstop-dbg1.R`): theta_2 - 10 reads 589.7186 against the fit's
589.0237, 0.69 lower. One and two Newton steps in the other parameters on
the fit's outer Hessian (`dev/nlstop-dbg2.R`) diverge: 712.8 and 2.7e11 at
-10, the same at -3 and -5. Only a refit from there reaches the higher
value, one more optimization per interior sd on every model with random
effects. The brief prefers a warning to such a search, and no cheap
warning exists; removed and filed.

## 4705: the cause is the quadratic model, not the Ritz-flat exception

`dev/nlstop-ridgeprior.R` (both review scripts' rows on one table; log
`ridgeprior-before.txt`), code-0 stops:

```
psd 100, k 5:   1.012e-3 short, no warning
psd 1000, k 20: 1.9e-4 short, no warning (below grad_tol: the rule)
psd 1000, k 50: 1.225e-3 short, no warning
psd 1000, k 80: 3.16e-3 short, gradient warning
```

On the stop psd 1000, k 50 (`dev/nlstop-dbg3.R`): grad_verdict() reads
gmax 0.0732 and headroom 4.23e-4, below grad_tol, and
headroom_ritz_flat() returns no column. So the reparameterized Newton
step under-reads the gain along the curved ridge by a factor of 2.9, and
the Ritz-flat exception is not reached. The prior's quadratic gain along
the ridge's tangent is the shortfall exactly ((k - 0.5)^2 / (2 psd^2)),
but the data-only Hessian's flattest direction carries the gradient
across the ridge: a gain read along it was 2.7e7 (`dev/nlstop-dbg4.R`).
A prior gain added in headroom_off_flat() does not reach these stops.
Removed and filed. dbg3 and dbg4 ran on the first lane build, whose
headroom_off_flat() had that gain; the path they read does not reach it.

## 2569: huber()

`dev/nlstop-huber.R` on rellib-r20 (`huber-before.txt`): `y ~ x + (1 |
g)`, 100 rows, 20 groups, t3 noise, seeds 1 to 20, against the best of 8
jittered starts:

```
re TRUE  scale 1:    rows 20, short 11, silent 6, max short 0.188
re TRUE  scale 1e3:  rows 20, short 11, silent 5, max short 0.272
re FALSE scale 1:    rows 20, short 0
re FALSE scale 1e3:  rows 20, short 0
```

Silent shorts: 0.0074 to 0.272. With `quadrature = TRUE`
(`dev/nlstop-huber-quad.R`, `huber-quad-before.txt`): scale 1, 0 short;
scale 1e3, 2 short (0.0274 and 6.7e-4), both silent. The Laplace
objective jumps where a residual crosses `k` (lane optrap, m6). Filed:
say so on a huber() fit with random effects under Laplace, or document
`quadrature = TRUE` as its route. Not changed, because a warning on every
such fit would fire on correct fits too, and the documentation is in
`R/families.R`, which this lane does not own.

## Cost

`dev/nlstop-cost.R`, one process, arms rotated per round, min of 3
rounds; "off" and "off2" replace nl_limit_message() by NULL and
prior_conflict_escape() by the identity (0.81.0's path), and off2 / off is
the control. nl: the scan's 45 bodies at seed 101; prior: the prior rows'
21 linear fits at seed 6201, 8 times.

```
cost1: nl: min s off 18.41 on 16.44 off2 18.53 | on/off 0.893 control 1.007
cost2: nl: min s off 21.60 on 21.28 off2 19.45 | on/off 0.985 control 0.900
cost2: prior: min s off 3.44 on 2.92 off2 3.01 | on/off 0.849 control 0.875
```

No cost is measurable against a control that moves by 10 to 13 percent.
By construction the limit check is 2 k gradients and 2 objective values
per nonlinear fit that nothing warned about, and the prior restart is one
optimization on a fit with a conflict and nothing otherwise. That last
sentence was wrong where it mattered: the review measured the restart at
2.1 times a logistic GLMM's fit and 2.5 times refit(), because these
rows were all gaussian `y ~ x`, whose restart is cheap. See "Fix round".

The scans' own wall times (`scansum-*.txt`: default fits 97.4 s before,
197.7 s after2, 311.8 s after4, a ratio of 2.0 and 3.2) are not a
measurement of the change. Each arm ran as separate processes at
different times beside eleven other lanes; the before arm ran 2
processes at a time, the after arms 3. The direct timings in one process
show no cost: the nl block above (on/off 0.893 and 0.985 against
controls of 1.007 and 0.900), the nonlinear test files ("Cost on the
nonlinear test files": 0.96 of base), and the review's own timing of the
check (`dev/reviews/nlstop-review-nlcost.R`: 0.000 to 0.010 s against
fits of 0.14 to 0.67 s). So the 3.2 is machine load.

## Tests

New `tests/testthat/test-fit-restart.R`: the 4850 rows (seeds 6201 and
6203) reach the body's value with the message; a prior the estimate
agrees with and a given start take no restart; a conflicting prior under
a poisson likelihood keeps its one mode with no message; the mm 55 row
warns `frmtmb_nl_limit`, sits at the linear limit to grad_tol, and every
one of 4 ulp draws warns; an identified mm body and a coefficient near 0
give no limit warning.

- rellib-r20 (`t-fit-restart-before2.txt`): pass 16, fail 8: the 4850
  test 6 (ll -2921.6 and -2938.7 against -2463.2 and -2457.8; slope
  49789.7 and 49718.0; no message), the mm test 2 (no limit warning; a
  silent ulp draw).
- lane build (`t-fit-restart-after2.txt`): pass 24, fail 0, err 0, warn 0;
  final build (`t-fit-restart-after4.txt`): the same.
- OpenBLAS 0.3.26 (LAPACK) failed 2 of that version: the mm 55 fit took
  the other stop, -206.27764, 0.0065 below the limit, where the SE check
  warns and the limit check does not run. The test was rewritten to
  accept either warned stop and to call nl_limit_message() on the fit
  directly (it finds the ridge from either stop: the growth side raises
  the log-likelihood by 0.0065 there). Final test file: reference BLAS
  pass 25 (`t-fit-restart-after5.txt`), OpenBLAS 0.3.26 pass 24 (the
  at-limit branch not taken), rellib-r20 pass 12, fail 7, error 1
  (`t-fit-restart-before5.txt`: the 4850 test 6, "every stop warns" 1,
  and nl_limit_message() absent 1).

## OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack` (no
denial), then `dev/nlstop-ob.sh <ver>` (18 files: test-fit-restart.R,
test-nl-start.R, test-nl.R, test-nl-family.R, test-nlf.R, test-map.R, the
10 test-prior-*.R and 2 test-priors-*.R; 4 threads; each log's blas: line
names the copy). Final test files on both copies: 18 of 18, fail 0, err
0, warn 0 (0.3.32: fit-restart 24, map 16, nl-family 333, nl-start 59, nl
92, nlf 78, prior-compat 202, prior-location-dpar 55, prior-mv-resp 36,
prior-one-slot 90, prior-only-objective 10, prior-sd-coef 37,
prior-sd-scope 25, prior-sds 52, prior-support 150, priors-autocor-classes
182, priors-bounds-grcov 53; 0.3.26 the same counts). `dev/ciharden-out/`
deleted.

## Suites

Plain tier, every core and extension file, one per process, 3 at a time
(`dev/release/run-par.sh nlstop-plain-suite`, `FRMTMB_RELLIB` = the
lane's library, extensions from rellib-r20; summary
`dev/nlstop-suitesum.R`, log `plain-suitesum1.txt`). It ran on the build
before the conditional-SE gate:

```
files 511 with RESULT 511 with the lane's core on the lib: line 511
             pkg files  pass fail err skip warn
          frmtmb   318 23510    1   0  239    0
 frmtmb.coupling    11   551    0   0    5    0
      frmtmb.eam    33  1848    0   0    5    0
   frmtmb.latent    12   687    0   0    2    0
    frmtmb.learn    16   474    0   0   13    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    87  4050    0   0   18    0
  frmtmb.spatial     8   478    0   0    6    0
 frmtmb.survival    15   556    0   0    0    0
total pass 32703 fail 1 err 0 skip 289 warn 0
```

The one failure: test-nl-start.R, "a prior's curvature enters a
coefficient's unit", a limit warning on a prior-held coefficient. The gate
was added; on the final build test-nl-start.R reads pass 59, fail 0
(`t-nl-start-after4.txt`). The gate only removes warnings, and no other
test names the limit warning (grep of every test file for
`nl_limit` and its text: test-fit-restart.R, and a comment in
test-nl-family.R), so the other 510 files stand.

Gated tier, the 55 files with a skip above, every gate on (NOT_CRAN,
FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS, FRMTMB_FUZZ,
FRMTMB_SCALE_TESTS), final build (`gated-suitesum.txt`): 55 of 55, every
lib: line the lane's core; pass 4655, fail 0, err 0, skip 0, warn 0
(core 38 files 3588; coupling 6, eam 31, latent 6, learn 85, ode 2,
sample 744, spatial 193).

## Cost on the nonlinear test files

`dev/nlstop-nlcost.sh` (`nlcost.txt`): each file on rellib-r20 and on the
final build, alternating, 3 rounds, one process at a time; minimum wall
seconds per file:

```
test-nl-start.R   base 14.96  lane 15.69  (1.05)
test-nl.R         base  6.36  lane  6.56  (1.03)
test-nl-family.R  base 54.20  lane 51.47  (0.95)
test-nlf.R        base  5.00  lane  4.72  (0.94)
total             base 80.52  lane 77.44  (0.96)
```

Within the rounds' own spread (test-nl-start.R base 14.96 to 20.62).

## Sessions

Not the 27 full sessions: the plain tier alone took several hours at 3
processes on the loaded machine, and one core session runs that suite's
work in one process. The lane adds one test file and no registry, option
or global, so the session run is on that file and test-nl-start.R, with
the package's helpers and setup (`dev/release/run-session.R frmtmb
<order> <dir>`, logs `session-forward.txt`, `session-preload.txt`):
forward and preload each `files=2/2 pass=83 fail=0 err=0 skip=0 warn=0`,
0 LEAK lines, every lib: line the lane's core. A full session run stays
for the consolidation.

## R CMD check --as-cran

Core, `dev/nlstop-check.sh` (R_LIBS lane, rellib-r20, user;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; pandoc and TinyTeX on PATH), on
the final tree: `Status: 1 NOTE`, the HTML manual's math rendering (V8
unavailable); tests 45 min, OK. Log kept as `dev/nlstop-out/check.log`;
the check tree is deleted. No extension changed, so no extension check.

Installed code equals the tree (`dev/nlstop-installed.R`): 133 objects of
R/fit-restart.R, R/fit.R and R/par-template.R, 0 differ.

## Bump and floor

Minor (two BREAKING behavior changes in fits). No extension changes and
no extension floor moves: the extensions' suites pass over this core with
their rellib-r20 builds.

## Not done, and why

- 2584 and 4705: measured, the tried fixes removed, filed (sections above).
- huber() (2569): filed; a warning on every such fit would fire on correct
  fits, and the documentation is in R/families.R.
- 4931 (`frmtmb_control()` switch for the probe): `frmtmb_control()`'s
  signature is in R/fit.R, which four other lanes edit; left to avoid a
  conflict on those lines.
- 4842, 4912, 4793, 4856, 4897, 5091: not reached.
- The 4705 debug reads (`dev/nlstop-dbg3.R`, `-dbg4.R`) ran on the first
  lane build, not on rellib-r20; the path they read (grad_headroom()'s
  return before headroom_checked()'s exception) is the same code on both.

## Files touched

- `R/fit-restart.R` (new): prior_conflict_z, prior_conflict_escape(),
  prior_conflict_message(), nl_limit_factor, nl_limit_message().
- `R/fit.R`: one hunk in fit_assembled() after nl_basis_descent() (the
  prior restart and its message). No other hunk; the two tried edits were
  reverted to the base text.
- `R/par-template.R`: nl_fit_probed() (the limit check at its two returns),
  nl_local_message() (penalized value).
- `tests/testthat/test-fit-restart.R` (new).
- `NEWS.md`: a `# frmtmb (development version)` section, three bullets.
- `dev/test-backlog.md`: three CLOSED lines (4850, 4916, 4960) and the
  section "Filed in the round of 2026-11-11", "From lane nlstop".
- `dev/nlstop-*` scripts and this record. Logs in `dev/nlstop-out/` and
  `dev/release/nlstop-*` (gitignored).

## Fix round (review dev/reviews/2026-11-11-nlstop.md, verdict FIX)

Each finding was reproduced on the pre-fix lane build first (logs in
`dev/nlstop-out/fix/`, `repro-*.txt`), with the reviewer's scripts.

### MAJOR 1: the prior restart where it cannot win

Reproduced (`repro-prior-a4.txt`, `repro-refit.txt`): the 20000-row
bernoulli GLMM, on 22.73 and 21.35 s against off 11.39 and 11.29 s, the
restart ran on every "on" fit and never won (llp -5586.801892 in every
arm); the gaussian LMM a4g on 1.90 and 1.73 against 1.58 and 1.41;
refit() of the 5000-row GLMM, on 20.45 and 20.19 s against 8.80 and
8.56 s, 3 calls of the escape per round of 3 refits.

Fixed, two gates:
- the likelihood must have a free dispersion in the outer vector (a
  `betad` entry: sigma, phi, shape, nu, ...); a bernoulli, binomial or
  poisson fit, whose likelihood has no scale to inflate, takes no
  restart;
- the user's own cold start only: the call site in fit_assembled() runs
  it when `announce_start` (frm()'s own call) and the template was not
  handed in by a caller (`cold_tpl <- is.null(template) ||
  prefit_seeded`, taken before make_start(), a hunk apart from
  ordstart's). refit(), influence(), frm_allfit(), anova()'s ML refit,
  the autoscale pre-fit and the autoscale choice's two fits never run it.
  This is ordstart's rule (`is.null(start) && is.null(template)`) with
  the pre-fit excluded as well.

Remeasured (`dev/nlstop-fix-cost.R`, `fix/fix-cost.txt`; one process,
arms off, on, off2 rotated, min of 3 rounds; off and off2 replace the
escape by the identity):

```
a4 (bernoulli GLMM, 20000 rows):  off 10.77 on 11.36 off2 10.64 | on/off 1.055 control 0.988
a4g (gaussian LMM, 20000 rows):   off  1.18 on  1.66 off2  1.21 | on/off 1.407 control 1.025
refit (bernoulli GLMM, 5000 x 3): off  7.75 on  7.72 off2  7.74 | on/off 0.996 control 0.999
```

The escape is called once on a4 and returns at the scale gate; on refit
it is not called (0 calls in every round); llp identical in every arm.
The gaussian LMM still pays one restart (1.4 times, 0.48 s here): it has
a free sigma, its conflict is real, and the restart did not win there.
That is the cost of the case the restart exists for; it is paid only on
a cold fit with a conflict.

The 4850 rows after the fix (`fix/prior-after5.txt`): 63 rows, 0 with
|gap| > 1e-4, 3 messages; the 20 ulp draws of seed 6201
(`fix/4850-ulp-after5.txt`): every draw at -2463.231227, the message on
20.

### MAJOR 2: the centered intercept

Reproduced (`dev/reviews/nlstop-review-intercept.R`,
`fix/repro-intercept.txt`): normal(0, 1e4) on class Intercept, default
-6668.727576 at code 0 with no warning or message, against -3055.084765
from the prior's location; brms 2.23.0's program from both points: lp
-6668.727576 and -3055.084765, rc 0 at each.

Fixed: prior_conflict_escape() reads an entry with an `offset` (the
centered intercept's prior is on b_Intercept + w' b): the centered value
`par[j] + sum(w * par[offset idx]) + shift` is compared with the
location, and the restart sets b_Intercept so that the centered value
is the location at the restart's start (after any plain coefficient in
conflict has moved). The message names "the centered Intercept".

After (`fix/intercept-after5.txt`, same script):

```
psd 1e3   default -3054.801403  (prior's mode from the default, unchanged)
psd 3e3   default -3055.739862  (unchanged)
psd 1e4   default -3055.084765  frmtmb_prior_multistart  (was -6668.727576)
psd 2e4   default -2915.698066  (data's mode, higher than -3049.084034)
psd 3e4   default -2220.942194  (data's mode, higher than -3035.359476)
psd 1e5   default -1716.063314  (unchanged)
```

brms's optimizer from the fixed fit's point reads lp -3055.084765, rc
0. The test asserts only the 1e4 case.

### MINOR 1, 3: filed

`dev/test-backlog.md`, section "Filed in the round of 2026-11-11", from
lane nlstop: the two limit-check false negatives (reproduced,
`fix/repro-fndebug.txt`: `a * (exp(b * x) - 1)`, a 1.40046e+04, b
1.43847e-04, side gains -923.8686 and -0.0199623 against 0.003285 along
(1000, 1/1000); Hill, side 1 NaN, gate 9.138484), and the corrected text
of binomial mm 11 and 44 (reproduced, `fix/repro-pole.txt`: -623.288529
at K -0.182532 against -618.966489 at K -6.0152; -626.898840 at K
-0.110348 against -625.980024 at K -6.45313) with the mm 55 user starts
at K = -2.5 (-207.650719 at K -1.00267, -207.001511 at K -2.91724, code
0, no warning).

### MINOR 2: the scan's 3.2 is machine load

Section "Cost" now says so, with the direct timings beside it.

### MINOR 4: the test file

`test-fit-restart.R` now goes through `allow_warnings()`: each fit's
messages are recorded and muffled, the warnings a test names are allowed
(the two ridge warnings in the mm test), and every other warning reaches
testthat. Two tests added: the centered intercept (psd 1e4 only), and
"the restart runs on a cold start with a free scale only" (a mocked
counter on prior_conflict_escape(): 0 calls on refit(), 1 on a cold
bernoulli fit, and with optimize_obj() mocked to stop, the escape
returns the bernoulli fit's optimum untouched).

- rellib-r20 (`fix/t-fit-restart-before6.txt`): pass 8, fail 9, error 2.
  Behavioral: the 4850 test 6, the intercept test 2 (-6668.7 against
  -3055.1, no message), the mm test 1 (no warning). The two errors are
  the missing functions (the gate test and nl_limit_message()); the gate
  test's behavioral failure on the pre-fix lane build is the
  reproduction above (3 escape calls per refit round, the restart run on
  the bernoulli GLMM).
- lane (`fix/t-fit-restart-after6.txt`): pass 29, fail 0, err 0, warn 0.

### MINOR 5: NEWS

The provenance sentence is gone. The first bullet now states the two
gates and the centered intercept; a Documentation bullet names the huber
note.

### huber()

`?huber` (the `huber()` details of frmtmb-families.Rd, R/families.R) now
says that with group-level terms the Laplace approximation changes in
steps where a residual crosses `k`, so a fit can stop on a lower piece at
code 0 with no warning, and that `quadrature = TRUE` (scalar group-level
terms) or several starts is the route. No warning. The numbers stay
here: 11 of 20 fits short at each response scale, 6 and 5 silent, 0.0073
to 0.272 (review `huber-1.txt`, `huber-1e3.txt`; lane `huber-before.txt`
0.0074 to 0.272); under quadrature 0 short at scale 1, 2 at 1e3
(`huber-quad-before.txt`). The backlog item stays open, with that note.

### Runs on the fix build

Core sessions on the review's 7 files (test-fit-restart.R, test-map.R,
test-nl-family.R, test-nl-start.R, test-nl.R, test-nlf.R,
test-prior-only-objective.R, with the helpers and setup.R;
`fix/session7-*.txt`): forward, reverse and preload each `files=7/7
pass=617 fail=0 err=0 skip=0 warn=0`, 0 LEAK lines, the lane's core on
every lib: line.

Full core suite, one file per process, 3 at a time
(`dev/release/run-par.sh nlstop-fix-core-suite`, `fix/core-suitesum.txt`):
318 of 318 with a RESULT line and the lane's core on the lib: line; pass
23516, fail 0, err 0, skip 239, warn 0.

Gated tier, the 55 files that skip in the plain tier, every gate on
(`fix/gated-suitesum.txt`): 55 of 55; pass 4655, fail 0, err 0, skip 0,
warn 0.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`, no
denial; `dev/nlstop-ob.sh 0.3.26`, the 18 nonlinear and prior files, each
blas: line naming the copy): 18 of 18, fail 0, err 0, warn 0
(test-fit-restart.R 28, the rest as before). `dev/ciharden-out/` deleted.

R CMD check --as-cran of core on the fix tree (`dev/nlstop-check.sh`, log
`fix/check.log`): `Status: 1 NOTE`, the HTML manual's math rendering (V8
unavailable); tests 38 min, OK. The check tree is deleted.

Installed code equals the tree (`dev/nlstop-installed.R`, now over
R/fit-restart.R, R/fit.R, R/par-template.R and R/families.R, sourced
against the namespace): 399 objects, 1 differs, `family_registry`, an
environment filled at load.

Not run in the fix round: the extension suites (the fix only removes
restarts and adds documentation, and the brief named the core suite),
the 450-row, fresh-row and bigx tables (no prior and no change to the
limit check reaches them; the scale gate only skips prior restarts, and
none of those rows has a prior), and OpenBLAS 0.3.32.

Files touched in the fix round: R/fit-restart.R (prior_conflict_escape():
the scale gate, the centered intercept; docs), R/fit.R (a second hunk,
`cold_tpl` before make_start(), and the call site's gate),
R/families.R and man/frmtmb-families.Rd (the huber() note; R/families.R
is outside the lane's area, a 6-line roxygen hunk), NEWS.md,
dev/test-backlog.md, tests/testthat/test-fit-restart.R,
dev/nlstop-fix-cost.R, dev/nlstop-installed.R, this record.

## Short check (SHIP), S2 and the S1 note

S2: the scale gate read any `betad` entry, so zero_inflated_poisson()
took the restart on its zi coefficient alone. The gate is now
free_dispersion(frame): a non-constant `betad` linear predictor whose
dpar (trailing digits dropped) is not a probability (zi, hu, zoi, coi,
theta<k>) and has a coefficient that `dpar = value` does not fix. Checked
(`fix/zi-probe.R`): zero_inflated_poisson `y ~ x` reads FALSE (its
linpreds are beta mu and betad zi), gaussian TRUE. NEWS names the
probabilities as not a dispersion. Test added: "a zero-inflation
probability is not a dispersion: no restart" (a zip fit with the slope
more than 4 prior scales from 0; with optimize_obj() mocked to stop,
prior_conflict_escape() returns the fit's optimum and runs 0 restarts).

S1, for the consolidation (not changed here, by instruction): ordstart's
`starts_cold = is.null(start) && is.null(template)` and this lane's
`cold_tpl <- is.null(template) || prefit_seeded` (with `announce_start`)
are two spellings of the cold-start rule and differ on the autoscale
pre-fit-seeded fit. And autoscale_choose() replays the plain arm's
messages whichever arm wins, so a `frmtmb_prior_multistart` message can
describe a fit frm() does not return. The consolidation unifies the
rule.

Runs after S2: test-fit-restart.R on the lane (`fix/t-fit-restart-after7.txt`)
pass 33, fail 0, err 0, warn 0; on rellib-r20 (`fix/t-fit-restart-before7.txt`)
pass 9, fail 9, err 3 (the new test errs there on the missing internals;
on the pre-S2 lane build its zip fit passes the old `betad` gate, so the
mocked optimize_obj() is reached and the test fails, by construction from
the zi-probe linpreds, not run). R CMD check --as-cran of core
(`fix/check2.log`): Status: 1 NOTE (V8 math rendering), tests 35 min OK.
