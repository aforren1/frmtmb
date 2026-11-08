# Lane nestquad: exact quadrature over two nested random-effect blocks

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-nestquad`, library `wt-nestquad-lib`, base library
`rellib-r18`.

Status: fix round 1 done; see "Fix round 1" at the end, which
wins over the sections before it.

## What changed

- New `R/quadrature-nested.R`: the nested route. `quad_nest_plan()`
  detects nesting from `Z`, `quad_nest_anchor()` reads the joint Laplace
  Gaussian, `quad_nest_objective()` writes the iterated Gauss-Hermite
  integral on the tape, `quad_nest_verify()` checks it against the
  model's own objective, `quad_nest_guard()` predicts memory, and
  `quad_nest_fit()` runs the anchor and re-anchor loop.
- `R/fit.R`: `quad_check_blocks()` returns the nested plan instead of
  refusing; `fit_assembled()` sends a plan to `quad_nest_fit()` and takes
  the conditional modes from the Laplace inner solve (as for
  `importance =`); `frmtmb_control()` gains `quad_nodes` (default
  `c(11, 15)`) and `quad_max_gb` (default 4); `frm()` docs.
- `R/compat.R`: the row `quadrature` x `group:quadrature_blocks` (text
  only, status stays "conditional").
- Tests: `test-quadrature-nested.R` rewritten; in
  `test-quadrature-defects.R` the nested refusal test became an exactness
  test and the non-Gaussian families loop runs `(1 | ga/gb)` again;
  `helper-fuzz.R` swaps the `quad-one-block-per-row` known refusal for
  `quad-nested-memory` (the memory guard on a censored Beta spec).
- NEWS: three bullets under `# frmtmb (development version)`.

## Design

The route that `dev/fuzzfix-nested-quad.patch` saved runs TMBad's
`marginal_gk` once per block; its tape grows by 12 to 551 MB per row
(`dev/fuzzfix-findings.md`, "To file"). This lane writes the iterated
integral on the tape itself, as the importance correction does.

For outer level a (effect u) and its inner levels c (effects v):

    m_a = int N(u) prod_c [int N(v) prod_{i in c} p_i dv] du

Both integrals use Gauss-Hermite rules centered and scaled by the joint
Laplace Gaussian at an anchor: u = uhat + s x_k, and given u,
v = vhat - r (u - uhat) + t y_j, with t^-2 = H_cc, r = H_ca / H_cc,
s^-2 = H_aa - sum_c H_ca^2 / H_cc. The nodes are constants on the tape,
so a row's predictor at a node is its fixed part plus a constant; every
row is repeated once per node pair and the family density is called on
the stacked vector (one block per inner node, the short dimension).
Re-anchoring at each optimum until the fresh tape is stationary
(gradient below `grad_tol`) makes the reported value that of a rule
anchored at the estimate.

Supported: exactly two blocks on a row, each a dimension-1 `us`, `diag`
or `homdiag` block, one nested in the other (every inner level inside
one outer level, read off `Z`), both on the same predictors, one
response, a rowwise family, no `nl` body or `cs()`, no `map` on `b`.
Refused by name: crossed blocks ("are crossed ... does not factor"),
three or more blocks, blocks on two predictors (`(1 | g)` in `mu` and
`sigma`), Student-t blocks, several responses, a structured family, a
nonlinear body, `cs()`, a mapped `b`, a multi-membership row, and a
predicted memory above `quad_max_gb`. A by-split `gr(g, by = f)` and
`mvbf()` with one block per response stay on the one-block
Gauss-Kronrod route, unchanged.

## Accuracy

References are written from the model alone (`dev/nestquad-ref.R`):
nested `stats::integrate()` with the integrand centered on its mode, and
non-adaptive Gauss-Hermite on the effects' own scale. The rule's value
is the nested objective re-anchored at the stated parameters
(`dev/nestquad-nodes.R`).

### Nodes (weighted Bernoulli, `dev/nestquad-nodes.R bern asym`)

Data: the test design, 12 ga x 2 gb x 3 rows, weights runif(0.5, 2),
seed 20261020; parameters at the fit's optimum. Reference converged:
`REFCONV bern gh80 -46.232110358425 gh120 -46.232110388498 diff
3.007e-08`. Log `dev/nestquad-nodes-bern.log`:

    REF bern ref -46.2321103885 laplace -47.0845448125 err_laplace 8.524e-01
    NODES bern k  1 err 8.524e-01
    NODES bern k  3 err 2.918e-01
    NODES bern k  5 err 4.389e-02
    NODES bern k  7 err 5.976e-03
    NODES bern k  9 err 8.142e-04
    NODES bern k 11 err 1.298e-04
    NODES bern k 15 err 8.528e-06
    NODES bern k 21 err 2.500e-07
    NODES bern k 31 err 6.421e-10
    ASYM bern ko 31 ki 11 err 1.264e-04
    ASYM bern ko 11 ki 31 err 3.289e-07
    ASYM bern ko 11 ki 15 err 6.968e-06
    ASYM bern ko 11 ki 21 err 1.062e-06
    ASYM bern ko 15 ki 21 err 2.835e-07
    ASYM bern ko  9 ki 15 err 1.155e-05

The inner level carries the error (ko 31 ki 11 is no better than 11 x
11; ko 11 ki 31 is 400 times better). The default `c(11, 15)` costs 165
stacked rows per row and gives 7.0e-6 here.

The unweighted fit of the same design (`dev/nestquad-proto.R`, 11 x 11):
`quad -41.6876706264 ref_int -41.6876657218 ref_gh80 -41.6876657218
lap_at -42.3300740252`, error 4.9e-6 against 0.642.

### Poisson, larger clusters (`dev/nestquad-nodes.R pois`)

8 ga x 3 gb x 8 rows, seed 20261107; reference by nested integrate()
(rel.tol 1e-9 and 1e-11 agree to the last digit). Log
`dev/nestquad-nodes-pois-linear.log`:

    REF pois ref -347.3535206657 laplace -347.3995995435 err_laplace 4.608e-02
    NODES pois k  5 err 3.841e-04
    NODES pois k  7 err 5.583e-06
    NODES pois k 11 err 7.510e-08
    NODES pois k 15 err 5.760e-10
    NODES pois k 21 err 1.251e-12

### Gaussian against its closed form (`dev/nestquad-nodes.R gaus`)

10 ga x 3 gb x 4 rows, seed 20261108; reference the multivariate normal
density with V = sd_a^2 Z_a Z_a' + sd_b^2 Z_b Z_b' + sigma^2 I. Every k
from 1 to 31: error 2.8e-14 to 4.3e-14 (the joint posterior is
Gaussian, so one node is exact).

### Inner sd near zero against the one-block route (`dev/nestquad-nodes.R zero`)

Bernoulli with no inner effect, seed 20261020, one-block GK fit `fs`;
nested objective at fs's beta and sd_a with sd_inner fixed, 11 x 11
nodes:

    ZERO sd_inner 1e-01 nested-gk -1.723e-02 gk-ref -1.015e-09
    ZERO sd_inner 1e-02 nested-gk -1.720e-04
    ZERO sd_inner 1e-03 nested-gk -1.720e-06
    ZERO sd_inner 1e-04 nested-gk -1.631e-08

The gap closes as sd_inner^2, to the one-block value, which is itself
1e-9 from a one-block integrate() reference.

### Standard errors (`dev/nestquad-se.R`, log `dev/nestquad-se.log`)

Weighted Bernoulli test design; sdreport of the fit against the inverse
numerical Hessian (step 1e-3) of the 80-node reference:

    SE par beta   fit 0.857664 ref 0.857711 ratio 0.99994
    SE par beta   fit 0.456183 ref 0.456194 ratio 0.99997
    SE par theta  fit 0.489925 ref 0.489940 ratio 0.99997
    SE par theta  fit 0.469382 ref 0.469467 ratio 0.99982

### Variants measured and not kept

1. Inner rule centered on each (inner level, outer node)'s own
   conditional mode and scaled by its curvature, found by Newton steps
   on doubles (`dev/nestquad-adapt-version.R.txt`, local). Poisson
   (`dev/nestquad-nodes-pois-adapt.log`): k 11 err 5.0e-8, k 15
   6.2e-10, the same as the linear rule. Bernoulli
   (`dev/nestquad-variants-*.log`, NQ_VARIANT=both|mode|scale|max,
   errors at (11,15), (21,21), (31,31)):

       linear 7.0e-06 2.5e-07 6.5e-10
       both   1.1e-04 8.4e-06 1.1e-07
       mode   1.6e-05 2.5e-07 1.2e-08
       scale  1.4e-04 1.1e-05 2.2e-07
       max    4.4e-05 6.9e-06 5.9e-08

   The adapted scale is what hurts. Kept the linear rule.
2. Blocks on two predictors. `(1 | ga)` in mu and sigma, gaussian, test
   design (`dev/nestquad-musigma.R`, logs `dev/nestquad-musigma-*.log`),
   reference by nested integrate() -87.5895224294, Laplace error 1.1:
   linear rule, mu outer: k 11 err 6.9e-3, 21 6.8e-4, 31 1.4e-4, 50
   1.8e-5; sigma outer: 11 6.2e-3, 31 1.3e-4; sigma outer with the
   adapted inner rule: (11,11) 1.1e-4, (15,15) 1.7e-6, (21,21) 7.4e-8.
   The sigma effect changes the width of the mu integrand from node to
   node, which one width per inner level cannot follow. Supporting it
   needs the adapted inner rule for this case only and a block order;
   not done. Such models are refused by name (they were refused before).

## Memory and time

Peak working set added by the fit (ps::ps_memory_info, one process per
configuration), `dev/nestquad-cost.R` via `dev/nestquad-cost-grid.sh`,
seed 20261106, 4 inner levels of 5 rows per outer level unless stated,
default nodes 11 x 15. `pred` is the guard's prediction at the time of
the run (tape bytes x 1, before the factor below); `tapeMB` the final
tape's bytes by the same formula. Log `dev/nestquad-cost3.log`:

    fam        n     time   deltaMB  pred   tapeMB
    bernoulli  1000   6.50s   345     22.5    23.6
    bernoulli  4000  28.83s   896     89.9    94.5
    bernoulli  8000  62.95s  1683    180     188.9
    poisson    1000   5.14s   455     35.8    36.9
    poisson    4000  17.08s  1273    147     151.5
    cumulative 4000  51.93s  1722    178     182.7
    gaussian   4000  12.28s  1631    168     172.6
    beta       1000  24.00s  1189    119     120.6
    beta       4000 108.25s  4452    478     482.4
    poisson    4000 (2 x 20 rows)  18.14s 1131 127 128.1
    bernoulli  4000 (20 x 1 row)   39.07s 1710 195 217.3
    bernoulli  4000 21 x 21        89.00s 2183 239 250.8

Memory is linear in rows and in the node product. The probe prediction
of the tape is 1 to 10 percent low (195 against 217.3 on the
one-row-per-inner-level design). The fit adds 7.9 to 14.6 times the
tape's bytes. For comparison, the marginal_gk route of the fuzzfix
patch: 1356 MB at 72 Bernoulli rows (19 MB per row); here 0.22 MB per
row at 4000 rows.

### Censored Beta, and the guard's two factors

The fuzz spec family=Beta aterm=cens re=nested (100 rows, seed
20535438; `dev/nestquad-betacens.R`) has 837 tape values per stacked
row against 5 to 25 for the families above, and taping it takes far
more beside the tape. Logs `dev/nestquad-betacens-*.log`:

    BETACENS nodes 5x5   time 40.2s  deltaMB 1753  tape 2093010 values
    BETACENS nodes 11x15 time 107.9s deltaMB 10579 tape 13812106 values
    (after releasing the old tape before each re-anchor)
    BETACENS nodes 5x5   time 41.3s  deltaMB 1570
    BETACENS nodes 11x15 time 112.0s deltaMB 9995

By phase (`dev/nestquad-betacens-phase.R`, 5 x 5): the peak rose from
518 to 1934 MB inside one `RTMB::MakeADFun()` while the steady size rose
by 85 MB; a `gc()` between inner-node blocks during taping left it at
1844 MB, so the transient is in the C++ taping, not R garbage. On
Bernoulli (`dev/nestquad-memphase.R 200 11 11`) the same phase was 406
to 784 MB.

Releasing the previous tape before a re-anchor (and a `gc()` before
each tape) took Bernoulli 4000 rows from 896 to 787 MB
(`dev/nestquad-cost4.log`).

The guard therefore uses two factors on the predicted tape bytes, plus
128 MB: 11 when the probes show at most 50 tape values per stacked row
(measured peak beyond the fixed part 7.3 to 9.2 times the final tape,
twelve fits), 25 above (measured 23 to 24, censored Beta). Prediction
against the measured peak on the final build (`dev/nestquad-cost6.log`,
`dev/nestquad-betacens-5c.log`, fuzz subset):

    fam                 n    pred MB  peak MB
    bernoulli         4000    1117     794
    bernoulli 20 x 1  4000    2270    1405
    poisson           4000    1744    1098
    gaussian          4000    1976    1415
    cumulative        4000    2087    1485
    beta              4000    5.26 GB refused (4452 MB measured earlier)
    censored beta 5x5  100    1713    1570
    censored beta      100    10.3 GB refused at 4 GB (9995 MB measured)

Every prediction is above the peak, by 1.06 to 1.62 times (10.3 GB is
10547 MB against 9995). The "beta 4000" refusal is correct: the fit
took 4452 MB when the limit did not apply.

## Tests

One file per process, NOT_CRAN, lane runner `dev/nestquad-run-tests.R`.

Seen to fail on rellib-r18 (logs `dev/nestquad-final-base-*.log`), and
on the lane (`dev/nestquad-all.log`):

    base test-quadrature-nested.R  pass=11 fail=5 err=5 skip=0 warn=0
    lane test-quadrature-nested.R  pass=37 fail=0 err=0 skip=0 warn=0
    base test-quadrature-defects.R pass=28 fail=0 err=2 skip=0 warn=0
    lane test-quadrature-defects.R pass=55 fail=0 err=0 skip=0 warn=0

Every base failure is the 0.73.0 refusal ("integrates one scalar random
effect per row"), raised by `frm()` on a model the lane fits or refuses
with another message: behavioral, not a missing symbol. The tests that
call internals (convergence, sd near zero, guard) error on that same
refusal before reaching them.

Under OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`: "Rblas: 74 to OpenBLAS", "Rlapack: 616 to OpenBLAS"), logs
`dev/nestquad-ob-*.log`:

    test-quadrature-nested.R  pass=37 fail=0 err=0 skip=0 warn=0
    test-quadrature-defects.R pass=55 fail=0 err=0 skip=0 warn=0

Whole suites, one file per process, 3 at a time (`dev/nestquad-run-par.sh
all dev/nestquad-all.jobs 3`, log `dev/nestquad-all.log`, "ALL ran 499
of 499"); every log's `lib:` line reads the lane's core from
`wt-nestquad-lib` and the extensions from rellib-r18 (unchanged):

    frmtmb          files=308 pass=22550 fail=0 err=0 skip=239 warn=0
    frmtmb.coupling files=11  pass=551   fail=0 err=0 skip=5   warn=0
    frmtmb.eam      files=33  pass=1848  fail=0 err=0 skip=5   warn=0
    frmtmb.latent   files=12  pass=683   fail=0 err=0 skip=2   warn=0
    frmtmb.learn    files=16  pass=474   fail=0 err=0 skip=13  warn=0
    frmtmb.ode      files=11  pass=549   fail=0 err=0 skip=1   warn=0
    frmtmb.sample   files=85  pass=3962  fail=0 err=0 skip=18  warn=0
    frmtmb.spatial  files=8   pass=478   fail=0 err=0 skip=6   warn=0
    frmtmb.survival files=15  pass=552   fail=0 err=0 skip=0   warn=0

Every skip is a gate (199 + 52 FRMTMB_BRMS_FIT_TESTS, 15 scale, 13
drmTMB, 1 fuzz, ...); the gated tiers were not run.

One session per suite, three orders (`dev/nestquad-run-sessions.sh
sessions 3`, log `dev/nestquad-sessions.log`, "SESSIONS ran 27 of 27"):
fail 0 and err 0 in all 27; all 46 LEAK lines carry `[with load of
...]`. Core reports files=295/308 because 13 files are gated whole
(the 13 files with pass=0 skip=0 in the per-process run);
frmtmb.sample 80/85 likewise.

Fuzz: the six specs of the gated plan (seed 20260901, 300 rows) with
mode=quadrature and re=nested (`dev/nestquad-fuzz-subset.R`, logs
`dev/nestquad-fuzz-{lane,base}.log`). Base: all six are refused
(`real_new` under this tree's helper, which no longer lists the
refusal). Lane: five fit and pass every invariant (student + weights,
cumulative + weights, binomial + trials, poisson + cens_chr, Gamma +
weights with a shape predictor); the sixth, Beta + cens, is refused by
the memory guard at 10.3 GB and matched by the new known refusal
`quad-nested-memory`.

## Composition

- `check_weights_unbounded()` runs before the route with `exact = TRUE`
  as for any quadrature fit; a one-row inner block, `(1 | ga) + (1 |
  obs)`, with weights above the rows is refused by it, and below them
  fits with the Laplace value (gaussian, exact either way); test "a
  one-row inner block composes with the weights check".
- frmtmb.sample's call to `check_weights_unbounded()` uses `route =
  "sample"` and does not read `quadrature`; nothing changed there, and
  its suite is green.
- The SE check skips quadrature fits as before (`fit$quadrature`); the
  standard errors come from `sdreport()` of the nested tape, which has
  no random effects, and agree with the reference Hessian
  ("Standard errors").
- `quad_check_blocks()` keeps its name and its one-block behavior; the
  by-split and one-block tests are unchanged and pass.
- `(1 + x || g)` is one `diag` block of dimension 2, not two blocks, and
  stays refused by the scalar check.

## Post-fit methods (`dev/nestquad-postfit.R`, log `dev/nestquad-postfit.log`)

RETRACTED in fix round 1. This section said `confint(method =
"profile")` "refits with a map through the nested route". It did not:
`tmbprofile()` read the fit's frozen tape, and the review measured an
interval of (0.520, 0.953) around an estimate of 1.171 (finding M1).
The check behind the claim was that the call returned. Profile and
root-search intervals are now refused on nested fits ("Fix round 1").

What did hold, on the unweighted test design: `confint()` (Wald),
`update(control = frmtmb_control(quad_nodes = 21))`, `predict(newdata
=)`, `simulate()` and `ranef()` return. logLik at the default
nodes -41.68766556, at 21 nodes -41.68766572. `predict()` warns, as on
a one-block quadrature fit, that the group effects are held at their
modes.

## Bump and floor

Minor bump (frmtmb 0.80.0 or whatever the round's version is): a model
class that was refused now fits, and `frmtmb_control()` has two new
arguments. No extension calls anything this lane added or changed, and
every extension suite passes against the lane's core with the
extensions of rellib-r18, so no extension needs a new floor.

## Check

`R CMD check --as-cran` of frmtmb, built with vignettes, R_LIBS the
lane library, rellib-r18 and the user library,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE` (`dev/nestquad-check.ps1`; logs
`dev/nestquad-00check.log`, `dev/nestquad-check-run.log`):

    Status: 1 NOTE

The NOTE is the known one: "checking HTML version of manual ... NOTE
Skipping checking math rendering: package 'V8' unavailable".

## Needs the user

Updated in fix round 1.

- M1: profile and root-search intervals are refused on nested fits. A
  correct one by nested refits costs 21 refits per parameter (23.5 s
  on 117 rows, against 0.91 s for a Laplace profile); building it is
  the user's decision ("Fix round 1", M1).
- The default nodes, `c(11, 15)`: 9e-11 to 7e-6 on seven designs with
  sds near 1, 5e-4 to 6e-4 with an inner sd near 3 or an outer one of
  3, 1.7e-2 on the near-separable design (warns). More nodes do not
  always help on the hard designs (21 x 21 is worse than 11 x 15 on
  two of them).
- The default memory limit, `quad_max_gb = 4`. With the per-family
  costs, the default refuses, by prediction, Bernoulli above about
  22000 rows, Poisson above about 14000, Student above about 1750,
  and a censored Beta above about 35 rows (linear in rows from
  `dev/nestquad-mem-pred.log`). An unmeasured family gets the largest
  measured cost.
- `quad_tape_size()` reads tape sizes through `TMB:::info()`, by
  `utils::getFromNamespace()`, because TMB exports no accessor;
  without it the guard uses the per-row costs (m6).

## Not done, and why

- Three or more nested levels: the iterated rule generalizes (one more
  loop per level), but the cost per row is the product of the node
  counts, and the brief asks for two levels first.
- Blocks on two predictors (see "Variants measured and not kept").
- Student-t blocks under the nested route: the per-level prior is read
  as Gaussian (`imp_block_prior()`).
- The gated tiers (FRMTMB_BRMS_FIT_TESTS, scale, drmTMB, fuzz in full):
  no gated test reads quadrature over two blocks; the fuzz specs that do
  were run as a subset.
- A profile of the tape for speed. Time is dominated by the optimizer's
  evaluations of a tape of n x 165 stacked rows (Bernoulli 4000 rows:
  0.22 s per fn and gr pair at 11 x 11, `dev/nestquad-memphase.R`); the
  fit took 28.8 to 40.9 s at 4000 Bernoulli rows and 108 s at 4000 Beta
  rows. Not optimized further.

## Defects found, not fixed

- None in other code. The 0.73.0 refusal is gone with the route.

## Fix round 1 (review `dev/reviews/2026-11-06-nestquad.md`)

Sections above this one describe round 0. Where they disagree with this
section, this section wins: the memory model of "Memory and time" and
"Censored Beta, and the guard's two factors" is replaced (M3 below),
and the profile claim of "Post-fit methods" is retracted there.

### M2: why the anchored rule exceeded the integral, and the loop

Cause. Round 0 froze the nodes on the effects' own scale and kept the
prior density in the integrand. At a frozen node u_k the prior density
N(u_k; 0, sd) grows like 1 / sd as sd shrinks while u_k stays near 0,
so the frozen tape grows without bound as an sd falls toward a node,
and the optimizer of a round walked there: the review's near-separable
Bernoulli (outer sd 0.001 at the end) reported +29.6 in round 1.

Fix 1, the standardized rule. The nodes are now frozen on z = u / sd:
the tape reads u = sd(theta) z, every weight of the rule (node weight,
phi(z) at the node, the scale) is a constant, and theta reaches the
integrand through the likelihood alone. At the anchor it is the same
rule; away from it no prior density can blow up. The rule's value at
its anchor is unchanged: the review's accuracy designs give the same
errors at the fit as round 0 (`dev/nestquad-fr1-acc-*.log` against
`dev/review-nestquad-out/acc-*.log`, e.g. unbal 2.410e-08 both, dummy
6.237e-04 against 6.1e-4, bigin 4.611e-04 against 4.6e-4). The lane's
own Bernoulli node table is unchanged to 4 digits
(`dev/nestquad-fr1-nodes-bern.log`: k 11 1.301e-04, k 15 8.563e-06,
k 21 2.516e-07, k 31 6.468e-10).

Fix 2, a trust region on the self-anchored objective. Each round
optimizes the frozen tape inside a box of the trust radius (in
parameter units) around the anchor, re-anchors at the round's optimum,
and accepts the round only if the re-anchored objective is no worse
than at the start (both are the rule anchored at its own point, so the
accepted values fall and the loop cannot cycle). An accepted round at
the box edge doubles the radius; a rejected one cuts it to a quarter of
its step; three rejections in a row end the loop. A loop that ends
without a stationary anchor warns with class `frmtmb_quad_unsettled`
and gives the gradient and the last move. On the near-separable design
(`dev/review-nestquad-acc.R sep`, MODE=anchor, log
`dev/nestquad-fr1-sep.log`), without the trust region the first round
stepped by 7114 and was rejected even at an eighth; with it:

    round accepted step 1, objective 9.8914406, grad 0.636
    round accepted step 2, objective 9.0047867, grad 1.19
    round accepted step 4, objective 8.3952863, grad 0.138
    then rejected steps of 8, 2, 0.5 (three in a row end the loop)
    WARN frmtmb_quad_unsettled ... gradient of 0.138 ...
    ANCHOR logLik -8.395286301939 re-anchored-at-par -8.395286301659

The value is the rule anchored at the estimate (difference 2.8e-10).
Against the review's integrate() reference at that estimate
(`dev/nestquad-fr1-sep-intref.log`): -8.3778083592, error 1.748e-02;
the rule's own error there is 1.7e-2 at 11 x 15, 6.0e-2 at 21 x 21,
5.8e-3 at 31 x 31, so no improvement smaller than the rule's error can
be seen, which is why the loop stops. Round 0 reported -8.3066 at 4
rounds and -8.4761 at 20 on the same design, silently.

Fix 3, the bound. `quad_nest_check_bound()` stops (class
`frmtmb_quad_bound`) when a discrete response's log-likelihood (the
objective less the prior) is above 1e-8 x rows. Not reached by any fit
after fixes 1 and 2; pinned by a direct test.

### M1: profile and root search

Refused by class `frmtmb_quad_nested_unsupported` in every route that
reads the tape away from the estimate: `confint(method = "profile")`,
`"uniroot"`, `profile()`, `hypothesis(method = "profile")`,
`conditional_effects(band = "profile")` (`quad_nested_refuse_profile()`,
called from `R/confint.R` and `R/conditional-effects.R`).

The correct route was measured, and not built
(`dev/nestquad-profcost.R`, log `dev/nestquad-profcost.log`, the
review's 117-row design, seed 20261107): a profile interval for x by
refitting through the nested route with x fixed by an offset,
`uniroot()` on 2 x drop = 3.84:

    PROF fit 1.78s logLik -63.86886708 x 1.171178 settled TRUE rounds 1
    PROF laplace tmbprofile points 76 time 0.91s
    PROF refit interval (0.7011, 1.7702) Wald (0.6463, 1.6960)
    PROF refits 21 total 23.5s mean 1.12s per refit
    PROF 2*drop at the endpoints 3.8403 3.8400

So a correct interval costs 21 nested refits per parameter, 23.5 s here
against 0.91 s for the Laplace profile; at 4000 Bernoulli rows (30 to 40
s per fit) about 12 minutes per parameter. Building it means keeping the
Laplace object and the fit's closures on the fit (lost on saveRDS) and
mapping one outer parameter in every component; refused instead, and
the decision is the user's. The test pins the refusal, not endpoints.

### M3: memory, per family

`dev/nestquad-mem.R` (lane copy of the review's script, seed 20261108,
20 rows per school in 4 classes of 5, nodes 11 x 15), peak working set
added by the fit, logs `dev/nestquad-mem1.log`, `-mem2.log`. "factor"
is (peak - 128 MB) / predicted tape bytes, "kB/row" is (peak - 128 MB)
per stacked row:

    poisson                    n  3000 peak    676 MB tape   97.60 MB factor   5.61 kB/row   1.13 values/row 7.7
    bernoulli                  n  3000 peak    527 MB tape   54.60 MB factor   7.30 kB/row   0.82 values/row 4.4
    student                    n  3000 peak   5270 MB tape  157.00 MB factor  32.75 kB/row  10.64 values/row 11.4
    gaussian                   n  3000 peak    904 MB tape  113.00 MB factor   6.87 kB/row   1.61 values/row 8.4
    exgaussian                 n   600 peak   5156 MB tape  604.00 MB factor   8.32 kB/row  52.01 values/row 203.4
    cumulative                 n  3000 peak   1080 MB tape  133.00 MB factor   7.16 kB/row   1.97 values/row 10.0
    Gamma                      n  3000 peak   1218 MB tape  140.00 MB factor   7.79 kB/row   2.25 values/row 10.4
    beta+cens                  n    40 peak   3622 MB tape   18.10 MB factor 193.04 kB/row 542.09 values/row 91.9
    beta                       n  3000 peak   3014 MB tape  346.00 MB factor   8.34 kB/row   5.97 values/row 24.4
    hurdle_poisson             n  3000 peak    989 MB tape  134.00 MB factor   6.43 kB/row   1.78 values/row 10.2
    negbinomial                n  3000 peak   2032 MB tape  251.00 MB factor   7.58 kB/row   3.94 values/row 18.0
    zero_inflated_negbinomial  n  3000 peak   1832 MB tape  216.00 MB factor   7.89 kB/row   3.52 values/row 15.5
    lognormal+cens             n  3000 peak   1041 MB tape  120.00 MB factor   7.61 kB/row   1.89 values/row 9.0
    weibull                    n  3000 peak   1734 MB tape  185.00 MB factor   8.68 kB/row   3.32 values/row 13.4
    skew_normal                n  3000 peak   1780 MB tape  194.00 MB factor   8.51 kB/row   3.42 values/row 14.4
    poisson                    n   600 peak    132 MB tape   18.90 MB factor   0.24 kB/row   0.05 values/row 7.5
    bernoulli                  n  1000 peak    207 MB tape   18.20 MB factor   4.35 kB/row   0.49 values/row 4.4
    poisson+cens               n  3000 peak  10768 MB tape  679.00 MB factor  15.67 kB/row  22.01 values/row 46.7
    gaussian                   n   600 peak    231 MB tape   22.70 MB factor   4.55 kB/row   1.07 values/row 8.4
    student                    n  1000 peak   1671 MB tape   52.20 MB factor  29.55 kB/row   9.57 values/row 11.4
    beta+cens                  n    20 peak   1861 MB tape   69.90 MB factor  24.79 kB/row 537.64 values/row 699.5

Two findings decide the model. The taping overhead depends on the
density's operations, not their count: Student needs 33 times its tape
with fewer values per row than negbinomial (7.6). A single RTMB
operation on 1e6 elements shows the same
(`dev/nestquad-transient-probe.R`): peak per tape value 86 bytes for
`dnorm`, 92 for `lgamma`, 88 for `plogis`, 218 for `dt`. And a censored
Beta costs about 540 kB per stacked row whatever its final tape: at 40
rows its tape was 18 MB, at 20 rows 70 MB (a different share of censored
rows), and the peak per stacked row was 542 and 538 kB. Its `pbeta()`
branch is taped for every row and then eliminated, so only the count of
stacked rows predicts it.

The guard now predicts 128 MB + max(factor x tape, kB/row x stacked
rows) with both numbers per family key (family, plus "+cens" when any
row is censored), each the larger of two runs times 1.3
(`quad_nest_mem_cost`). A family not listed gets the largest factor
(42.9), and when censored the largest censored cost per row (711 kB);
without the tape's size (no `TMB:::info()`), the per-row term alone,
with the largest uncensored cost per row (67.6 kB) for an unlisted
uncensored family. The probes now run on a placeholder anchor (nonzero
values, because RTMB folds a constant zero operand), before the Laplace
warm start (m8).

Prediction against the measured peak, a third run of every family
(`dev/nestquad-mem3.log`), the final constants applied to each run's
logged tape and stacked rows (`dev/nestquad-mem-pred.R`, log
`dev/nestquad-mem-pred.log`; arithmetic on the logged numbers):

    PRED poisson                    n  3000 measured    750 MB predicted    950 MB ratio 1.27 | without tape size    950 MB ratio 1.27
    PRED bernoulli                  n  3000 measured    507 MB predicted    660 MB ratio 1.30 | without tape size    660 MB ratio 1.30
    PRED student                    n  3000 measured   5314 MB predicted   6896 MB ratio 1.30 | without tape size   6896 MB ratio 1.30
    PRED gaussian                   n  3000 measured    936 MB predicted   1191 MB ratio 1.27 | without tape size   1191 MB ratio 1.27
    PRED exgaussian                 n   600 measured   5154 MB predicted   6664 MB ratio 1.29 | without tape size   6664 MB ratio 1.29
    PRED cumulative                 n  3000 measured   1053 MB predicted   1385 MB ratio 1.32 | without tape size   1385 MB ratio 1.32
    PRED Gamma                      n  3000 measured   1232 MB predicted   1570 MB ratio 1.27 | without tape size   1530 MB ratio 1.24
    PRED beta+cens                  n    40 measured   3661 MB predicted   4711 MB ratio 1.29 | without tape size   4711 MB ratio 1.29
    PRED beta                       n  3000 measured   3031 MB predicted   3899 MB ratio 1.29 | without tape size   3899 MB ratio 1.29
    PRED hurdle_poisson             n  3000 measured   1009 MB predicted   1288 MB ratio 1.28 | without tape size   1288 MB ratio 1.28
    PRED negbinomial                n  3000 measured   2043 MB predicted   2613 MB ratio 1.28 | without tape size   2593 MB ratio 1.27
    PRED zero_inflated_negbinomial  n  3000 measured   1833 MB predicted   2353 MB ratio 1.28 | without tape size   2352 MB ratio 1.28
    PRED lognormal+cens             n  3000 measured   1078 MB predicted   1385 MB ratio 1.28 | without tape size   1385 MB ratio 1.28
    PRED weibull                    n  3000 measured   1742 MB predicted   2255 MB ratio 1.29 | without tape size   2255 MB ratio 1.29
    PRED skew_normal                n  3000 measured   1749 MB predicted   2303 MB ratio 1.32 | without tape size   2303 MB ratio 1.32
    PRED bernoulli                  n  1000 measured    201 MB predicted    305 MB ratio 1.52 | without tape size    305 MB ratio 1.52
    PRED student                    n  1000 measured   1685 MB predicted   2384 MB ratio 1.41 | without tape size   2384 MB ratio 1.41
    PRED beta+cens                  n    20 measured   1838 MB predicted   2419 MB ratio 1.32 | without tape size   2419 MB ratio 1.32
    PRED poisson+cens               n  1500 measured   5457 MB predicted   8778 MB ratio 1.61 | without tape size   7041 MB ratio 1.29

No under-prediction on any measured family: 1.24 to 1.61, with and
without the tape's size. m6: without the tape's size the 72-row
Bernoulli test design is predicted at about 141 MB (128 MB + 72 x 165
stacked rows x 1.1 kB), and the test asserts it is under the default
limit. The fuzz Beta + cens spec is refused at 12.3 GB (measured 10.0
GB in round 0).

### Minor findings

- m1: `?frm` (`quadrature`) and NEWS give the measured range: 9e-11 to
  7e-6 on seven designs with sds near 1 (the review's accuracy script
  rerun on this build: unbal 2.4e-8, slope 2.4e-7, intslope 2.2e-6,
  binom 9.0e-11, negbin 2.7e-10, bigout at the fit 5.9e-8, the lane's
  Bernoulli 7e-6), 5e-4 to 6e-4 with an inner sd near 3 or an outer one
  of 3 (dummy 6.2e-4, bigin 4.6e-4, bigout at sds 3 / 0.1 5.7e-4), and
  1.7e-2 on the near-separable design, which warns. Default nodes
  unchanged.
- m2: the first test's bound is now `err_q / err_l < 5e-5` (measured
  8e-6 / 0.85), and a new test on the review's unbalanced design asserts
  `< 2e-7` (measured 2.4e-8 / 0.85). On the review's mutant m1 (inner
  rule on the inner marginal; `rev-nestquad-lib/m1`, log
  `dev/nestquad-fr1-t1-m1.log`) they fail at 4.0e-4 (8 times the bound)
  and 7.9e-7 (4 times).
- m3: the block-count refusal names its reason in order: a by-split
  term; blocks on different predictors (the grouped mixture); a row
  reading three or more blocks. `mvbf()` is refused before it as "for
  one response". The crossed message gives each direction its clause
  ("so tt does not nest in ga", "so ga does not nest in tt").
  `dev/nestquad-fr1-struct.log` (the review's struct script on this
  build) shows each.
- m4: `(1 + x || g)` is one `diag` block of two coefficients; the
  scalar refusal now says so and gives `(1 | g) + (0 + x | g)`, which
  fits. The message no longer says "intercepts only", and a
  dimension-1 block of another structure is named by its structure.
- m5: `ranef(condVar = TRUE)` on a nested fit stops with class
  `frmtmb_quad_nested_unsupported`.
- m7: `?frmtmb_control` gives the measured ratios (1.24 to 1.61).
- m8: the guard runs before the Laplace warm start; the test asserts no
  "quadrature warm start" message precedes the refusal.

### Merge notes

- arstep: its `opt <- gp_second_descent(...)` after the Laplace
  branch's `}))` in `fit_assembled()` goes inside the first branch,
  before this lane's `} else if (!is.null(qnest)) {`.
- bigx: it attaches its optimization basis to `obj$env` of the Laplace
  object; `quad_nest_fit()` optimizes tapes it builds itself
  (`RTMB::MakeADFun(qo$fn, otpl, ...)` in `freeze()`), so they carry no
  basis, and its trust region is a box in `par_units`. For a bigx-style
  design under nested quadrature, consolidation must attach the same
  basis to each tape `freeze()` returns (or pass it to the
  `optimize_obj()` call of each round, with the box in the basis's
  units), or the rounds run unscaled.
- loglik: under its split, a nested fit's `logLik()` is
  `-opt$objective - fit_log_prior(fit)`; `quad_nest_objective()` adds
  `neg_log_prior(pars)` as `nll` does, and `quad_nest_check_bound()`
  already reads the log-likelihood that way.

### Tests and check, fix round 1

One file per process, NOT_CRAN, `dev/nestquad-run-tests.R`; the whole
core suite rerun because `R/fit.R`, `R/confint.R`, `R/methods-fit.R` and
`R/conditional-effects.R` changed (`dev/nestquad-run-par.sh core`, log
`dev/nestquad-core.log`, "CORE ran 308 of 308", every `lib:` line the
lane library):

    frmtmb files=308 pass=22576 fail=0 err=0 skip=239 warn=0

The quadrature files, from that log (the review's list of 21):

    RESULT test-quadrature-nested.R pass=63 fail=0 err=0 skip=0 warn=0
    RESULT test-quadrature-defects.R pass=55 fail=0 err=0 skip=0 warn=0
    RESULT test-autocor.R pass=192 fail=0 err=0 skip=0 warn=0
    RESULT test-compat.R pass=759 fail=0 err=0 skip=0 warn=0
    RESULT test-covstruct-seam.R pass=131 fail=0 err=0 skip=0 warn=0
    RESULT test-dry-run.R pass=23 fail=0 err=0 skip=0 warn=0
    RESULT test-importance.R pass=344 fail=0 err=0 skip=0 warn=0
    RESULT test-latres.R pass=298 fail=0 err=0 skip=0 warn=0
    RESULT test-me.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT test-mo-factor-re.R pass=34 fail=0 err=0 skip=2 warn=0
    RESULT test-predfix.R pass=92 fail=0 err=0 skip=0 warn=0
    RESULT test-predict-re-uncertainty.R pass=71 fail=0 err=0 skip=0 warn=0
    RESULT test-ps.R pass=71 fail=0 err=0 skip=0 warn=0
    RESULT test-review-v25.R pass=60 fail=0 err=0 skip=0 warn=0
    RESULT test-se-exact-hessian.R pass=35 fail=0 err=0 skip=0 warn=0
    RESULT test-se-tier.R pass=95 fail=0 err=0 skip=0 warn=0
    RESULT test-student-share.R pass=28 fail=0 err=0 skip=0 warn=0
    RESULT test-tmb-examples.R pass=27 fail=0 err=0 skip=0 warn=0
    RESULT test-tre.R pass=102 fail=0 err=0 skip=0 warn=0
    RESULT test-v11.R pass=22 fail=0 err=0 skip=0 warn=0
    RESULT test-v19.R pass=23 fail=0 err=0 skip=0 warn=0

frmtmb.latent's two files that match the quadrature refusal text, with
the lane's core (`dev/nestquad-fr1-lca.log`, `-hmm.log`):

    RESULT test-lca.R pass=105 fail=0 err=0 skip=0 warn=0
    RESULT test-hmm.R pass=120 fail=0 err=0 skip=0 warn=0

Seen to fail: on the review's mutant m1 build, `test-quadrature-nested.R`
pass=37 fail=14 err=2 (`dev/nestquad-fr1-t1-m1.log`): the two accuracy
bounds (4.0e-4 and 7.9e-7 against 5e-5 and 2e-7), every new message,
the profile and condVar refusals, the loop's warning, and the bound
check (missing function). m1 carries round 0's loop and messages, so
those failures are the round-0 behavior. On rellib-r18 the file fails
on the 0.73.0 refusal throughout.

The fuzz subset on this build (`dev/nestquad-fr1-fuzz.log`): five of
the six nested quadrature specs pass every invariant; the Beta + cens
spec is the known refusal `quad-nested-memory`, at 12.3 GB.

`R CMD check --as-cran` of frmtmb on this build (`dev/nestquad-check.ps1`,
logs `dev/nestquad-fr1-00check.log`, `dev/nestquad-fr1-check-run.log`):

    Status: 1 NOTE

The NOTE is the known one: "Skipping checking math rendering: package
'V8' unavailable".

## Fix round 2 (short check SHIP; three minors)

1. `frmtmb_quad_unsettled` now gives the rule's own error at the
   estimate: the absolute log-likelihood difference between the fit's
   nodes and a larger rule (21 x 21, or 10 more per level), on doubles at
   the same anchor, also stored as `fit$quad_nested$rule_error`. `?frm`
   says more nodes are not a sure remedy there. On the near-separable
   design (`dev/nestquad-fr2-sep.log`): "the rule's own error there is
   of the order of 0.043 (the log-likelihood at 11 x 15 nodes against
   21 x 21)"; the error against the integrate() reference is 1.7e-2
   (fix round 1), so the estimate is of the right order and errs high.
   The near-separable test pins the message's content and the stored
   value.
2. Memory entries for binomial (trials), lognormal, von_mises and
   beta_binomial (Gamma, negbinomial and cumulative were listed), two
   runs each of 3000 rows and one binomial run of 1000
   (`dev/nestquad-mem4.log`; constants the larger run times 1.3).
   Predicted over measured (`dev/nestquad-mem-pred.R nestquad-mem4.log`,
   log `dev/nestquad-mem-pred4.log`):

       PRED lognormal      n 3000 measured  974 MB predicted 1240 MB ratio 1.27
       PRED von_mises      n 3000 measured 1129 MB predicted 1482 MB ratio 1.31
       PRED binomial       n 3000 measured  581 MB predicted  715 MB ratio 1.23
       PRED lognormal      n 3000 measured  949 MB predicted 1240 MB ratio 1.31
       PRED binomial       n 3000 measured  579 MB predicted  715 MB ratio 1.24
       PRED beta_binomial  n 3000 measured 4383 MB predicted 5677 MB ratio 1.30
       PRED binomial       n 1000 measured  207 MB predicted  324 MB ratio 1.56
       PRED von_mises      n 3000 measured 1154 MB predicted 1482 MB ratio 1.28
       PRED beta_binomial  n 3000 measured 4384 MB predicted 5677 MB ratio 1.29

   Without the tape's size: 1.22 to 1.55. Before the entries, the
   unlisted factor of 42.9 predicted binomial at 3092 MB for 581 MB
   measured (5.3 times). Over all 20 measured families the ratio is now
   1.22 to 1.61; `?frmtmb_control` says so.
3. `profile(fit)` without `parm` raises
   `frmtmb_quad_nested_unsupported` before `parm` is read; pinned.

Tests, one file per process, lane library (`dev/nestquad-fr2-t*.log`):

    RESULT test-quadrature-nested.R pass=69 fail=0 err=0 skip=0 warn=0
    RESULT test-quadrature-defects.R pass=55 fail=0 err=0 skip=0 warn=0
