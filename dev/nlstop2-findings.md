# Lane nlstop2: nonlinear short stops, second round

Round opened 2026-11-13, base ec8178d3 (frmtmb 0.82.0). Worktree
`frmtmb-wt-nlstop2`, branch `wt-r83-nlstop2`, library `wt-nlstop2-lib`,
base library `rellib-r21`. Nothing committed. Logs are in
`dev/nlstop2-out/` (gitignored); every number below is pasted from the log
named beside it. Before-arm runs set `NLSTOP2_ARM=before` (rellib-r21
alone); after-arm runs put `wt-nlstop2-lib` first. Each log's first line
names its library. Scripts: `dev/nlstop2-*.R`; common setup
`dev/nlstop2-common.R`; runner `dev/nlstop2-run.sh`; job queue
`dev/nlstop2-queue.sh`.

Status: DONE after the fix round (review `dev/reviews/2026-11-13-nlstop2.md`,
verdict FIX; section "Fix round" at the end). The huber() quadrature
route of the first pass is REVERTED by the user's decision; the sections
on it below are kept as the evidence for the documentation item.

## Items, in the brief's order

| line | item | verdict |
|---|---|---|
| 5079 | MM fits trapped between or beside poles, 4.32 and 0.92 short; user starts at K -2.5 silent | FIXED: a restart beyond the poles from the default start, a warning from a given start |
| 5093 | `frmtmb_nl_limit` misses `a * (exp(b x) - 1)` and the Hill overflow | explin FIXED (a Newton step across the ridge); Hill MISDIAGNOSED: no limit (h to infinity is 0.147 lower), but the fit is 1.8e-3 short of h doubled, silent; filed LOW |
| 5041 = 2617 | F11 seed 8, 0.14 short at the sd boundary | NOT FIXED: 3 of 20 seeds short; only a refit per smooth sd reaches it; filed; deferred by the user |
| 5053 = 5487 | ridge tilted by a weak prior, 1.2e-3 short | NOT FIXED: narrowed to `restarts = 0` or a custom optimizer; filed |
| 5065 = 2602 | huber() with a group term stops short in silence | DOCUMENTATION ONLY (user's decision): the quadrature route was built, measured and reverted; `?huber` gives the numbers, when quadrature helps, its cost and its limits |
| B3 | 5626, 5697, 5881, 5576, 5641, 5682, 3359, 2818, 2590, 2599, 3013, 3017, 3022 | 3359 rerun: still open (warned), filed; 2818 rerun: does not reproduce, CLOSED; the rest not reached |
| 5718 | sign-probe switch in `frmtmb_control()` | ANSWERED: built at the consolidation, not by the lane |

## What changed

- `R/fit-restart.R` (owned): `nl_limit_message()` reads each move along
  the flattest direction as it lands and after one Newton step across it.
  New: `nl_denoms()`, `nl_pole_walls()`, `nl_pole_pressed()`,
  `nl_pole_escape()`, `nl_pole_message()`, `nl_pole_warning()`,
  `nl_refit_held()`.
- `R/par-template.R` (owned): `nl_body_at()` takes `exprs =` (the values
  of parts of the body, read by `nl_pole_walls()`); `nl_fit_probed()`
  calls `nl_pole_escape()` on the fit it returns, on both of its paths.
- `R/fit.R`: unchanged. The first pass's huber hunk and `@param
  quadrature` lines are reverted; `git hash-object R/fit.R` equals
  `git rev-parse ec8178d3:R/fit.R` (b68a5bdb), and `man/frm.Rd` likewise
  (36acd839).
- `R/families.R`: the `huber()` details paragraph (the documentation
  item); `man/frmtmb-families.Rd` regenerated.
- `tests/testthat/test-fit-restart.R`: `mm_data()` draws a binomial
  response too; 5 new tests.
- `NEWS.md`: a `# frmtmb (development version)` section: `## Fitting`
  with 2 bullets (the pole restart, BREAKING; the limit check) and
  `## Documentation` with 1 bullet (`?huber`).
- `dev/test-backlog.md`: CLOSED lines under 2818, 5065 (documentation
  only), 5079 and 5093, and the section "Filed in the round of
  2026-11-13", "From lane nlstop2" (8 items).

## 5079: Michaelis-Menten fits beside the poles of the body

Reproduced on rellib-r21 (`dev/nlstop2-pole.R`, the review's
`nlstop-review-pole.R` with the lane's setup; `pole-before.txt`):

```
ROW binomial 11 default -623.288529 0 - -0.087112 -0.182532 POLE
ROW binomial 11 Vm-2_K-6 -618.966489 0 - -0.12533 -6.0152 -
ROW binomial 44 default -626.898840 0 - -0.090366 -0.110348 -
ROW binomial 44 Vm-2_K-6 -625.980024 0 - -0.18333 -6.45313 -
ROW bernoulli 55 Vm-2_K-2.5 -207.650719 0 - 0.015741 -1.00267 POLE
ROW bernoulli 55 Vm0.5_K-2.5 -207.001511 0 - -0.0044968 -2.91724 POLE
```

Cause. `K + xp` is 0 at K = -xp_i for every row: the body has a pole
there, and the likelihood cannot be continued across it, so the walls
-xp_i cut the line of K into pieces and the optimizer stays in the piece
it starts in. The sign probe starts K at -1, inside the band of walls
(-5 to -0.1). Where each default sits (`dev/nlstop2-polegeo.R`,
`polegeo-before.txt`; dist is the distance from K to the nearest wall,
gap the median gap between walls):

```
ROW binomial 11 -623.288529 0 - - -0.087112 -0.182532 -4.989 -0.132 0.00547 0.012
ROW binomial 33 -630.923617 0 - - -0.12411 -0.164165 -4.998 -0.1 0.00584 0.011
ROW binomial 44 -626.898840 0 - - -0.090366 -0.110348 -4.998 -0.114 0.00365 0.0115
```

In median gaps: held fits 0.32 to 0.53 from a wall (binomial 11, 33, 44),
the code-1 fits of bernoulli 11, 33, 44 and binomial 55 at 1.6e-14 to
5.2e-13 (on a wall), and the fits clear of the band 68 to 669 (gaussian,
poisson, gamma, cumulative, 20 fits).

Fix (`nl_pole_escape()`). For each denominator of a body (`/` and a
negative literal power) that is affine in one of its coefficients, the
walls are read from three evaluations of the denominator (`nl_pole_walls()`).
A band holds the fit when the denominator takes both signs over the rows
(a pole between rows of the data) or the fit is nearer a wall than the
walls' median gap (`nl_pole_pressed()`). Then the optimizer runs on the
fit's tape from past each end of the band the fit is not beyond, half the
band's width out, the rest of the model at the fit. A run is kept if it
gains more than `nl_probe_gain()` and ends clear of the band; frm() then
refits from it, and the refit must end clear too. From the default start
the better fit is returned with a message (`frmtmb_nl_pole_start`); from
a given start the fit is kept and warns (`frmtmb_nl_pole`).

The first build required code 0 of the run; bernoulli 55 from K = -2.5
then stayed silent, because the run to the linear limit ends at code 7
(`poledbg1.txt`: from K 2.311, objective 206.278, "singular convergence
(7)"). The run's code is not read now; the refit reports its own.

The first build also counted a fit between the first and last wall as
held. On `a / (1 + b * x)` with x of both signs the walls lie on both sides
of every b, and the fit sits in the wide piece between them, 7.5 to 20.0
median gaps from a wall, with the denominator positive on every row
(`dev/nlstop2-gaps.R`, `gaps-after2.txt`, seeds 1 to 5). The pieces of
`K + xp` are up to 10.4 median gaps wide (same log), so no gap count
separates the two; the sign of the denominator over the rows does.

After (`pole-after5.txt`, `polegeo-after2.txt`):

```
ROW binomial 11 -618.966489 0 - frmtmb_nl_pole_start -0.12533 -6.0152
ROW binomial 44 -625.980024 0 - frmtmb_nl_pole_start -0.18333 -6.45313
ROW binomial 55 -627.468943 0 - frmtmb_nl_pole_start -0.12401 -6.04249
ROW bernoulli 55 Vm-2_K-2.5 -207.650719 0 frmtmb_nl_pole 0.015741 -1.00267 POLE
ROW bernoulli 55 Vm0.5_K-2.5 -207.001511 0 frmtmb_nl_pole -0.0044968 -2.91724 POLE
```

binomial 55 was 2.05 short at code 1 with warnings; it now reaches the
fresh rows' best, -627.468943. binomial 33 (inside the band, silent) is
not moved: no run beats it, and it is 9.3 above the best of 17 starts of
lane nlstop's reference (-640.258237). bernoulli 11, 33, 44 stop on a wall
at code 1 with warnings on both builds; their runs end on a wall again.

Rounding (`dev/nlstop2-ulp.R`, 20 draws `xp * (1 + 2e-15 k)`):
`ulp-pole-after1.txt` (first build): seed 11, 20 of 20 at K -6.0152,
2.0e-9 from the start at K = -6, message on 20; seed 44, 20 of 20 at K
-6.4531 within 6.6e-12, message on 20. `ulp-user-after1.txt`: 20 of 20
warn `frmtmb_nl_pole` at K -1.0027. On the final build of the first pass
(`ulp-pole-after2.txt`, `ulp-user-after2.txt`): 40 of 40
`frmtmb_nl_pole_start`, 20 of 20 `frmtmb_nl_pole`. The fix round does not
touch this code.

## 5093: the limit check's two misses

explin (`dev/nlstop2-fndebug.R`, `fndebug-before.txt`): the flattest
direction (-0.707, 0.707) moved by 1000 loses 0.0200 (side -1), the exact
ridge (1000, 1/1000) gains 0.003285. Fix: each move is also read after one
Newton step across the move's direction, on the curvature at the fit, and
the lower objective is taken. After (`fndebug-after5.txt`): "explin ll
-707.2658 code 0 warn frmtmb_nl_limit". Rounding (`ulp-explin-after1.txt`,
20 draws of x): 3 stop at code 0 with the limit warning, 17 at code 1 with
the optimizer's warnings; every draw warns.

Hill on step data: misdiagnosed. The profile in h, every other parameter
optimized with h held (`dev/nlstop2-hilldbg.R`, `hilldbg3.txt`):

```
profile h x 1.5: gain -0.00261905 code 0
profile h x 2: gain 0.00175847 code 0
profile h x 4: gain -0.146924 code 0
profile h x 8: gain -Inf code 0
```

h alone times 2, 4: -0.144, -0.147. The step limit h to infinity is
0.147 below the fit, so no limit warning is right. But the fit is not a
maximum: the profile with h doubled gains 1.76e-3, above `grad_tol`, and
the fit gives no warning. That short stop is filed LOW (corrected in the
fix round; the first pass said "the fit is a maximum"). A build that moved by smaller factors where
the body overflows, and along a direction rounded to (0, 0, 1), was tried
and removed: on this fit it said nothing either, and no case needs it.

The review's limit battery (`dev/nlstop2-limit.R`, 1680 fits; the final
build is in "Acceptance sets"). First build (`limit-after/*.txt`) against the review's
logs on 0.81.0's lane build: A, B, C-prior 0 `frmtmb_nl_limit` on 1440
fits on both; C-start 4 on the lane (fn_explin n 200 seed 5, n 1000 seeds
2 and 10, fn_mmlin n 200 seed 5), 0 in the review's log; those 4 rows
were silent there.

## 5041: F11, the sd of a smooth of sigma at its boundary

`dev/nlstop2-f11.R 1:20` on rellib-r21 (`f11-seeds-before.txt`): seeds 5,
8 and 14 are short by 0.532335, 0.141838 and 0.589451, code 0, no
warning; each time a start of theta_2 (the sd of `s(x2)` in sigma's
formula) at -10 reaches the higher value and stays at -10. The other 17
seeds read gap 0. Lane nlstop's probes (objective at log sd - 10 with
the rest held: 0.69 lower; Newton steps: diverge) stand. Not fixed: the
only test found is a refit per smooth sd (about one fit, 0.14 to 0.88 s
on these rows), the cost the brief declined for 2584. Filed with the
3-of-20 count.

## 5053: the weak-prior ridge

`dev/nlstop2-ridgedbg.R` (`ridgedbg-before.txt`, `ridgedbg-before2.txt`),
the four stops of `dev/nlstop2-ridgeprior.R`:

```
psd 100 k 5: short 0.001012 gmax 0.00543 pred 0.0003513 warn -
  default optimizer from the stop: short 0.0002959 code 0 warn frmtmb_conv_gradient
  nlminb from the stop, restarts = 0: short 0.0005174 code 0 warn frmtmb_conv_gradient
psd 1000 k 20: short 0.0001901 gmax 0.0293 pred 1.46e-05 warn -
  default optimizer from the stop: short 0.0001592 code 0 warn -
  nlminb from the stop, restarts = 0: short 0.0001725 code 0 warn -
psd 1000 k 50: short 0.001225 gmax 0.0732 pred 0.0004232 warn -
  default optimizer from the stop: short 0.0007029 code 0 warn -
  nlminb from the stop, restarts = 0: short 0.0009309 code 0 warn -
psd 1000 k 80: short 0.00316 gmax 0.117 pred 0.001817 warn frmtmb_conv_gradient
  default optimizer from the stop: short 0.0006945 code 0 warn -
  nlminb from the stop, restarts = 0: short 0.002014 code 0 warn -
```

Under the defaults every stop ends within `grad_tol`. The silent stops
above `grad_tol` come from the custom optimizer of the reproducer, which
moves the fit along the ridge after nlminb, and from `restarts = 0` (psd
1000 k 80: 2.0e-3 short). Newton steps from the stop do not see the gain:
the first chord step loses 40.8, 0.159, 66.2 and 588 and the next are
NaN; full Newton steps lose 1154, 1.05e-6, 5.58e-6 and 1.34e-5. A check
that sees it would have to run the optimizer, which is what `restarts`
does. Not fixed; filed.

## 5065: huber() with a group term

Detector looked for (`dev/nlstop2-huberkink.R`, `huberkink-before.txt`):
the Laplace marginal jumps where a residual at the modes crosses k, so is
a short stop at a residual on k? No: at the 11 silent short stops the
residual nearest k is 2.9e-4 to 1.1e-2 from it, like the fits that are not
short; only the code-1 stops (warned) sit on k (2e-15 to 5e-10).

Rows of lane nlstop, best of 8 jittered starts on the same route:

```
Laplace (huberkink-before.txt):   scale 1 rows 20 short 11 silent 6 max 0.188
                                  scale 1000 rows 20 short 11 silent 5 max 0.272
quadrature (huberquad-before.txt): scale 1 rows 20 short 0
                                   scale 1000 rows 20 short 0
lane default (huberkink-after1.txt): scale 1 rows 20 short 0 silent 0
                                     scale 1000 rows 20 short 0 silent 0
```

(Lane nlstop's quadrature run on 0.81.0 had 2 short at scale 1e3; on
0.82.0 none.) The first pass built `huber_quad_route()`, which sent such
a model to `quadrature = TRUE` when the call did not give `quadrature`
(rounding, `ulp-huber-after1.txt`: 20 of 20 by quadrature, the best of 4
jittered starts 2.1e-9 above; cost in "Cost"). The review then measured
where quadrature fails on huber() (`dev/reviews/nlstop2-review-out/
huber-*.txt`, `huber2-*.txt`): with `sigma ~ (1 | g)`, "could not
marginalize" on 6 of 10 seeds at each scale; with `huber(link = "log")`
and `(1 | g)`, on 4 of 5; on a quadrature fit `predict()`, `fitted()`,
`residuals()`, `frm_linpred(se.fit = TRUE)` and `avg_slopes()` warn that
group effects are held at their modes, and `vcov_cluster()` refuses.
The user decided against changing a default for a special case: the
route is reverted, and `?huber` states the shortfall, when quadrature
helps, its cost and these limits.

## Part B cluster 3

- 3359, rerun (`dev/nlstop2-rerun.R 3359`, `rerun3359-before.txt`):
  degree 8 seeds 1 to 3, logLik minus the constrained lm -2644, -779.4,
  -2373, code 0, warned (`frmtmb_ill_design`, `frmtmb_conv_gradient`).
  Still open; the fix is a start, lane optstart's code. Filed as a rerun.
- 2818, rerun (`dev/nlstop2-rerun.R 2818`, `rerun2818-before.txt`): does
  not reproduce. "2818: ll -14964.054453 code 0 t 1656.6 warn -" against
  "mgcv::gamm ML logLik -14964.054453". CLOSED in the backlog. The fit
  took 1656.6 s in one process beside the other lanes; filed as a LOW
  cost item, not measured on a quiet machine.
- Not reached, for time: 5626, 5697, 5881, 5576, 5641, 5682, 2590, 2599,
  3013, 3017, 3022. 2590 is in `R/autoscale.R` (lane slopescale's) and
  5641 in the profile code of `R/confint.R` (lane selimit2's area).

## 5718: question for the user (answered: built at the consolidation)

`options(frmtmb.nl_sign_probe = FALSE)` turns off the sign probe; nothing
in `frmtmb_control()` does, so the choice does not travel with the fit
(`update()`, `refit()`, a saved fit). This round adds a second restart
from the default start, the pole restart, which no switch turns off.

Recommendation: one entry, `frmtmb_control(nl_starts = TRUE)`, that
turns off both the sign probe and the pole restart (`FALSE`), with the
option kept as its default (`getOption("frmtmb.nl_sign_probe", TRUE)`)
for one release and then removed. Not built: `frmtmb_control()`'s
signature is in `R/fit.R`, which lane slopescale edits at lines 2732 to
2764, and the brief says to ask.

## Tests

`tests/testthat/test-fit-restart.R`, first pass, 6 new tests: the
bending ridge (explin), the poles (binomial seeds 11 and 44, 4 ulp draws
each), a given start between the poles, the absent case (gaussian seed 55
clear of its walls, and `a / (1 + b * x)` with walls on both sides), and
the huber route (removed in the fix round; final counts there).

- rellib-r21 (`t-fit-restart-before1.txt`): pass 40, fail 18, err 2,
  warn 2. Behavioral: no `frmtmb_nl_limit` on explin; binomial 11 and 44
  4.322 and 0.919 below the start at K = -6, K -0.2 and -0.1, on every
  ulp draw; no `frmtmb_nl_pole` from K = -2.5; huber `quadrature` FALSE
  and a jittered start 0.030 higher, with two warnings escaping (code 8).
  The two errors are missing functions (`nl_pole_walls()` in the absent
  case, `huber_quad_route()`).
- lane, final build (`t-fit-restart-after2.txt`): pass 71, fail 0, err 0,
  skip 0, warn 0 (first pass, with the huber test).

Other test files that fit huber(): test-families.R pass 252,
test-numerical-robustness.R 723, test-simulate-density.R 620,
test-update-pool.R 21, each fail 0, err 0, warn 0
(`t-*-after1.txt`).

## Acceptance sets

Final build, both arms run by this lane (`jobs-acc1.txt` before,
`jobs-acc2.txt` after; queue `dev/nlstop2-queue.sh`, 2 processes),
summary `dev/nlstop2-accsum.R after2 before` (`accsum-after2.txt`). The
fresh rows' reference is the best of 16 or 17 started fits from lane
nlstop's before arm (same data, `frmtmb-wt-nlstop/dev/nlstop-out/`) and
of both arms' defaults; the scan's is lane nlstop's best of 8 starts.

```
== fresh rows
rows 425, default error after 90, before 90
after: short > 1e-4 8 (silent 2), short > 1 4 (silent 1)
before: short > 1e-4 11 (silent 4), short > 1 6 (silent 2)
after - before: higher > 1e-4 3, lower < -1e-4 0, within 332
frmtmb_nl_limit after 2 before 1; frmtmb_nl_pole after 0
== scan rows
rows 450, default NA after 0 before 0
after: short > 1e-4 0 (silent 0), warned rows 100
before: short > 1e-4 0 (silent 0), warned rows 100
after - before: higher > 1e-4 0, lower < -1e-4 0, within 450
warning classes differ on 0 rows
== limit battery
A after2: rows 480 err 120 warned 0 nl_limit 0 nl_pole 0
A before: rows 480 err 120 warned 0 nl_limit 0 nl_pole 0
B after2: rows 480 err 360 warned 0 nl_limit 0 nl_pole 0
B before: rows 480 err 360 warned 0 nl_limit 0 nl_pole 0
Cprior after2: rows 480 err 0 warned 0 nl_limit 0 nl_pole 0
Cprior before: rows 480 err 0 warned 0 nl_limit 0 nl_pole 0
Cstart after2: rows 240 err 0 warned 54 nl_limit 4 nl_pole 0
Cstart before: rows 240 err 0 warned 50 nl_limit 0 nl_pole 0
```

The 3 fresh rows that moved are binomial mm 11, 44 and 55, each now at
the reference (-618.9665, -625.9800, -627.4689). The 2 silent rows left
are cumulative pow 22 (1.03, no denominator; filed again) and bernoulli
sig4 55 (2.8e-4, below grad_tol). The new `frmtmb_nl_limit` on the fresh
rows is binomial sat_z 44, at the reference value -566.708119 with b1
406.9, b2 0.002542: the linear limit of `b1 * (1 - exp(-b2 * x))`, which
the move now reads after the step across. The 4 new warnings of the
battery are the fn_ rows it exists for (fn_explin n 200 seed 5, n 1000
seeds 2 and 10, fn_mmlin n 200 seed 5); no logLik of the battery moved
by more than 1e-6. The bigx grid fits no nonlinear body and no huber()
family, so neither change reaches it (frm() calls `nl_fit_probed()` only
for a body, and `huber_quad_route()` returns FALSE for any other
family); not rerun.

The wall times (fresh 1072.7 s after, 864.1 s before; scan 457.8 and
379.0) come from runs at different times beside other lanes and are not
a measurement of the change; see "Cost".

## Cost

`dev/nlstop2-cost.R`, one process, arms rotated per round, minimum of 3
rounds, "off2" the control that must read 1.

Pole restart. Wall time does not resolve it (`cost-pole.txt`):

```
pole clear: min s on 10.99 off 13.21 off2 11.22 | on/off 0.832 control 0.849
pole pressed: min s on 2.61 off 0.64 off2 1.58 | on/off 4.078 control 2.469
```

The controls move by 15 percent and by a factor of 2.5, so the evaluation
counts of every optimize_obj() run of the fit decide it
(`cost-polecount.txt`, load-independent):

```
polecount binomial mm 11: evaluations on 439 off 331 ratio 1.326
polecount binomial mm 33: evaluations on 502 off 402 ratio 1.249
polecount binomial mm 44: evaluations on 409 off 377 ratio 1.085
gaussian mm 11 to 55, poisson mm 11 to 55: ratio 1.000 on all 10
```

A fit clear of its poles pays 3 evaluations of each denominator per
coefficient and no optimization; a held fit pays 8.5 to 33 percent more
evaluations, plus the standard errors of the refit where the restart
wins (binomial 11 and 44).

Limit check (`cost-limit.txt`, 20 calls on each of 15 fits per block):
on 0.39 s, off (0.82.0's function) 0.37, off2 0.41; on/off 1.054, control
1.108; per call 0.0013 s against a fit's 0.4873 s. At 100 calls per
block (`cost-limit2.txt`): on 2.69 s, off 2.19, off2 2.08; on/off 1.228,
control 0.950; per call 0.0018 s against 0.0015 s, a fit 0.7620 s. About
0.3 ms more per call (2 gradients and 2 objective values), 0.04 percent
of a fit.

huber() route, reverted (`cost-huber.txt`, rerun on an unedited script as
`cost-huber2.txt`); these numbers are the cost `?huber` now gives for
`quadrature = TRUE`:

```
huber 100 rows x 10 seeds: default 4.57 laplace 0.67 laplace2 0.62 | default/laplace 6.82 control 0.925
huber 2000 rows x 3 seeds: default 43.16 laplace 1.11 laplace2 1.01 | default/laplace 38.88 control 0.910
huber 100 rows x 10 seeds: default 4.97 laplace 0.58 laplace2 0.63 | default/laplace 8.57 control 1.086
huber 2000 rows x 3 seeds: default 44.25 laplace 0.91 laplace2 1.01 | default/laplace 48.63 control 1.110
```

At the larger size (`dev/nlstop2-huber-big.R 20`, `huber-big-20.txt`,
100 groups of 20 rows, seeds 1 to 5): every Laplace fit stopped at code 1
with the false-convergence and gradient warnings, 4 of 5 short of the
best of 8 jittered Laplace starts (0.036, 0.153, 0.099, 0.050); every
quadrature fit at code 0 with no warning, at the best of its 2 jittered
starts (1.8e-7 at most), 13.6 to 16.1 s against 0.36 to 1.78 s.

## User decision: the huber() route (answered 2026-11-13: revert)

The first pass asked whether to keep the route (7 to 9 times the fit at
100 rows, 39 to 49 times at 2000; 0 of 45 fits short) and recommended
keeping it. The user decided against changing a default for a special
case. The route is reverted in the fix round and 5065 is a documentation
item.

## OpenBLAS, suites, sessions, check

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`, no
denial; `dev/nlstop2-ob.sh 0.3.26 2`, each log's blas: line naming the
copy, lib: line the lane's core): test-fit-restart.R 69, test-nl-start.R
59, test-nl.R 92, test-nl-family.R 333, test-nlf.R 78, test-families.R
252, test-numerical-robustness.R 723, test-simulate-density.R 620,
test-update-pool.R 21; fail 0, err 0, skip 0, warn 0 on each. (69 against
71 on reference BLAS: two branches OpenBLAS does not take, the explin
test's code-0 branch, removed in the fix round, and the bernoulli mm 55
test's value check at the limit ending, described in "Fix round".) OpenBLAS 0.3.32
with LAPACK: the same 9 files, fail 0, err 0, warn 0 (test-fit-restart.R
69). Fix-round runs are in "Fix round".

Sessions (`dev/nlstop2-sessions.sh`: the 10 files above plus test-map.R,
with the helpers and setup.R, `dev/release/run-session.R`):

```
SESSION frmtmb forward files=10/10 pass=2265 fail=0 err=0 skip=0 warn=0
SESSION frmtmb reverse files=10/10 pass=2265 fail=0 err=0 skip=0 warn=0
SESSION frmtmb preload files=10/10 pass=2265 fail=0 err=0 skip=0 warn=0
```

0 LEAK lines; every lib: line the lane's core, extensions from rellib-r21.

Plain suite, every core and extension file, one per process, 3 at a time
(`dev/release/run-par.sh nlstop2-plain-suite dev/nlstop2-out/suite-plain.jobs
3`, `FRMTMB_RELLIB` the lane's library, extensions from rellib-r21; log
`dev/release/nlstop2-plain-suite.log`, "ran 529 of 529"; every file's
lib: line names the lane's core):

```
              pkg files   pass fail err skip warn
           frmtmb   330  24424    2   0  240    0
  frmtmb.coupling    11    551    0   0    5    0
       frmtmb.eam    33   1851    0   0    5    0
    frmtmb.latent    13    697    0   0    2    0
     frmtmb.learn    16    474    0   0   13    0
       frmtmb.ode    11    549    0   0    1    0
    frmtmb.sample    91   4238    0   0   18    0
  frmtmb.spatial      8    529    0   0    6    0
  frmtmb.survival    16    575    0   0    0    0
total files 529 pass 33888 fail 2 err 0 skip 290 warn 0
```

The 2 failures were test-bracket-access.R: `huber_quad_route()` read
`r$family$family` with `$`, a hazard container. Fixed (`[[`); the file
then reads pass 36, fail 0 (`t-test-bracket-access-after3.txt`), and
test-fit-restart.R pass 71, fail 0 (`t-test-fit-restart-after3.txt`).
The fix changes no value, so the other 528 files stand.

Gated tier, the 55 files with a skip above, every gate on (NOT_CRAN,
FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS, FRMTMB_FUZZ,
FRMTMB_SCALE_TESTS; `dev/release/run-par.sh nlstop2-gated-suite
dev/nlstop2-out/suite-gated.jobs 1`, log
`dev/release/nlstop2-gated-suite.log`): "NLSTOP2-GATED-SUITE ran 55 of
55"; summed from the RESULT lines, files 55 pass 4711 fail 0 err 0 skip 0
warn 0. None of the 55 files mentions huber, so the revert does not
reach them; not rerun.

R CMD check --as-cran of core, first pass, with the huber route
(`dev/nlstop2-check.sh`, log kept as `dev/nlstop2-out/check1.log`):
"Status: 2 NOTEs": examples timing ("checking examples ... [169s] NOTE",
load) and the HTML manual (V8 unavailable). The fix-round check is in
"Fix round". `dev/nlstop2-check/` and `dev/ciharden-out/` were deleted
after reading.

## Bump and floor

Minor: one BREAKING behavior change in fits (the pole restart). No
extension changes; no extension floor moves (the extension suites pass
over this core with their rellib-r21 builds, above).

## Fix round (review dev/reviews/2026-11-13-nlstop2.md, verdict FIX)

User decisions of 2026-11-13: revert the huber() quadrature route (5065
becomes documentation); 5718 is built at the consolidation; 5041 stays
filed open.

### F5: the huber() route, reverted

- `R/fit.R` restored to ec8178d3: `git hash-object R/fit.R` =
  `git rev-parse ec8178d3:R/fit.R` = b68a5bdbd9b9811dd9afd2df4ce1d7e6dcf6779b.
  The lane has no other line in `R/fit.R`. `man/frm.Rd` regenerated equals
  ec8178d3's (36acd8390c6d9ee9714f3d13e86f69f9125c4adb).
- `R/fit-restart.R`: `huber_quad_route()` and its roxygen block removed;
  nothing else called it.
- `tests/testthat/test-fit-restart.R`: `huber_data()`, the huber test and
  the header's "huber()'s Laplace marginal" removed.
- `NEWS.md`: the BREAKING huber() bullet replaced by a `## Documentation`
  bullet.
- `R/families.R`, `huber()` `@details`, the documentation item: Laplace
  with group-level terms can stop short (11 of 20 at each of two response
  scales, by up to 0.188 and 0.272; 6 and 5 silent); `quadrature = TRUE`
  reached the maximum on all 40 fits, at about 7 to 9 times the Laplace
  fit (100 rows, 20 groups) and 39 to 49 times (2000 rows, 100 groups);
  use it for a scalar group-level term in the mean with the identity
  link; with `sigma ~ (1 | g)` it stopped with "could not marginalize" on
  6 of 10 data sets, with `huber(link = "log")` on 4 of 5; on a quadrature
  fit `predict()`, `fitted()`, `residuals()`, `frm_linpred(se.fit =
  TRUE)` and `marginaleffects::avg_slopes()` warn that group effects are
  held at their modes, and `vcov_cluster()` refuses. The last three
  facts are the review's measurements (`dev/reviews/nlstop2-review-out/
  huber-*.txt`, `huber2-*.txt`). `man/frmtmb-families.Rd` regenerated.
- `dev/test-backlog.md`: the 5065 CLOSED line now says documentation
  only, by the user's decision.

### F1: the explin test

The test asserted the limit warning only when the fit ended at code 0;
on the review's draw k = 3 the fit ends at code 0 with the gradient
warning, and `limit()` is skipped there because another warning was
given. The test now asserts, on every ending, a warning and
`expect_false(is.null(nl_limit_message(fit)))`.

`dev/nlstop2-explintest.R` (the new test body on y * (1 + 2e-15 k), k 0
to 11):
- lane, reference BLAS (`explintest-after.txt`): pass 12 of 12, k 3 "code
  0 warn frmtmb_conv_gradient limit_message TRUE test-body PASS".
- rellib-r21 (`explintest-before.txt`): pass 0 of 12 (no warning on the
  code-0 draws, and `nl_limit_message()` NULL on every draw).
- OpenBLAS with LAPACK, `OPENBLAS_CORETYPE=Haswell` (`dev/nlstop2-ob2.sh`;
  copies rebuilt with `dev/ciharden-openblas.sh`, no denial): 0.3.26
  explin draws pass 12 of 12 (`ob2-0.3.26-explintest.txt`);
  test-fit-restart.R pass 58, fail 0, err 0, warn 0 under 0.3.26 and
  0.3.32 (`ob2-*-fit-restart.txt`; 58 against 59 on reference BLAS: the
  one assertion is test-fit-restart.R line 203, inside
  `if (has_class(got$conds, "frmtmb_nl_limit"))` in the bernoulli mm 55
  test, unchanged from ec8178d3. It checks the value at the limit ending
  (17 of 20 draws on reference BLAS); OpenBLAS takes the other valid
  ending, 0.0065 below with the standard-error warning; both endings get
  the same unconditional assertions). 0.3.32 explin draws
  (`ob2-0.3.32-explintest.txt`): pass 12 of 12 (1 at code 0 with the limit
  warning, 11 at code 1 with the optimizer's warnings).

### F2, F3, F4: backlog

- F2: "cumulative pow 22" filed (LOW, 1.03 short, code 0, silent, no
  denominator).
- F3: the 5093 CLOSED line and the 5093 section now say the Hill fit has
  no limit (h to infinity is 0.147 lower) and stops 1.76e-3 short of h
  doubled, silent; filed LOW.
- F4: the review's `a / (x + b)` band case filed HIGH (591 to 646 short,
  silent; the escape starts only past the ends of the band). Not fixed.

### Fix-round runs (lane build, reference BLAS)

test-fit-restart.R pass 59 (71 less the 12 huber assertions), fail 0;
test-bracket-access.R 36; test-families.R 252; test-numerical-robustness.R
723; test-simulate-density.R 620; test-update-pool.R 21; each fail 0, err
0, warn 0 (`t-*-after4.txt`). rellib-r21 on the final test file
(`t-fit-restart-before3.txt`): pass 35, fail 17, err 1.

R CMD check --as-cran of core on the fix-round tree
(`dev/nlstop2-check.sh`, log kept as `dev/nlstop2-out/check2.log`):
"Status: 2 NOTEs", examples timing ("checking examples ... [153s] NOTE",
machine load) and the HTML manual (V8 unavailable); tests "[66m] OK".
`dev/nlstop2-check/` and `dev/ciharden-out/` deleted after reading.

Short check of the fix round (`dev/reviews/2026-11-13-nlstop2-short.md`):
SHIP after record edits. Edits: the two OpenBLAS count explanations
above; "6 of 10 data sets at each of two scales" in `?huber`; "on all 40
measured fits" in NEWS; `man/frmtmb-families.Rd` regenerated.
