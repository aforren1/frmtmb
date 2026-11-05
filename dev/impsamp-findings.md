# Lane impsamp: `importance =` converges on the models it exists for

Round opened 2026-10-29, base d5d0056d (frmtmb 0.78.0). Scope: D3 of
`dev/equiv-findings.md`. Status: fix in the worktree; see the last sections
for what ran. Logs are local, under `dev/impsamp-out/` (gitignored); every
number below was generated into one of them and is quoted with its script.

## 1. Diagnosis

The weights were never the problem. At p124's data (seed 124), 200 draws:
ESS min 0.998, sd of the log weights 0.03 to 0.07, MCSE 0.004 to 0.007
(`dev/impsamp-diag2.R`, `dev/impsamp-out/diag2.txt`, draw seeds 1:20, base
build). The draws are common random numbers: one seed per fit, the same z
every round, so the objective is smooth. The proposal is the Laplace
Gaussian (`spHess(random = TRUE)` and the inner mode) with antithetic pairs.
The outer optimizer sees a smooth deterministic function. None of those is
the cause.

The cause is the frozen proposal's CURVATURE. `u_i = u_hat(a) + L(a) z_i` does
not move while the optimizer does, so a step d in an intercept moves each
group's integrand by d under draws that stay put. The second derivative of
`-log mean exp(...)` along it is, per group,
`I_g (1 - mean z^2) + marginal information`, with `I_g` the conditional
information (43 to 53 here, against 3.7 for the prior precision). That is a
sample statistic many times larger than the marginal curvature it has to
resolve. The iteration "minimize on the frozen tape, re-anchor" is then a
Newton iteration on an accurate gradient with a noisy Hessian `Hn`, and its
growth factor is `rho = max |1 - eig(Hn^-1 H_exact)|`.

Measured at the Laplace optimum, base build, 20 draw seeds:

| draws | `H[1,1]` exact | sd of `H[1,1]` over seeds | min `H[1,1]` | seeds with rho > 1 | median rho |
|---|---|---|---|---|---|
| 200 | 34.328 | 21.2 | -3.411 | 5 of 20 | 0.472 |
| 512 | 34.328 | 12.9 | 4.839 | 4 of 20 | 0.212 |
| 1000 | 34.328 | 11.5 | 11.34 | 2 of 20 | 0.248 |
| 4000 | 34.328 | 3.9 | 24.99 | 0 of 20 | 0.084 |

The theory formula `sum_g I_g (1 - mean z_g^2) + H_exact` correlates 0.9998 to
0.9999 with the measured `H[1,1]`, and the antithetic pairs double the noise of
`mean z^2`, because a pair shares its `z^2`: `I sqrt(4 / N) sqrt(G)` predicts
21, 13 and 9.4 at 200, 512 and 1000 draws.

The default seed is an unlucky one, and the measured growth matches: seed 1
has rho 5.84, 5.47 and 2.06 at 200, 512 and 1000 draws, and the verbose rounds
(`dev/impsamp-diag1.R`, `diag1-base.txt`) move by a factor of 4 to 5 per round
at 200 and 512 and by 2.1 at 1000: 0.00108, 0.00227, 0.00477, 0.01. So the 1000
in the survey did not fit because 1000 is enough. It fit because five rounds
were too few to run away, and it ended with a capped warning.

p069 is the same mechanism in a harder place: the correlation sits on its
boundary (chart value 335), the prior is nearly singular, and the frozen draws
are hypersensitive to the correlation's chart: a round reached 318.4 on its
tape, 2 nats below the honest 320.3, by aligning that narrow direction with
the draws.

**The refusal's text was wrong.** It blamed an unidentified covariance. The
design its own test pins as such (test-importance.R, 40 groups of 4 binary
rows, 2 x 2 covariance, seed 11) now settles on the exact maximum at every
draw count from 4 up: -102.15312 at 2000 draws against an independent
reference of -102.15295 (`dev/impsamp-hard.R`, `hard-v4.txt`).

Also found on the way, separately (section 5): the numeric path read the wrong
model under a constant dpar.

## 2. The fix

`R/importance.R`. The proposal follows the optimizer's theta to first order
around the anchor; the standard normals stay fixed. For theta:

- `H(theta) = Hlik(a) * exp(rate' dtheta) + P(theta)`, with the prior's
  precision `P` exact (from the registry, on the tape) and the likelihood's
  curvature scaled at a rate measured at the anchor.
- `m(theta) = uhat(a) - H(theta)^-1 [ (P(theta) - P(a)) uhat(a) + cross dtheta
  ]`, one Newton step from the anchor's mode.
- `u_i(theta) = m(theta) + R(theta)^-1 z_i`, `log|L|` on the tape.

`cross` and `rate` are central differences at the anchor of the AD gradient
and the sparse Hessian (`imp_anchor_slopes()`), and nothing has to be exact:
a proposal chosen before the draws are read leaves the estimate unbiased at
each theta, so an error costs efficiency and never validity. Three design
points, each measured:

- **Mode only is not enough** (build v1, `diag1-v1.txt`): p124 settled, p069
  still walked the correlation chart by 20 to 80 per round and capped, its
  rounds reaching 318.3 on the tape against an honest 320.3.
- **The prior's precision on the tape** (build v2, grid arm v2) settled p069
  and took the grid from 162 rises to 1, a gaussian slope at 64 draws, and left
  a fixed point with a Monte Carlo bias in the dispersion on gaussian fits:
  -0.8564 against -0.8788 in log sigma at 200 draws.
- **The curvature's rate must be the TOTAL derivative along the moved mode**.
  A step in a log-link intercept moves the mode by minus that step, so the
  predictor at the mode and the curvature there stay put. Taking the partial
  derivative (build v3) halved the tape's `H[1,1]` on a Poisson slope, 4.1
  against 8.1, and left 2 of 648 rising and 5 capped. With the total derivative
  (v4, final) the gaussian fixed point is the exact optimum at any draw count,
  and the tape's Hessian is the exact one to about 1 percent.

Why a scale for the curvature and not a matrix: it stays positive semi-definite
at any theta, so no optimizer step can leave the Gaussian undefined. A negative
direction of `Hlik` at the anchor (a likelihood that is not log-concave) is cut
off by `imp_psd_part()`, and a log-concave likelihood keeps the Laplace Hessian
to the last bit.

After the fix, at p124, 20 draw seeds (`diag2-v4.txt`):

| draws | sd of `H[1,1]` over seeds | seeds with rho > 1 | median rho |
|---|---|---|---|
| 200 | 0.004 | 0 of 20 | 0.000 |
| 512 | 0.003 | 0 of 20 | 0.000 |
| 1000 | 0.002 | 0 of 20 | 0.000 |

## 3. The grid, before and after

`dev/impsamp-grid.R` (`dev/impsamp-agh.R` holds the cells, the data and the
reference; `dev/impsamp-grid-summary.R`; logs `grid-summary-final.txt`). Six
families (binomial, poisson, bernoulli, negbinomial, gaussian, Gamma), three
shapes (`(1 | g)`, `(1 + x | g)`, `(1 | g) + (0 + x | g)`), three sizes (10
groups of 20, 40 of 5, 30 of 3), data seeds 1:3, draws 64, 200, 512, 2000, draw
seed 1: 648 fits, same cells and seeds in each arm. Base is rellib-r16, after
is the worktree.

Two crossed blocks are not in the grid: `importance =` refuses blocks over
different factors, by name and with its reason (the integral does not split
per group). That refusal stands; see section 8.

The reference is an own adaptive Gauss-Hermite written in R (25 nodes for
one dimension, 11 per dimension for two, `optim()` on top), sharing no code
with the package. Against `quadrature = TRUE` on the 54 intercept cells it
differs by at most 9.78e-05 (the optim tolerance) at the maximum, and at a fit's
own estimates to all ten digits printed by `dev/impsamp-agh-check.R`.

| | base | after |
|---|---|---|
| fits ok | 486 of 648 | 648 of 648 |
| stops "ROSE" | 162 | 0 |
| capped at 5 rounds | 204 | 0 |
| median rounds (paired ok fits) | 4 | 1 |
| median abs(logLik - exact max) (paired) | 0.01202 | 0.00664 |
| median max abs(beta - exact) (paired) | 0.00410 | 0.00021 |
| median max abs(log sd - exact) (paired) | 0.01714 | 0.00103 |
| gaussian: max abs(logLik - exact) | 0.0894 | 9.9e-10 |
| total fit seconds | 5372 (rerun; 4509 first run) | 2680 |
| ok to not ok | | 0 |

After, by draw count, median and max of abs(logLik - exact max), against
Laplace's median gap of 0.09546 over the same cells:

| draws | median | max | median reported MCSE | min ESS |
|---|---|---|---|---|
| 64 | 0.018949 | 0.28518 | 0.0547 | 0.268 |
| 200 | 0.012314 | 0.06939 | 0.0306 | 0.458 |
| 512 | 0.004751 | 0.06857 | 0.0181 | 0.587 |
| 2000 | 0.001712 | 0.02828 | 0.0095 | 0.742 |

The fit approaches the exact value at the rate the draws buy. Rounds on the 648
(`docnumbers.txt`): one in 398, two in 237, three in 12, four in 1. The rise
the guard reads runs from -0.365 to +0.0784 and never passes 0.15 of its
tolerance. The worst group's ESS across the fits: 0.27 at the minimum, median
0.97.

19 of the 648 carry a warning, none about the correction: 17 say a standard
error of a theta is not available at a boundary, and 2 are `nlminb`'s "false
convergence".

## 4. The Monte Carlo error the fit reports

`dev/impsamp-cover.R`, `cover-v4.txt`: 40 draw seeds, exact value from the
reference at the fit's own estimates.

| design | draws | mean err | sd err | median MCSE | sd / MCSE | cover 95 |
|---|---|---|---|---|---|---|
| p124 binomial | 64 | +0.0003 | 0.0032 | 0.0080 | 0.40 | 1.00 |
| p124 binomial | 200 | +0.0003 | 0.0017 | 0.0044 | 0.38 | 1.00 |
| p124 binomial | 512 | +0.0002 | 0.0011 | 0.0030 | 0.35 | 1.00 |
| hard bernoulli (sd 2, 40 x 3) | 64 | +0.0046 | 0.2318 | 0.2026 | 1.14 | 0.95 |
| hard bernoulli | 200 | +0.0269 | 0.1636 | 0.1172 | 1.40 | 0.97 |
| hard bernoulli | 512 | +0.0178 | 0.1080 | 0.0790 | 1.37 | 1.00 |
| poisson slope 40 x 5 | 64 | -0.0039 | 0.0811 | 0.2074 | 0.39 | 1.00 |
| poisson slope 40 x 5 | 200 | +0.0086 | 0.0839 | 0.1202 | 0.70 | 1.00 |
| poisson slope 40 x 5 | 512 | +0.0004 | 0.0440 | 0.0759 | 0.58 | 1.00 |

`imp_ess()`'s MCSE counts draws as independent. On an ordinary design it is
conservative, by a factor of 1.4 to 2.9, because the pairs cancel the odd part
of the weights. An MCSE that treats a pair as the unit is close on ordinary
designs (sd / MCSE 0.86 to 1.70) but covers 0.88 to 1.00 (median 0.92) where
the iid one covers 0.95 to 1.00. I kept the iid one: it is the safe side, and
it is the number every existing reader knows, and `vignette("diagnostics")`
now says it reads high on an ordinary design. The hard row's tail is real: see
section 6.

## 5. Constant dpar: the numeric path read another model

`imp_par_list()` recycled the free `betad` entries over the entries a map
fixes, so with `bf(y ~ x + (1 | g), nu = 5)` the offsets, the ESS and the MCSE
read `nu` from `sigma`'s value (`dev/impsamp-map.R`, seed 31, 200 draws):

| | betad the numeric path reads | numeric objective | tape objective | ESS min | MCSE |
|---|---|---|---|---|---|
| base | -0.7125, -0.7125 | 135.2524 | 126.0120 | 0.503 | 0.102 |
| after | -0.7126, 1.3863 | 126.0107 | 126.0107 | 0.980 | 0.0149 |

The fit's own objective was right. The defect dates from the first version.
`imp_par_list()` now takes `map`; the callers pass it.

## 6. Gaussian or t

`dev/impsamp-proposal.R`, `proposal.txt`: at the estimates of a 2000-draw fit,
with the reference's own modes and Hessians, 512 antithetic draws, 30 seeds,
error of the summed `log m_g` against adaptive quadrature.

| design | Gaussian sd | t16 sd | t4 sd | 0.25 t4 + 0.75 Gaussian sd | Gaussian worst ESS |
|---|---|---|---|---|---|
| p124 binomial | 0.0029 | 0.0499 | 0.1373 | 0.0521 | 0.994 |
| hard bernoulli | 0.6683 | 0.0359 | 0.4372 | 0.1001 | 0.171 |
| bernoulli slope, 30 x 2 | 0.5724 | 0.2613 | 0.4541 | 0.2222 | 0.033 |
| poisson slope, 30 x 3 | 0.1205 | 0.2214 | 0.5458 | 0.1923 | 0.760 |
| gaussian, 30 x 3 | 0.0000 | 0.1548 | 0.4147 | 0.1602 | 1.000 |

Kept the Gaussian. A t, or a mixture with one, is 9 to 47 times worse on the
ordinary designs, because the chi-squared scale of a draw is not paired away,
and it gives up the gaussian identity (weights exactly equal). The Gaussian
is the one with a heavy tail on the hard binary designs, where a group's
integrand has the prior's tail and the proposal is narrower: there its error
sd is a rare seed's, and the ESS warning fires on that seed (0.171, 0.033). A
draw-count-free fix would be to switch to a mixture only when the ESS asks;
filed in section 8, not done.

## 7. Refusal and guard

The stop that said ROSE is replaced by `imp_refuse_draws()`. It fires on the
same condition (the corrected objective at the final anchor exceeds the
objective at the Laplace estimates by more than `max(3 MCSE, 3 MCSE, 0.1)`),
and runs the same iteration at up to four doublings, up to 8192 draws and 2e6
rows times draws. A count counts as settled only if it neither rose nor capped.
The message says what was tried and no more (section 15, M1, M2).

Where it must not fire: it fires on none of the 648 fits and none of the
review's off-grid fits (section 15). Where it fires on its own: `bf(y ~ x + (x |
g))`, bernoulli, 20 groups of 4, a slope near separation, at 2 to 64 draws: 12
of 60 fits at the five draw counts (`dev/impsamp-natural.R`, `natural.txt`). I
first wrote that I found no natural case. That was wrong: the review found one
through B1, and this one needs no B1. After fix round 2 it is gone (section 17).

## 8. Not done, and defects filed

- **Crossed factors under `importance =`** stay refused (the brief's "two
  crossed blocks" shape), by name: the integral does not split per group, and
  that is a different estimator, not this fix. The grid measures two blocks
  over ONE factor instead.
- **A heavy-tail guard for hard binary designs** (section 6). A mixture proposal
  used only when the first anchor's ESS is low would remove the rare-seed error
  at no cost to ordinary designs. Measured tail: sd 0.67 at 512 draws. It is a
  feature and I did not build it.
- **`x + 1e5`, poisson**: false convergence, NaN SEs (`scales.txt`). The
  tape's Hessian has a condition number of 9e18 there, past double precision;
  `+ 1e4` is 1.0000. It warns, and the Laplace fit, which holds its Hessian in
  the design's QR coordinates, reads finite SEs. Not fixed.
- **The MCSE** could treat pairs as the unit (section 4); left as is, with the
  reason.
- `extensions/frmtmb.learn/dev/learn-importance.R` reads `fixef()` and
  `VarCorr()` in shapes older than 0.71.0 and stops; `dev/impsamp-learn.R` is
  a corrected copy, and I left the script alone.

## 9. frmtmb.learn, rerun

`?frmtmb.learn`, `vignette("learning")` and the compat row quoted numbers that
the frozen proposal produced, so they are measured again (`dev/impsamp-learn.R`,
`learn-importance.txt`, seeds 4301 to 4306, 40 subjects, 100 draws, truth 0.5):

| trials | seed | sd Laplace | sd corrected | log shift | mcse | min ESS | secs |
|---|---|---|---|---|---|---|---|
| 100 | 4301 | 0.0004 | 0.0004 | +0.000 | 0.000 | 1.000 | 12.7 |
| 100 | 4302 | 0.4170 | 0.4304 | +0.032 | 0.038 | 0.954 | 15.8 |
| 100 | 4303 | 0.3807 | 0.3931 | +0.032 | 0.034 | 0.967 | 16.8 |
| 100 | 4304 | 0.4543 | 0.4732 | +0.041 | 0.056 | 0.906 | 21.5 |
| 100 | 4305 | 0.3798 | 0.3987 | +0.048 | 0.038 | 0.959 | 25.7 |
| 100 | 4306 | 0.3078 | 0.3154 | +0.024 | 0.024 | 0.988 | 23.5 |
| 20 | 4301 | 0.3707 | 0.4984 | +0.296 | 0.044 | 0.964 | 3.8 |

The other five at 20 trials stay at 0.0001 with ESS 1.000. Six of six
complete, where the old record has two of six refused and five walking by
0.3645 a round. The old shift of +0.11 to +0.18 at 100 trials was noise in a
frozen proposal, not Laplace error. The new one holds as draws grow
(`dev/impsamp-learn-ref.R`, `learn-ref.txt`): seed 4302, 0.4304, 0.4279,
0.4272 at 100, 400, 1600 draws; seed 4305, 0.3987, 0.3958, 0.3950; seed 4301
at 20 trials, 0.4984, 0.4720, 0.4693. Quadrature refuses this family, so no
exact reference exists here. That vignette claim changes meaning: at 100
trials the downward bias is NOT Laplace error. `dev/learn2-findings.md` keeps
its record with a note at the top.

## 10. Equivalence pairs and compat rows

`FRMTMB_EQUIV_LIB=... bash dev/equiv-run.sh 2 dev/impsamp-equiv-out`, p124 and
p069 remeasured, the other 170 copied from the r16 baseline. `dev/impsamp-out/
equiv-verdicts.txt`:

| pair | | baseline | this build |
|---|---|---|---|
| p124 binomial x importance | measured verdict | defect (frm stops) | **equivalent**, L1+L4 |
| | L1 constant, gradient | none | 8.0e-13, 3.1e-09 |
| | L4 fitted-mean max z (Laplace control) | none | 0.13 (0.13) |
| p069 us x importance | measured verdict | defect (frm stops) | **defect** by the flow, L1+L4 |
| | L1 constant, gradient | none | -1.9e-10, 1.5e-10 |
| | L4 max z (Laplace control) | none | 0.627 (0.627) |

**p069: agree, with weak evidence, by the review's ruling** (coordinator's
decision; `dev/reviews/2026-11-02-impsamp.md`, "p069: the verdict"). frm()
fits, L1 is an identity, and L4 form (ii) reads a max z of 0.627 over its 0.5
bound, the SAME 0.627 the Laplace control reads. Over 11 seeds of this design
the correction's reading is the control's within 0.022, so the route has no
power on this pair in either direction. My first word for the cause, a
boundary, was too strong: 3 of 6 boundary seeds read above 0.5, and 0 of 5 off
it, so the boundary is associated with it and neither necessary nor
sufficient. The evidence is the exact marginal: the reference's maximum is
-320.26470, and the fits read -320.26810, -320.26509 and -320.26381 at 200,
512 and 2000 draws (MCSE 0.0211, 0.0138, 0.0072). The flow's measured verdict
stays "defect" until consolidation applies the ruling.

**The adjudication table holds both pairs as defects** (`dev/equiv-data/
adjudication.tsv` lines 11 and 12, which I must not edit), so results.tsv's
final `verdict` column still says defect for both. Drop both rows at
consolidation (and give p069 its reason), or the record rerun will keep them.

Compat rows, from these measurements:

| row | status | change |
|---|---|---|
| `binomial` x `importance` | works | NEW explicit row with the measured note; was inherited from `kind:family` x `kind:mode`, same status |
| `group:importance_blocks` x `importance` (reads for `us`) | conditional | status kept (the condition, one grouping factor, is real); a measured sentence added to the note |
| frmtmb.learn `<family>` x `importance` | works | status kept; the note's measurements replaced |

So no status changes: the proposals `p124 works -> broken` and `p069
conditional -> broken` in `compat-proposals.tsv` are overtaken and should not
be applied.

## 11. Profile and run time

`dev/impsamp-prof2.R`, both builds at once, `prof2-before.txt`,
`prof2-after.txt`, minimum over blocks over 1.2 s, two passes; control (fixed
arithmetic) 3.62 to 4.49 s before and 3.75 to 4.49 s after.

| case | N | gr before (s) | gr after (s) | tape before (s) | tape after (s) | R peak before (Mb) | R peak after (Mb) |
|---|---|---|---|---|---|---|---|
| binomial, `(1 \| g)`, 200 rows | 200 | 0.0197 to 0.0226 | 0.0235 to 0.0308 | 0.11 to 0.17 | 0.14 to 0.17 | 182 to 185 | 184 to 191 |
| | 1000 | 0.121 to 0.124 | 0.122 to 0.151 | 0.70 to 0.88 | 0.72 to 0.91 | 212 to 228 | 211 to 230 |
| poisson, `(1 + x \| g)`, 200 rows | 200 | 0.0061 to 0.0081 | 0.0068 to 0.0069 | 0.08 to 0.09 | 0.10 to 0.12 | 181 | 196 |
| | 1000 | 0.033 to 0.036 | 0.030 to 0.035 | 0.42 to 0.54 | 0.50 to 0.75 | 215 to 217 | 221 to 240 |
| poisson slope, 1000 groups of 5 | 64 | 0.063 to 0.068 | 0.069 to 0.078 | 0.89 to 1.25 | 1.45 to 1.58 | 294 to 295 | 278 to 284 |

The load-independent count: tape operators on the binomial, 173835 to 178134
at N = 200 and 866655 to 886934 at N = 1000, +2.5 and +2.3 percent. These small
designs read about the same per evaluation. That does not carry to larger ones
(section 15, m2). A whole fit, poisson slope,
512 draws, Rprof with memory: 3.27 s and 206.7 Mb over 5 rounds before; 1.16 s
and 141.3 Mb over 1 round after. In the before arm the optimizer is 61 percent
of the time and the freezes 25 percent; in the after arm the freezes are 41
percent, the SE check 35 percent (a Hessian it now gets to compute), and the
evaluations 37.

`dev/impsamp-testfile.R`, both builds at once, `times-before.tsv`,
`times-after.tsv`, the 36 tests of test-importance.R common to both arms:
229.5 s before, 189.8 s after. The whole file: 325.2 s (38 tests) before,
336.1 s (46 tests) after. The other files that use `importance =` are in
section 12.

## 9b. The core vignettes' numbers, rerun

`dev/impsamp-rl.R`, `rl.txt`, seeds 4201 to 4206, 100 draws, the correlated
`(1 | p | id)` model of `vignette("reinforcement-learning")`:

| | 40 x 100 trials | 40 x 20 trials |
|---|---|---|
| replicates completed | 6 of 6 (was 4 of 6) | 6 of 6 (was 5 of 6) |
| mean shift, `log sd(alpha)` | +0.121 (was +0.357) | +0.451 (was +0.480) |
| largest single shift | 0.207 (was 0.545) | 0.914 (was 1.386) |
| Monte Carlo standard error | 0.11 to 0.15 | 0.02 to 0.83 |
| smallest ESS per draw | 0.69 | 0.03 |

I did not have the old table's seeds, so these are a new draw of the same
design, recorded here. Dataset 4206 at 20 trials, whose text said the
correction refused at 400 draws: -0.2879 at 100 draws (MCSE 0.162, ESS 0.63)
and -0.6501 at 400 (MCSE 0.036, ESS 0.96), a sd of 0.52 against a truth of
0.5, Laplace -1.2024. One fit of the twelve at 100 draws is capped at 5 rounds
(seed 4204, 20 trials, ESS 0.11), and that fit warns. Two of the six at 20
trials are below the ESS floor at 100 draws.

Exams_2 (`dev/impsamp-exams.R`, `exams.txt`, integrate() over each subject):
the vignette said the Laplace log-likelihood is off by 0.013, more than its
noise. That was the frozen proposal. Laplace -44.09545, corrected -44.09566
(MCSE 0.00068), exact at the corrected estimates -44.09563, at Laplace's
-44.09563. The true error is 0.0002.

The other files that use `importance =`, each pair at once
(`dev/impsamp-times.sh`, `times-others.txt`), seconds before and after; a file
whose test the lane rewrote runs its original on the base build:
test-grad-verdict 151.9, 42.6; test-se-exact-hessian 6.9, 5.8; test-rl-example
15.0, 13.8; frmtmb.learn test-rlddm-ndt 18.3, 14.0; test-session 4.3, 3.9;
frmtmb.spatial test-car 10.4, 10.3. No file regresses.

## 12. Tests

New, in `tests/testthat/test-importance.R`, seen to fail on rellib-r16
(`dev/impsamp-out/t-importance-newtests-on-r16.txt`, RESULT pass=226 fail=14
error=6 warn=4 on the new file):

| test | failure on rellib-r16 |
|---|---|
| the textbook binomial settles at the draw counts users pick | stops: "ROSE from 418.11055 ... to 418.27218 after 5 rounds" |
| a correlated slope at a boundary settles at the draw counts | stops: "ROSE from 320.26599 ... to 320.38537 after 2 rounds" |
| the objective's curvature along an intercept is not a sample | ratio off by 0.85, at seeds 1 and 14 (behavioral) |
| a gaussian fit stays at the exact optimum at any draw count | stops: "ROSE from 108.66084 ... to 109.33497 after 5 rounds" |
| a collapsed component settles where Laplace put it | capped, "moved by the same amount, 0.909" |
| too few draws are named with their ESS and the count that works | old message: fails 5 expectations on its text |
| no count that works is not reported as one | old message |
| a constant dpar's numeric path reads the model the tape reads | weak form on the new test (the third argument is new); the behavioral failure is `dev/impsamp-map.R`, numeric 135.25 against tape 126.01 |
| a negative likelihood curvature is cut off, a positive one kept | new helper, missing symbol (weak form); behavior reached by every fit |

Rewritten because the fix removes what they constructed, with what each
asserts now:

| test | was | is |
|---|---|---|
| the effective sample size diagnostic separates the two regimes | half a log sd of displacement drops ESS | a real design whose group holds 0.019 and warns; displacement no longer degrades it, because the proposal follows |
| a diverging iteration refuses instead of reporting nonsense | expected "ROSE from" | renamed "a weakly identified covariance settles on the exact value": settles, within 3 MCSE of GLMMadaptive |
| a profile bound in the degenerate regime warns | theta_1's profile warned | renamed "a profile bound warns only where the weights fail": silent and within 0.01 of the interval of quadrature's; the warning asserted at a bound past where it covers (ESS 0.069) and absent at the real bounds |
| (1 + x \| g) plus a sigma block: gradient and coverage hold | numDeriv defaults | a smaller Richardson step; both defaults and AD measured (`tune.txt`): defaults carry 3e-04 of their own error |
| a correction that is still shrinking still asks for rounds | cap of 2 | cap of 1: moves 0.0844, 0.000266 |
| the rise tolerance separates drift from divergence | numbers of the old rises | the floor, and a real fit's rise against its tolerance |
| an importance fit over grad_tol does not print the clean line (test-grad-verdict.R) | a draw count that landed over `grad_tol`; it skips on this build | `grad_tol = 1e-9`: the same assertion without the skip |
| the importance correction runs on a learn fit (frmtmb.learn test-factorization.R) | warned of the same step | settles |

Test count of test-importance.R: 38 tests before, 55 after (`tests-old.txt`,
`tests-new.txt`). 19 new, 2 renamed, 2 dropped with their branch ("a stalled
correction is not told to raise the round count", "imp_stalled() separates the
two measured regimes"; the branch is gone, section 15, m5). 4 of the rewrites
above changed again in fix round 1 (m3, m4).

## 13. Overlaps

`R/importance.R` is mine alone. `NEWS.md`: a new top section. `R/compat.R`: a
note and a new line under the importance group. `R/fit.R`: `@param` text only,
with `man/frm.Rd` and `man/frmtmb_control.Rd`. Lane latres adds a covstruct
`latres` with `levels_coupled = TRUE` (read-only look at its worktree):
`imp_crosslevel()` reads that flag and refuses it, and its test pins the
refusal, so the two do not meet. If latres ever drops the flag or registers a
block with a level-wise density that is not Gaussian, `imp_block_prior()`
would read a wrong precision, as before this round.

## 14. What ran, on the final build

One file per process, on this library, plain tier. Core 286 files: pass 20534,
fail 0, error 0, skip 206, warn 0 (test-importance 314, test-grad-verdict 74).
learn 429 pass, 13 skip; spatial 478, 6; sample 3767, 15; latent 681, 2; eam
1838, 5; each fail 0, error 0, warn 0. The skips are the gated tier. coupling,
ode and survival did not rerun in this round: they do not call `importance =`,
and they passed before it (551, 549, 552) on the build that added the compat
row.

Gated files that use `importance =`, with `FRMTMB_BRMS_FIT_TESTS=true`
(`dev/impsamp-final-b.sh`, `gated/`), all fail 0, error 0, skip 0, warn 0:
test-me 71, test-rl-example 102, test-importance 314, test-bcm-simple 13,
frmtmb.spatial test-brms-likelihood-car 67, frmtmb.sample test-loo 146.

One session (`sessions-impsamp.log`), forward, reverse, preload: core pass
20534, 20684, 20777; learn 429; sample 3767; spatial 478; latent 681. All fail
0, err 0, warn 0, 15 of 15 sessions, 37 LEAK lines, every one tagged "[with
load of ...]". Core 276 of 286 and sample 73 of 78 files: the rest source to
nothing outside the gated tier, as before.

OpenBLAS with LAPACK, 4 threads (`dev/ciharden-scan.sh`, list `ob-list.txt`),
on the final code. 0.3.26 and 0.3.32 each 10 of 10, each fail 0, err 0, warn
0: test-importance 314, test-grad-verdict 74, test-se-exact-hessian 35,
test-rl-example 102, test-covstruct-seam 131, test-factorization 34,
test-rlddm-ndt 53, test-session 22, test-surface 146, test-car 132. The R
copies are deleted.

The grid on the final build (`grid-summary-final.txt`): 648 of 648, 0 stops, 0
capped, rounds as in section 3, gaussian max gap 9.9e-10. Against base,
paired: median abs logLik error 0.01202 to 0.00650, beta 0.00410 to 0.00021,
log sd 0.01714 to 0.00102. Total seconds 4722 under this load, against 5372
and 4509 on base and 2680 for the same code earlier: time on this machine
reads load, which is why NEWS gives a range.

The review's off-grid designs (`dev/impsamp-offgrid.sh`, its scripts and
designs, `offgrid-summary.txt`): 116 fits of 28 designs, against the exact
optimum: largest max z 0.148, Newton decrement at most 0.153, SE ratio 0.90 to
1.08, Laplace's max z 3.63. The same as the review's on these designs. The
0.90 is its noted `nb_mushape` seed 2. 6 fits warned. The refusal fired on
none. The review's other 9 fits were reruns of these designs.

The year table (`dev/impsamp-year.R`, `year-after.log`), `y ~ year + (1 | g)`,
`year = round(2000 + 5 z)`, 15 groups of 10, 300 draws, data seeds 1:6,
against the centered column:

| poisson, seed | base | review's build | this build: rounds, SE ratio, slope shift (SE) |
|---|---|---|---|
| 1 | capped | capped, 1.47 | 2, 1.00, 0.0002 |
| 2 | stop | capped, 1.23 | 2, 1.00, 0.0005 |
| 3 | stop | capped, 1.22 | 2, 1.00, -0.0002 |
| 4 | stop | stop | 2, 1.00, 0.0009 |
| 5 | stop | capped, 1.47 | 2, 1.00, 0.0001 |
| 6 | stop | stop | 2, 1.00, -0.0005 |

Bernoulli, all six: rounds 2, no warning, SE ratio 1.00, slope shift at most
0.005 SE, fresh gradient 0.003 to 0.014. The review's build read 1.02 to 1.09,
silently.

Extreme forms (`dev/review-impsamp-scale-se.R`, `scale-se.txt`), one data
seed, SE ratio slope and theta 1.000 to 1.0002 on every row:

| column | family | result | the review's build |
|---|---|---|---|
| `1e4 x` | poisson | rounds 2, false convergence warning | silent, gradient 19.3 |
| `1e5 x` | poisson | rounds 2, false convergence warning, fit moved, SEs right | the Laplace optimum, NaN SEs |
| `x + 2000` | bernoulli | rounds 2, silent | silent, SE ratio 1.135 |
| `x + 2000` | poisson | rounds 2, silent | capped, moved 55.4 |
| `100 x + 1e4` | both | rounds 2, silent | stop |

Equivalence (`dev/impsamp-verdicts.R`, `equiv-verdicts.txt`), against the r16
baseline: p124 measured verdict defect to **equivalent**, L1+L4 (constant
8.0e-13, max z 0.13, control 0.13). p069 frm() fits, L1 identity (-1.9e-10),
L4 0.627, control 0.627, measured verdict stays defect by the flow; by the
review's ruling it is agree with weak evidence (section 10). No other verdict
moves. `dev/equiv-data/adjudication .tsv` still holds both rows as defects:
consolidation drops them.

## 15. Fix round 1: what the review changed

Review `dev/reviews/2026-11-02-impsamp.md`. Numbers are from the final build
unless a line says otherwise.

**B1, the slopes' step.** `imp_outer_steps()` gives each outer parameter a
step in units of the predictor it moves: `1e-4 u / max |X_j|`, with `u` the
response's unit on an identity link, for `beta` and for `betad` through its
own column. A parameter with no column (theta, a family's own) keeps `1e-4
max(1, |par|)`. Nlpars cannot reopen it: `importance =` refuses a nonlinear
predictor. An offset moves no column, and `dev/impsamp-scales.R` reads SE
ratio 1.0000 on `offset(log 1000)`. The fix is in the step only: the review's
mutant `fix_x` covered `beta`; `betad` is covered here and has its own test
(`beta_only` fails it).

`moved < imp_move_tol` stays in coefficient units, so a column at `1e5 x`
stops on the move earlier than `x` does and a column at `1e-4 x` later. All
read 1.0000 on the SE (`scales.txt`), so this changes rounds and not
estimates.

**M1, the refusal's text.** `settled <- !rose && !capped`. The message names a
count only if that fit settled. Otherwise it lists the counts it tried, and
says "No count above k was tried", with the bound's reason if the bound
stopped it. It no longer says "more draws alone are not the remedy" at all:
nothing it runs shows that. Natural firing after B1 (`dev/impsamp-natural.R`,
7 designs, data seeds 1:12, draws 2, 4, 8, 16, 64, 420 fits): 12 fire, all on
a slope near separation (`bf(y ~ x + (x | g))`, bernoulli, 20 groups of 4),
and 96 are capped and warn, at those small counts. One message names a count
that settled (data seed 3, 8 draws: "settled at importance = 128", with a
worst group of 0.012, so its own number says to read the fit with care). The
rest list what they tried. I do not claim these are the only natural firings;
these are the ones I looked for.

**M2, the refusal's cost.** `imp_refuse_cells` is 2e6, `iterate()` is in
`tryCatch()` and an allocation failure becomes part of the message (test: a
mock that stops with `std::bad_alloc` at 64 draws). The comment reads 1.1 KB
per cell. `dev/impsamp-refuse.R`, `refuse.txt`, the wrapper displacing every
round, under load:

| rows, draws | the fit | the refusal | peak |
|---|---|---|---|
| 244, 512: the worst the bound allows, 1024 to 8192 | 8.5 s | 129.9 s, 15.3 times | 3.3 GB |
| 2400, 512: the review's case, no ladder now | 23.2 s | 39.4 s, 1.7 times | 3.3 GB (from the row above) |
| 2400, 64 | 6.5 s | 84.5 s, 13.0 times | |
| 20000, 500: no ladder | 155.5 s | 296.4 s, 1.9 times | 16.5 GB, the fit's own |

The review measured 796.7 s, 33.8 times, and 22.1 GB on the 2400 rows.

**M3, tests with power.** New tests, each run against the mutants below. The
objective off the anchor against an in-test integral at 4 displaced points;
the tape's `H[1,1]` on the poisson slope; the rounds' accuracy against
quadrature on a design where Laplace is far from the optimum (intercept -8.9
against -2.5, sd 4, 40 groups of 3, data seed 1003: the default at z 0.10, one
round at z 1.32, at 200 draws). `no_cut` is `abs()` in place of the cut of a
negative likelihood curvature at q = 1. The cut is reachable: Student rows
that disagree make the curvature negative at the mode (`scan-cut.txt`: at `nu
= 2`, 30 rows, 1 to 2 cut groups on 2 of 6 data seeds to 3 of 6). The new test
asserts a cut group's proposal is the prior's.

Mutants, `dev/impsamp-mutants.R`, `dev/impsamp-mutrun.sh`, each in its own
library outside the shared ones, with a stop if the mutant's library did not
load, on the file then in the worktree, baseline pass 314, fail 0:

| mutant | edit | result |
|---|---|---|
| `partial_rate` | the curvature's partial derivative | fail 2: curvature on a slope; rounds' accuracy |
| `abs_step` | the old step | fail 11, error 1: center or scale (8); dispersion's column (2); natural firing |
| `beta_only` | a column's scale for beta only | fail 3: dispersion's column |
| `half_jacobian` | half the log-determinant's move | off-anchor test fail 4 of 4 (run alone) |
| `no_jacobian` | a frozen log-determinant | off-anchor test fail 4 of 4 (run alone) |
| `one_round`, `grad_always` | stop after round 1 | fail 2, error 1: rounds' accuracy (z 1.32 over 0.5); natural firing |
| `no_rate` | curvature frozen | fail 3, error 2: gaussian exact optimum; natural firing |
| `no_cross` | mode ignores the likelihood's slope | fail 9, error 3 |
| `no_prior_mode` | the prior's move out of the mode | fail 8: 4 tests |
| `frozen` | the old estimator | fail 12, error 5 |
| `map_recycle` | the old `imp_par_list()` | fail 3: constant dpar |
| `no_refusal` | refusal removed | fail 5, error 5: 5 tests |
| `no_cut` | `abs()` for the cut | fail 1, error 1: its test |
| `capped_settles` | a capped run counts as settled | fail 3: its test |
| `no_prior` | prior out of the objective | fail 5: its test |

The whole-file runs of the two Jacobian mutants hit the 25 minute limit, since
every fit climbs the ladder. The named test is what I have for those two.

**M4, the MCSE.** Formula kept; the vignette and `?frm` now say it can read
low as well as high, with the Student numbers, and give something true to do:
refit with another seed or twice the draws and read how far the log-likelihood
moves. Could a cheap warning catch the Student case? I think not. The weights
there look ordinary: min ESS 0.81 to 0.87, in the range of designs where the
number reads high. The one statistic in the weights that separates even from
odd is the pair formula, and on the review's run it moves the Student number
by 1.08 (0.00711 against 0.00657, sd of the error 0.00929) and covers 0.88
against 0.82, so it does not reach the error. A warning that fires would need
a second fit at another seed, which doubles the cost and is not cheap. I did
not build one.

**m1, rounding.** `grad < grad_tol` and `moved < imp_move_tol` decide the round
count, and the textbook model at 200 draws reads a fresh gradient of 8.69e-4
against 1e-3. My claim that no branch is decided by rounding was wrong and is
gone. No test depends on the side: the textbook test takes 1 or 2 rounds, and no
other asserts a round count. OpenBLAS 0.3.26 and 0.3.32 are both clean on this
code (section 14).

**m2, the cost.** `dev/impsamp-perf.sh`, `perf.txt`, both builds interleaved,
2 reps. The machine is shared with other lanes, so this is not an unloaded
run: the control reads 0.67 to 9.22 s. Poisson, 64 draws:

| design | gradient before | gradient after | tape before | tape after |
|---|---|---|---|---|
| 5000 groups, q = 1 | 0.100, 0.114 | 0.116, 0.123 | 6.4, 6.7 | 5.9, 8.8 |
| 5000 groups, q = 2 | 0.126, 0.128 | 0.143, 0.158 | 7.0, 8.7 | 10.8, 12.4 |
| 5000 groups, q = 3 | 0.153, 0.156 | 0.217, 0.270 | 6.8, 8.2 | 11.8, 13.6 |
| 1000 groups, q = 2, 500 draws | 0.196, 0.222 | 0.262, 0.268 | 11.1, 12.0 | 13.7, 15.4 |
| 2000 groups, q = 2, p = 40 | 0.061, 0.062 | 0.075, 0.077 | 2.7, 3.5 | 3.6, 4.5 |

So a gradient costs 1.0 to 1.2 times at q = 1, 1.1 to 1.25 at q = 2 and 1.4 to
1.8 at q = 3; a tape up to 2.0 times. Whole fits where both arms take 2
rounds: q = 2, 92.6, 97.1 s against 91.7, 93.3 s; q = 1, 58.2, 64.9 against
59.4, 61.0; p = 40, 60.1, 67.1 against 52.2, 59.9 (up to 15 percent). So the
time NEWS gives for the grid (12 to 50 percent, by load: 4509 and 5372 s
before, 2680 and 4722 s after) is a fact about the grid, where the old fit
used 4 rounds.

The freeze at p = 40 reads 1.3 to 1.6 s against 0.28 to 0.36 (3.6 to 5.6
times; the review read 7). I took the cheap part: one-sided differences for a
fixed effect, and a map from `@x` built once instead of a triplet conversion
per Hessian. What is left is the p gradient and p Hessian sweeps themselves:
Rprof gives 56 percent to the gradient and 22 to the Hessian. A shortcut
through the predictor would need the family's second derivative, which no
structured family gives. Not done.

**m3, power restored.** The MCSE assertion: at 3 log sds off the probe's
anchor the worst group falls under the floor and the MCSE exceeds 5 times the
anchor's. The weakly identified test compares estimates with GLMMadaptive, as
a ratio to their SEs. The boundary test reads an in-test integral at each
fit's own estimates. The rise is read from a real fit against its tolerance.

**m4, absolute tolerances.** The three, and test-factorization's, are ratios
to an SE or to the fit's MCSE.

**m5, the stalled branch.** Deleted, with `imp_stalled()`, `imp_stall_tol` and
its text in `?frmtmb_control`. Nothing I ran reached it after B1, including
the 96 capped fits of `natural.txt`, whose moves are not equal.

**m6, like with like.** Both vignettes now say what the rerun at more draws
reads, instead of comparing a log sd's shift with an MCSE in nats: RL, seed
4201, -0.8703 under Laplace, then -0.7068, -0.7509, -0.7464 at 100, 400, 1600
draws (and seed 4205, -1.4841, then -1.2767, -1.3112, -1.3139); learning, seed
4302, 0.4170, then 0.4304, 0.4279, 0.4272 (`rl2.txt`, `learn-ref.txt`).

**m7, a prior.** Honored, since it was cheap: `neg_log_prior_fn()` already
returns a function of the outer parameters, so the objective adds it and the
pin expects it. A prior reads no group coefficient. Test against quadrature
with the same prior, class `b` and class `sd`; `no_prior` fails it. The row
stays works, and is now true.

**NITs.** The two "(quadrature = TRUE)" messages name `importance =` too.
`imp_verify()`'s per-group difference now tolerates the rounding of its sums:
the uncentered random slope no longer stops there
(`dev/review-impsamp-year3.R`, data seed 1; log `year3.txt`), and the pin test
still catches a per-group error. Stale words are gone (`imp_fit_proposal()`,
the header, `confint.R`). `test-importance.R`'s blank lines are back: 70
restored, the diff is 585 added, 106 removed before this round's additions. "A
proposal chosen before the draws are read" now says it holds at a given
anchor. "0.02 of 500" is now a range over ten data sets, 0.02 to 1
(`rl2.txt`).

Filed, found this round:

- `dev/impsamp-scales.R`: Student at `y x 1e6` stops the pin on both builds, and
  `x + 1e5` poisson gives NaN SEs with a false convergence warning (section 8).
- My first `integrate()` in a test read 0.04 nats wrong, because its default
  absolute tolerance accepts anything on an integrand of size e^-40.
  `imp_log_integral()` centers on the peak. The grid's reference is
  Gauss-Hermite and was not affected.

## 16. R CMD check on the final build

`--as-cran`, `dev/impsamp-rcheck.sh`: frmtmb Status 2 NOTEs (examples timing
under load, `frm_hazard_reads` 5.7 s and `dharma_residuals` 8.9 s, not touched
here; V8 for math rendering), tests FAIL 0, WARN 0, SKIP 234, PASS 20287.
frmtmb.learn Status OK, tests FAIL 0, WARN 0, SKIP 15, PASS 421.
frmtmb.spatial's only change is a comment in a test that names a renamed
function; see the Status below.

Bump: minor for frmtmb (BREAKING NEWS bullets, a removed warning, a new
refusal). The frmtmb.learn floor is the frmtmb that carries this fix, because
test-factorization.R asserts the settled fit. frmtmb.spatial needs no floor
change.

frmtmb.spatial Status 1 NOTE (V8 for math rendering), tests FAIL 0, WARN 0,
SKIP 7, PASS 410.

## 17. Fix round 2

The coordinator's fix round 2, from "Short check of fix round 1" in
`dev/reviews/2026-11-02-impsamp.md` (logs in `dev/review-impsamp-out/fr1/`).
Numbers are from the final build unless a line says otherwise.

**F1, a step that read a parameter's own size: fixed, for every parameter.**
The review's `fix_tau` handled `tau_raw`. I took the rule out instead, since
the rule is the defect: `imp_outer_steps()` no longer reads any parameter's
size. A parameter with no column steps by `imp_cross_step`. The review's
`fix_abs` is this change, and it noted that it moves the natural-firing test
off its seed; see below. `dev/impsamp-ordfix.R` (the review's script),
`ordfix.txt`, `y ~ w + (1 | g)`, `w = 5 z + m`, 300 draws, slope SE as a ratio
to the centered column:

| m, data seed | fix round 1 | rellib-r16 | now |
|---|---|---|---|
| 1e4, 5 | 4.56, no warning | 0.996 | 1.0001 |
| 1e4, 3 | 1.053, no warning | 1.0007 | 0.9996 |
| 1e4, 6 | 1.16, false convergence | 0.998 | 1.0001 |
| 3e4, 5 | NaN, capped | 0.996 | 1.0002 |
| 1e5, 5 | NaN, capped | | 1.0008 |

At 3e4 on seed 6, and at 1e5 on seeds 3 and 6, the optimizer reports false
convergence and the SE is NaN, with the warning; log sd and logLik equal the
centered fit's. That is the class of `x + 1e5` poisson, filed below. At 1e5 on
seed 8 there is no warning, and the fit reads log sd -0.772 against -0.763,
slope SE 0.9931, logLik 0.0088 low: small and silent, at a center 2e4 times the
column's scale. I file it with that class.

Each class of parameter with no column, what it held and what it now reads
(`dev/impsamp-classes.R`, `classes.txt`):

| class | largest value seen | can its size reflect a column? | reads now |
|---|---|---|---|
| `tau_raw`, cumulative | 6141 at m = 3e4 | yes, the first threshold takes the center | SE 1.000, dll 0 against centered |
| `tau_raw`, sratio, cratio, acat | 5581, every threshold | yes, all of them | SE 0.9998 to 1.000 |
| `tau_raw`, `thres(gr = h)` | 2057, one per level | yes | SE 1.000 |
| `theta`, a log sd | 12.3 at `z x 1e-6`, -15.4 at `z x 1e6` | it moves the proposal in log units, not by the column | SE 0.9987 to 1.0003 |
| `theta`, a correlation's chart | 119 at the limit | no | beta and log sd SE 0.997 to 1.000 against 2000 draws |
| `zeta`, a `mo()` simplex's chart | -0.02 | no | SE 0.999 to 1.002 against 2000 draws |
| a family's `extra_pars` | not measured | possible: a custom family may add one to its location | the rule is the same, so size cannot enter |

The last row is by construction. One thing I did not settle: a family's own
parameter whose natural scale is far below 1e-4 takes too large a step. That
is true of the old rule as well, and it costs efficiency and not validity.

**The natural firing is gone, and I do not call it impossible.**
`dev/impsamp-natural.R`, the same 420 fits: 0 fire and 106 are capped, against
12 and 96 on fix round 1. The 12 were all on the slope near separation, and
that design is now 57 of 60 capped, each with the capped warning. The only
change that can have moved them is the step: taking size out for `theta` too.
So only the test's wrapper reached the refusal at this reading. Corrected by the
review's short check: a wider search finds 8 natural firings in 1260 fits, 5 at
the default rounds (`dev/review-impsamp-fr2-natural.R`), and they are chaotic.
That leaves its condition untested on a real fit, which is a weakness of the
tests and I report it as one. The test that pinned the firing now asserts what
must hold whichever way the design lands: the fit stops with the refusal, or
it is capped and warns.

**F2, a precision that rounded: fixed.** `imp_prior_scales()` searches, off
the tape and once per anchor, for a `t` per coefficient at which the block's
density at `t e_j` lies about one nat below its value at zero, and
`imp_block_prior()` reads there. Per coefficient rather than per block, so a
slope and an intercept on different scales both read true. Test: the reading
at sd 1e-3, 1 and 1e8 within 1e-8; a fit on gaussian and on Student at `y x
1e6` and `1e7`. `classes.txt`, data seeds 1:4, against the fit in units with
the Jacobian, worst of the four:

| family | 1e3 | 1e6 | 1e7 |
|---|---|---|---|
| gaussian, abs dll, SE | 7e-10, 1.0000 | 6e-10, 1.0000 | 5e-9, 1.0000 |
| Student `nu = 4`, abs dll, SE | 1e-8, 1.0000 | 2e-9, 1.0000 | 1e-8, 1.0000 |

It is off the filed list.

**F3, the rounds' start: fixed.** `importance_fit()` calls `sd_up_escape()` on
the warm start, through a small fit-shaped list. Same 40 seeds and the
review's script (`dev/impsamp-escape.R`, `escape-after.txt`, against the
review's `escape.txt`):

| design, 40 data seeds, fits below Laplace by more than noise | before | after |
|---|---|---|
| ordinal, `w = 5 z + 2000` | 15 | **0** |
| ordinal, standard column | 0 | 0 |
| bern_small | 0 | 0 |
| pois_small | 0 | 0 |
| gauss_small | 0 | 0 |
| unbal_bern | 0 | 0 |
| wtbin | 0 | 0 |
| t_outlier | 1 | 1, the same fit |
| smallpois | 26 | 26 |

The 6 standard-column designs the review counts (240 fits) read 0 before and 0
after. The other two rows are not collapses: Laplace itself reads high there,
and the log sds agree (t_outlier's seed 4: -0.315 and -0.302). On
`imp_ord_data(8, 2000)` the fit now reads log sd -0.763, the centered fit's,
where it read -13.49. I took the escape alone. The Laplace route's other steps
(`escape_stationary()`, `mo_search()`, `newton_polish()`) answer other stops,
and nothing here measured a need for them.

**F4, the named count: fixed.** A count is named only if that fit settled and
its worst group is at or above `importance_ess`. "or more" is gone. The test
uses the wrapper at a floor no group clears; mutant `low_ess_named` fails it.

**F5, NEWS: done.** About 40 percent less where the old fit took 3 to 5
rounds, beside median rounds 4 to 1, as the review reads it. The
per-evaluation costs stay.

**F6: `u` stays, with a test; the wording is changed.** I removed `u` as asked
and measured it. On large units the review is right: no effect. On small units
it has one, which the review's `y x 1e2` to `1e5` did not read: Student at `y
x 1e-6` without `u` stops on 2 of 4 data seeds (seeds 3 and 4) and reads SE
1.0073 and 1.0015 on the other two, since a step of 1e-4 is 100 of that
response's scale. With `u` it reads 1.0000. A new test holds it, and
`no_units` fails it. "this check stops at 100 draws" to a user at 500 now
reads "No count above 500 was tried: on 20000 rows it would tape more than
2e+06 rows times draws".

Filed, not fixed:

- **MAJOR, for core, not this lane.** `y ~ (x + 1e6) + (1 | g)`, poisson: the
  Laplace fit's slope is off by 5.4, 8.9 and 9.3 SE on data seeds 1 to 3, with
  no warning on 2 of 3, and importance reads the same. The base build's
  importance fit reads -5.3 on seed 1. At `x + 1e7` the column drops as rank
  deficient (the review's `own2-after.txt`).
- `x + 1e5` poisson stays filed (section 8). `logLik` holding the prior on all
  three routes is consistent, and the coordinator raises it with the user.

Notes on the tables above, from the final logs. The `theta` log sd row gives the
SE of that log sd. On the same fits the slope reads 1.001 and 0.9986, the other
log sd and the correlation's chart 0.9941 to 1.0012, and the `z x 1e6` fit warns
once. The correlation row's chart itself has no SE at the limit, on either side.

The time in F5 is the coordinator's wording. This grid does not reproduce it:
over the 341 paired fits at 3 to 5 old rounds the total is 33 percent less, 16,
31 and 41 percent at 3, 4 and 5 rounds, and the median of the per-fit ratios is
23 percent less. The whole grid is 26 percent less (5372 s to 3978 s). This run
shared the machine with other lanes and the base arm did not run beside it, so
seconds read load. The coordinator decides whether NEWS keeps the sentence.

Tests against the base (`dev/impsamp-r16check.R`, `r16check.txt`). "the rounds
start where the Laplace route ends" fails on rellib-r16, fail 2, error 0: F3
is the base's defect. The base stops at the other new tests too (threshold
fail 2, error 1; sd far from one error 1; small units error 1; named count
fail 1). Those read this round's internals or wording, or they pin defects fix
round 1 brought in, so they are not evidence about the base.

New mutants on the final code (`dev/impsamp-mutrun.sh`), each fails only its
test:

| mutant | what it does | result (baseline 343) |
|---|---|---|
| `own_size` | step by the parameter's own size | fail 7, threshold test |
| `t_one` | read the prior at `t = 1` | error 1, sd far from one |
| `no_escape` | no `sd_up_escape()` | fail 2, rounds' start |
| `low_ess_named` | name a count below the floor | fail 3, named count |
| `no_units` | no response unit in the step | error 1, small units |

The mutants of fix round 1 did not rerun. Their count of the natural-firing test
is out of date, and the `no_refusal` site in `dev/impsamp-mutants.R` no longer
matches. No real fit now reaches the refusal, so `no_refusal` would read only
through the wrapper tests.

## 18. What ran on the final build of fix round 2

The grid (`grid-summary-final.txt`): 648 of 648 fit, 0 stops, 0 capped, against
162 and 204 on the base; median rounds 4 to 1, largest 4; median error 0.0041 to
0.00021 in beta and 0.017 to 0.0010 in log sd. Warnings on 22 fits, as listed
there.

Year table (`year-after.log`), data seeds 1:6, poisson and bernoulli, against
the centered column: largest slope z 0.0049, SE ratio 1 to the digits shown,
largest abs dll 1.2e-5, no warnings. `scale-se.txt`: slope 0.99999, theta 1.

Off-grid (`offgrid-summary.txt`): 116 fits of 28 designs, largest max z 0.148,
Newton decrement at most 0.153, Laplace's max z 3.63, 6 warned. The SE ratios
are 0.90 to 1.08 but for two fits, both as in the review's own logs: `null2`
seed 3, 7246 (review 7290), the correlation's chart on a flat likelihood, with
the false convergence warning, the other four at 0.98 to 1.03; and `year` seed
2, 0.867 to 18.27 (review 0.867 to 18.63), with no warning. That reference was
the reviewer's own wrong difference, 1e-3 on a raw column near 2000: against the
exact Hessian on the centered column, by the delta rule, seed 2 reads SE 1.0011
(dev/review-impsamp-fr2-year.R), and no warning is due. Section 14's "0.90 to 1.08" left these two out and
did not say so.

Core, one file per process: 286 of 286 files, pass 20563, fail 0, error 0,
skip 206, warn 0. learn 15 of 15, pass 429, skip 13. spatial 8 of 8, pass 478,
skip 6. sample 78 of 78, pass 3767, skip 15. latent 12 of 12, pass 681, skip
2. eam 32 of 32, pass 1838, skip 5. All fail 0, error 0, warn 0.

Gated, `FRMTMB_BRMS_FIT_TESTS=true`: test-me 71, test-rl-example 102,
test-importance 343, test-bcm-simple 13, frmtmb.spatial test-brms-likelihood-car
67, frmtmb.sample test-loo 146. All fail 0, error 0, skip 0, warn 0.

One session (`sessions-impsamp.log`), forward, reverse, preload: core pass
20563, 20713, 20806; learn 429; sample 3767; spatial 478; latent 681. All fail
0, err 0, warn 0, 15 of 15 sessions, 37 LEAK lines, all 37 tagged "[with load
of ...]". Core 276 of 286 and sample 73 of 78 files, as before.

OpenBLAS with LAPACK, 4 threads, list `ob-list.txt`: 0.3.26 and 0.3.32 each 10
of 10, with identical counts: covstruct-seam 131, grad-verdict 74, importance
343, rl-example 102, se-exact-hessian 35, factorization 34, rlddm-ndt 53,
session 22, surface 146, car 132. All fail 0, err 0, skip 0, warn 0.

`R CMD check --as-cran`: frmtmb Status 2 NOTEs (`frm_hazard_reads` 5.59 s under
load, not touched here; V8 for math rendering), tests FAIL 0, WARN 0, SKIP 234,
PASS 20316. frmtmb.learn Status OK, tests FAIL 0, WARN 0, SKIP 15, PASS 421.
frmtmb.spatial Status 1 NOTE (V8), tests FAIL 0, WARN 0, SKIP 7, PASS 410.

p124 and p069 (`dev/impsamp-verdicts.R`): both fit. p124 reads `equivalent`, L1
const 8e-13, L4 z 0.0487 and 0.13 against control 0.0488 and 0.13. p069 reads
`defect` on L4, z 0.0872 and 0.627 against control 0.0885 and 0.627, so the fit
sits at its control's reading, and L1 const is 2e-10. The reviewer's ruling
stands: agree with weak evidence, the compat row stays conditional, and
consolidation drops both adjudication rows.
