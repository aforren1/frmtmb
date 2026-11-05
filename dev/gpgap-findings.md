# Lane gpgap: the exact gp() L1 gap against brms 2.23.0

Round opened 2026-10-29, base d5d0056d (frmtmb 0.78.0), with fix round 1 of the
review `dev/reviews/2026-11-02-gpgap.md`. The numbers are in section 9, generated
from the logs by `dev/gpgap-summary.R`, and the last section before them, 10, is what the
installed build reads. Logs: `dev/gpgap-out/`,
`dev/gpgap-equiv-out/` (local, gitignored).

## 1. Result

The 94 to 1127 nats were the INSTRUMENT. frmtmb sorts a term's distinct
positions (`posdf[do.call(order, posdf), ]`, `R/frame.R`). brms keeps them
in order of first appearance (`Xgp[!duplicated(Jgp), ]`,
`brms:::.data_gp`). The translator `stan_pars_from_fit()` read
`z = L^-1 f` in frmtmb's order and gave it to brms, which multiplies by
its factor in brms's order, so each latent value landed at another
position. The suite never saw this: the exact row of
`test-brms-likelihood.R` has `x = sort(runif(n, 0, 10))`, and
`test-brms-agreement.R` compares only the partitions of the two indices
("Level numbering is an internal detail"). Fixed in the translator.

Reordered, on brms's own program at frmtmb's estimate, the gaps are:

| pair | before | brms's order |
|---|---|---|
| p022 | -1126.63 | -5.117 |
| p023 | -527.37 | -2.215, net of the free prior constant 0.671 |
| p094 | -94.43 | -0.0282 |
| p109 | -180.04 | -0.485 |
| p136 | -151.65 | -0.531 |

The rest is the one line where the models differ, the nugget: with brms's
line `cov[n, n] += 1e-12` rewritten to frmtmb's
`cov[n, n] += 1e-6 * square(sdgp)`, the constant is at most 3.8e-9 in
absolute value at the estimate and 6.9e-9 at the two stepped points on all
five pairs, and the gradient in `zgp` is at most 1.9e-9. So scaling, kernel,
Cholesky map, prior of z and likelihood are one model in the two packages,
apart from that line.

The review's two "unexplained" numbers close:
- -83.81 at 1e-9 and -83.74 at 1e-12 reproduce to the digit through the
  old order, and read -0.0203 and -0.0197 in brms's order
  (`dev/gpgap-review83.R`).
- `dev/equiv-gp-gap.R`'s B 911.2 against C 1013.5 reads B 1013.508674
  against C 1013.508670, B - C = 4.1e-6, once the script reads brms's
  order. Double precision was not the obstacle; the order was. The
  extended-precision run says the same (section 3).

## 2. Verdicts (M1; this replaces my first round's reading)

An edited program that matches shows that frmtmb equals the edit. That
attributes the gap. It is not agreement with brms's program, which is what L1
is defined on. So:

- p022, p023, p094 and p109 are **divergence by decision**, `agree = FALSE`.
  The reason is in NEWS. The instrument reads L1 on brms's own line, which fails
  the tolerance, and reports the rewrite as the attribution.
- p136 stays **equivalent at L4 only, weak**, as at the baseline.
- `results.tsv` gives these rows no "every level agrees". The reason and
  `summary.md` state the rewrite, brms's own-line constants and the coincidence of
  the maxima at ordinary noise, each measured by the row
  (`eq_gp_brms_max_gap()`).
- The `spec$l1_decision` branch is gone; nothing used it.

My first round's headline of 84.9 was wrong, and it goes. The review's table
(`dev/review-gpgap-rates.R`, survey summary run unchanged, only the adjudication
differs):

| reading | agree | no defect | weak evidence out | measured only |
|---|---|---|---|---|
| baseline (0.78.0, cannot measure) | 83.6 | 84.6 | 80.5 | 92.7 |
| my first round: 4 rows dropped, edited L1 (withdrawn) | 84.9 | 85.9 | 82.9 | 92.8 |
| **divergence on 4, p136 L4: this lane's result** | **83.6** | **85.9** | **80.5** | 91.3 |
| divergence on all 5 | 82.5 | 85.9 | 80.5 | 90.2 |
| closed form as L4 on p022, p023, p094; p109 divergence; p136 L4 | 84.6 | 85.9 | 81.5 | 92.5 |

I reproduced the bold row from the flow's own verdicts: 83.6, 85.9, 80.5, 91.3
(`dev/gpgap-out/summ/`, a scratch copy without the 4 adjudication rows).

`dev/equiv-data/adjudication.tsv`, not mine to edit, still forces rows p022, p023,
p094 and p109 to "cannot measure" and quotes the withdrawn gap, so
`dev/gpgap-equiv-out/summary.md` reads 83.6 and 84.6 until consolidation drops or
replaces them.

## 3. Round-off, separated (`dev/gpgap-mpfr.R`, 256 bits, Rmpfr)

Rmpfr is in the user library; nothing was installed. On p094 and on p023
without its prior, in brms's order:

- The instrument's reading on brms's own program is correct in double:
  p094 -0.02819900 exact, -0.02819902 through Stan, difference 1.8e-8. The
  gap comes from the factors: the same z through the two factors moves
  the field by up to 1.1e-2 (p094) and 2.2e-1 (p023).
- The gap at the same FIELD is a different and larger number, and it is
  the nugget too: C = 1013.5066 exact, 1013.5087 in R's double (p094),
  520.43 and 520.60 (p023). It lives in the kernel's near-null directions:
  sum z_b^2 is 281.8 and sum z_f^2 is 6.4. Stan rebuilds f from the exact
  z_b: its likelihood there equals ll(f) to the 6 decimals printed on p094,
  and to 1.6e-4 on p023.
- The last row of `dev/gpgap-review83.R` puts brms's exact nugget in
  frmtmb's kernel: -4.5e-3 at a false-convergence fit, in double, on a
  kernel conditioned about 1e13. I did not attribute this residual
  further: the rewritten-line identity settles the model without it.

So the gap is not round-off on both sides. The answer to brief item 4 is
no, and the equiv specs need no better-conditioned data.

## 4. The warning (B1, M2, M4, M5, m1)

The first form compared sigma alone and threw on a `sigma` with no coefficient.
It is replaced. For each exact block, `gp_position_info()` reads the likelihood's
curvature in each latent value off the tape (the joint Hessian's diagonal less the
block's prior), the inverse of the working noise at that position. Then `r` is the
nugget's variance over that noise, and the rule is

    90th percentile of r >= sqrt(2 / m), m the block's positions

The arithmetic: with one row per position, `sigma_frm^2 = sigma_brms^2 - nu sdgp^2`,
so `log sigma` moves by `r / 2` against a standard error of about `1 / sqrt(2 m)`;
that is the review's `1e-3 (n / 2)^(1/4)` in sigma over sdgp. Why the 90th percentile
and not the median: a quarter of the positions at low noise moves sdgp from 2.36 to
1.44 (`sigma low on share 0.250`, gap 2.40 nats) with a median r of 2e-5, and one
low row in 60 moves nothing.

It cannot abort a fit (`tryCatch` around the reading; a test mocks the failure),
it reads an offset under sigma, and it has one rule for `se()`, `sigma ~ (1 | g)`,
`sigma ~ s(z)`, mixtures, nl bodies, a Poisson at a large mean, a binomial at many
trials, and replicates, because all of them are the curvature. It costs 0.65 percent
of a 200-position fit (`dev/gpgap-msgcost.R`, control 1.09).

Calibration: `dev/gpgap-calib.R`, 281 designs (the review's floor, rep, silent and
trend designs on its seeds, and mine), each fitted and read against the same tape
non-centered at brms's line, from 3 starts, and at a relative 1e-12 where y's units
are not near 1 or the field is steep (M7). "Moved" is any outer parameter at 1
standard error of the reference, or 0.5 nats either way. Generated in section 9:

- **150 of 152 moved rows fire.** The two misses: Poisson log mean 10 (0.52 SE, 0.88
  nats), and `sigma ~ (1 | g)` at noise 0.1, where frm() stopped with false
  convergence for a reason of its own: r there is 3e-4.
- **5 of 128 unmoved rows fire.** Four sit at 0.89 to 0.99 SE. The fifth, `rep 30x10
  r0.01`, reads 0.26 SE, with frm() 0.16 nats above the reference.
- The first form had 42 of 101 warned rows with a bias below 1 SE (review, m1). On
  one row per position here: 3 of 97. One design is unread (gap below -2 nats, a
  reference that stopped short).

The row the review could not read: binomial at 1e5 trials on its seed 5601, frm()
0.85 nats above its reference. The reference did not stop short: from three starts at
absolute and relative 1e-12 it reaches the same -444.8539 with sdgp 2.00; at 1e-6 it
reads -444.0001 with sdgp 2.55 (`dev/gpgap-binom1e5.R`). My first explanation, that the
nugget's white noise pays for itself there, was one-sided, and I withdraw it. On the
review's seeds 6001 to 6012 (`dev/review-gpgap-binom.R`, importance-sampled at each
model's own maximum) the gap of 1e-6 less 1e-12 runs from -0.77 to +0.54 with mean
+0.04, and sdgp moves both ways. It is a seed-sized fluctuation at r near 0.14, against
a rule of 0.18 at that m, so the warning's silence is right. frm() above the reference
is not by itself a failed reference.

The message says what moves: an estimated sigma on unreplicated data; sdgp and the
length scale where positions repeat, sigma is known, or the family has none. It does
not promise `k =`.

Not covered: `importance =` and `integrate =`, which keep their own objective. NEWS
and `?"frmtmb-gp"` say so.

## 5. The sampling route (M3)

`frm_sample()` was non-centered and still read 1e-6, whose stated reason does not
apply there. It now reads a relative 1e-12 on its non-centered exact blocks:

- `gp_block_nugget(blk)` in core reads `blk$gp_nugget`, falling back on `gp_nugget`;
  `gp_corr()`, the cross-covariance and all four kriging readers use it. The sampler
  sets the field on its copy of the fit (`gp_sample_nugget()`), so tape, draws and
  kriging read one model. The centered route and `laplace = TRUE` keep 1e-6.
- The chain's start comes from the fit's own factor. Through the smaller nugget
  the mode's field has a z the size of the near-null directions (sum z^2 of 1464
  against 7 on p023's data).
- `gp_chol_num()`: a numeric factor that takes the smallest jitter that works, at
  most 1e-6, so a draw the tape accepted cannot lose a finished chain in LAPACK, and
  a prediction cannot stop.
- The formula route holds: `frm_sample(bf, data)` is what the runs below use.
  `test-reparam.R`'s identity (centered objective = non-centered + log |L|) holds at
  each block's own nugget, unchanged, and its new rows read the run's tape against
  the closed form: on the base build they fail (4 failures, 2 errors).

L1 against brms's own line at ordinary noise, the sampler's tape at a typical z, on
the five pairs (`dev/gpgap-sampler-l1.R`), now in the factor's order (section 10):
2.2e-3, 9.4e-2, -2.5e-3, 2.7e-2, 3.5e-2 nats, where `frm()`'s line reads -0.028 to -5.1. It is not an identity: relative
against absolute. The same tape at 1e-6 against brms rewritten reads at most 1.1e-13,
so nothing else in the tape moved.

`laplace = TRUE` also holds an exact block by its factor, for the nugget. It stayed
centered at 1e-6 before, where it reads the same floor. It is an addition to the
coordinator's list, and I made it for a reason I measured: **it is the route that
mixes at true sigma 1e-3**, and the integral is exact for a gaussian response.
`check_laplace()` compares against frm()'s fit and keeps 1e-6.

Converged runs against brms (`dev/gpgap-sampler.R`, seeds as named, 4 chains;
brms's default priors on both sides; q05, q50, q95 of sigma):

| data | route | sigma | sigma Rhat | sdgp q50 | wall |
|---|---|---|---|---|---|
| true 1e-3, n 15, seed 6601 | `frm_sample(laplace = TRUE)`, base (1e-6) | 7.0e-5, 7.4e-4, 2.8e-3 | 1.000 | 1.49 | 129 s |
| | **`frm_sample(laplace = TRUE)`, now (1e-12)** | **5.3e-4, 1.03e-3, 2.3e-3** | **1.000** | 2.17 | 168 s |
| | brms drawn, 20000 iterations, tree depth 13 | 5.2e-4, 1.03e-3, 2.4e-3 | 1.031 | 2.17 | 3189 s |
| | field integrated at brms's line (`dev/gpgap-marg15.R`) | 5.2e-4, 1.04e-3, 2.4e-3 | 1.000 | 2.12 | |
| true 1e-3, n 40, seed 6601 | `laplace = TRUE`, base | 2.3e-5, 3.0e-4, 8.7e-4 | 1.001 | 1.19 | 195 s |
| | **`laplace = TRUE`, now** | **8.7e-4, 1.08e-3, 1.37e-3** | **1.001** | 2.62 | 386 s |
| | field integrated at brms's line (review, `marg.txt`) | 8.7e-4, 1.07e-3, 1.36e-3 | 1.00 | 2.74 | |
| y x 1000, noise 100, n 60, seed 7701 | **`frm_sample()`, drawn, now** | **90.6, 106.5, 128** | **1.001** | 1192 | 708 s |
| | `laplace = TRUE`, now | 90.4, 106.9, 127.6 | 1.004 | 1175 | 645 s |
| | brms, drawn, its own line | 96, 260, 642 | Inf (sdgp 10.0, 186 divergent) | 490 | 2021 s |

**Superseded by section 10:** the factor's order was the cause, and the installed build
reads differently. These are the runs over SORTED positions, kept as the record.
At n = 15, seed 6601, 4 chains of 20000:

| drawn run | sigma q05, q50, q95 | sigma Rhat | sdgp Rhat | at max tree depth | wall |
|---|---|---|---|---|---|
| brms, diag, delta 0.99, depth 13 | 5.2e-4, 1.03e-3, 2.4e-3 | 1.031 | 1.154 | 39867 of 40000 | 3189 s |
| `frm_sample()`, diag, delta 0.99, depth 13 | 5.0e-4, 6.5e-4, 1.8e-3 | 2.461 | 2.978 | 39996 of 40000 | 14347 s |
| `frm_sample()`, dense, delta 0.95, depth 12 | 3.9e-4, 4.3e-4, 1.5e-3 | 2.834 | 3.088 | 39999 of 40000 | 7769 s |

The target is the row of section 5's table: q50 1.03e-3. So the coordinator's standard is
met at true sigma 1e-3 by `laplace = TRUE` (Rhat 1.000, 168 s), and on y x 1000 by the
drawn route (Rhat 1.001), and NOT by the drawn route at true sigma 1e-3. Its median
sigma there is no estimate. Of the two candidates I filed, the review measured the
second: the sort. The first was wrong for these runs, which started at random.

## 6. For the user

**Non-centering `frm()` is a separate round, and yours to decide.** Until then
`frm()` keeps its centered 1e-6 with the repaired warning. The facts for it:

1. **It is not only noise-free data** (M4, the review's table, closed form at each
   nugget, 100 rows, x on 0 to 10, true sigma 0.3; `dev/review-gpgap-silent.R`, seeds
   5701 and 5702):

| y | frm() sigma, each seed | brms's sigma, each seed | sdgp frm, brms, seed 5701 | ll gap, each seed |
|---|---|---|---|---|
| `0.2 x` | 0.340, 0.283 | 0.340, 0.283 | 0.74, 0.74 | 0, 0 |
| `x` | 0.340, 0.281 | 0.340, 0.283 | 7.5, 7.5 | 0, 0 |
| `5 x` | 0.330, 0.0001 | 0.340, 0.283 | 82, 82 | 0, 0.07 |
| `x^2` | 0.0001, 0.0001 | 0.344, 0.283 | 347, 813 | -0.65, -0.64 |
| `exp(0.5 x)` | 0.0001, 0.0001 | 0.352, 0.286 | 356, 664 | 0.40, 0.25 |
| `exp(0.8 x)` | 0.0000, 0.0001 | 0.352, 0.272 | 657, 19730 | 67.7, 75.7 |
| `100 plogis(x - 5)` | 0.353, 0.286 | 0.355, 0.288 | 33, 33 | 0, 0 |

   The decision of 2026-10-06 rested on sigma 0.3 to 0.4 on signals of range about
   2. `y ~ gp(x)` on a curve of range 100 with sigma 0.3 reads sigma 3 orders low,
   and `exp(0.8 x)` is 68 to 76 nats from brms's maximum. The warning fires on these
   rows, and stays quiet on `x` and `plogis`, where the ll gap is 0.

2. **brms's absolute line is not a target** (M7). It has no scale invariance. On y
   times 1000 (`dev/review-gpgap-ncp.R`, `ncp2.R`, seeds 911 to 913) a non-centered
   tape at absolute 1e-12 reaches relative convergence on 0 of 18, against 9 of 12 at
   relative 1e-12. brms itself: its kernel does not factor at frm()'s estimates
   (kappa 2.3e18), and its run reads Rhat 3.77, 267 divergent, sigma 258 against a true
   100 (`dev/review-gpgap-brmsunits.R`, seed 7701), and I read the same again in
   section 5. In small units it floors: y times 1e-3, true sigma 3e-7, brms's closed
   form reads 5.68e-11 (seed 122). So my first round's "brms has no such floor" was
   false, and it is gone. **Any later `frm()` route must be relative.** It does not
   make L1 a plain identity either: relative is not brms's line.

3. What the separate round would cost, by count, not built: 25 reads of
   `estimates$b` in core, 27 of `jointPrecision`, `par.random` or `diag.cov.random`,
   15 of a block's `vcov()`, every start built from a fit, and kriging that reads z.
   A fit took 4 to 9 times the wall clock (14 s against 3 s at 200 rows here, 10.8 s
   against 1.2 s in the review). On the low-noise design the non-centered tape at
   brms's line reaches the closed form in 5 of 6 rows.

**`gp(x, k = )` is not the remedy by itself** (M6). It has no nugget and reaches
brms's numbers when k is large enough for the field. The review's case,
`dev/review-gpgap-approx.R`, seed 8102, 200 rows, `sin(8 x)` at noise 1e-3: k = 30
reads sigma 0.70 against 9.7e-4, 720 times high, with `conv 0` and no warning. I read
why (`dev/gpgap-kstart.R`): the k = 30 basis stops at frequency 6.51 and the signal
is at 8, and from k = 60's scales the k = 30 tape reaches
-214.2, no higher than the -213.6 it found. So it is the approximation, not its convergence, and no
convergence defect is filed. The absence of a check is filed in `dev/test-backlog.md`.
The message says k must be large enough, and to compare several.

## 7. Scope and bump (m8)

Minor, for the new condition class and the sampler's change. frmtmb.sample's floor
rises: it sets a block field that an older core ignores, which would sample the old
nugget in silence. No other extension reads a changed symbol.

## 8. What changed

First round: `tests/testthat/helper-brms.R` and the spatial copy
(`brms_gp_order()`, `lscale_` under `scale = FALSE`), `dev/equiv-gp-gap.R`.

This round:

- `R/covstruct.R`: `gp_block_nugget()`, `gp_position_info()`, `gp_nugget_message()`
  replaced, `gp_chol_num()`, `chol_L`'s numeric route. `R/predict.R`: the kriging
  readers. `R/fit.R`: the warning. `R/gp-by.R`: help topic `frmtmb-gp`
  (alias `frmtmb_gp_nugget`), since gp() has no topic of its own in frmtmb.
  `R/conditions.R`, `R/sampling-api.R`, `_pkgdown.yml`, `man/`.
- `extensions/frmtmb.sample/R/sample.R`: `gp_sampler_nugget`, `gp_sample_nugget()`,
  `ncp_start_pars()`; its NEWS. `test-reparam.R`: three rows.
- `tests/testthat/test-gp-nugget-floor.R`: rewritten: class by name, the boundary
  each side, every route of M2, B1's rows, a failure to read the tape, share and
  one-row low noise, replicates, the block's nugget, the numeric factor.
- `tests/testthat/test-brms-lik-gp-exact.R`: the unsorted approximate row (m4), the
  sampler's factor against brms's own line on a well-conditioned kernel.
- `.github/workflows/brms-likelihood.yaml`, `ubuntu-next.yaml`: a step for the exact
  file and the cache key's md5 list (m3).
- `dev/equiv-flow.R`, `equiv-lib.R`, `equiv-specs.R`, `equiv-summary.R`: section 2.
- `vignettes/brms-migration.Rmd`, `frmtmb.Rmd`, `compatibility.Rmd` (m7);
  `dev/brms-likelihood-tests.md`, `dev/equiv-findings.md` (m5: lines 277 and 545, 2.5,
  2.6, section 4), `dev/test-backlog.md`. m6: NEWS says "usually".
- Nits: the dead `se()` branch and the `min()` went with the old condition; the test's
  absolute tolerances are ratios; logs name their R copy (extSoftVersion is empty on
  Windows); 0.3.32 ran; check Status lines come from the logs.
- `b[i]`'s order: documented in `?"frmtmb-gp"`, filed in the backlog. The fix touches
  the block's `cnms`, which every label reader shares, so it is not cheap.
- m2: mutants of the rule in section 9, all killed but the control.

## 10. Fix round 2 (review, "Short check of fix round 1")

Read section 5's table as the record of the factor over sorted positions. This section
is what the installed build reads. Runs: `dev/gpgap-sampler.R`, seed 6601 or 7701,
4 chains, one process per run; logs `dev/gpgap-out/samp/r2-*.txt`, each generated into
section 9.

### F1: the factor's order

`gp_sample_nugget()` gives each exact block it moves a `gp_order`, from
`gp_farthest_first()` in frmtmb.sample; `chol_L` factors `K[p, p]` and returns
`Lp[order(p), ]`, and the kriging uses the same order. The field stays sorted, so draw
names and the field's positions do not move. frm()'s own blocks carry no order.

**The drawn route, on this build, as measured:**

| data | route | sigma q05, q50, q95 | Rhat sigma, sdgp, lscale | bulk ESS sigma | divergent | at max depth | wall |
|---|---|---|---|---|---|---|---|
| true 1e-3, n 15, 12000, delta 0.95, depth 12 | `frm_sample()` | 5.9e-4, 1.05e-3, 2.32e-3 | **1.025**, 1.069, 1.061 | 191 | 508 | 23420 of 24000 | 3474 s |
| the same data, fix round 1, 20000, delta 0.99, depth 13, sorted | `frm_sample()` | 5.0e-4, 6.5e-4, 1.8e-3 | 2.461, 2.978, 2.230 | 5 | 0 | 39996 of 40000 | 14347 s |
| the same data, 20000, delta 0.99, depth 13 | brms | 5.2e-4, 1.03e-3, 2.38e-3 | 1.031, 1.154, 1.112 | 152 | 0 | 39867 of 40000 | 3189 s |
| field integrated at brms's line | Stan, 4 chains of 4000 | 5.2e-4, 1.04e-3, 2.41e-3 | 1.000, 1.001 | | 0 | | |
| true 1e-3, n 40, 3000, delta 0.95, depth 12 | `frm_sample()` | 8.8e-4, 1.10e-3, 1.40e-3 | **1.072**, 1.144, 1.174 | 55 | 0 | 6000 of 6000 | 3550 s |
| field integrated at brms's line, n 40 (review) | Stan | 8.7e-4, 1.07e-3, 1.36e-3 | 1.00 | | | | |
| y x 1000, n 60, 2000, delta 0.9, depth 12 | `frm_sample()` | 90.3, 106.3, 126.5 | **1.002**, 1.007, 1.005 | 3816 | 44 | 0 | 341 s |
| the same data | brms, its own line | 96, 260, 642 | Inf, 10.0, 3.9 | 4 | 186 | 3814 of 4000 | 2021 s |

**The shortfall, plainly:** at true sigma 1e-3 the drawn route does NOT reach your Rhat
1.01 on sigma, at n 15 (1.025 after 58 minutes) or at n 40 (1.072 after 59 minutes), and
sdgp and lscale read worse. The medians are now right: 1.05e-3 against 1.04e-3 at n 15,
1.10e-3 against 1.07e-3 at n 40, where the sorted factor read 6.5e-4 at n 15. The n 15
run is level with brms (Rhat 1.025 against 1.031, ESS 191 against 152, at about the same
clock), with 508 divergent at delta 0.95 where brms at 0.99 has none, so I would not call
it better than brms: that comparison has two different deltas and depths. y times 1000
converges (1.002), as it did before. Rhat 1.01 at 1e-3 would take longer chains than I
ran; I did not run them.

Other routes through the order (`dev/gpgap-order-check.R`, section 9, and tests):

- `test-reparam.R`'s identity at each block's own nugget:
  centered minus non-centered minus log |L| reads at most 2.8e-3 of a log
  |L| of 150 to 1341 at the nugget of 1e-12, the rounding of the centered density at that
  condition, and at nugget 1e-6 the test holds it to 1e-8 relative. log |L| is a
  determinant now, not a sum over a diagonal.
- `laplace = TRUE`: n 15 sigma 1.03e-3, Rhat 1.001; n 40 sigma 1.07e-3, Rhat 1.001; y x
  1000 sigma 106.7, Rhat 1.001. Unchanged against the target.
- Kriging at new positions, and newdata at observed rows against no newdata: difference 0.
  A block's order moves no number the kriging reports (test-gp-nugget-floor.R).
- `by =`, `gr = FALSE` with rows tied, tied positions, two terms, anisotropic: each block
  has its own permutation, and the map runs.
- `zgp`: the stanfit's z is in that order's basis. It was never brms's `zgp` coordinate by
  coordinate (brms's order is first appearance, frmtmb's was sorted). The draws never
  named z. Stated in `?frm_sample` and NEWS.
- L1 on the five pairs, against brms's own line: 2.19e-3, 9.44e-2, -2.49e-3, 2.65e-2,
  3.47e-2, against 2.2e-3, 9.5e-2, -2.9e-3, 2.6e-2, 3.5e-2 over sorted positions. The 1e-6
  tape against the rewritten program: at most 1.1e-13. The verdicts do not read this route.

### F2: `frm_sample(fit)` from a collapsed fit

The default init is `"random"` when the run moves an exact block's nugget, and the call
says so; `init = "last.par.best"` is honored; a fit with no exact gp keeps the mode.
Same data, `frm()` sigma 5.38e-8, 4 chains of 3000, delta 0.95, depth 12 (`r2-low-fit-15`):
sigma 5.4e-4, 1.01e-3, 2.33e-3, Rhat 1.067, sdgp 1.035, lscale 1.079, 18 divergent, 1097 s.
It read 6.9e-8 in the review. `frm_sample(fit, laplace = TRUE)` from the same fit: 1.03e-3,
Rhat 1.001, 41 s.

### F3

`?"frmtmb-gp"`, frmtmb.sample's NEWS and the warning's message now say that `laplace =
TRUE` serves the outer parameters only and refuses `posterior_epred()`,
`posterior_linpred()`, `posterior_predict()`, `fitted()`, `predict()`, `log_lik()`, `loo()`
and `conditional_effects()`, and that the drawn route is the one for the field. I did not
re-measure the refusals: they are the review's (`lap.txt`), and F1 and F2 do not touch them.

### F4

NEWS states the clock. Measured on this build (`dev/review-gpgap-lapcost.R`, seed 6701, 2
chains of 300, ordinary noise, one run of each, side by side): 25 s against 16 s, 233 s
against 64 s, 1233 s against 388 s at n 40, 100, 200, so **1.6, 3.6 and 3.2 times**. The
review read 1.7, 2.2 and 3.8 on the factor over sorted positions. I put my numbers in NEWS,
since they are of this build; one run each, on a loaded machine, so read them as about 2 to
4 times.

### F5

The binomial wording (help and NEWS), the withdrawn explanation in section 4, the workflow
comments ("FIVE", "All five") and the step's position in `ubuntu-next.yaml` are fixed.
The backlog item is rewritten with these numbers.

### Mutants

`dev/gpgap-mut2.R` on `test-reparam.R`: no_order, no_nugget, first_last, start_medoid,
one_dim, no_start_fix, no_random all killed, the control survives (section 9). I reran
the warning's mutants at fix round 1; the rule did not change here, and the message's text
did, in a part they do not read.

### For consolidation

- `R/covstruct.R`: lane ncgp replaces `chol_L` and moves `gp_nugget`, so the hunks need a
  hand merge, and ncgp's sorted factor should take this order. `gp_farthest_first()` lives
  in frmtmb.sample; if core needs it, it moves.
- `R/predict.R`, kriging; `extensions/frmtmb.sample/R/sample.R`.
- `dev/equiv-data/adjudication.tsv`: section 2.

## 9. Numbers

Generated by `dev/gpgap-summary.R` from the logs of the final build; long lines are logs.

### perm (dev/gpgap-perm.R, base build, before the fix; p094 printed to the terminal: A -94.434083, A1 -0.028199)

```
GPPERM gp_pred() x mvbf: positions 197, sdgp 0.9369, rho 0.137; eig<1e-6 91, eig<1e-12 0
  sorted x already: FALSE
  A  translator as is        -1126.629990
  A1 z through frm kernel, brms order  -5.116873
  kernel off-diagonal, brms scale and order: 3.33e-16
GPPERM gp_pred() x prior: positions 196, sdgp 1.023, rho 0.2018; eig<1e-6 91, eig<1e-12 0
  sorted x already: FALSE
  A  translator as is        -527.367504
  A1 z through frm kernel, brms order  -1.544724
  kernel off-diagonal, brms scale and order: 3.33e-16
GPPERM student x gp: positions 195, sdgp 0.7913, rho 0.1514; eig<1e-6 91, eig<1e-12 0
  sorted x already: FALSE
  A  translator as is        -180.040681
  A1 z through frm kernel, brms order  -0.485363
  kernel off-diagonal, brms scale and order: 2.22e-16
GPPERM gp x nl: positions 45, sdgp 0.7672, rho 0.2697; eig<1e-6 15, eig<1e-12 0
  sorted x already: FALSE
  A  translator as is        -151.650789
  A1 z through frm kernel, brms order  -0.530673
  kernel off-diagonal, brms scale and order: 1.67e-16
```

### nugget line (dev/gpgap-nugget.R)

```
GPNUG lognormal x gp: ours -224.5175, tol 1e-6*|ours| 0.000225, cond(K+nug) 1.23e+08
  brms's program, 1e-12        constant -0.028199
  frmtmb's nugget in that line constant 2.351e-09
  gradient in zgp at frmtmb's nugget, max abs 3.320e-10
GPNUG gp_pred() x mvbf: ours -567.6701, tol 1e-6*|ours| 0.000568, cond(K+nug) 1.09e+08
  brms's program, 1e-12        constant -5.116873
  frmtmb's nugget in that line constant -2.021e-10
  gradient in zgp at frmtmb's nugget, max abs 5.336e-10
GPNUG gp_pred() x prior: ours -267.2522, tol 1e-6*|ours| 0.000267, cond(K+nug) 1.2e+08
  brms's program, 1e-12        constant -1.544724
  frmtmb's nugget in that line constant 6.707e-01
  gradient in zgp at frmtmb's nugget, max abs 1.191e-09
GPNUG student x gp: ours -306.0649, tol 1e-6*|ours| 0.000306, cond(K+nug) 1.02e+08
  brms's program, 1e-12        constant -0.485363
  frmtmb's nugget in that line constant -3.804e-09
  gradient in zgp at frmtmb's nugget, max abs 1.875e-09
GPNUG gp x nl: ours -170.7576, tol 1e-6*|ours| 0.000171, cond(K+nug) 2.72e+07
  brms's program, 1e-12        constant -0.530673
  frmtmb's nugget in that line constant 1.670e-09
  gradient in zgp at frmtmb's nugget, max abs 7.095e-10
```

### extended precision (dev/gpgap-mpfr.R)

```
GPMPFR gp_pred() x prior
  A1 exact -2.19831887, through Stan in double -2.19830659, difference -1.228e-05
  max |L_brms z - f| 2.229e-01 (sd of the field's residual scale: sigma 0.3714)
  C exact 520.431568, in R's double 520.602229
  sum z_b^2 exact 1463.927201, double 1463.513005; sum z_f^2 exact 6.678584, double 6.678584
  ll at Stan's L z: z_b from R double -82.356254; z_b exact -82.356136; ll(f) -82.356292
  smallest pivot of L_brms 1.021e-06, of L_frm 1.046e-03
  R double factor of K_brms against the exact one: max 1.663e-05
total 82 s
GPMPFR lognormal x gp
  A1 exact -0.02819900, through Stan in double -0.02819902, difference 1.832e-08
  max |L_brms z - f| 1.115e-02 (sd of the field's residual scale: sigma 0.2885)
  C exact 1013.506613, in R's double 1013.508670
  sum z_b^2 exact 281.806813, double 281.809741; sum z_f^2 exact 6.436307, double 6.436307
  ll at Stan's L z: z_b from R double -41.187390; z_b exact -41.187394; ll(f) -41.187394
  smallest pivot of L_brms 1.017e-06, of L_frm 6.036e-04
  R double factor of K_brms against the exact one: max 7.387e-07
total 86 s
```

### the review's numbers (dev/gpgap-review83.R)

```
GP83 relative nugget 1e-06: conv 0, old order -94.4341, brms order -0.028199
GP83 relative nugget 1e-09: conv 1, old order -83.8141, brms order -0.020255
GP83 relative nugget 1e-12: conv 1, old order -83.7425, brms order -0.019696
GP83 brms's nugget in frmtmb's kernel: conv 1, old order -83.9427, brms order -4.464e-03
```

### identities (dev/gpgap-identity.R, seed 20261102)

```
GPID gaussian gp(x)              old -2203.993074  brms -1.39779e-01  frmnug  2.43972e-10  grad 4.23e-10  ours -34.664  cond 5.8e+06 conv 0
GPID lognormal gp(x)             old -2203.994142  brms -1.39777e-01  frmnug  9.16600e-12  grad 1.35e-10  ours 69.496  cond 5.8e+06 conv 0
GPID student gp(x)               old  -369.800004  brms -1.89638e-02  frmnug -1.97666e-10  grad 5.13e-11  ours -50.245  cond 5.4e+06 conv 0
GPID poisson gp(x)               old  -136.554073  brms -3.48125e-04  frmnug  1.13687e-12  grad 2.23e-10  ours -144.600  cond 2.3e+06 conv 0
GPID gp(x, z) iso                old -3287.937971  brms -2.35802e-02  frmnug -1.23862e-10  grad 7.56e-11  ours -82.285  cond 1.8e+07 conv 0
GPID gp(x, z, iso = FALSE)       old -3380.064430  brms -2.52857e-02  frmnug -4.71942e-11  grad 9.58e-11  ours -82.394  cond 1.8e+07 conv 0
GPID gp(x, by = f)               old  -753.269032  brms -3.43833e-03  frmnug  8.87113e-11  grad 8.31e-11  ours -48.317  cond 5.5e+06,2.1e+06 conv 0
GPID gp(x, by = f, cmc = FALSE)  old -4191.640219  brms -4.88766e-02  frmnug -1.80165e-10  grad 2.38e-10  ours -48.607  cond 5.7e+06,4.9e+06 conv 0
GPID gp(x, by = w)               old -2283.401850  brms -1.07687e-01  frmnug  4.54605e-11  grad 1.56e-10  ours -32.181  cond 5.3e+06 conv 0
GPID gp(xu, gr = FALSE)          old   -21.516021  brms -2.15160e+01  frmnug  7.07089e-10  grad 1.46e-09  ours -96.839  cond 3.7e+07 conv 0
GPID gp(x, scale = FALSE)        old -2203.993074  brms -1.39779e-01  frmnug  2.43958e-10  grad 5.96e-11  ours -34.664  cond 5.8e+06 conv 0
GPID sigma ~ gp(x)               old   -98.600312  brms -1.00321e-05  frmnug  0.00000e+00  grad 7.02e-12  ours -45.058  cond 45 conv 0
GPID gp(x) + gp(z)               old -4235.489513  brms -2.48263e-01  frmnug -1.12266e-11  grad 4.85e-10  ours -46.358  cond 5.9e+06,6.6e+06 conv 0
GPID nl b1 ~ gp(z)               old  -357.407767  brms -1.06596e-02  frmnug -7.95382e-11  grad 5.21e-11  ours -32.252  cond 5.4e+06 conv 0
GPID mvbf gp(x), gp(z)           old -3003.813868  brms -1.50576e-01  frmnug -4.31726e-11  grad 4.44e-10  ours -66.119  cond 5.8e+06,5.9e+06 conv 0
```

### nugget variants (dev/gpgap-nugvar.R)

```
GPVAR lognormal x gp 
  rel 1e-6       conv 0 (relative convergence (4)) gr 7.49e-05, ll -56.515960, sdgp 0.59595, rho 0.96364 (3.0s)
  rel 1e-9       conv 1 (false convergence (8)) gr 4.19e-03, ll -56.515951, sdgp 0.59596, rho 0.96364 (3.6s)
  rel 1e-12      conv 1 (false convergence (8)) gr 6.79e-01, ll -56.515808, sdgp 0.58270, rho 0.95539 (7.3s)
  abs 1e-12      conv 1 (false convergence (8)) gr 1.46e-01, ll -56.515727, sdgp 0.59860, rho 0.96549 (3.9s)
  abs 1e-12 ncp, default conv 0 (relative convergence (4)) gr 3.91e-05, ll -56.515953, sdgp 0.59595, rho 0.96364 (14.0s)
  abs 1e-12 ncp, at fit  conv 0 (relative convergence (4)) gr 1.30e-04, ll -56.515953, sdgp 0.59594, rho 0.96363 (13.6s)
  closed form rel 1e-06: ll -47.275861, sdgp 0.59595, rho 0.96364, sigma 0.28849, conv 1
  closed form abs 1e-12: ll -47.275853, sdgp 0.59595, rho 0.96364, sigma 0.28849, conv 1
  brms's model read at frmtmb's maximum: ll -47.275853 (loss 8.431e-10)
GPVAR gp_pred() x mvbf 
  rel 1e-6       conv 0 (relative convergence (4)) gr 3.95e-06, ll -401.829337, sdgp 0.93695, rho 0.78683 (5.0s)
  rel 1e-9       conv 1 (false convergence (8)) gr 3.46e-03, ll -401.829336, sdgp 0.93688, rho 0.78681 (6.6s)
  rel 1e-12      conv 1 (false convergence (8)) gr 1.30e-01, ll -401.829733, sdgp 0.92453, rho 0.78129 (7.1s)
  abs 1e-12      conv 1 (false convergence (8)) gr 3.44e-01, ll -401.829106, sdgp 0.92107, rho 0.77919 (7.5s)
  abs 1e-12 ncp, default conv 0 (relative convergence (4)) gr 1.64e-04, ll -401.829337, sdgp 0.93694, rho 0.78683 (17.9s)
  abs 1e-12 ncp, at fit  conv 0 (relative convergence (4)) gr 2.18e-04, ll -401.829337, sdgp 0.93695, rho 0.78683 (14.3s)
  closed form rel 1e-06: ll -121.790942, sdgp 0.93695, rho 0.78683, sigma 0.41233, conv 0
  closed form abs 1e-12: ll -121.790942, sdgp 0.93695, rho 0.78683, sigma 0.41233, conv 0
  brms's model read at frmtmb's maximum: ll -121.790942 (loss 1.269e-09)
GPVAR gp_pred() x prior 
  rel 1e-6       conv 0 (relative convergence (4)) gr 1.09e-04, ll -99.621529, sdgp 1.03020, rho 0.94166 (3.1s)
  rel 1e-9       conv 1 (false convergence (8)) gr 5.06e-03, ll -99.621532, sdgp 1.02999, rho 0.94159 (3.9s)
  rel 1e-12      conv 1 (false convergence (8)) gr 9.47e-01, ll -99.620006, sdgp 1.02250, rho 0.93794 (5.1s)
  abs 1e-12      conv 1 (false convergence (8)) gr 7.91e-01, ll -99.625768, sdgp 1.07317, rho 0.95876 (5.9s)
  abs 1e-12 ncp, default conv 0 (relative convergence (4)) gr 1.27e-06, ll -99.621534, sdgp 1.03019, rho 0.94166 (14.4s)
  abs 1e-12 ncp, at fit  conv 0 (relative convergence (4)) gr 3.89e-06, ll -99.621534, sdgp 1.03019, rho 0.94166 (13.7s)
  closed form rel 1e-06: ll -99.621529, sdgp 1.03020, rho 0.94166, sigma 0.37141, conv 0
  closed form abs 1e-12: ll -99.621534, sdgp 1.03019, rho 0.94166, sigma 0.37141, conv 0
  brms's model read at frmtmb's maximum: ll -99.621534 (loss 2.963e-09)
GPVAR student x gp 
  rel 1e-6       conv 0 (relative convergence (4)) gr 7.47e-06, ll -140.833313, sdgp 0.79129, rho 0.88608 (3.5s)
  rel 1e-9       conv 1 (false convergence (8)) gr 7.07e-04, ll -140.833301, sdgp 0.79124, rho 0.88605 (4.0s)
  rel 1e-12      conv 1 (false convergence (8)) gr 5.38e-01, ll -141.133590, sdgp 0.62296, rho 0.57670 (5.5s)
  abs 1e-12      conv 1 (false convergence (8)) gr 8.90e-02, ll -141.133721, sdgp 0.63250, rho 0.58657 (5.5s)
  abs 1e-12 ncp, default conv 0 (relative convergence (4)) gr 1.53e-05, ll -140.833303, sdgp 0.79129, rho 0.88608 (17.0s)
  abs 1e-12 ncp, at fit  conv 0 (relative convergence (4)) gr 3.56e-05, ll -140.833303, sdgp 0.79130, rho 0.88608 (17.0s)
GPVAR gp x nl 
  rel 1e-6       conv 0 (relative convergence (4)) gr 3.28e-05, ll -138.502817, sdgp 0.76723, rho 1.48322 (1.8s)
  rel 1e-9       conv 0 (relative convergence (4)) gr 2.38e-04, ll -138.502795, sdgp 0.76725, rho 1.48324 (0.2s)
  rel 1e-12      conv 1 (false convergence (8)) gr 1.50e-01, ll -138.502853, sdgp 0.75345, rho 1.47453 (0.6s)
  abs 1e-12      conv 1 (false convergence (8)) gr 8.11e-03, ll -138.502410, sdgp 0.76533, rho 1.48115 (0.5s)
  abs 1e-12 ncp, default conv 0 (relative convergence (4)) gr 5.61e-06, ll -138.502795, sdgp 0.76724, rho 1.48323 (0.4s)
  abs 1e-12 ncp, at fit  conv 0 (relative convergence (4)) gr 1.65e-04, ll -138.502795, sdgp 0.76724, rho 1.48323 (0.4s)
```

### low noise (dev/gpgap-lownoise.R, seed 20261102)

```
GPLOW true sigma 1.0e-01 | closed form frm nugget: sigma 1.0877e-01 sdgp 1.0747 rho 0.9506 ll 25.2288 | brms nugget: sigma 1.0878e-01 sdgp 1.0747 rho 0.9506 ll 25.2288 | brms model at frm max loses 1.213e-07 | frm(): sigma 1.0877e-01 sdgp 1.0747 conv 0
  ncp at brms nugget: conv 0 (relative convergence) gr 1.1e-06 ll 25.2288 sdgp 1.0747 rho 0.9506 
  gp(x, k = 30): sigma 1.0878e-01 sdgp 1.0747 ll 25.2288 conv 0
GPLOW true sigma 1.0e-02 | closed form frm nugget: sigma 1.0769e-02 sdgp 2.0003 rho 1.3039 ll 141.8048 | brms nugget: sigma 1.0954e-02 sdgp 2.0003 rho 1.3039 ll 141.8048 | brms model at frm max loses 1.461e-02 | frm(): sigma 1.0769e-02 sdgp 2.0003 conv 0
  ncp at brms nugget: conv 0 (relative convergence) gr 3.9e-06 ll 141.8048 sdgp 2.0003 rho 1.3039 
  gp(x, k = 30): sigma 1.0954e-02 sdgp 2.0003 ll 141.8048 conv 0
GPLOW true sigma 3.0e-03 | closed form frm nugget: sigma 2.5035e-03 sdgp 2.1935 rho 1.3586 ll 201.1157 | brms nugget: sigma 3.3285e-03 sdgp 2.1934 rho 1.3586 ll 201.1157 | brms model at frm max loses 4.831e+00 | frm(): sigma 2.5035e-03 sdgp 2.1934 conv 0
  ncp at brms nugget: conv 0 (relative convergence) gr 3.0e-05 ll 201.1157 sdgp 2.1935 rho 1.3586 
  gp(x, k = 30): sigma 3.3285e-03 sdgp 2.1935 ll 201.1157 conv 0
GPLOW true sigma 1.0e-03 | closed form frm nugget: sigma 9.0168e-08 sdgp 1.1908 rho 1.1823 ll 252.3122 | brms nugget: sigma 1.1166e-03 sdgp 2.4066 rho 1.3892 ll 254.0843 | brms model at frm max loses 2.392e+07 | frm(): sigma 1.2725e-07 sdgp 1.1908 conv 0
  ncp at brms nugget: conv 0 (relative convergence) gr 3.0e-05 ll 254.0843 sdgp 2.4068 rho 1.3892 
  gp(x, k = 30): sigma 1.1166e-03 sdgp 2.4067 ll 254.0842 conv 0
GPLOW true sigma 3.0e-04 | closed form frm nugget: sigma 1.7977e-08 sdgp 0.5054 rho 0.9179 ll 291.1191 | brms nugget: sigma 3.3125e-04 sdgp 4.2198 rho 1.5783 ll 311.7284 | brms model at frm max loses 1.919e+06 | frm(): sigma 2.6944e-08 sdgp 0.5054 conv 0
  ncp at brms nugget: conv 1 (false convergence (8) gr 2.7e-03 ll 311.7289 sdgp 4.1823 rho 1.5743 
  gp(x, k = 30): sigma 3.3211e-04 sdgp 4.1211 ll 311.7218 conv 0
GPLOW true sigma 1.0e-04 | closed form frm nugget: sigma 1.2142e-08 sdgp 0.3910 rho 0.8408 ll 301.4596 | brms nugget: sigma 1.0714e-04 sdgp 4.8711 rho 1.6516 ll 364.5176 | brms model at frm max loses 2.036e+05 | frm(): sigma 1.7596e-08 sdgp 0.3909 conv 0
  ncp at brms nugget: conv 0 (relative convergence) gr 8.0e-05 ll 364.6165 sdgp 5.6685 rho 1.6870 
  gp(x, k = 30): sigma 1.0948e-04 sdgp 5.6755 ll 364.5729 conv 0
```

### equiv pairs against the baseline

| pair | baseline verdict (level) | L1 const | now (level) | L1 const | at the 2 points | max grad | brms's own line | s |
|---|---|---|---|---|---|---|---|---|
| p022 | cannot measure (NA) | -1127 | divergence by decision (NA) | -5.12 | NA | 5.3e-10 | -5.117 | 40 |
| p023 | cannot measure (NA) | -527.4 | divergence by decision (NA) | -1.54 | NA | 1.2e-09 | -2.215 | 32 |
| p094 | cannot measure (NA) | -94.43 | divergence by decision (NA) | -0.0282 | NA | 3.3e-10 | -0.0282 | 32 |
| p109 | cannot measure (NA) | -180 | divergence by decision (NA) | -0.485 | NA | 1.9e-09 | -0.4854 | 36 |
| p136 | equivalent (L4) | -151.7 | equivalent (L4) | -0.531 | NA | 7.1e-10 | -0.5307 | 363 |

### the warning's calibration (dev/gpgap-calib.R, dev/gpgap-calib-sum.R)

```
rows 281 unreadable 1 errors 0 
old       fires   0 | of 152 moved: hit   0 miss 152 | of 128 unmoved: false   0
r_med     fires 153 | of 152 moved: hit 148 miss   4 | of 128 unmoved: false   5
r_mean    fires 155 | of 152 moved: hit 149 miss   3 | of 128 unmoved: false   6
r_q90     fires 155 | of 152 moved: hit 150 miss   2 | of 128 unmoved: false   5
scan      fires 154 | of 152 moved: hit 149 miss   3 | of 128 unmoved: false   5
fix_0.1   fires 165 | of 152 moved: hit 149 miss   3 | of 128 unmoved: false  16
fix_0.25  fires 143 | of 152 moved: hit 141 miss  11 | of 128 unmoved: false   2
the rule is r_q90
MISS:
       poisson log mean 10  59    0  FALSE 5.42e-02 1.48e-01 3.09e-01 1.62e+00 1.98e+00   0.8800   0.515  0.515       0 4.9e-05 FALSE  TRUE FALSE
     sigma~(1|g) noise 0.1  59    1  FALSE 1.13e-04 3.21e-04 7.18e-04 1.03e+00 1.03e+00   7.2000  17.100  0.002       0 7.5e-09  TRUE  TRUE FALSE
FALSE ALARM:
    floor n200 r0.007 s101 200    0  FALSE 1.06e-01 1.06e-01 1.06e+00 2.22e+00 2.22e+00   0.0000   0.977  0.000       0 4.7e-05 FALSE FALSE  TRUE
           rep 30x10 r0.01  30    0  FALSE 3.88e-01 3.88e-01 1.50e+00 2.10e+00 2.28e+00  -0.1640   0.264  0.221       0 1.3e-05 FALSE FALSE  TRUE
    scale A1000 n60 r0.005  60    0  FALSE 2.19e-01 2.19e-01 1.20e+00 1.79e+03 1.79e+03   0.0000   0.924  0.000       0 2.1e-04 FALSE FALSE  TRUE
    scale A1e-3 n60 r0.005  60    0  FALSE 1.99e-01 1.99e-01 1.09e+00 2.21e-03 2.21e-03   0.0000   0.897  0.000       0 2.5e-05 FALSE FALSE  TRUE
     sigma ~ z noise 0.003  59    0  FALSE 5.78e-01 1.35e+00 2.92e+00 1.92e+00 2.28e+00   0.3860   0.994  0.350       0 1.2e-04 FALSE FALSE  TRUE
one row per position (floor, scale, xs10): 97 of 168 fire, 3 of those below 1 SE
```

### binomial at 1e5 trials (dev/gpgap-binom1e5.R)

```
frm(): sdgp 2.5531 rho 1.3857 ll -444.0001
default start, 1e-12 abs           ll -444.8539 sdgp 2.0007 rho 1.3174 conv 0 gr 1.4e-04
from frm(), 1e-12 abs              ll -444.8539 sdgp 2.0007 rho 1.3174 conv 0 gr 4.2e-05
sd 2 rho 1, 1e-12 abs              ll -444.8539 sdgp 2.0007 rho 1.3174 conv 0 gr 4.5e-05
default start, 1e-12 rel           ll -444.8539 sdgp 2.0007 rho 1.3174 conv 0 gr 1.4e-04
from frm(), 1e-12 rel              ll -444.8539 sdgp 2.0007 rho 1.3174 conv 0 gr 5.9e-05
sd 2 rho 1, 1e-12 rel              ll -444.8539 sdgp 2.0007 rho 1.3174 conv 0 gr 2.1e-05
default start, 1e-09 rel           ll -444.8530 sdgp 2.0013 rho 1.3175 conv 0 gr 1.4e-04
from frm(), 1e-09 rel              ll -444.8530 sdgp 2.0013 rho 1.3175 conv 0 gr 5.3e-05
sd 2 rho 1, 1e-09 rel              ll -444.8530 sdgp 2.0013 rho 1.3175 conv 0 gr 3.8e-05
default start, 1e-06 rel           ll -444.0001 sdgp 2.5530 rho 1.3857 conv 0 gr 6.9e-06
from frm(), 1e-06 rel              ll -444.0001 sdgp 2.5531 rho 1.3857 conv 0 gr 2.2e-05
sd 2 rho 1, 1e-06 rel              ll -444.0001 sdgp 2.5530 rho 1.3857 conv 0 gr 8.6e-05
z of the residual at frm(): sd 0.946 (1 if binomial alone)
nugget variance 6.52e-06, median working noise 4.61e-05
```

### mutants of the warning (dev/gpgap-mut.R)

```
MUT class_drop     pass=56 fail=1 error=1 warn=0 KILLED
MUT control        pass=65 fail=0 error=0 warn=0 SURVIVED
MUT first_only     pass=63 fail=0 error=1 warn=0 KILLED
MUT msg_k_promise  pass=64 fail=1 error=0 warn=0 KILLED
MUT msg_no_moves   pass=64 fail=1 error=0 warn=0 KILLED
MUT no_prior       pass=47 fail=14 error=1 warn=12 KILLED
MUT no_trycatch    pass=63 fail=0 error=1 warn=0 KILLED
MUT no_warning     pass=56 fail=1 error=1 warn=0 KILLED
MUT nug_x10        pass=62 fail=3 error=0 warn=3 KILLED
MUT one_block_loop pass=61 fail=0 error=1 warn=0 KILLED
MUT q_max          pass=64 fail=1 error=0 warn=1 KILLED
MUT q_med          pass=64 fail=1 error=0 warn=0 KILLED
MUT sd1            pass=59 fail=2 error=1 warn=0 KILLED
MUT thr_double     pass=64 fail=1 error=0 warn=0 KILLED
MUT thr_fixed      pass=62 fail=3 error=0 warn=3 KILLED
MUT thr_half       pass=64 fail=1 error=0 warn=1 KILLED
MUT thr_x30        pass=55 fail=5 error=2 warn=0 KILLED
MUT thr_x5         pass=59 fail=2 error=1 warn=0 KILLED
```

### the sampling route (dev/gpgap-sampler.R, dev/gpgap-sampler-l1.R, dev/gpgap-marg15.R)

```
GPSAMPL1 lognormal x gp: sdgp 0.5959, ours -224.548, tol 0.00022 | sampler tape against brms's own line 2.187e-03 | tape at 1e-6 against brms rewritten -2.842e-14
GPSAMPL1 gp_pred() x mvbf: sdgp 0.9369, ours -572.882, tol 0.00057 | sampler tape against brms's own line 9.498e-02 | tape at 1e-6 against brms rewritten 1.137e-13
GPSAMPL1 gp_pred() x prior: sdgp 1.03, ours -268.003, tol 0.00027 | sampler tape against brms's own line -2.864e-03 | tape at 1e-6 against brms rewritten 5.684e-14
GPSAMPL1 student x gp: sdgp 0.7913, ours -306.577, tol 0.00031 | sampler tape against brms's own line 2.649e-02 | tape at 1e-6 against brms rewritten 1.137e-13
GPSAMPL1 gp x nl: sdgp 0.7672, ours -171.323, tol 0.00017 | sampler tape against brms's own line 3.473e-02 | tape at 1e-6 against brms rewritten 8.527e-14
MARG n 15 s 1e-03 frm_sample before, rel 1e-6  sigma 6.431e-05 0.0007348 0.002796 (rhat 1.000) sdgp 0.8714 1.488 2.623 (rhat 1.000) divergent 0
MARG n 15 s 1e-03 frm_sample now, rel 1e-12    sigma 0.0005382 0.001028 0.002379 (rhat 1.002) sdgp 1.091 2.132 4.194 (rhat 1.001) divergent 0
MARG n 15 s 1e-03 brms, abs 1e-12              sigma 0.0005197 0.001036 0.002407 (rhat 1.000) sdgp 1.09 2.121 4.242 (rhat 1.001) divergent 0
-- round 2, the installed build, farthest-first factor
GPSAMPL1 lognormal x gp: sdgp 0.5959, ours -224.548, tol 0.00022 | sampler tape against brms's own line 2.188e-03 | tape at 1e-6 against brms rewritten -2.842e-14
GPSAMPL1 gp_pred() x mvbf: sdgp 0.9369, ours -572.881, tol 0.00057 | sampler tape against brms's own line 9.442e-02 | tape at 1e-6 against brms rewritten 1.137e-13
GPSAMPL1 gp_pred() x prior: sdgp 1.03, ours -268.003, tol 0.00027 | sampler tape against brms's own line -2.492e-03 | tape at 1e-6 against brms rewritten 5.684e-14
GPSAMPL1 student x gp: sdgp 0.7913, ours -306.577, tol 0.00031 | sampler tape against brms's own line 2.650e-02 | tape at 1e-6 against brms rewritten 1.137e-13
GPSAMPL1 gp x nl: sdgp 0.7672, ours -171.323, tol 0.00017 | sampler tape against brms's own line 3.471e-02 | tape at 1e-6 against brms rewritten 8.527e-14
GPSAMP low frm_sample n 15 seed 6601 sigma      q05 0.0005925 q50 0.001048 q95 0.002317 rhat 1.025 ess 191
GPSAMP low frm_sample n 15 seed 6601 sdgp_gpx   q05 1.101 q50 1.997 q95 3.979 rhat 1.069 ess 47
GPSAMP low frm_sample n 15 seed 6601 lscale_gpx q05 0.2107 q50 0.2502 q95 0.2998 rhat 1.061 ess 65
GPSAMP low frm_sample n 15 seed 6601 iter 12000 chains 4 delta 0.950 td 12 metric diag_e divergent 508 treedepth 23420 wall 3474
GPSAMP low frm_sample n 40 seed 6601 sigma      q05 0.0008769 q50 0.0011 q95 0.001395 rhat 1.072 ess 55
GPSAMP low frm_sample n 40 seed 6601 sdgp_gpx   q05 1.415 q50 2.13 q95 3.855 rhat 1.144 ess 24
GPSAMP low frm_sample n 40 seed 6601 lscale_gpx q05 0.2037 q50 0.2364 q95 0.2642 rhat 1.174 ess 18
GPSAMP low frm_sample n 40 seed 6601 iter 3000 chains 4 delta 0.950 td 12 metric diag_e divergent 0 treedepth 6000 wall 3550
frm() sigma 5.38e-08
GPSAMP low fit n 15 seed 6601 sigma      q05 0.0005448 q50 0.001014 q95 0.002329 rhat 1.067 ess 59
GPSAMP low fit n 15 seed 6601 sdgp_gpx   q05 1.199 q50 1.931 q95 3.117 rhat 1.035 ess 53
GPSAMP low fit n 15 seed 6601 lscale_gpx q05 0.2104 q50 0.2432 q95 0.2797 rhat 1.079 ess 44
GPSAMP low fit n 15 seed 6601 iter 3000 chains 4 delta 0.950 td 12 metric diag_e divergent 18 treedepth 5962 wall 1097
frm() sigma 5.38e-08
GPSAMP low fitlap n 15 seed 6601 sigma      q05 0.0005238 q50 0.001025 q95 0.002347 rhat 1.001 ess 2220
GPSAMP low fitlap n 15 seed 6601 sdgp_gpx   q05 1.12 q50 2.131 q95 4.273 rhat 1.002 ess 1496
GPSAMP low fitlap n 15 seed 6601 lscale_gpx q05 0.2018 q50 0.252 q95 0.3023 rhat 1.002 ess 1581
GPSAMP low fitlap n 15 seed 6601 iter 2000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 41
GPSAMP low laplace n 15 seed 6601 sigma      q05 0.0005179 q50 0.001026 q95 0.002415 rhat 1.001 ess 4277
GPSAMP low laplace n 15 seed 6601 sdgp_gpx   q05 1.107 q50 2.148 q95 4.227 rhat 1.001 ess 2982
GPSAMP low laplace n 15 seed 6601 lscale_gpx q05 0.2022 q50 0.2528 q95 0.3031 rhat 1.001 ess 3075
GPSAMP low laplace n 15 seed 6601 iter 4000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 56
GPSAMP low laplace n 40 seed 6601 sigma      q05 0.0008715 q50 0.001072 q95 0.001366 rhat 1.001 ess 4812
GPSAMP low laplace n 40 seed 6601 sdgp_gpx   q05 1.439 q50 2.599 q95 5.106 rhat 1.001 ess 3197
GPSAMP low laplace n 40 seed 6601 lscale_gpx q05 0.2078 q50 0.2448 q95 0.2824 rhat 1.001 ess 3246
GPSAMP low laplace n 40 seed 6601 iter 4000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 143
GPSAMP units frm_sample n 60 seed 7701 sigma      q05 90.32 q50 106.3 q95 126.5 rhat 1.002 ess 3816
GPSAMP units frm_sample n 60 seed 7701 sdgp_gpx   q05 635.7 q50 1183 q95 2523 rhat 1.007 ess 672
GPSAMP units frm_sample n 60 seed 7701 lscale_gpx q05 0.09189 q50 0.1561 q95 0.2257 rhat 1.005 ess 517
GPSAMP units frm_sample n 60 seed 7701 iter 2000 chains 4 delta 0.900 td 12 metric diag_e divergent 44 treedepth 0 wall 341
GPSAMP units laplace n 60 seed 7701 sigma      q05 89.7 q50 106.7 q95 128.3 rhat 1.001 ess 2438
GPSAMP units laplace n 60 seed 7701 sdgp_gpx   q05 623.9 q50 1184 q95 2792 rhat 1.003 ess 1081
GPSAMP units laplace n 60 seed 7701 lscale_gpx q05 0.09108 q50 0.157 q95 0.23 rhat 1.004 ess 1170
GPSAMP units laplace n 60 seed 7701 iter 2000 chains 4 delta 0.900 td 10 metric diag_e divergent 0 treedepth 0 wall 488
GPORD plain          blocks 1 positions 60 | centered - noncentered - log|L| 6.484e-04 of 618.8 | order a permutation TRUE
GPORD two_terms      blocks 2 positions 60,60 | centered - noncentered - log|L| -2.771e-03 of 1341.0 | order a permutation TRUE
GPORD aniso          blocks 1 positions 60 | centered - noncentered - log|L| 7.366e-07 of 184.1 | order a permutation TRUE
GPORD by_f           blocks 2 positions 24,36 | centered - noncentered - log|L| 3.859e-04 of 424.7 | order a permutation TRUE
GPORD rows_tied      blocks 1 positions 23 | centered - noncentered - log|L| 1.867e-04 of 153.7 | order a permutation TRUE
GPORD gr_false_tied  blocks 1 positions 60 | centered - noncentered - log|L| -4.363e-04 of 669.4 | order a permutation TRUE
GPORD poisson        blocks 1 positions 60 | centered - noncentered - log|L| 8.435e-04 of 542.9 | order a permutation TRUE
GPORD sigma_gp       blocks 1 positions 60 | centered - noncentered - log|L| -1.278e-03 of 1021.0 | order a permutation TRUE
GPORD draws: 60 b columns, first b[1] last b[60]; cor with frm()'s field 0.9999
GPORD kriging at new positions: means 0.373 0.352 -1.065 against sin(2 x) 0.203 0.356 -1.011; sd 0.240 0.176 0.155
GPORD epred at observed rows, newdata against none: max diff 0.00e+00
lapcost, now, then the centered route at 1e-6 (dev/review-gpgap-lapcost.R, one run):
LCOST rellib-r16 n 40: 25 s, sigma median 0.2557 sdgp median 1.009
LCOST rellib-r16 n 100: 233 s, sigma median 0.3232 sdgp median 1.080
LCOST rellib-r16 n 200: 1233 s, sigma median 0.3342 sdgp median 1.056
LCOST rellib-r16 n 40: 16 s, sigma median 0.2600 sdgp median 1.005
LCOST rellib-r16 n 100: 64 s, sigma median 0.3202 sdgp median 1.010
LCOST rellib-r16 n 200: 388 s, sigma median 0.3338 sdgp median 1.011
MUT2 control        pass=423 fail=0 error=0 SURVIVED
MUT2 first_last     pass=406 fail=17 error=0 KILLED
MUT2 no_nugget      pass=416 fail=7 error=0 KILLED
MUT2 no_order       pass=416 fail=7 error=0 KILLED
MUT2 no_random      pass=421 fail=2 error=0 KILLED
MUT2 no_start_fix   pass=422 fail=1 error=0 KILLED
MUT2 one_dim        pass=422 fail=1 error=0 KILLED
MUT2 start_medoid   pass=422 fail=1 error=0 KILLED
-- before: the factor over sorted positions (fix round 1)
GPSAMP low laplace n 15 seed 6601 sigma      q05 0.0005326 q50 0.001034 q95 0.002313 rhat 1.000 ess 4183
GPSAMP low laplace n 15 seed 6601 sdgp_gpx   q05 1.114 q50 2.165 q95 4.289 rhat 1.001 ess 2847
GPSAMP low laplace n 15 seed 6601 lscale_gpx q05 0.2025 q50 0.2524 q95 0.3033 rhat 1.001 ess 2968
GPSAMP low laplace n 15 seed 6601 iter 4000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 168
GPSAMP low baselap n 15 seed 6601 sigma      q05 7.01e-05 q50 0.0007429 q95 0.002763 rhat 1.000 ess 4304
GPSAMP low baselap n 15 seed 6601 sdgp_gpx   q05 0.8534 q50 1.491 q95 2.591 rhat 1.001 ess 2407
GPSAMP low baselap n 15 seed 6601 lscale_gpx q05 0.1829 q50 0.2323 q95 0.2772 rhat 1.001 ess 2491
GPSAMP low baselap n 15 seed 6601 iter 4000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 129
GPSAMP low laplace n 40 seed 6601 sigma      q05 0.0008723 q50 0.001076 q95 0.001366 rhat 1.001 ess 4898
GPSAMP low laplace n 40 seed 6601 sdgp_gpx   q05 1.447 q50 2.617 q95 5.084 rhat 1.001 ess 2987
GPSAMP low laplace n 40 seed 6601 lscale_gpx q05 0.2085 q50 0.2454 q95 0.2832 rhat 1.001 ess 2879
GPSAMP low laplace n 40 seed 6601 iter 4000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 386
GPSAMP low baselap n 40 seed 6601 sigma      q05 2.276e-05 q50 0.0003041 q95 0.0008679 rhat 1.001 ess 2738
GPSAMP low baselap n 40 seed 6601 sdgp_gpx   q05 0.9337 q50 1.187 q95 1.529 rhat 1.001 ess 3932
GPSAMP low baselap n 40 seed 6601 lscale_gpx q05 0.1793 q50 0.2041 q95 0.2266 rhat 1.000 ess 3987
GPSAMP low baselap n 40 seed 6601 iter 4000 chains 4 delta 0.950 td 10 metric diag_e divergent 0 treedepth 0 wall 195
GPSAMP units laplace n 60 seed 7701 sigma      q05 90.38 q50 106.9 q95 127.6 rhat 1.004 ess 1990
GPSAMP units laplace n 60 seed 7701 sdgp_gpx   q05 629.7 q50 1175 q95 2703 rhat 1.001 ess 1509
GPSAMP units laplace n 60 seed 7701 lscale_gpx q05 0.09323 q50 0.1549 q95 0.2279 rhat 1.001 ess 1420
GPSAMP units laplace n 60 seed 7701 iter 2000 chains 4 delta 0.900 td 10 metric diag_e divergent 0 treedepth 0 wall 645
GPSAMP low frm_sample n 15 seed 6601 sigma      q05 0.0004953 q50 0.0006546 q95 0.001792 rhat 2.461 ess 5
GPSAMP low frm_sample n 15 seed 6601 sdgp_gpx   q05 2.315 q50 2.916 q95 4.622 rhat 2.978 ess 5
GPSAMP low frm_sample n 15 seed 6601 lscale_gpx q05 0.255 q50 0.2834 q95 0.305 rhat 2.230 ess 5
GPSAMP low frm_sample n 15 seed 6601 iter 20000 chains 4 delta 0.990 td 13 metric diag_e divergent 0 treedepth 39996 wall 14347
GPSAMP low frm_sample n 15 seed 6601 sigma      q05 0.0003908 q50 0.000428 q95 0.0015 rhat 2.834 ess 5
GPSAMP low frm_sample n 15 seed 6601 sdgp_gpx   q05 1.623 q50 3.422 q95 5.774 rhat 3.088 ess 5
GPSAMP low frm_sample n 15 seed 6601 lscale_gpx q05 0.2527 q50 0.2768 q95 0.3152 rhat 2.824 ess 5
GPSAMP low frm_sample n 15 seed 6601 iter 20000 chains 4 delta 0.950 td 12 metric dense_e divergent 0 treedepth 39999 wall 7769
GPSAMP low brms n 15 seed 6601 sigma      q05 0.0005181 q50 0.001032 q95 0.002381 rhat 1.031 ess 152
GPSAMP low brms n 15 seed 6601 sdgp_gpx   q05 1.167 q50 2.171 q95 4.229 rhat 1.154 ess 20
GPSAMP low brms n 15 seed 6601 lscale_gpx q05 0.2006 q50 0.2539 q95 0.3053 rhat 1.112 ess 28
GPSAMP low brms n 15 seed 6601 iter 20000 chains 4 delta 0.990 td 13 metric diag_e divergent 0 treedepth 39867 wall 3189
GPSAMP units frm_sample n 60 seed 7701 sigma      q05 90.59 q50 106.5 q95 128 rhat 1.001 ess 3110
GPSAMP units frm_sample n 60 seed 7701 sdgp_gpx   q05 633.1 q50 1192 q95 2638 rhat 1.003 ess 584
GPSAMP units frm_sample n 60 seed 7701 lscale_gpx q05 0.09416 q50 0.1566 q95 0.2273 rhat 1.006 ess 400
GPSAMP units frm_sample n 60 seed 7701 iter 2000 chains 4 delta 0.900 td 12 metric diag_e divergent 56 treedepth 0 wall 708
GPSAMP units brms n 60 seed 7701 sigma      q05 96.09 q50 260.2 q95 641.8 rhat Inf ess 4
GPSAMP units brms n 60 seed 7701 sdgp_gpx   q05 196.8 q50 490.5 q95 564.9 rhat 10.012 ess 4
GPSAMP units brms n 60 seed 7701 lscale_gpx q05 0.04017 q50 0.04148 q95 0.04261 rhat 3.890 ess 4
GPSAMP units brms n 60 seed 7701 iter 2000 chains 4 delta 0.900 td 12 metric diag_e divergent 186 treedepth 3814 wall 2021
GPSAMP low frm_sample n 10 seed 6601 sigma      q05 0.0008107 q50 0.003229 q95 0.02227 rhat 1.153 ess 19
GPSAMP low frm_sample n 10 seed 6601 sdgp_gpx   q05 0.7467 q50 1.559 q95 3.983 rhat 1.104 ess 49
GPSAMP low frm_sample n 10 seed 6601 lscale_gpx q05 0.1407 q50 0.2188 q95 0.2936 rhat 1.087 ess 40
GPSAMP low frm_sample n 10 seed 6601 iter 2000 chains 4 delta 0.990 td 13 metric dense_e divergent 407 treedepth 1981 wall 698
GPSAMP low brms n 10 seed 6601 sigma      q05 0.0007933 q50 0.003334 q95 0.02548 rhat 1.059 ess 57
GPSAMP low brms n 10 seed 6601 sdgp_gpx   q05 0.7077 q50 1.313 q95 2.517 rhat 1.066 ess 108
GPSAMP low brms n 10 seed 6601 lscale_gpx q05 0.1242 q50 0.2126 q95 0.2722 rhat 1.065 ess 73
GPSAMP low brms n 10 seed 6601 iter 2000 chains 4 delta 0.990 td 13 metric dense_e divergent 211 treedepth 2503 wall 253
GPSAMP low frm_sample n 15 seed 6601 sigma      q05 0.0008001 q50 0.001368 q95 0.003044 rhat 1.058 ess 33
GPSAMP low frm_sample n 15 seed 6601 sdgp_gpx   q05 2.275 q50 2.968 q95 4.357 rhat 1.912 ess 3
GPSAMP low frm_sample n 15 seed 6601 lscale_gpx q05 0.2359 q50 0.2638 q95 0.304 rhat 1.906 ess 3
GPSAMP low frm_sample n 15 seed 6601 iter 1200 chains 2 delta 0.950 td 12 metric dense_e divergent 70 treedepth 1126 wall 381
GPSAMP low brms n 15 seed 6601 sigma      q05 0.0005203 q50 0.0005275 q95 0.002053 rhat 2.210 ess 3
GPSAMP low brms n 15 seed 6601 sdgp_gpx   q05 1.352 q50 2.935 q95 2.948 rhat 2.130 ess 3
GPSAMP low brms n 15 seed 6601 lscale_gpx q05 0.2033 q50 0.278 q95 0.2782 rhat 2.258 ess 3
GPSAMP low brms n 15 seed 6601 iter 1200 chains 2 delta 0.950 td 12 metric dense_e divergent 597 treedepth 584 wall 196
```

### k and the rest

```
KS k  15 default start: sigma 0.6984 sdgp 0.1204 rho 1.413 ll -213.606 conv 0 msg relative convergence ( | basis top frequency 3.25
KS k  30 default start: sigma 0.6984 sdgp 0.1204 rho 1.413 ll -213.606 conv 0 msg relative convergence ( | basis top frequency 6.51
KS k  60 default start: sigma 0.0009714 sdgp 6.119 rho 0.001149 ll 889.811 conv 0 msg relative convergence ( | basis top frequency 13.01
KS k 100 default start: sigma 0.0009691 sdgp 6.853 rho 0.4375 ll 885.250 conv 0 msg relative convergence ( | basis top frequency 21.68
KS k  15 from k=60's scales: ll -214.094 conv 0 msg relative convergence (, sdgp 2.17 rho 0.0006846 | default ll -213.606
KS k  30 from k=60's scales: ll -214.211 conv 0 msg relative convergence (, sdgp 1.704 rho 0.0006066 | default ll -213.606
GPCOST fit 1.070 s, message 0.0069 s (0.65 percent of the fit), control 1.09
```

### test files

- frmtmb (core-gated): 288 files, pass 23142, fail 0, error 0, skip 1, warn 0. lib: C:/Users/adf44/source/r/wt-gpgap-lib
- frmtmb.coupling (ext-gated): 11 files, pass 551, fail 0, error 0, skip 5, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.eam (ext-gated): 32 files, pass 1844, fail 0, error 0, skip 3, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.latent (ext-gated): 12 files, pass 681, fail 0, error 0, skip 2, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.learn (ext-gated): 15 files, pass 501, fail 0, error 0, skip 2, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.ode (ext-gated): 11 files, pass 549, fail 0, error 0, skip 1, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.sample (ext-gated): 78 files, pass 4169, fail 0, error 0, skip 1, warn 0. lib: C:/Users/adf44/source/r/wt-gpgap-lib; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.spatial (ext-gated): 8 files, pass 561, fail 0, error 0, skip 0, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 
- frmtmb.survival (ext-gated): 15 files, pass 552, fail 0, error 0, skip 0, warn 0. lib: C:/Users/adf44/source/r/rellib-r16; core lib: C:/Users/adf44/source/r/wt-gpgap-lib 

### OpenBLAS, LAPACK routed (dev/ciharden-openblas.sh); each log names its BLAS

- 0.3.26: `RESULT test-brms-lik-gp-exact.R pass=43 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.26-lapack | t`
- 0.3.26: `RESULT test-gp-by.R pass=86 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.26-lapack | t`
- 0.3.26: `RESULT test-gp-multidim.R pass=67 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.26-lapack | t`
- 0.3.26: `RESULT test-gp-nugget-floor.R pass=69 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.26-lapack | t`
- 0.3.26: `RESULT test-reparam.R pass=423 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.26-lapack | t`
- 0.3.32: `RESULT test-brms-lik-gp-exact.R pass=43 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.32-lapack | t`
- 0.3.32: `RESULT test-gp-by.R pass=86 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.32-lapack | t`
- 0.3.32: `RESULT test-gp-multidim.R pass=67 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.32-lapack | t`
- 0.3.32: `RESULT test-gp-nugget-floor.R pass=69 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.32-lapack | t`
- 0.3.32: `RESULT test-reparam.R pass=423 fail=0 error=0 skip=0 warn=0`; `blas: R home C:/Users/adf44/source/r/frmtmb-wt-gpgap/dev/ciharden-out/Rob0.3.32-lapack | t`

### one session (ungated, nine packages, three orders)

```
SESSION frmtmb forward files=278/288 pass=20518 fail=0 err=0 skip=211 warn=0
SESSION frmtmb reverse files=278/288 pass=20668 fail=0 err=0 skip=211 warn=0
SESSION frmtmb preload files=278/288 pass=20761 fail=0 err=0 skip=213 warn=0
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
SESSION frmtmb.sample forward files=73/78 pass=3817 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample reverse files=73/78 pass=3817 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample preload files=73/78 pass=3817 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
```

LEAK lines 46, untagged 0.

`R CMD check --as-cran` frmtmb: `Status: 1 NOTE; * checking HTML version of manual ... [24s] NOTE;`

`R CMD check --as-cran` frmtmb.spatial: `Status: 1 NOTE; * checking HTML version of manual ... NOTE;`

`R CMD check --as-cran` frmtmb.sample: `Status: 1 NOTE; * checking HTML version of manual ... NOTE;`

Seen to fail on the base build: `test-gp-nugget-floor.R` first form on rellib-r16 `pass=4 fail=2 error=3`; `test-brms-lik-gp-exact.R` on the old translator `pass=23 fail=10 error=1`, gr = FALSE passing as it must; `test-reparam.R` on rellib-r16's frmtmb.sample `pass=376 fail=4 error=2`. Mutants of the new F1 and F2 rows: above.

The core skip is test-curve-scale.R's scale gate; drmTMB and fuzz gates were on. Extension skips are `test-scale.R`'s gate in six packages. Headline on a scratch copy with the 4 adjudication rows out: 83.6, 85.9, 80.5, 91.3.
