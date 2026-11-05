# The release fixes' NEWS text, and sentences split or corrected, for the
# 0.79.0 NEWS sections. dev/rel079-news-merge.R applies each pair (lane text,
# new text) to the top section of the file named, before it rewraps the
# section; a run of whitespace in the lane text matches any whitespace, so a
# pair does not depend on where the lane wrapped its lines. Each pair must
# match once. The numbers are the reviews' (dev/round-20261104.md).
list(
  "NEWS.md" = list(
    # impsamp RF1
    c(r"---[from 4 to 1, and the time is about 40 percent less where the old fit took 3 to 5 rounds.]---",
      r"---[from 4 to 1, and the total time is about a third less where the old fit took 3 to 5 rounds.]---"),
    # predread F1
    c(r"---[On a `cumulative()` row at `sdx = 1.2`, category 1 read 0 to 0.419 and now reads 0.049 to 0.462, where brms's quantiles are 0.048 and 0.469.]---",
      r"---[On a new `me()` row they are the draw's own quantiles. On a `cumulative()` row at `sdx = 1.2`, category 1 read 0 to 0.419 and now reads 0.048 to 0.464, and category 3 reads 0.121 to 0.251, where brms's quantiles are 0.048 to 0.469 and 0.135 to 0.244.]---"),
    # predread F3 and predread x ncgp
    c(r"---[A slope, contrast or test of a coefficient is right, -0.491 with a standard error of 0.0485.]---",
      r"---[A test of a coefficient is right, and so is a slope or contrast on an identity link, -0.491 with a standard error of 0.0485. On a log link a slope in another term reads 0.067 against 0.312.]---"),
    c(r"---[0.46 on every row of a `gp()`, against 0.089 to 0.880.]---",
      r"---[0.31 on every row of a `gp()`, against 0.087 to 0.695.]---"),
    c(r"---[Under `re_formula = NA` the two agree, and no warning comes.]---",
      r"---[Under `re_formula = NA` the two agree, and no warning comes. Over 168 fits on four families with `(1 | g)` or `(1 + x | g)`, `avg_predictions()` was too wide by more than 25 percent on 77 percent of them, median 2.03 times, and `avg_slopes()` on 52 percent, median 1.30 times. A gaussian with no group-level slope reads 0.86 to 1.15.]---"),
    # refuse RF1, RF2
    c(r"---[when that approximation is off by more than 1 in `logLik()`, by `integrate()` over each read value at the estimate. The warning names the remedy: `frm_sample()`, or deleting the rows if the reader carries little about the value.]---",
      r"---[when that approximation moves an estimate by 0.5 standard errors or more. Three inner Hessians over every read value give Tierney-Kadane's term, a gate that costs 0.1 percent of a fit, and at a term of at most 1 nothing else runs. Past the gate, `integrate()` over 10 read values at the estimate and at one standard error above it, per parameter of the read response, gives the move. Over 144 fits against an exact maximum it warned on 75 of 76 shifts of 1 standard error or more and on 3 of 39 below 0.5, and the largest silent shift was 1.01. The warning names the parameter, the size and the remedy: `frm_sample()`, or deleting the rows if the reader carries little about the value. `frmtmb_control(check_mi = )` sets what the check does.]---"),
    c(r"---[The check costs about a third of the time of `frm()` on such a fit; past 30 read values it scales the total from an even sample of them.]---",
      r"---[Past the gate the check costs about 37 percent of the time of `frm()` on such a fit.]---"),
    c(r"---[`von_mises` loses about 1 per value at a mean near the window's end, where the value has two modes: -1.22 at mean 3.14 and kappa 20, -0.98 at mean 3.05 and kappa 2, -0.068 at mean 0, with a reader of slope 0.2.]---",
      r"---[`von_mises` loses -1 to -4.2 per value at a mean near the window's end, where the value has two modes, and estimates are off by up to 17 standard errors there. So `frm()` warns by name, with no integral, when 10 percent or more of read rows have `pi - abs(mu) < 2 / sqrt(kappa)`. At mean 0 it loses -0.068 per value.]---"),
    # latres RF1
    c(r"---[where the Laplace approximation of a latent residual is biased: binary or ordinal rows, at most 3 trials, a mean share of the latent below 0.18, or a family with an estimated dispersion that trades off against the latent.]---",
      r"---[where the Laplace approximation of a latent residual is biased: ordinal rows, or a first-order Laplace correction that moves log `sderr` or `ar[1]` by more than 0.15, or cannot be read. Over 84 designs and 3130 fits it warned on 0.864 of the wrong fits and on 0.157 of the right ones, 0.14 to 0.17 on each set of designs. Above 3000 rows it reads the rows' share of the latent, as the correction's dense matrix is the cost that limit avoids.]---"),
    # ncgp R1
    c(r"---[from the best cell of each of 4 length scales and the best 3 cells overall,]---",
      r"---[from the best cell of each of 6 length scales and the best 3 cells overall,]---"),
    c(r"---[Gaussian, on 384 rows of 16 fresh curves: 5 short by more than 1e-3 nats, worst 9.4, conv 0 and no warning, where 0.78.0 read 57 short of 192, worst 613.]---",
      r"---[Gaussian, on 744 rows of fresh curves: 1 short by more than 1e-3 nats, worst 0.56, conv 0 and no warning, where 0.78.0 read 40 short and 62 short on the two tables of 192 it was read on.]---"),
    # ncgp R2: the off-gaussian gaps are named under Limits, by the decision that
    # this ships with them
    c(r"---[* What it does not reach. Low-count Poisson read 3 of 6 short, worst 1.8, and bernoulli 2 of 6, worst 1.4, conv 0 and no warning. 0.78.0 reads 3 and 2 of 6 there, worst 2.6 and 1.4, not always on the same rows: the field's own information is low. A nonlinear body takes the scan, and stops 98.5 nats short on a row 0.78.0 stops 98.5 short on as well. No start is certain: a marginal can have more than one maximum.]---",
      r"---[* No start is certain: a marginal can have more than one maximum.]---"),
    c(r"---[* `sigma` near `1e-6 * sdgp`, where a nugget of `1e-12` reaches the noise, is the new limit,]---",
      r"---[* Off the gaussian the start reads a working response, and it is not certain: on 150 fresh rows 13 stopped short by more than 1e-3 nats with no warning, against 51 on 0.78.0 (worst 420). Most are under 2 nats at low information, on bernoulli, Poisson and binomial, and some of those 0.78.0 reached. On a step-shaped field the gaps are larger: beside `ar()` 26 nats where 0.78.0 reached the maximum, with `student()` 11 and 6 where 0.78.0 was 1.6 and 2.8 short. A nonlinear body is 98.5 short on both builds on a trend plus a wiggle. A `start` for `theta` is how to check a fit that matters.

* `sigma` near `1e-6 * sdgp`, where a nugget of `1e-12` reaches the noise, is the new limit,]---"),
    # ncgp R4
    c(r"---[on none above `3e-6 * sdgp`.]---",
      r"---[on none above `3e-6 * sdgp`; on 3 of 79 and 10 of 67 rows it fired with `sigma` within 0.8 to 1.25 of the truth, where the noise is at the nugget.]---")
  ),
  "extensions/frmtmb.sample/NEWS.md" = list(
    # the range of the two measures, the lane's and the reviewer's
    c(r"---[takes 3.5 times 0.78.0's wall clock at 300 positions, 570 s against 165 s, from the dense gradient. Memory holds, 1.91 GB against 1.92 GB.]---",
      r"---[takes 3.2 to 3.5 times 0.78.0's wall clock at 300 positions, 531 to 570 s against 165 s, from the dense gradient. Memory holds, 1.91 to 1.95 GB against 1.92 GB.]---")
  ),
  "extensions/frmtmb.learn/NEWS.md" = list(
    c(r"---[Needs the frmtmb that carries the fix, for these numbers.]---",
      r"---[Needs frmtmb 0.79.0 for these numbers.]---")
  )
)
