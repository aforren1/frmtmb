# EAM follow-ups lane

Worktree `C:/Users/adf44/source/r/frmtmb-wt-eam-followups`, branch
`wt-eam-followups`, base 68d6782 (frmtmb 0.52.0, frmtmb.eam 0.3.0).
Written as the work happens. Every number below came from a script under
the `ea-` prefix in the shared scratchpad and can be re-run.

## Environment

Private library
`.../scratchpad/ea-lib`, built fresh:

* `frmtmb` 0.52.0 and `frmtmb.eam` 0.3.0, both `R CMD INSTALL`ed from
  this worktree with `--library=`.
* `Rmpfr` 1.1-2, `statmod` 1.5.2, `WienR` 0.3-17, `rtdists` 0.11-6,
  `RWiener` 1.3-3, `numDeriv` 2016.8-1.1, `testthat` 3.3.2, CRAN Windows
  binaries.
* The user library is on `R_LIBS` as a read fallback only.

R 4.6.1. EMC2 is NOT installed and is not on CRAN's Windows binary index
for 4.6 at the time of writing; see item 1 for what that costs.

## The three items

1. Across-trial variability (`sv`, `sz`, `st0`) for `wiener_gng()`.
2. `lccdf` for `rdm()` and `wiener_gng()`, so `cens()` and `trunc()`
   work through core 0.52.0's slot.
3. Adopt core's any-of `required_aterms` in `wiener()` and `gddm()`.

## Progress log

* Library built, core and eam installed from the worktree.

## Reference machinery (item 1)

`scratchpad/ea-work/ea-ref.R`. Two independent 200-bit routes to the
no-go probability, plus Gauss-Legendre nodes refined by Newton in
`mpfr` so the reference is not capped at 1e-16 by its own node
positions.

Validated against the six constants `test-rdm-gng.R` already pins:

| t | v | a | w | 200-bit no-go | abs(A-B)/A | shipped blend rel |
|---|---|---|---|---|---|---|
| 1.25 | 1.0 | 1.4 | 0.45 | 0.26333325052669940 | 2.4e-60 | 2.0e-15 |
| 20 | 2.0 | 1.0 | 0.50 | 0.11920292202211756 | 1.3e-59 | 2.7e-16 |
| 2.5 | 5.0 | 4.0 | 0.75 | 1.1883283058299342e-13 | 2.7e-48 | 1.9e-15 |
| 0.05 | 1.0 | 1.4 | 0.45 | 0.99878687191466810 | 6.2e-61 | 3.7e-15 |
| 3.0 | 1.5 | 2.0 | 0.90 | 0.0021602702829594682 | 3.4e-58 | 1.7e-15 |
| 2.5 | 5.0 | 4.0 | 0.90 | 8.2014416716939389e-16 | 1.7e-45 | 8.8e-13 |

Every digit the review reports is reproduced. `A` is the large-time
eigenfunction route, `B` the small-time image sum; they share no
algebra.

## Item 1: the cost wall, and what it forces

The first design gave the no-go branch the SAME node sets the density
uses, which is a three-dimensional product where the density's is
two-dimensional. Measured, on 500 rows with a constant deadline:

| variability | wall | outcome |
|---|---|---|
| none | - | logLik -237.336779089870, bit-identical to 0.3.0 |
| sv | 3.1 s | fits |
| sz | 2.6 s | fits |
| st | 9.0 s | fits |
| sv+sz | 32.2 s | fits |
| sv+sz+st | 460 s | `std::bad_alloc` |

`sv+sz+st` is 11 x 7 x 21 = 1617 evaluations of a two-series
distribution function per row, against the density's 7 x 21 = 147
evaluations of a two-series density. The tape does not fit in memory.
So the no-go branch needs its own node counts, measured on its own
integrand, and the sweep below is what sets them.

## Item 2: censoring, measured

`rdm()` declares BOTH `lccdf` and `lcdf`; `wiener_gng()` declares
`lccdf` only.

The race's survivor is the product over accumulators of the same
start-point integral `rdm_law$lsurv()` already forms for every loser of
every observed trial, so the seam is three lines
(`R/rdm.R`, `rdm_lccdf()`). The distribution function is
`-expm1(log S)`, floored, so it does not return exactly zero where the
survivor is within a rounding of one.

Every identity is against the likelihood written out by hand at the
fitted parameters, through the finalized family's own links:

| identity | rows | frmtmb | by hand | rel |
|---|---|---|---|---|
| `rdm(3)` `cens()`, right only | 400, 100 censored | -7.15755969740973 | -7.15755969740975 | 2.2e-15 |
| `rdm(2)` `cens()`, all four codes | 300 | -170.821242566795 | -170.821242566795 | 5.0e-16 |
| `rdm(2)` `trunc(lb = 0.30)` | 1146 kept | 606.388731053518 | 606.388731053518 | 9.4e-16 |
| `wiener_gng()` no-go branch vs `cens()` | 400 | -109.40944838876700 | -109.40944838876700 | **0, exactly** |

Two further facts, both measured:

* The `vint()` winner on a CENSORED race row is genuinely ignored:
  moving every censored row's winner from accumulator 1 to accumulator
  3 changes the log likelihood by **exactly zero**. It is still
  required, because a declaration cannot be conditional on a censoring
  code.
* `wiener_gng()` refuses left censoring, interval censoring and
  `trunc()` by name, through core's own message, because the family
  declares no `lcdf`. That is deliberate: this likelihood is a
  defective density plus a point mass, and a window normalizer on the
  response scale would renormalize the density while saying nothing
  about the mass.

The go/no-go identity is the interesting one. A right-censored go/no-go
trial and a no-go trial are the same statement about the same process,
so the two spellings do not merely agree to a tolerance, they are the
same arithmetic and differ by no bits at all.

## Item 1: what each of the three integrals costs, measured

Worst relative error over six parameter settings spanning deadlines
0.8 to 2.5, drifts 0.5 to 3, boundary separations 0.8 to 2.5 and start
points 0.4 to 0.8, against the 200-bit reference. One variability
parameter live at a time.

### sz, Gauss-Legendre over the uniform start point

| nodes | sz=0.1 | sz=0.3 | sz=0.5 |
|---|---|---|---|
| 3 | 4.5e-06 | 2.1e-03 | n/a |
| 5 | 1.8e-11 | 6.4e-07 | n/a |
| 7 | 5.1e-15 | 4.2e-11 | n/a |
| 9 | 4.0e-15 | 2.1e-15 | n/a |

The `sz = 0.5` column is not a failure of the quadrature. One of the six
settings has `bias = 0.8`, where a width of 0.5 puts the top of the
uniform range at 1.05, outside the boundaries, and the reference is
undefined there. See the clamp below.

### st, Gauss-Legendre over the uniform non-decision time

| nodes | st=0.05 | st=0.15 | st=0.3 |
|---|---|---|---|
| 3 | 3.4e-11 | 2.5e-08 | 1.6e-06 |
| 5 | 8.6e-15 | 3.5e-13 | 4.0e-10 |
| 7 | 6.3e-15 | 1.2e-14 | 6.1e-14 |
| 9 | 4.2e-15 | 4.0e-15 | 2.7e-15 |
| 21 | 2.0e-15 | 5.3e-15 | 8.5e-16 |

Seven nodes is the knee, where the DENSITY's own `st` integral needs 21.
The difference is real and is the point: the density's range is cut at
the response time and its integrand turns on sharply at the cut; the
probability has no cut, because a non-decision time past the deadline
leaves no window and the answer there is exactly one.

### sv, Gauss-Hermite over the normal drift

| nodes | sv=0.3 | sv=0.8 | sv=1.5 |
|---|---|---|---|
| 5 | 8.3e-05 | 3.0e-02 | 1.6e-01 |
| 7 | 9.0e-07 | 8.9e-04 | 1.8e-02 |
| 9 | 7.6e-09 | 7.1e-04 | 2.0e-02 |
| 11 | 9.5e-12 | 3.8e-06 | 1.0e-02 |
| 15 | 8.4e-15 | 1.3e-06 | 1.4e-03 |
| 21 | 9.5e-15 | 3.3e-09 | 3.5e-05 |
| 31 | 1.2e-14 | 1.2e-12 | 6.6e-07 |
| 41 | 1.4e-14 | 2.8e-14 | 9.1e-08 |

This is the one integral of the three that does not saturate, and it is
the honest cost of the family. The drift enters the distribution
function through the eigenvalues, so there is no closed form, and the
integrand is a transition in the drift whose width shrinks as the
boundary separation grows: the node count a given accuracy needs rises
with the product of the separation and `sv`.

### A splitting idea that does NOT work, measured before it was believed

The obvious remedy is to split the large-time route into the
gambler's-ruin probability, which is two exponentials and therefore
cheap enough for a fine rule, and the correction series, which is
thirty terms and therefore wants a coarse one. Measured, at the same
four settings:

| nodes | whole | `P_lower` alone | correction alone |
|---|---|---|---|
| 11 | 3.8e-06 | 1.1e-03 | 4.7e-04 |
| 21 | 3.3e-09 | 1.5e-05 | 6.2e-06 |
| 31 | 1.2e-12 | 4.8e-07 | 2.1e-07 |
| 41 | 2.9e-14 | 3.2e-08 | 1.4e-08 |

(`sv = 0.8`; the `sv = 1.5` table has the same shape.) Each half
converges two to six orders WORSE than the sum of the two. The
quadrature errors of the two halves cancel, so splitting them destroys
the cancellation. The idea is not merely no help, it is a regression,
and it is recorded here so that the next person does not spend the
afternoon on it.

### The out-of-boundary start point: a NaN, not a barrier

`wiener()` documents that a wide `sz` at a biased start can push the
uniform range past a boundary, "where the density is near zero and the
likelihood is a barrier rather than a cliff". That is true of the
density. It is NOT true of the no-go probability. At `t = 1.25`, drift
1, separation 1.4:

| w | no-go log prob | density log |
|---|---|---|
| 0.95 | -4.441 | -5.988 |
| 1.00 | -41.04 | -40.81 |
| 1.05 | **NaN** | -6.057 |
| -0.05 | **+0.1412** | -5.089 |

Above one the gambler's-ruin term takes `log1p(-w)` of a negative and
returns `NaN`; below zero it returns a positive log probability, which
is a probability above one. A `NaN` is not a barrier, it is the end of
the tape, and one node of one row takes the whole fit with it. Measured
end to end: `bias = 0.80` with `sz = 0.50` returns `NaN` from the
variability integral, `sz = 0.30` returns -2.78694.

`ddm_wclamp()` holds the start point in `(1e-9, 1 - 1e-9)` inside the
no-go integral only. The clamped value is the correct LIMIT rather than
a fudge: a start point at the upper boundary is a trial that has
already finished at the go boundary, so its no-go probability is zero,
and one at the lower boundary is a trial that finished at the other, so
its no-go probability is one. The go branch is left alone, because it
is `wiener()`'s density shared bit for bit and clamping it would make
one model score differently under two family names.

### The image sum was three times longer than it needs to be

`ddm_cdf_ks` was 12, meaning 25 image terms and 50 `pnorm` calls in
every evaluation of the no-go probability. The blend gives the
small-time route a non-zero weight only below `u = 0.197`, where `tanh`
has not yet saturated, and at that `u` the `j = 2` term already carries
`exp(-2 * 2 * 2.5 / 0.197)`, which is 9e-23.

Measured over 840 rows spanning `t` in {0.05 .. 15}, `v` in {-2 .. 5},
`a` in {0.8 .. 4} and `w` in {0.25 .. 0.9}:

| ddm_cdf_ks | rows identical to ks=12 | max abs log difference | max gradient difference |
|---|---|---|---|
| 8 | 840/840 | 0 | - |
| 6 | 840/840 | 0 | - |
| 4 | 840/840 | 0 | **0** |
| 3 | 840/840 | 0 | - |
| 2 | 840/840 | 0 | - |

Values AND gradients are bit-identical, the gradients checked by
`numDeriv` in all four arguments on every row. `ddm_cdf_ks` is now 4,
which keeps a margin of three terms over `|j| = 1`, the largest index
that contributes anything. `ddm_cdf_kl` is NOT reducible the same way:
its terms decay like `exp(-k^2 pi^2 u / 2)` and the blend hands over to
it from `u = 0.01`, where `k = 30` is only `exp(-44)`.

It is a pure cost change and it is worth 2.3x. The same 500-row
three-variability fit, before and after, on the 175-point grid: 285.3 s
against 123.8 s, and the log likelihood is -226.737921 either way.

### The three dimensions do not compound

Joint error against the same rule at 61 x 21 x 21, four settings:

| grid | pts | sv=0.3 joint / sv alone | sv=0.8 joint / sv alone | sv=1.5 joint / sv alone |
|---|---|---|---|---|
| 7, 5, 5 | 175 | 8.2e-07 / 9.0e-07 | 1.7e-03 / 8.9e-04 | 2.6e-02 / 1.8e-02 |
| 11, 5, 5 | 275 | 3.7e-10 / 9.5e-12 | 2.8e-05 / 3.8e-06 | 1.0e-02 / 1.0e-02 |
| 11, 7, 7 | 539 | 3.7e-12 / 9.5e-12 | 2.8e-05 / 3.8e-06 | 1.0e-02 / 1.0e-02 |
| 15, 7, 7 | 735 | 8.8e-15 / 1.4e-14 | 1.5e-06 / 1.3e-06 | 8.7e-04 / 1.4e-03 |
| 21, 7, 7 | 1029 | 9.6e-15 / 1.5e-14 | 5.3e-09 / 3.3e-09 | 8.2e-05 / 3.5e-05 |
| 31, 9, 9 | 2511 | 1.4e-14 / 6.4e-15 | 5.5e-13 / 1.2e-12 | 2.1e-06 / 6.6e-07 |

The joint error tracks the `sv` error alone to within a factor of
three everywhere, so there is no interaction penalty. At `sz` and `st`
of 5 the uniforms become the binding term (row two, first column: 3.7e-10
against 9.5e-12); at 7 they never do. Hence
`nogo_nodes = c(sv = 15L, sz = 7L, st = 7L)`: seven so the uniforms are
not the constraint, and `sv` as the only knob worth turning.

### The cost, on 500 rows with all three live

| nogo_nodes | grid | wall |
|---|---|---|
| 7, 5, 5 | 175 | 123.8 s |
| 11, 5, 5 | 275 | 184.5 s |
| 11, 7, 7 | 539 | 609.9 s |

## Item 1: the identities, measured

### The go branch IS wiener()'s density

Required to hold to 1e-12. It holds to ZERO. On 223 all-go rows, the
two families' `lpdf` slots evaluated at the same distributional
parameters:

| variability | max abs difference | rows bit-identical |
|---|---|---|
| sv | 0 | 223/223 |
| sz | 0 | 223/223 |
| st | 0 | 223/223 |
| sv+sz | 0 | 223/223 |
| sv+sz+st | 0 | 223/223 |

Not a tolerance but an identity: `gng_lpdf_var()` calls
`ddm_lpdf_var()` with `up = 1`, the same node sets built from the same
`nodes` argument, and the same `delta`. Making the two `delta`s agree
is why `wiener_gng()` grew a `gng_family(cfg, ub, delta)` builder: at
0.3.0 it floored the decision time at a flat 1e-12 while `wiener()`
used `1e-9 * min(y)`, and the identity cannot hold across two margins.

### The two branches still add to one

Go density integrated from its own support edge to the deadline, plus
the no-go probability, at drift 1.0, separation 1.4, ndt 0.25, bias
0.5, deadline 1.5:

| variability | go mass | no-go | 1 - total |
|---|---|---|---|
| none | 0.7775817880134 | 0.2224182119866 | -4.4e-16 |
| sv | 0.7586018631478 | 0.2413981368522 | 0 |
| sz | 0.7744621349195 | 0.2255378650805 | 0 |
| st | 0.7774883284769 | 0.2225116715235 | -4.5e-13 |
| sv+sz | 0.7559421783347 | 0.2440578216653 | 0 |
| sv+sz+st | 0.7558594019322 | 0.2441405980628 | 4.9e-12 |

**CORRECTED in the punch round: these are a BEST case, not a bound.**
Every row above is inside the boundaries and at a moderate width. The
reviewer's wider grid finds -9.9e-08 at `st` = 0.45 and -8.1e-05 at
`sv` = 2.0 with the same default nodes, and percents once the `sz`
range leaves the boundaries. See the punch round below.

Under `st` the go density's SUPPORT starts at `ndt - st / 2`, not at
`ndt`. An integral started at `ndt` misses that mass and reports a
defect of 8.1e-04 that is the integration limit's, not the model's.

### The AD gradient of the new integral

`RTMB::MakeTape()$jacobian` against `numDeriv::grad` in all seven
arguments, four settings spanning drift -0.5 to 3:

| v, a, w | max rel |
|---|---|
| 1.0, 1.4, 0.50 | 4.0e-10 |
| 2.5, 1.0, 0.60 | 8.0e-11 |
| -0.5, 2.0, 0.35 | 7.8e-10 |
| 3.0, 2.5, 0.80 | 2.2e-10 |

That is `numDeriv`'s own accuracy, so the tape records what the
function computes.

### Simulation from the process, 40000 trials each

`wiener_gng_simulate()` draws per-trial parameters, runs the diffusion
through RWiener and records a response only if the go boundary was
reached before the deadline. It evaluates neither the density nor the
probability.

| variability | observed no-go | predicted | z | chi-square (9 df) |
|---|---|---|---|---|
| none | 0.22503 | 0.22242 | 1.25 | 9.25 |
| sv | 0.24140 | 0.24140 | 0.00 | 7.44 |
| sz | 0.22567 | 0.22554 | 0.07 | 15.95 |
| st | 0.22323 | 0.22251 | 0.34 | 14.13 |
| sv+sz+st | 0.24063 | 0.24414 | -1.64 | 10.35 |

The chi-square is on ten equal-count cells of the go response times
against the density's own integral over each. The 0.999 critical value
at 9 df is 27.88; the worst statistic is 15.95.

### Against EMC2, whose DDMGNG does carry all three

It does. `EMC2::DDMGNG()$p_types` is `v a sv t0 st0 s Z SZ`, identical
to its own `DDM()`. So the comparison is available and was made.

**One convention differs, and it is worth stating because it is not a
defect on either side.** EMC2's `st0`, through `WienR`, is a uniform on
`[t0, t0 + st0]`; this package's `st`, following brms, is centred on
`ndt`. With that shift applied the two agree; without it they disagree
by 14 to 45 percent, in the DENSITY as much as in the probability,
which is what identifies it as a convention rather than a bug in the
new code.

No-go probability, `1 - EMC2:::pDDM(precision = 1e-12)`, 108 grid
points spanning deadlines 0.8 to 2.5, separations 0.8 to 2.5, drifts -1
to 3 and start points 0.4 to 0.6:

| sv, sz, st | max rel |
|---|---|
| 0, 0, 0 | 9.4e-14 |
| 0.5, 0, 0 | 2.4e-13 |
| 0, 0.2, 0 | 3.9e-13 |
| 0, 0, 0.1 (shifted) | 1.5e-13 |

Go density against `EMC2:::dDDM`, 200 response times: 1.2e-15 at every
combination without `st0`, and 1.6e-15 with `st0` once shifted.

### sv is weakly identified in a go/no-go design, and that is the design

Recovery on 3000 trials, drift 1.2, separation 1.4, true `sv` 0.6:

| family | what it sees | fitted sv |
|---|---|---|
| `wiener_gng()` | one boundary | **0.001** |
| `wiener()` | both boundaries | 0.612 |

This is not an optimizer failure and not a defect in the quadrature.
Profiled at the true values of the other parameters the likelihood does
peak at the truth (-1028.476 at `sv` = 0.6 against -1040.056 at
`sv` = 0.001), but the fit reaches -1027.928 at `sv` = 0.0005 with the
drift at 1.09 and the separation at 1.33. It found a genuinely higher
point: with only one boundary observed, drift variability trades
against the drift and the separation along a ridge. `wiener()` on the
same generative parameters, seeing both boundaries, recovers 0.612.

The family documents this rather than pretending otherwise.

## Item 3: the adoption, line by line

`wiener()`: `required_aterms = list(c("dec", "vint1"))` at
`R/wiener-family.R:382`, and the `is.null(up)` branch of
`ddm_check_response()` deleted. Net **-5** lines in the file (16 added,
21 removed), of which the deleted branch is 18 lines: 13 of refusal, 3
of comment, 2 of brace. The roxygen paragraph calling this "the one
hand-rolled check left in this package" is retired with it.

`gddm()`: the boundary group prepended to `req` at `R/gddm.R:1035`, so
`required_aterms` at `R/gddm.R:1040` is now a list, and the
`is.null(up)` branch of `gd_check_response()` deleted. Net **-14**
lines (9 added, 23 removed); the deleted branch is 18 lines, of which
12 are the refusal and 5 the comment explaining why `required_aterms`
could not express it. The review predicted "about twelve lines
including the multi-line refusal" and put the branch at `R/gddm.R:1085`;
both are right.

`lba()` and `rdm()` gain nothing and were not touched: both already read
`required_aterms = "vint1"`.

**What still cannot be declared, and why.** `gddm()`'s CONDITION index
stays hand-rolled. The boundary arrives as `dec` or `vint1`, and which
slot carries the condition MOVES with that choice: alongside `dec()` it
is `vint1`, inside `vint(upper, cond)` it is `vint2`. The requirement is
therefore a disjunction of conjunctions, and core's seam expresses a
conjunction of disjunctions. The comment at the declaration says so.

### The refusal texts that changed

| before | after |
|---|---|
| `wiener: the decision indicator is missing. ... frm(bf(rt \| dec(decision) ~ x), family = wiener(), ...) ... vint(upper) carries the same thing ...` | `wiener: the density needs one of `dec` or `vint1`, which nothing on this response supplies. Write the addition term: rt \| dec(<column>) ~ ...` |
| `gddm: the decision indicator is missing. ...` | `gddm: the density needs one of `dec` or `vint1`, ... rt \| dec(<column>) ~ ...` |

The cost is real and is the review's: the new sentence names the term
VALUES rather than the spellings, so it no longer writes
`dec(decision)` or mentions `vint` at all beyond the value `vint1`. It
is a fair trade only because the alternative is named in the
`one of ... or ...` clause.

### The pins that changed, and why each one had to

| file:line | was | now |
|---|---|---|
| `test-family.R:109-117` | `"decision indicator is missing"`, `"dec\(decision\)"`, plus a comment saying frmtmb has no such seam | `"the density needs one of ..."`, `"rt \| dec(<column>) ~"` fixed, plus `expect_identical()` on the declaration itself |
| `test-defects.R:154-161` | `"decision indicator is missing"`, `"vint\(upper\)"` | `"the density needs one of ..."`, `"vint1"` |
| `test-gddm-family.R:32-33` | `expect_identical(gddm()[["required_aterms"]], character(0))`, `expect_setequal(..., "vreal1")` | `expect_identical(..., list(c("dec", "vint1")))`, `expect_identical(..., list(c("dec", "vint1"), "vreal1"))` |
| `test-gddm-family.R:156-168` | `"decision indicator is missing"`, `"vint(upper, cond)"` fixed | `"the density needs one of ..."`, `"rt \| dec(<column>) ~"` fixed |

**CORRECTED in the punch round.** "Four pins changed" undercounts. It
is EIGHT assertions changed plus ONE added, over three files, and two
of the eight are not text: `test-gddm-family.R` changes
`expect_identical(gddm()[["required_aterms"]], character(0))` to
`list(c("dec","vint1"))` and `expect_setequal(..., "vreal1")` to
`expect_identical(..., list(c("dec","vint1"), "vreal1"))`, which pin a
DATA STRUCTURE rather than a message. The breakdown: `test-defects.R`
2 regex; `test-family.R` 2 regex plus a new
`expect_identical(wiener()[["required_aterms"]], ...)`;
`test-gddm-family.R` 2 regex plus the 2 structural ones.

The review named the first two files' four expectations and predicted
`test-gddm-family.R:32-33` as "a separate pin that adoption of `gddm()`
would break, which the lane does not mention at all". Both were
correct. `test-gddm-family.R:158` the review called "unaffected by a
wiener-only adoption", which is true, and it does break once `gddm()`
adopts too, which this lane did.

Two comment blocks that had become false went with the expectations:
`test-family.R:109-111` and `test-defects.R:154-156`, both of which
said frmtmb had no way for a family to declare a required addition term.

**NEWS.** `extensions/frmtmb.eam/NEWS.md:429` under `# frmtmb.eam
0.2.0` says "One hand-rolled check remains, and is not
`required_aterms`'s fault". It is now false as a statement about the
package. It is NOT edited: a released section records what was true at
that release, and rewriting it would misrepresent 0.2.0. The
development entry retires it by name instead, so a reader searching for
"hand-rolled" finds both and the order tells the story.

## Cross-lane

`R/zzz.R` was edited, but only the four compat ROWS in
`ddm_compat_rules()` and `rdm_gng_compat_rules()` (the `rdm`/`cens()`,
`rdm`/`trunc()`, `wiener_gng`/`cens()` and `wiener_gng`/`trunc()`
entries). The `.onLoad()` body, the `frmtmb_register_aterm()` call and
the feature declarations are untouched, which is where the protocol
lane's `accepts_aterms` work would land. `R/ddm-shared.R` was not
touched at all.

The `required_aterms =` argument was added on its own line inside the
`frmtmb::custom_family()` calls in `R/wiener-family.R` and `R/gddm.R`.
If the protocol lane adds `accepts_aterms =` to the same calls the two
edits are adjacent but distinct lines.

## Verification

Private library `.../scratchpad/ea-lib`: `frmtmb` 0.52.0 and
`frmtmb.eam` from this worktree, plus `Rmpfr` 1.1-2, `statmod` 1.5.2,
`WienR` 0.3-17, `rtdists` 0.11-6, `RWiener` 1.3-3, `numDeriv`
2016.8-1.1, `testthat` 3.3.2 and `EMC2` 3.5.0 as CRAN Windows
binaries. Explicit `--library=` throughout.

### The seven named files, one process each

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-rdm-gng.R | 221 | 0 | 0 | 0 |
| test-defects.R | 59 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 |
| test-variability.R | 140 | 0 | 0 | 0 |
| test-gddm-family.R | 77 | 0 | 0 | 0 |
| test-family.R | 35 | 0 | 0 | 0 |

`test-rdm-gng.R` was 169 before this lane and is 221: eleven new
`test_that()` blocks, 38 to 49.

### The whole suite, one process, NOT_CRAN=true

**1339 pass, 0 fail, 0 error, 1 skip.** 0.3.0 was 1284 with 1 skip, so
+55 and the same single skip (`test-sampling.R`, unchanged).

### Core's test-compat.R against this eam

**273 pass, 0 fail, 0 error, 0 skip**, run with this lane's
`frmtmb.eam` loaded so its rewritten rows and its aterm registration
are in the vocabulary.

### R CMD check --as-cran

`_R_CHECK_CRAN_INCOMING_=false`, core in `ea-lib`, pandoc from
`RSTUDIO_PANDOC` on PATH. **Status: 1 ERROR, 1 WARNING, 2 NOTEs**, and
all four are the machine rather than the package:

| finding | cause |
|---|---|
| ERROR, PDF manual without index | `pdflatex is not available`; confirmed absent from PATH |
| WARNING, PDF version of manual | the same missing LaTeX |
| NOTE, HTML manual | `Skipping checking math rendering: package 'V8' unavailable` |
| NOTE, non-standard things | `frmtmb.eam-manual.tex`, the leftover of the failed PDF build |

Everything that tests the package passed: `checking tests ... [14m] OK`,
`checking re-building of vignette outputs ... [32s] OK`, and every Rd,
code and dependency check OK.

### roxygen

`roxygenise()` run to a fixed point: a second run leaves every `man/*.Rd`
and `NAMESPACE` byte-identical (md5 compared). Three Rd files change,
`rdm.Rd`, `wiener_gng.Rd` and `wiener_gng_simulate.Rd`. One roxygen
warning remains and is pre-existing: `gddm.R:169` cannot resolve a link
to `gd_tri_df`, an internal.

## What was left out, and why

* **`sv` for `wiener_gng()` is shipped but weakly identified**, and the
  documentation says so rather than the family refusing it. Refusing
  would have been defensible; it is the wrong call because `sv` is
  identified when a go/no-go block sits beside a two-choice condition,
  and because a family that silently lacks a parameter its sibling has
  is worse than one that has it with a caution.
* **No better quadrature for the drift.** Gauss-Hermite reaches 1e-2 at
  `sv = 1.5` with 11 nodes. The obvious remedy, splitting the
  gambler's-ruin probability off from the correction series, was
  measured and is a REGRESSION: each half converges two to six orders
  worse than the sum, because their quadrature errors cancel. Recorded
  above so the next person does not repeat it.
* **The go branch is not clamped** at a start point outside the
  boundaries, only the no-go branch is. Clamping it would change
  `wiener()`'s shipped density, which is a decision for that family and
  for a lane that owns its pins.
* **`ddm_cdf_kl` is unchanged at 30.** Unlike `ddm_cdf_ks` it is not
  reducible: its terms decay like `exp(-k^2 pi^2 u / 2)` and the blend
  hands over from `u = 0.01`, where `k = 30` is only `exp(-44)`.
* **A defect found and NOT fixed.** `wiener_gng()`'s starting value for
  the boundary separation is `log(1.4)`, at `R/wiener-gng.R:489`.
  `init_dpars` takes RESPONSE-scale values and frmtmb applies the link
  itself, as `rdm_inits()` spells out in a comment and as `wiener()`
  and `lba()` both do, so this starts the separation at 0.336 rather
  than at 1.4. **CORRECTED in the punch round: the line is
  `R/wiener-gng.R:489`, not `:566`, and the fix shipped rather than
  being deferred. See the punch round below.**
* **`nogo_nodes` defaults are a cost compromise, not an accuracy one.**
  `c(sv = 15L, sz = 7L, st = 7L)` is 735 grid points; the measured wall
  for 539 points on 500 rows is 610 s. A user who needs 1e-9 at a wide
  `sv` must raise `sv` and pay for it, and `?wiener_gng` carries the
  table to decide with.

## Touched files

| file | + | - |
|---|---|---|
| `extensions/frmtmb.eam/NEWS.md` | 121 | 0 |
| `extensions/frmtmb.eam/R/gddm.R` | 9 | 23 |
| `extensions/frmtmb.eam/R/rdm.R` | 64 | 0 |
| `extensions/frmtmb.eam/R/wiener-cdf.R` | 167 | 4 |
| `extensions/frmtmb.eam/R/wiener-family.R` | 16 | 21 |
| `extensions/frmtmb.eam/R/wiener-gng.R` | 417 | 44 |
| `extensions/frmtmb.eam/R/zzz.R` | 7 | 7 |
| `extensions/frmtmb.eam/man/rdm.Rd` | 32 | 0 |
| `extensions/frmtmb.eam/man/wiener_gng.Rd` | 185 | 18 |
| `extensions/frmtmb.eam/man/wiener_gng_simulate.Rd` | 8 | 1 |
| `extensions/frmtmb.eam/tests/testthat/test-defects.R` | 11 | 5 |
| `extensions/frmtmb.eam/tests/testthat/test-family.R` | 13 | 6 |
| `extensions/frmtmb.eam/tests/testthat/test-gddm-family.R` | 15 | 9 |
| `extensions/frmtmb.eam/tests/testthat/test-rdm-gng.R` | 293 | 2 |
| `dev/eam-followups-findings.md` | new | |

Nothing outside `extensions/frmtmb.eam/` except this findings file.
No commit, no staging, no branch moved. The main checkout was read
only. No `.log` file anywhere in the worktree: every log lives under
the `ea-` prefix in the scratchpad.

# Punch round, 2026-09-06

Against `dev/reviews/2026-09-06-eam-followups.md`, verdict PUNCH: five
items, four documentation and one word of code.

## 1. The vignette said the opposite of what shipped (F2)

`vignettes/ddm.Rmd`, the "What it does not do" subsection, said
`wiener_gng()` "offers none of them" and argued that shipping a subset
would be wrong. It was written before this lane and this lane did not
look at it, which is the defect: the vignette ships in the tarball and
renders into the site, so it was the first document a user would read
and it described a feature that no longer matched the code.

Replaced with two sections, the old heading kept for what really
remains:

* `vignettes/ddm.Rmd:788` "Across-trial variability" - what shipped,
  the go-branch bit-identity, why `sv` needs a quadrature the density
  does not, and both node arguments with `nogo_nodes =` named and its
  defaults explained.
* `vignettes/ddm.Rmd:828` "What a go/no-go design cannot tell you" -
  the `sv` weak-identification caution, which previously existed only
  in the Rd, plus the out-of-boundary `sz` limit and the EMC2 `st0`
  convention.
* `vignettes/ddm.Rmd:863` "What it does not do" now holds only the
  refusing mean, which is what the heading is still true of.

The vignette rebuilds; `R CMD check` reports
`checking re-building of vignette outputs ... OK`.

**Not rebuilt: `docs/`.** The pkgdown site under `docs/frmtmb.eam/` is
a repository-level artifact outside `extensions/frmtmb.eam/`, which is
what the brief assigns to this lane. Regenerating it would rewrite
pages owned by every other extension. `docs/frmtmb.eam/articles/ddm.md`
is therefore stale against this change and needs a site rebuild at
consolidation.

## 2. The mass identity was overstated (F1)

Three places said the pair sums to one without qualification, and the
lane's own clamp is what makes that false outside the boundaries.

* `R/wiener-gng.R:141` (shipped as `man/wiener_gng.Rd`) now says the
  identity holds WHILE THE START-POINT RANGE STAYS INSIDE THE
  BOUNDARIES, gives the in-range quadrature sizes (3.9e-14 at
  `sv` = 0.6, 8.1e-05 at `sv` = 2.0, 2.0e-10 at `st` = 0.20, 9.9e-08 at
  `st` = 0.45) instead of "exactly one", carries the reviewer's table
  of the out-of-range deficit (0.970639 at `sz` = 0.5 / `bias` = 0.85
  and 0.843370 at `sz` = 0.9 / `bias` = 0.90, so up to 16 percent), and
  states the barrier result that makes it tolerable: the interior peak
  beats the best exterior point by 189 log units and the surface is
  monotone across the crossing.
* `R/wiener-cdf.R`, the `ddm_wclamp_dev()` comment, no longer says "the
  correct LIMIT rather than a fudge" flat. It says the clamped value is
  the correct limit FOR THE BRANCH IT IS APPLIED TO, and that it does
  not repair the pair.
* `NEWS.md`, the clamp bullet, the same correction with the same
  numbers.

The reviewer's own revision of F1 from MEDIUM to LOW is carried into
the text rather than dropped: the barrier measurement is quoted,
because it is the reason the defect is a wording problem and not a
numerical one.

## 3. The bs init defect, fixed (F6)

`R/wiener-gng.R:547` is now `bs = function(y, aterms) 1.4`, with a
comment naming the contract (`?frmtmb_family` step 4, restated at this
package's `R/rdm.R:414`) and recording what the wrong start cost.

**It changes no answer, measured.** The two starts on the ten
replicates the recovery test uses, same data, one start each:

| replicate | logLik at 1.4 | logLik at log(1.4) | difference | bs difference |
|---|---|---|---|---|
| 1 | -844.651 | -844.651 | 0 | 0 |
| 2 | -913.391 | -913.391 | 0 | 0 |
| 3 | -861.591 | -861.591 | 0 | 0 |
| 4 | -838.073 | -838.073 | 0 | 0 |
| 5 | -831.327 | -831.327 | 0 | 0 |
| 6 | -799.648 | -799.648 | 0 | 0 |
| 7 | -802.358 | -802.358 | 0 | 0 |
| 8 | -874.485 | -874.485 | 0 | 0 |
| 9 | -812.651 | -812.651 | 0 | 0 |
| 10 | -848.305 | -848.305 | 0 | 0 |

Every one the same optimum to optimizer tolerance (the review measured
logLik differing by 5e-12 and 1e-11 and `bs` by 1e-8 and 6e-8 on replicates 1 and 2; the zeros above are three-decimal rounding), which is why **no pin in `test-rdm-gng.R` needed changing**: the file is 221 pass, 0 fail before and after.

One thing worth recording. The run after the fix emits
`Large maximum absolute gradient at the optimum (0.00184)` at
`test-rdm-gng.R:733`. The fix introduced it: over the same ten
replicates the new start warns once (replicate 4) and the old start
not at all (the review counted both; the first version of this
section said one warning each, which was wrong). One replicate in
ten at a gradient of 0.00184 is not worth tuning a start over.

## 4. The rdm trunc() row now names its range (F8)

`R/zzz.R`, the `rdm`/`trunc()` row, moves from `works` to
`conditional` and carries the reviewer's table. Core forms the
left-truncation normalizer as `1 - F(lb)` on the probability scale, so
the declared `lccdf` is not used for it and the cancellation returns at
the bound: exact at `log S(lb) = -0.406`, 5.0e-11 at -14.0, 1.1e-07 at
-21.5, 8.7e-05 at -28.7 and 20 percent wrong at -35.8. The row points
at core's own statement of the limitation at `R/families.R:278-285`,
which calls the fix a windowed log-difference slot core does not yet
have. Both `trunc()` directions are named, the right-truncated fit at
4.1e-16 included.

## 5. The lane's own report, corrected (F3)

All three corrections are made in place above, each marked
**CORRECTED in the punch round** so the original claim and its
correction sit together: the `bs` line is `:489` (now `:547` after the
comment), the pin change is eight assertions plus one addition over
three files with two structural rather than textual, and the 4.9e-12
mass figure is a best case rather than a bound.

## The optional items

**F4, the censoring refusal, addressed in documentation rather than in
code, and here is why.** Core's guard is at `R/frame.R:1309` and
`valid_y` is called at `R/frame.R:1406`, so core refuses before any
family-supplied check runs; there is no seam at which `wiener_gng()`
can substitute a message. The only mechanism that would work is
declaring an `lcdf` that stops, and that is worse than the problem: it
would make frame assembly ADMIT left censoring, interval censoring and
`trunc()`, then fail later from inside the objective, and it would arm
a stop in a slot any other machinery may call. So the Rd's Censoring
section now quotes core's generic message verbatim, says it does not
name this family, and tells the reader to read it as a decision rather
than an omission. The claim "refused by name" is gone.

**F5, the dropped guidance, measured and largely not dropped.** Both
pieces still reach a user at the point of failure:

| what a user does | what they are told now |
|---|---|
| `gddm`, `vint(cond, upper)`, the wrong order | "the condition index must be a positive integer ... It is **the second value of vint()** here, and gddm_conditions() builds one." |
| `wiener`, `vint()` given 1/2 rather than 0/1 | "dec() coerces a factor or a character vector for you, **taking its second level as the upper boundary**; vint() does not, so recode it with as.integer(...)." |

Both survive in checks this lane did not delete, and both are in the
reference pages as well (`man/gddm.Rd:74` shows the `vint(upper, cond)`
ordering, `man/gddm.Rd:77` the factor convention). The only case with
no such guidance is a model supplying NO addition term at all, where
the declared refusal names the terms and one spelling. That is the
trade the brief asked for. No code change.

**F7, the sv attribution, softened.** The Rd no longer rests on the
contrast alone. It now records that a grid search on two-boundary data
from the same generative parameters finds the same corner, `sv` near
zero beating the truth by 2.4 log units against 1.8 on the go/no-go
data, and says that the ridge is a property of the drift-diffusion
likelihood rather than only of this design, with the go/no-go design's
contribution being the loss of the boundary-proportion information.
The practical advice is unchanged.

## Punch round verification

Re-run after every change above. Nothing moved.

### The seven named files, one process each

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-rdm-gng.R | 221 | 0 | 0 | 0 |
| test-defects.R | 59 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 |
| test-variability.R | 140 | 0 | 0 | 0 |
| test-gddm-family.R | 77 | 0 | 0 | 0 |
| test-family.R | 35 | 0 | 0 | 0 |

Identical to the pre-punch counts, which is the point: the `bs` fix
changes no answer and the rest is prose.

### Whole suite, one process, NOT_CRAN=true, counts by name

| file | passed | failed | error | skipped |
|---|---|---|---|---|
| test-brms-parity.R | 13 | 0 | 0 | 0 |
| test-defects.R | 59 | 0 | 0 | 0 |
| test-density.R | 132 | 0 | 0 | 0 |
| test-family.R | 35 | 0 | 0 | 0 |
| test-gddm-family.R | 77 | 0 | 0 | 0 |
| test-gddm-gradients.R | 22 | 0 | 0 | 0 |
| test-gddm-recovery.R | 28 | 0 | 0 | 0 |
| test-gddm-reference.R | 166 | 0 | 0 | 0 |
| test-gddm-solver.R | 94 | 0 | 0 | 0 |
| test-lba.R | 109 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 |
| test-moments.R | 29 | 0 | 0 | 0 |
| test-rdm-gng.R | 221 | 0 | 0 | 0 |
| test-sampling.R | 94 | 0 | 0 | 1 |
| test-simulate-density.R | 69 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 |
| test-variability.R | 140 | 0 | 0 | 0 |
| **total** | **1339** | **0** | **0** | **1** |

The one skip is `test-sampling.R`, "frm_sample runs a short chain on a
wiener model", the same skip 0.3.0 had.

### Core's test-compat.R with the extension loaded

**273 pass, 0 fail, 0 error, 0 skip.**

### roxygen

Idempotent: two consecutive `roxygenise()` runs leave every `man/*.Rd`
and `NAMESPACE` byte-identical by md5. The one pre-existing warning
remains, `gddm.R:169` cannot resolve `gd_tri_df`, at a line this lane
does not touch.

### R CMD check --as-cran, with the manual actually building

`_R_CHECK_CRAN_INCOMING_=false`, core in `ea-lib`, no `--library=` on
check, `RSTUDIO_PANDOC` and
`C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows` prepended to PATH.

**Status: 1 NOTE**, down from 1 ERROR, 1 WARNING and 2 NOTEs before
TinyTeX was on PATH. The remaining NOTE is
`checking HTML version of manual ... Skipping checking math rendering:
package 'V8' unavailable`, which is the machine and not the package.

What the LaTeX toolchain settles, now that it runs:

* `checking PDF version of manual ... OK` - the earlier ERROR and
  WARNING were `pdflatex is not available` and nothing else. The Rd is
  clean.
* `checking for non-standard things in the check directory ... OK` -
  the earlier NOTE was the `frmtmb.eam-manual.tex` left by the failed
  PDF build.
* `checking tests ... [558s] OK`
* `checking re-building of vignette outputs ... [34s] OK`, which is the
  rewritten `ddm.Rmd`.
