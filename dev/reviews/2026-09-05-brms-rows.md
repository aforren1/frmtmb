# Review: lane wt-brms-rows (brms log-density identity, second round)

Reviewer notes, written incrementally. Worktree
`C:/Users/adf44/source/r/frmtmb-wt-brms-rows`, branch `wt-brms-rows`,
branched at 564e185. Main checkout at 5dfdd84, clean, untouched.

Scope confirmed: four modified files, no untracked files, `R/` untouched.

```
NEWS.md                               |  14 +
dev/brms-likelihood-tests.md          | 417 ++++++++++++++++++++-
tests/testthat/helper-brms.R          | 554 +++++++++++++++++++++++++++-
tests/testthat/test-brms-likelihood.R | 666 ++++++++++++++++++++++++++++++++++
```

## 1. Tier run, warm, one process

Reproduced. Private library `scratchpad/rbr-lib` with the worktree's
core installed; `FRMTMB_STAN_CACHE` pointed read-only at the lane's
`dev/stan-cache`.

| quantity | claimed | measured |
| --- | --- | --- |
| tests | 31 | 31 |
| assertions | 347 | 347 |
| failed / errors / skipped | 0 | 0 / 0 / 0 |
| wall, warm | 51.0 s | 49.1 s |
| programs compiled | 0 (36 restored) | 0 |

The cache directory listing, `ls -la` with full ISO timestamps, is
byte-identical before and after the run: 36 `.rds` files, no new file,
no mtime moved. Nothing was recompiled.

**Correction to the write-up.** The gate is not one variable. A run with
`FRMTMB_BRMS_FIT_TESTS=true` alone skips all 31 blocks with
`Reason: On CRAN`, because `skip_unless_brms_fit()` calls
`skip_unless_brms()`, which calls `testthat::skip_on_cran()`
(`tests/testthat/helper-brms.R:10`). `NOT_CRAN=true` is also required
outside `R CMD check`. The tier file's own header comment
(`tests/testthat/test-brms-likelihood.R:31-32`) documents only
`FRMTMB_BRMS_FIT_TESTS`; the plan doc's results table does say
"`FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`". The header is the
one a reader hits first and it is incomplete. See punch list.

## 2. The two CAR normalizers

Derived from the compiled program text and brms's standata, without
reading `brms_car_const()`. Both are **exact closed forms in the data
alone**, not fitted: neither depends on `sdcar`, `car`, or the field.

**escar.** brms's `sparse_car_lpdf` returns

```
0.5 * (Nloc*log(tau) + sum(log1m(car*eigenMcar))
       - tau*(phi'D phi - car*(phi'W phi)))
```

A proper CAR is `phi ~ N(0, (tau (D - car W))^-1)`, whose log density is
`-Nloc/2 log(2pi) + 0.5(Nloc log tau + logdet(D - car W)) - 0.5 tau q`.
Since `D - car W = D^(1/2) (I - car M) D^(1/2)` with
`M = D^(-1/2) W D^(-1/2)`, and brms ships the eigenvalues of that `M` as
`eigenMcar`, `logdet(D - car W) = sum(log Nneigh) + sum(log1m(car*eigenMcar))`.
Everything cancels term for term and what brms drops is

```
Nloc/2 log(2pi) - 0.5 sum(log Nneigh)
```

Both steps checked numerically: `max|eigenMcar - eig(M)| = 2.9e-15`, and
the log-determinant identity holds to 3.6e-15 at car = 0.1, 0.5, 0.9.
Value on the row's 4x4 rook lattice (Nloc 16, Nedges 24, degrees
4x2 + 8x3 + 4x4): **6.149684293242652**. Claim 6.149684293, delta
2.4e-10, which is the rounding of the claim.

**icar.** brms's target on `zcar` is
`-0.5*dot_self(zcar[edges1]-zcar[edges2]) + normal_lpdf(sum(zcar)|0, 0.001*Nloc)`,
i.e. `-0.5 z'Kz - 0.5 log(2pi) - log(s)` with `K = L + 11'/s^2`,
`s = 0.001*Nloc`, `L` the graph Laplacian of the edge list. A proper
`N(0, K^-1)` is `-Nloc/2 log(2pi) + 0.5 logdet K - 0.5 z'Kz`, so brms
drops

```
0.5(Nloc-1) log(2pi) - log(s) - 0.5 logdet K
```

Value **5.253269633695128**. Claim 5.253269634, delta -3.0e-10, again
rounding. The edge list `edges1`/`edges2` reproduces the adjacency
matrix exactly, so the constant is built from brms's own graph.

Both agree with `brms_car_const()`'s formulas line for line. The
derivation is sound and the numbers are right.

**One fragility, not a defect.** `brms_car_const()` hardcodes
`s <- 0.001 * n` rather than reading `con_sd` off the fit. brms's 0.001
is a literal in its generated Stan; frmtmb's is a default. The test does
pin it (`expect_identical(brms_car_block(fit2)$aux_car$con_sd, 0.001)`,
`test-brms-likelihood.R:1036`), so a change to frmtmb's default fails
loudly rather than silently shifting the constant. Acceptable as written.

## 3. esicar

Reproduced on the installed core alone, on the row's 4x4 lattice, all
four types frmtmb's `car()` accepts:

| type | logLik | df | theta | sd(car) |
| --- | --- | --- | --- | --- |
| escar | -92.292141924828968 | 5 | 0.37572404271851972, 0.29272851263966632 | 1.456045, car 0.572664 |
| esicar | -90.316306819585037 | 4 | 0.34600842212903699 | 1.413415 |
| icar | -90.316306819585037 | 4 | 0.34600842212903699 | 1.413415 |
| bym2 | -90.316306831849971 | 5 | -0.029771575329913876, 20.744440327292085868 | 0.9706672, rhocar 1.000000 |

esicar and icar agree **bit for bit**: logLik gap 0, theta gap 0, beta
gap 0. Every field of the block's `aux_car` is equal except the `type`
string itself: `n`, `W`, `L`, `deg`, `con_sd`, `kappa0`, `Sgrp`, `K`,
`ldet_K`, `n_comp` all identical. The lane's claim is exact.

**Is the sum-to-zero constraint implemented at all? Yes, softly, and for
both spellings identically.** `car_aux()` (`R/covstruct.R:1204`) has no
`esicar` branch: `escar` returns early at `R/covstruct.R:1209-1230`, and
everything else falls into the shared intrinsic path at
`R/covstruct.R:1231-1245`, which builds

```
K <- L + t(Sgrp) %*% Diagonal(n_comp, kappa0) %*% Sgrp
kappa0 <- 1 / (con_sd * nj)^2
```

That rank-`c` update is the constraint. It is soft, never hard. There is
no code path anywhere that declares `Nloc - 1` free values.

**Classification: documented choice, with a brms-compatibility gap the
user should know about.**

The design note is explicit (`R/covstruct.R:916-947`):

> We adopt brms's remedy - a soft sum-to-zero constraint whose precision
> rides on tau, as it does in brms's non-centered zcar parameterization
> ... The component sums are then pinned at an sd of `con_sd n_j sdcar`
> rather than exactly zero, so the fit approaches the hard-constrained
> (brms esicar) likelihood as `con_sd -> 0`; `esicar` selects the same
> density.

with a measured bias table: at the 1e-3 default the log-likelihood is
off the hard-constrained reference by 4.7e-4 (3.6e-5 relative in
`sdcar`), and tightening costs optimizer robustness (nlminb false
convergence 0 times in 25 refits at 1e-3, once at 1e-4, 6 times at
1e-5). The user-facing text says the same thing
(`vignettes/frmtmb.Rmd:150-151`, `:178-183`):

> `"icar"` and its alias `"esicar"` (intrinsic)

> An intrinsic CAR is improper ... so `icar`, `esicar` and `bym2` add
> brms's soft sum-to-zero constraint, whose scale is `con_sd`.

So this is deliberate, argued, and quantified. The gap is one of
labeling, not of arithmetic: the same vignette section opens with

> `car()` fits a conditional autoregressive field over an adjacency
> matrix, with brms's spelling and **all four of its types**

and frmtmb has three distinct densities under four spellings, while in
brms `esicar` is not an alias for `icar` - brms declares `Nloc - 1` free
values and sets `rcar[Nloc] = -sum(zcar)`, which the lane's test asserts
directly (`test-brms-likelihood.R:1093`). A user porting a brms `esicar`
model gets a different (softly constrained) likelihood under the same
call. Worth one sentence in the vignette; see the defect list.

Incidental, outside the lane's scope: the row's `bym2` fit puts
`rhocar` on the boundary at 1.0 (theta2 = 20.7) and lands on icar's
log-likelihood to 1.2e-8. Not a finding about this lane, but the row's
data is fully spatial, so the mixture is not being exercised.

## 4. The exact-GP nugget

**The 216.6 nats reproduce exactly.** On the row's design (seed 23, 80
points of `sin(x)` on [0, 10]), at frmtmb's estimates
`sdgp = 0.91365746`, `lscale = 1.73170388` (data scale),
`dmax = 9.768956`, so brms's scale is 0.177266:

| quantity | measured | doc |
| --- | --- | --- |
| joint gap, brms - frmtmb | 216.6043 | 216.6 |
| logdet frmtmb | -967.11 | -967.1 |
| logdet brms | -1826.56 | -1826.6 |
| quadratic form frmtmb | 7.108 | 7.1 |
| quadratic form brms | 433.32 | 433.3 |
| max off-diagonal kernel gap | 0 (same expression) | 3.9e-16 |

Every number is right. Two prose nits, no more: the sentence "216.6
nats, of which 859.4 is log-determinant" is loose, since 859.4 is the
raw log-determinant difference and its *contribution* is half that,
429.72, offset by -213.11 from the quadratic form (429.72 - 213.11 =
216.61). And the eigenvalue count I measure with the test's own
expression, `eigen(kf / sdgp^2)$values < 1e-6`, is **30**, not 29. The
test only asserts `> 10`, so nothing depends on it.

**What the nugget does to a FIT.** `gp_corr()` was replaced in the
installed namespace in a scratch process by `assignInNamespace`; the
tree was not touched. Two designs, nugget at 1e-6 (frmtmb's) and 1e-12
(brms's).

Row 10 design, fitted lscale ~1.73, 30 of 80 correlation eigenvalues
below 1e-6:

| nugget | logLik | sdgp | lscale | sigma |
| --- | --- | --- | --- | --- |
| 1e-6 | -31.809111359 | 0.91365746 | 1.73170388 | 0.3038748 |
| 1e-12 | -31.808525538 | 0.91162184 | 1.72739776 | 0.3039232 |
| shift | +5.86e-4 | -0.223% | -0.249% | +0.016% |

Long-length-scale design (`2 sin(x/5) + 0.3x`, fitted lscale ~12.2, 75
of 80 eigenvalues below 1e-6 and 72 below 1e-12):

| nugget | logLik | sdgp | lscale | sigma |
| --- | --- | --- | --- | --- |
| 1e-6 | -32.993657016 | 4.80882840 | 12.19501382 | 0.3223246 |
| 1e-12 | -32.994143942 | 4.63557344 | 11.96783216 | 0.3225587 |
| shift | -4.87e-4 | -3.603% | -1.863% | +0.073% |

**Answer: only the joint density. The nugget does not measurably bias
the fit.** The largest parameter move is 3.6% of `sdgp` on the design
built to be hostile, and it moves the log-likelihood by 4.9e-4 nats
against the 1.92 nats of a one-parameter 95% confidence boundary: about
0.03% of a confidence interval in likelihood units. `sigma`, the
parameter anyone would actually report, moves in the fourth decimal.
This is the lane's claim and it holds.

**A result the lane did not have, and it strengthens the write-up.**
At 1e-12 the fit does not converge. On *both* designs the 1e-12 run
raises

```
Optimizer did not report convergence: false convergence (8)
Large maximum absolute gradient at the optimum (0.117)   # row 10
Large maximum absolute gradient at the optimum (0.157)   # long
```

while the 1e-6 run raises neither. So part of the small shift tabulated
above is optimizer failure rather than bias, and the true bias is
smaller still. More to the point, the doc's own sentence

> that is a real problem: at 1e-12 the factorization of an 80-point
> kernel is already marginal

is not a hedge, it is a measured fact: frmtmb cannot fit at brms's
constant on either design. "Match brms's constant" is therefore not one
of the three options the maintainer has; it is the one option that
demonstrably does not work as the code stands. The doc should say so
rather than list it first among equals.

**Classification: documented defect, low severity, but the documentation
half is missing.** The nugget is a white-noise term with SD
`1e-3 * sdgp` that is part of the model and is not in any user-facing
document. `R/covstruct.R:1411-1413` explains it to a maintainer ("The
nugget on the diagonal keeps the factorization stable at long
lengthscales") and nothing tells a user that `gp()` is not the squared
exponential the vignette advertises. The arithmetic consequence is
negligible; the documentation gap is real.

## 5. The length-scale scale "inconsistency" - the write-up is wrong

The *fact* is right and the translator rule that follows from it is
right. The *characterization* is not, and it should be corrected before
this lands.

**What is true.** frmtmb parameterizes the two GP forms internally on
different scales. `gp_corr()` (`R/covstruct.R:1415-1422`) reads raw
squared differences out of `blk$aux_D2`, so the exact form's `theta` is
a data-scale length-scale. The Hilbert-space form divides the
coordinates by `dmax` at frame construction
(`R/frame.R:1763-1772`, `gp_max_dist()` reproducing
`brms:::.data_gp` bit for bit), so its `theta` is on brms's unit-maximum-
distance scale. Measured on the row's design, `dmax = 9.768956024`:

| form | internal `exp(theta[2])` | block `gp_dmax` |
| --- | --- | --- |
| `gp(x)` | 1.731703876 | NULL |
| `gp(x, k = 25)` | 0.177265455 | 9.768956024 |

The translator's `lscale_<i> = exp(theta[2]) / dmax_<i>` for the exact
form is therefore correct and necessary, and the row's 4e-16 kernel
agreement confirms it.

**What is false.** The doc
(`dev/brms-likelihood-tests.md:719-723`) says:

> The inconsistency is therefore internal as well: `gp(x)` and
> `gp(x, k = 10)` **report** length-scales on two different scales in the
> same package, and only the second is comparable with brms's posterior
> for the same term.

Both clauses fail. `confint_varcorr()` converts the Hilbert-space
length-scale back to data units, deliberately and exactly
(`R/confint.R:498-510`):

```r
# hsgp estimates the lengthscale on brms's rescaled inputs, but the
# reported range belongs in data units. The scale factor is a data
# constant, so the shift on the log scale is exact and the se rides
# through unchanged. The exact gp keeps the raw scale (dmax NULL).
log_dmax <- log(bk[["gp_dmax"]] %||% 1)
...
add(term_j, "range", t0[1 + j] + log_dmax, se_t[1 + j], bk)
```

Measured, same data, same seed:

| form | reported `range(gp)` | reported / internal |
| --- | --- | --- |
| `gp(x)` | 1.731703876 | 1 |
| `gp(x, k = 25)` | 1.731698435 | 9.768956024 = dmax |

The two agree to 3e-6 relative. Nothing user-visible is inconsistent,
and the vignette's promise (`vignettes/frmtmb.Rmd:133`)

> Reported lengthscales stay in data units.

is exactly what the code implements, for both forms. The second clause
is wrong too: brms rescales its coordinates for the exact form as well
(the row's own assertion that `Xgp_1` has unit maximum squared distance
is the proof), so brms's `lscale` is unit-scale for *both* its forms and
neither of frmtmb's *reported* ranges is directly comparable with it
without multiplying by `dmax`. What is directly comparable with brms is
frmtmb's *internal* hsgp theta, which is not a reported quantity.

**Classification: neither a defect nor a documentation gap in frmtmb. A
write-up error in this lane's own doc.** The docs promise data units and
deliver them. The one-line replacement is on the punch list.

## 6. sigma under ar(cov = TRUE)

Both halves confirmed.

**brms's sigma is the innovation SD.** From the generated program for
the row's model:

```stan
matrix cholesky_cor_ar1(real ar, int nrows) {
  ...
  return cholesky_decompose(mat ./ (1 - ar^2));
}
...
Lcortime = cholesky_cor_ar1(ar[1], max_nobs_tg);
target += normal_time_hom_lpdf(Y | mu, sigma, Lcortime, nobs_tg, begin_tg, end_tg);
```

and `normal_time_hom_lpdf` forms `Cov = multiply_lower_tri_self_transpose(sigma * chol_cor)`.
So the covariance is `sigma^2 * rho^|i-j| / (1 - rho^2)`: the marginal
variance is `sigma^2 / (1 - rho^2)` and `sigma` multiplies the
innovation. brms's own comment above the function calls it "the cholesky
factor of an AR1 correlation matrix", which it is not - a correlation
matrix has a unit diagonal and this one has `1/(1 - ar^2)`. The lane is
right and brms's naming is what misleads.

**frmtmb's sigma is the marginal SD.** `R/autocor.R:16-30`:

> R is UNIT-DIAGONAL here, so `sigma` is the MARGINAL residual standard
> deviation, as it is in nlme and as it is everywhere else in this
> package (`sigma()`, pearson residuals, `se()`, every distributional
> model). brms diverges: its `cholesky_cor_ar1(ar, n)` returns
> `chol(rho^|i-j| / (1 - rho^2))` ... so under `ar()`/`ma()`/`arma()`
> brms's `sigma` is the INNOVATION sd.

and the exported roxygen carries it to the user, `R/autocor.R:161-169`,
`@section Divergence from brms`:

> brms parameterizes `ar()`, `ma()` and `arma()` by the INNOVATION
> standard deviation (its `cholesky_cor_ar1()` divides by `1 - ar^2`),
> while `cosy()` and `unstr()` use the marginal one. Here every
> structure uses the marginal `sigma` ... the scales relate by
> `sigma_marginal = sigma_innovation / sqrt(1 - phi^2)` for AR(1).

The test asserts both ends: the Stan line by `expect_match`
(`test-brms-likelihood.R:894`) and the conversion to 1e-12
(`test-brms-likelihood.R:896-898`). Sound.

**One accuracy point.** The doc heads this
"### Finding: brms's sigma under ar(cov = TRUE) is the innovation SD"
and reads as a discovery. It is not: frmtmb documented this divergence
before the lane started, in a source comment and in an exported
`@section`. What the lane genuinely adds is the first *numeric*
verification of that documented claim (5.2834 nats and gradient 58.17
without the factor), which is worth more than a discovery would be. The
heading should say so; see the punch list.

## 7. Merge overlap with wt-mo-terms

Both lanes branch at the same commit, 564e185, so this is a clean
three-way merge. I ran it: base from `git show 564e185:<path>`, the two
sides from the two worktrees, `git merge-file -L rows -L base -L mo`.

### Result

| file | textual conflicts | note |
| --- | --- | --- |
| `tests/testthat/test-brms-likelihood.R` | **0** | merges clean, parses |
| `tests/testthat/helper-brms.R` | **0** | merges clean, parses |
| `NEWS.md` | 1 | both lanes opened the same new section |
| `dev/brms-likelihood-tests.md` | 1 | both lanes rewrote the Status paragraph |

Why the two R files do not conflict:

- `test-brms-likelihood.R`. The rows lane is a **pure append**: one hunk,
  `@@ -439,0 +440,666 @@`, no deletions and no context change. The mo
  lane touches only base lines 23-31 (the header) and 172-195 (row 3).
  Disjoint.
- `helper-brms.R`. The mo lane has two hunks, `@@ -199,0 +200,17 @@`
  (`brms_mo_terms_of()`) and `@@ -434,3 +451,21 @@` (the `simo_<j>`
  branch). The rows lane's nearest hunks are `@@ -188,13 +188,17 @@` and
  `@@ -452 +471,166 @@`. The first pair share context but not a changed
  line - rows edits base 191 and 194-197, mo inserts after base 199 -
  and base 198-199 are untouched by both, which is enough separation.
  The second pair are 16 lines apart. Both lanes' rules survive: the
  merged `stan_pars_from_fit()` carries mo's `^simo_` branch and rows's
  `mi_map`, `us_chol_L`, `autocor_natural`, `brms_lrescor` and
  `brms_car_const` branches as separate arms of the same else-if chain.

### The merged tier actually runs

I built the merged pair into a scratch test directory and ran it against
the **rows** lane's installed core, warm cache, one process:

```
blocks: 31   passed: 349   failed: 3   errors: 0   skipped: 0
only failing block: "row 3: mo(inc) * z is the same model in both packages"
```

All three failures are inside that one block and are exactly the
expected dependency: it asserts two simplexes and a
`const = 2 * lgamma(3)` identity, which needs the mo lane's `R/frame.R`
change. Against a core that still shares a simplex per variable it
misses by design (gradient 1132). **Every one of the rows lane's 30
other blocks passes unchanged under the merged helper.** Assertions go
347 -> 349, the two added by mo's rewritten row 3.

So: merge the R files with git and do not hand-edit them. Install the mo
lane's core and the block goes green.

### The semantic conflict git will not show you

This is the important part. `test-brms-likelihood.R` merges **clean and
wrong**. The merged file's header takes the mo lane's text:

> That list is now empty. It held one entry, row 3's `mo(inc) * z`,
> until every mo() TERM was given its own simplex; **every row below is
> an identity.**

while the merged file below it carries three `EXEMPTION` blocks at lines
780, 934 and 1048, the last covering two shapes: **four divergences
across three rows.** Both lanes are individually right and the merge is
false.

The rows lane anticipated this and left the header alone on purpose
(`dev/brms-likelihood-tests.md:16-18`: "That comment belongs to another
lane this round and was left alone deliberately; it needs to say four").
Good discipline, but it means the mo lane's text wins silently.

### Reconciled text

**`tests/testthat/test-brms-likelihood.R`, header, replacing merged
lines 25-30 (the mo lane's paragraph):**

```
# There is no known-divergence list for numeric mismatches. The only
# admissible exemption is a design choice frmtmb states in its own
# source, recorded here with its reason, and every one is recorded by
# asserting the structural difference rather than by skipping the row.
# That list held one entry, row 3's `mo(inc) * z`, until every mo()
# TERM was given its own simplex; that row is now an identity. Four
# entries remain, over three rows: the exact `gp()` nugget (row 10a),
# brms's `ar(cov = FALSE)` likelihood (row 18d), and the esicar and
# bym2 CAR parameterizations (row 19c).
```

**Same file, add to the opt-in note at merged lines 32-33**, which
currently names only one variable and is wrong on its own terms (see
section 1):

```
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
```

**`NEWS.md`.** Not a real conflict: both lanes opened a
`# frmtmb (development version)` section under one heading and git could
not order the bullets. Keep **both** bullet sets under the single
heading, mo's first because a BEHAVIOR CHANGE leads, then rows's. Then
fix the two bullets that now contradict:

- mo's last bullet ends "and the tier's exemption list is empty."
  Replace with: "and the tier's exemption list no longer holds a `mo()`
  entry."
- rows's bullet says "Three divergences the tier found are recorded and
  asserted rather than skipped". Replace "Three divergences" with "Four
  divergences, over three rows," to agree with the doc and the header.
  (This inconsistency is internal to the rows lane, not caused by the
  merge: its NEWS says three, its doc says four, its coverage table
  lists three EXEMPT rows. All three countings are defensible; pick one.)

**`dev/brms-likelihood-tests.md`.** The conflict is the Status
paragraph, which both lanes rewrote wholesale. Take the rows lane's
version, which is written for the post-merge world, and fold in mo's
resolution. Reconciled:

```
Status: implemented in two rounds, on branches `wt-brmslp` and
`wt-brms-rows`, 2026-09-05. The harness, the translator with its
round-trip unit tests and the CI workflow came from the first round,
with rows 1, 2, 3, 5, 7, 12, 13, 14, 15, 16, 17, 20 and 21. The second
round added rows 4, 6, 8, 9, 10, 11, 18 and 19, which is every
remaining row of the matrix. The first round's one exemption, row 3's
`mo(inc) * z`, was closed on branch `wt-mo-terms`: frmtmb now builds
one simplex per mo() TERM and row 3b asserts the identity with the flat
Dirichlet admitted once per simplex. Four real divergences between the
packages remain, over three rows, and each is asserted rather than
skipped. Results, deviations, findings and the remaining work are under
"Implementation log" at the end of this document, and the second round
has its own section there.
```

Drop the rows lane's two-sentence note about the header comment - it is
discharged once the header is fixed.

Then three stale lines survive the clean merge and need a touch:

- **merged line 504**, in "Coverage against the plan's matrix": "Row 3's
  `mo(inc) * z` spelling is the documented exemption above, and it is
  asserted structurally rather than skipped." This is a live status
  statement and is now false. Replace with: "Row 3's `mo(inc) * z`
  spelling was the first round's one exemption and is now an identity;
  see the RESOLVED note above."
- **merged line 444**, "`3a` and `C3` are the one documented exemption",
  and **merged line 308**, "This is the first and only entry on the
  exemption list the plan keeps." Both sit inside first-round history,
  and the mo lane's RESOLVED note at merged line 259 already says to
  read that section as the record of what was found. Leave them, or
  prefix line 308's paragraph with "At the time of the first round,".
  Not blocking.
- **merged line 816**, "it is the only one of the four exemptions that a
  user cannot walk into by accident" - already correct at four. No
  change.

The mo lane's other doc hunk (`@@ -254,2 +259,8 @@`, the RESOLVED note)
does not overlap anything the rows lane touched and needs no attention.

## Defect list for the user

Three items the lane surfaced, with my classification. None blocks this
lane; all three are about frmtmb, not about the tests.

### D1. `esicar` is an alias for `icar` and brms's is not - DOCUMENTED CHOICE, labeling gap

frmtmb's `car(type = "esicar")` returns a fit identical to
`type = "icar"` bit for bit (verified: logLik, theta and beta gaps all
exactly 0, every `aux_car` field equal but the type string). The
sum-to-zero constraint is implemented, softly, as a rank-`c` update
`K = L + t(Sgrp) Diagonal(kappa0) Sgrp` with
`kappa0 = 1/(con_sd * nj)^2` (`R/covstruct.R:1231-1245`); there is no
hard-constrained code path. This is deliberate, argued and quantified
in the source (`R/covstruct.R:916-947`) and stated in the vignette
(`vignettes/frmtmb.Rmd:150-151`, "`\"icar\"` and its alias
`\"esicar\"`"). The measured bias against a hard constraint at the 1e-3
default is 4.7e-4 in log-likelihood, and tightening costs optimizer
robustness.

The gap is that brms's `esicar` is *not* an alias for its `icar` - brms
declares `Nloc - 1` free values and sets `rcar[Nloc] = -sum(zcar)` -
while the same vignette section claims `car()` supports "brms's spelling
and all four of its types". A user porting a brms `esicar` model gets a
different likelihood under the same call, silently.

**Recommended:** one sentence in `vignettes/frmtmb.Rmd` after line 183:
"brms's `esicar` imposes the constraint exactly, on `Nloc - 1` free
values; here `esicar` and `icar` are the same soft-constrained model, so
a brms `esicar` fit and a frmtmb one are close but not identical."

### D2. The exact-GP nugget is part of the model and is undocumented - DOCUMENTATION GAP

`gp_corr()` (`R/covstruct.R:1421`) adds `diag(1e-6, n)` to the
*correlation*, i.e. a white-noise term of SD `1e-3 * sdgp` on the
covariance. brms adds an absolute `1e-12`. Nothing user-facing says the
kernel is not the squared exponential the vignette advertises.

**Severity: low on the arithmetic, real on the documentation.** I
measured the effect on fits, which the lane had not: it is negligible.
Largest parameter move across two designs is 3.6% of `sdgp` on a
deliberately hostile long-length-scale design, worth 4.9e-4 nats against
the 1.92 nats of a 95% confidence boundary; `sigma` moves in the fourth
decimal. The joint density of a fixed field moves by 216.6 nats, but
that is not a quantity anyone reports.

**And the obvious remedy does not work.** At 1e-12 the fit fails to
converge on *both* designs ("false convergence (8)", max abs gradient
0.117 and 0.157), where 1e-6 converges clean. So "match brms's constant"
is not a live option as the code stands. The right fix is to document
the nugget, not to shrink it.

### D3. The length-scale "inconsistency" - NOT A DEFECT

The lane reports that `gp(x)` and `gp(x, k =)` report length-scales on
different scales. They do not. `confint_varcorr()` converts the
Hilbert-space length-scale back to data units exactly
(`R/confint.R:498-510`); measured on one design, `gp(x)` reports
`range(gp) = 1.731703876` and `gp(x, k = 25)` reports `1.731698435`,
agreeing to 3e-6 relative. The vignette's promise (`frmtmb.Rmd:133`,
"Reported lengthscales stay in data units") is what the code delivers.
Only the internal `theta` differs in scale, which is an implementation
detail the translator correctly absorbs. **The lane's own doc needs
correcting, not frmtmb.**

## Punch list

Ordered: must-fix first, then should-fix, then nits. File and line
refer to the **rows worktree** unless the item says "merged", which
means the file after merging with wt-mo-terms.

### Must fix before the merge lands

**P1. The merged header will claim an empty exemption list while four
exemptions sit below it.** `tests/testthat/test-brms-likelihood.R`,
merged lines 25-30. `test-brms-likelihood.R` merges with **zero** git
conflicts and the mo lane's text wins silently: "That list is now empty
... every row below is an identity", above `EXEMPTION` blocks at merged
780, 934 and 1048. Reconciled text in section 7. This is the single
highest-value item in the review: git will not warn anyone.

**P2. `dev/brms-likelihood-tests.md:719-723` states something false
about frmtmb.** The claim that `gp(x)` and `gp(x, k = 10)` "report
length-scales on two different scales in the same package, and only the
second is comparable with brms's posterior" is wrong on both clauses
(section 5, D3). Replacement:

```
The inconsistency is internal only. `confint_varcorr()` shifts the
Hilbert-space length-scale back to data units by `log(gp_dmax)`
(R/confint.R:498-510), so both forms REPORT `range(gp)` in data units -
1.731704 and 1.731698 on the row's design - and the vignette's promise
that "reported lengthscales stay in data units" holds. What differs is
the internal theta, which the translator absorbs. brms rescales for both
of its forms too, so neither of frmtmb's reported ranges is directly
comparable with a brms `lscale` without multiplying by `dmax`.
```

### Should fix

**P3. Count the exemptions the same way in all three places.**
`NEWS.md` says "Three divergences", `dev/brms-likelihood-tests.md:11`
says four, the coverage table (`:898-905`) lists three EXEMPT rows. All
defensible, none consistent. Use "four divergences over three rows"
everywhere, including the P1 header.

**P4. `tests/testthat/helper-brms.R:14-15`** carries the same incomplete
gate note I fixed in the tier file:

```r
# Numeric tier: Stan compilation, minutes per model. Opt in with
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
```

`skip_unless_brms_fit()` calls `skip_unless_brms()` calls
`skip_on_cran()`, so `NOT_CRAN = "true"` is needed too. One line. I left
it alone to keep my edit to a single hunk; see "Edits made".

**P5. `dev/brms-likelihood-tests.md:752`**, "### Finding: brms's sigma
under ar(cov = TRUE) is the innovation SD". Not a finding - frmtmb
documented this divergence before the lane began, in `R/autocor.R:16-30`
and in the exported `@section Divergence from brms`
(`R/autocor.R:161-169`). Retitle to "### Confirmed: brms's sigma under
ar(cov = TRUE) is the innovation SD" and open with one sentence saying
the lane verifies a documented claim numerically for the first time
(5.2834 nats, gradient 58.17). That is a stronger result than a
discovery and the section currently undersells itself.

**P6. Merged `dev/brms-likelihood-tests.md:504`** (coverage section)
still says row 3's `mo(inc) * z` "is the documented exemption above, and
it is asserted structurally rather than skipped". False after the merge.
Replacement in section 7.

**P7. `dev/brms-likelihood-tests.md:663-702`, the gp nugget section,
should record that brms's constant is not reachable.** The section lists
three options for the maintainer, "to match brms's constant" first. I
measured it: at 1e-12 the fit does not converge on either of two designs
(false convergence 8, max abs gradient 0.117 and 0.157). Add the
measurement and drop that option to last. Numbers in section 4.

### Nits

**P8. `dev/brms-likelihood-tests.md:684`**, "216.6 nats, of which 859.4
is log-determinant" - 859.4 exceeds the total it is a part of. It is the
raw log-determinant difference; the contribution is 429.72, offset by
-213.11 from the quadratic form. Reword or add the halves.

**P9. `dev/brms-likelihood-tests.md:686`** says 29 of 80 eigenvalues sit
below the nugget. Measured with the test's own expression
(`eigen(kf / sdgp^2)$values < 1e-6`) it is **30**. The test asserts
`> 10` so nothing depends on it.

**P10. `brms_car_const()` (`helper-brms.R:973-982`) hardcodes
`s <- 0.001 * n`** instead of reading `con_sd` off the fit. Guarded by
`expect_identical(..., 0.001)` at `test-brms-likelihood.R:1036`, so it
fails loudly rather than drifting. Leave as is; noted so the next reader
does not have to re-derive it.

## 8. Scope and the three suites

**Scope.** `git diff --name-only 564e185` returns exactly the four
claimed files and nothing else:

```
NEWS.md
dev/brms-likelihood-tests.md
tests/testthat/helper-brms.R
tests/testthat/test-brms-likelihood.R
```

Untracked: `dev/review-brms-rows.md`, this file, which I created.
`R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `inst/` and `src/` are
untouched, confirmed by an explicit path-limited diff. The lane's claim
that `R/` belongs to other lanes this round is honored.

**Main checkout.** `C:/Users/adf44/source/r/frmtmb` is at
5dfdd84deb2c64937bd14ca95ba6bd2859929140 with a clean working tree, the
same as at the start of the review. I did not touch it.

**Suites, one file per process, against the worktree's installed core.**

| suite | blocks | assertions | failed | errors | skipped | s |
| --- | --- | --- | --- | --- | --- | --- |
| test-brms-agreement.R | 19 | 187 | 0 | 0 | 0 | 277.1 |
| test-message-uniqueness.R | 1 | 6 | 0 | 0 | 0 | 5.8 |
| test-bracket-access.R | 3 | 7 | 0 | 0 | 0 | 2.0 |

All three match the lane's claimed counts exactly (187, 6, 7).

**Tier re-run after my edit:** 31 tests, 347 assertions, 0 failures,
0 skipped, 67.6 s. Still green.

## A correction to my own section 1, and a note on the cache

The tier file needs **36** programs, as claimed. I briefly measured 37
and chased it down: the extra program is the two-simplex
`y ~ mo(inc) * z` model, and it was compiled by **my** merged-tier
experiment in section 7, not by this lane's tier. Timestamps settle it:
the lane's warm run finished 09:48:10 with 36 in the cache, the 37th
program was written at 10:02:56, and the merged run finished 10:03:50.
The program's Stan text confirms it - it declares `simo_2`, `Xmo_2` and
`Jmo[2]`, which only the mo lane's row 3b produces.

I deleted that one file and re-verified the directory listing is
byte-identical to the state I found: 36 `.rds` plus
`makevars-cxx17.mk`. **The lane's cache is as I found it.**

Worth carrying into the merge: the post-merge tier needs **37**
programs, because mo's row 3b compiles a shape the rows lane's row 3
only ever generated as text. Budget one more compile on the first CI
run after the merge.

## Edits made to the worktree

**One edit, one hunk.**

`tests/testthat/test-brms-likelihood.R`, lines 31-33. Before:

```r
# Stan compiles here. The whole file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
```

After:

```r
# Stan compiles here. The whole file is opt-in, and skip_unless_brms()
# calls skip_on_cran(), so outside R CMD check BOTH are needed:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
```

Why: following the comment as written runs the file to 31 silent skips
with `Reason: On CRAN`, which is how I lost the first run. `.github/
workflows/brms-likelihood.yaml` sets both (lines 49-51) and even asserts
both (lines 121-122), so CI is right and only the in-file guidance was
incomplete.

Verified after the edit: the tier is still 31/347/0, and the three-way
merge against wt-mo-terms still reports **zero** conflicts in both R
files, with my two-line change and the mo lane's header rewrite both
landing. I deliberately did *not* fix the exemption-count sentence three
lines above, because the mo lane rewrites exactly those lines and an
edit there would turn a clean merge into a conflicting one for no gain.
That fix belongs in the merge; the text is P1.

Nothing else in the worktree was modified. No commits, no staging, no
branch operations.

## Verdict: GO WITH FIXES

The work is sound and unusually well evidenced. Everything I could check
independently held up:

- the tier reproduces at 31/347/0 warm with no recompilation;
- both CAR normalizers are exact closed forms and I re-derived them from
  the Stan text without reading the helper, matching to 1e-10 (the
  rounding of the stated figures);
- the esicar/icar identity is bit-exact and is a documented design
  choice, not an accident;
- the 216.6-nat GP gap reproduces to four figures, and all four of the
  doc's decomposition numbers are right;
- the ar(cov = TRUE) innovation-SD finding is confirmed on both sides;
- the three suites match their claimed counts exactly;
- `R/` is untouched and main is clean and unmoved.

The lane is also honest where it matters: it records missing accessors
rather than reaching for `R/`, it asserts divergences instead of
skipping them, and its `brms_basis_map()` refuses a non-orthogonal
change of basis rather than letting a genuinely different model pass. It
even left the header comment alone on purpose to protect another lane's
merge. That is the right instinct and it is why this merges cleanly.

Two things keep it from a plain GO.

**P1 is the real one.** `test-brms-likelihood.R` merges with the mo lane
at zero git conflicts and the result is false: the header will claim an
empty exemption list directly above four asserted exemptions. Nothing in
the tooling will catch it. Fix it as part of the merge, not after.

**P2** is a statement in the lane's own doc that is wrong about frmtmb.
`confint_varcorr()` puts both GP length-scales in data units and I
measured them agreeing to 3e-6; the doc says they are reported on
different scales. Left standing it would send a maintainer to fix a bug
that does not exist.

P3 through P7 are worth doing before this lands but block nothing.
P8-P10 are nits.

For the maintainer, separately from this lane: D1 (esicar is an alias
here and is not in brms - one vignette sentence) and D2 (the exact-GP
nugget is an undocumented part of the model, and brms's 1e-12 is not a
usable alternative because the fit stops converging there).

---

# Punch re-check, 2026-09-05

Second pass over the same worktree after the lane's punch round. Same
rules: private library, absolute paths, the lane's cache reused
read-only, no commits. **P1 was deliberately skipped by the lane and
stays with the coordinator for merge time**; I did not touch it and it
is excluded from the verdict below. I made **no edits** this round.

Scope unchanged and still correct: the same four tracked files, `R/`,
`man/`, `NAMESPACE` and `DESCRIPTION` untouched (path-limited diff
returns 0 files), main at 5dfdd84 clean and unmoved, cache at 36
programs. The doc grew 417 to 510 inserted lines, the helper 554 to 560,
the tier file 666 to 671; NEWS stayed at 14.

## Item by item, against the diff

| item | claim | verdict |
| --- | --- | --- |
| P2 | doc retitled, confint.R text, draft named wrong | done |
| P3 | four divergences over three rows, 3 sites | done |
| P4 | helper 14-17, doc 118-122, my edit intact | done |
| P5 | doc retitled "Confirmed" | done |
| P6 | doc names wt-mo-terms | done |
| P7 | refit measurement + non-convergence, option last | done |
| P8/P9 | halves and eigenvalue count | done, recomputed |
| P10 | hardcoded 0.001 recorded | done |
| D1 | esicar reclassified | done, beyond what was asked |

**P2** (`dev/brms-likelihood-tests.md:734`). Retitled to "Finding: the
two GP forms carry different INTERNAL length-scales". The body at
:751-759 carries the `R/confint.R:498-510` correction, the two reported
values 1.731704 and 1.731698, and the vignette citation. :761-765 names
the earlier draft wrong on both clauses and says why it was kept rather
than deleted: "the wrong version named frmtmb as the inconsistent
party". That is the right call, the error is on the record instead of
quietly gone.

The three code comments the lane says stand do stand. I read all three:
`test-brms-likelihood.R:786` ("once the length-scale is put on brms's
scale", about the kernel comparison), `:813` ("frmtmb's length-scale is
dmax times brms's. That part IS translatable and the rule does it"), and
`helper-brms.R:565-569` ("An exact gp() measures distances on the DATA
scale in frmtmb and on brms's unit-maximum-distance scale in Stan ...
The HSGP form rescales its inputs the same way brms does and needs no
correction"). All three are about the internal theta and the
translator's `dmax` division. None claims anything about what is
reported. Correct to leave them.

**No residual stale claim anywhere.** I grepped NEWS and the doc for
"reported", "different scales", "lengthscale", "lscale" and
"range(gp)". The only hit inside the corrected region is :762, the
sentence naming the earlier draft wrong. Everything else is pre-0.22
NEWS history, parameter-name lists, or the corrected text.
Independently corroborating: `tests/testthat/test-gp-multidim.R:221`
already said the lengthscales "are estimated on the rescaled inputs but
reported in" data units, which predates this lane and is what the
corrected section now agrees with.

**P3.** All three sites carry "four ... over three rows":
`NEWS.md:10-11` ("Four divergences over three rows are recorded and
asserted rather than skipped"), `dev/brms-likelihood-tests.md:8-10`, and
the coverage table. The table site is better than I asked for: :961-963
adds an explicit reconciliation, "Three EXEMPT rows carrying four
divergences: row 19c holds two, esicar and bym2, which is why the count
of divergences and the count of rows that stand aside are not the same
number." That closes the ambiguity rather than merely picking a number.

**P4.** `helper-brms.R:14-17` names both variables and says why ("This
calls skip_unless_brms(), which calls skip_on_cran(), so outside R CMD
check BOTH variables are needed to opt in"). `doc:118-122` matches and
states the failure mode ("or the whole file skips silently"). My edit at
`test-brms-likelihood.R:31-33` is intact and unmodified.

**P5.** `doc:794` retitled "### Confirmed: brms's sigma under
ar(cov = TRUE) is the innovation SD", opening "Not a discovery" and
citing `R/autocor.R:16-30` and the exported `@section` at `:161-169`
down to the conversion formula, then claiming the right credit: "What
the lane adds is the first NUMERIC verification of that documented
claim, which is worth more than a discovery". The `cholesky_cor_ar1`
misnomer is called out at :811-814, matching what I read in the
generated program.

**P6.** `doc:484-486`: "Row 3's `mo(inc) * z` spelling was the first
round's one exemption and is now an identity, closed on branch
`wt-mo-terms`; see the RESOLVED note above." Correct.

**P7.** `doc:698-706` carries the refit measurement (0.22 and 3.60
percent on `sdgp`, 0.25 and 1.86 on the length-scale, `sigma` in the
fourth decimal, log-likelihood under 6e-4 nats against 1.92) and
:709-717 the non-convergence on both designs with gradients 0.117 and
0.157. :722-726 lists "match brms's constant" last and says why: "as the
code stands it is the one that demonstrably does not work". All of that
matches my measurements in section 4.

**P10.** `doc:1053-1054` records the hardcoded `0.001 * Nloc` and its
reason. The `expect_identical` guard is still at
`test-brms-likelihood.R:1036`.

**D1, esicar.** Reclassified more thoroughly than the punch list asked.
`doc:857-874` keeps the divergence, correctly, since the densities
really are functions of different arguments, but adds "The soft
constraint is a documented choice, not an omission", citing
`R/covstruct.R:916-947`, the bias table, and
`vignettes/frmtmb.Rmd:178-183`. `doc:886-893` adds the bit-for-bit
result with my numbers reproduced exactly (log-likelihood
-90.316306819585, `sdcar` 1.413415, "every field of `aux_car` is equal
except the `type` string") and the `car_aux()` citations. `doc:895-904`
states the labeling contradiction: the vignette's "all four of its
types" opening against three distinct densities under four spellings,
and "A user porting a brms `esicar` model gets a different, softly
constrained likelihood under the same call and is told nothing at the
call site."

## Recomputed, as asked

**Attribution correction.** The coordinator's note says "confirm the
eigenvalue count (30 vs your 29)". That is the wrong way round: **I
measured 30**, the lane's original draft said 29, and the punch round
adopted 30. There is no disagreement between us to resolve.

Recomputed from a fresh fit on the row's design:

| quantity | doc now says | I measure |
| --- | --- | --- |
| eigenvalues < 1e-6 | 30 of 80 | 30 |
| eigenvalues < 1e-12 | none | 0, min 1e-06 |
| raw logdet frmtmb | -967.1 | -967.1114 |
| raw logdet brms | -1826.6 | -1826.5580 |
| raw logdet gap | 859.4 | 859.4467 |
| logdet CONTRIBUTION | +429.72 | +429.7233 |
| raw quadform frmtmb | 7.1 | 7.1089 |
| raw quadform brms | 433.3 | 433.3232 |
| quadform CONTRIBUTION | -213.11 | -213.1072 |
| total | 216.6 | 216.6162 |

`+429.7233 - 213.1072 = +216.6161`, so the two halves now add up to the
total they are halves of, which was the whole of P8. Confirmed.

On P9, one clarification worth having on the record. The count depends
on which matrix "the correlation" means, and the doc's is the right one:
frmtmb's correlation is `exp(-Q) + diag(1e-6)` by `gp_corr()`'s own
definition, which is also the test's expression `eigen(kf / sdgp^2)`,
and on that matrix the counts are 30 below 1e-6 and 0 below 1e-12. On
the bare `exp(-Q)` they would be 66 and 60. The doc's sentence pairs "30
below 1e-6" with "none below 1e-12", which is true only of the nuggeted
matrix, so it is self-consistent as written.

## Runs, this pass

Private library `scratchpad/rbr-lib` with the worktree's core, the
lane's `dev/stan-cache` reused read-only via `FRMTMB_STAN_CACHE`, one
`test_file()` per process, `FRMTMB_BRMS_FIT_TESTS=true` and
`NOT_CRAN=true`.

| run | blocks | assertions | failed | errors | skipped | s |
| --- | --- | --- | --- | --- | --- | --- |
| test-brms-likelihood.R | 31 | 347 | 0 | 0 | 0 | 42.2 |
| test-brms-agreement.R | 19 | 187 | 0 | 0 | 0 | 217.5 |
| test-message-uniqueness.R | 1 | 6 | 0 | 0 | 0 | 6.6 |
| test-bracket-access.R | 3 | 7 | 0 | 0 | 0 | 2.5 |

Every count matches the lane's report and my first pass. The tier's
wall time differs from the lane's 71 s only through machine load; the
counts are identical.

**Cache untouched.** `ls -la` with full ISO timestamps before and after
the tier run is byte-identical: 36 `.rds` plus `makevars-cxx17.mk`, no
new file, no mtime moved. Nothing recompiled. I wrote nothing into the
lane's cache this pass.

## Merge, re-verified after the punch round

The punch round rewrote the doc's Status paragraph, which is exactly
where the mo lane also edits, so I re-ran the three-way merge from base
564e185 against `C:/Users/adf44/source/r/frmtmb-wt-mo-terms`:

| file | conflicts |
| --- | --- |
| `tests/testthat/test-brms-likelihood.R` | 0 |
| `tests/testthat/helper-brms.R` | 0 |
| `NEWS.md` | 1 |
| `dev/brms-likelihood-tests.md` | 1 |

Unchanged from my first pass. The two R files still merge clean; the two
markdown conflicts are the same two the recipe in section 7 covers. The
P4 helper edit at `helper-brms.R:14-17` is far from either mo hunk
(base 199 and base 434) and does not disturb the merge.

## Residual items

Nothing blocking. Four items, all for merge time, none a defect in the
lane's work.

**R1 (coordinator, merge time). P1 itself.** Untouched by design.
`tests/testthat/test-brms-likelihood.R:25-30` will take the mo lane's
"That list is now empty ... every row below is an identity" above three
`EXEMPTION` blocks at `:781`, `:935` and `:1049` carrying four
divergences. Reconciled text is in section 7 of this file. Git reports
zero conflicts here, so nothing will flag it.

**R2 (coordinator, merge time). Delete the P1 handoff note.**
`dev/brms-likelihood-tests.md:16-18` reads "The file-header comment in
`tests/testthat/test-brms-likelihood.R` still says that list has one
entry. That comment belongs to another lane this round and was left
alone deliberately; it needs to say four." Correct to keep while P1 is
open; it should go in the same commit that applies P1, or the doc will
describe a state that no longer exists.

**R3 (coordinator, merge time). The Status paragraph still does not say
who closed row 3.** `dev/brms-likelihood-tests.md:3-14` now says "Four
real divergences ... over three rows ... so the exemption list has
exactly those four", which is right, but it never says the first round's
one exemption was closed on `wt-mo-terms`. The doc says so elsewhere
(`:484-486`, and the RESOLVED note), so this is only about the paragraph
a reader hits first. The reconciled Status text in section 7 folds it
in; it is still the text to use.

**R4 (optional, lane).** `tests/testthat/test-brms-likelihood.R:1049-1060`,
the esicar/bym2 EXEMPTION comment, still frames the soft constraint as a
bare divergence. Everything it says is true, and the doc now carries the
"documented choice, not an omission" classification with its citations,
but a reader of the test alone gets the pre-punch framing. One clause
pointing at `R/covstruct.R:916-947` would close it. The file header
already says "See dev/brms-likelihood-tests.md", so this is a nicety.

## Updated verdict: GO

Upgraded from GO WITH FIXES.

Every punch item I raised is addressed at the cited lines, and three of
them are addressed better than I asked: P3 added a reconciliation
sentence explaining why three rows carry four divergences rather than
just picking a number; P2 kept the wrong draft on the record with a note
saying why, instead of deleting it; and D1 went past reclassification to
name the vignette's own internal contradiction. The lane also added an
honest provenance note at `dev/brms-likelihood-tests.md:607-614`
recording what the review corrected, including that the review supplied
the non-convergence measurement.

The two numbers the coordinator asked me to recompute both confirm, and
the eigenvalue count was never in dispute: 30 was my measurement and the
lane adopted it.

No claim anywhere in NEWS or the doc still says the GP length-scales are
reported on different scales, and the three code comments that mention
scale are all about the internal theta and are correct as they stand.

All four suites are green at the claimed counts, the cache is untouched
at 36, `R/` is untouched, and main is clean and unmoved at 5dfdd84.

The only thing standing between this lane and a merge is P1, which is
the coordinator's by agreement, plus R2 and R3, which are one deletion
and one paragraph in the same commit.
