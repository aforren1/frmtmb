# Review: wt-car-jacobian (exact delta method through the esicar centering)

Reviewer lane. Base `68d6782`; branch `wt-car-jacobian`; the work under
review is the uncommitted diff against that commit plus untracked files.
`main` at `b131fe1` has an `R/` tree identical to `68d6782`
(`git diff 68d6782 HEAD -- R/` is empty), so the base commit is the
correct baseline for every "unchanged" claim and the comparison library
built from it is a `main` core.

Status: COMPLETE. Verdict at the end: MERGE, with four follow-ups,
none blocking. (Written incrementally; the machine restarted at 22:35 on
2026-09-05 and this file was restarted from the beginning.)

## (a) Diff hygiene

`git -C ... status --porcelain` and `git diff --stat 68d6782`:

    M NEWS.md                        |  28 ++-
    M R/covstruct.R                  |  54 +++
    M R/methods-fit.R                |  19 +-
    M R/predict.R                    |  49 ++-
    M tests/testthat/test-car-spde.R | 104 ++++--
    M vignettes/frmtmb.Rmd           |  21 +-
    ?? dev/car-jacobian-findings.md

Exactly the six tracked files claimed plus the one untracked findings
document. No `NAMESPACE`, no `man/`, no `dev/*.log`, no `docs/`, no
stray scratch files. Clean.

## (b) The mathematics

### The Jacobian is P, derived independently

`car_center()` (`R/covstruct.R:1309`) is
`c_i = b_i - (1/n_{j(i)}) sum_{k in comp j(i)} b_k`, i.e. `c = P b` with
`P = I - sum_j (1/n_j) s_j s_j'` and `s_j` the 0/1 indicator of
component `j`. The map is LINEAR in `b`, so `dc/db = P` with no
higher-order term, and the second derivative of `c` in `b` is exactly
zero (this matters below: it makes the Hessian argument family-free).

`car_center_jacobian()` (`R/covstruct.R:1332`) emits, per component of
size `k`, `as.vector(diag(1, k) - 1/k)` at `i = rep(mem, times = k)`,
`j = rep(mem, each = k)`. Column-major `as.vector` puts entry `(r, c)`
at index `(c-1)k + r`, so the triplet is `(mem[r], mem[c], d_rc - 1/k)`.
Correct, and orientation-safe anyway since `P` is symmetric.

Singletons: `k == 1L` is `next`, so no entries. That is right, not a
shortcut. For a singleton component `{i}`, `P` row `i` is
`e_i - (1/1) e_i = 0` and column `i` likewise, so the row and column are
structurally empty; `c_i = b_i - b_i` is the constant zero, and a
constant's conditional SD is zero, not "small".

Verified numerically rather than by reading: on the singleton graph
(components `{1..4}`, `{5..8}`, `{9}`), `car_center_jacobian()` densified
equals a `P` I built here from the memberships to `max|diff| = 0`
(bitwise), and the `c_idx x b_idx` block of `rr_jacobians(fit)$Jb`
equals it to `max|diff| = 0`, at every `con_sd` tried.

### The closed form, and the structure of V it needs

The claim `sqrt(diag(V) - con_sd^2) == sqrt(diag(P V P'))`
(`car_center_condsd()`, `R/covstruct.R:1364`) needs exactly this: the
inert coordinates `m_j = s_j' b / n_j` must be uncorrelated with
everything else in `V` and have variance exactly `con_sd^2`. Then
`V = cov(f) + con_sd^2 sum_j s_j s_j'` on the orthogonal split
`b = f + sum_j m_j s_j`, `P` annihilates the second term, `P f = f`, and
the diagonals differ by `con_sd^2` at every level of every component,
`n_j`-free. A singleton has `cov(f)_ii = 0`, so it lands on exactly
zero.

The code guarantees that structure, and more strongly than the lane
claims. From the esicar `nll` (`R/covstruct.R:1424`) the block penalty is
`0.5 (tau b'Lb + sum_j kappa0_j (s_j'b)^2)` with
`kappa0_j = 1/(con_sd n_j)^2`, and `L s_j = 0` because `car_components()`
(`:1201`) is an undirected BFS over the pattern and `car()` refuses an
asymmetric `M` (`:1050`), so `s_j` is in the Laplacian's null space per
component by construction. The data enters only through `c = P b`
(`expand_b()`, `:1836`), and `P s_j = 0`. Hence for the FULL joint
Hessian `H` at ANY parameter value,

    H s_j = tau L s_j + P0 s_j + P'(d2 data/dc2)P s_j
          = 0 + kappa0_j n_j s_j + 0

so `s_j` is an exact eigenvector with eigenvalue `1/(con_sd^2 n_j)`, and
`var(m_j) = s_j' H^-1 s_j / n_j^2 = con_sd^2` exactly. Three consequences
the lane understates:

1. It is structural, not asymptotic. It does not need the optimizer to
   be at the mode, and `kappa0` is tau-free (`car_esicar_prec()`, `:1374`)
   so it does not need theta to be converged either.
2. It holds for any family, because `c` is linear in `b` so the data
   term's `b`-Hessian is `P' M P` whatever `M` is.
3. It holds for the joint covariance too, not only `cov.random`:
   `H_{m,theta} = 0` because `con_sd` is data and the `m` quadratic
   carries no tau.

Measured on the singleton graph, at `con_sd` in `{1e-1, 1e-3, 1e-4}`,
from `frm_joint_cov()` rebuilt here:

    var(m_j) - con_sd^2         <= 4.0e-19  (all j, all con_sd)
    max |cov(m_j, m_l)|, j != l <= 1.2e-20
    max |cov(m_j, f)|           <= 2.2e-18
    max |cov(m_j, beta/theta)|  <= 4.3e-19
    diag(P V P') - (diag(V) - con_sd^2) <= 6.9e-18

So the closed form is exact to solver roundoff on the joint covariance,
which is the harder case.

### Can the sqrt argument go negative?

Only by solver roundoff, and only where `cov(f)_ii` is itself zero, i.e.
at a singleton. Everywhere else `diag(V) - con_sd^2 = cov(f)_ii` is a
genuine positive variance and `con_sd^2` is subtracted from a quantity
that structurally exceeds it. The lane's worry about the optimizer
stopping off the mode does NOT apply: per the eigenvector argument
above, `var(m_j) = con_sd^2` holds at any parameter value, so the
convergence tolerance of the esicar objective (the lane's 1e-11) has no
route into this subtraction.

Measured: `min(diag.cov.random - con_sd^2)` over all nine coordinates was
`0` EXACTLY at all three `con_sd`, the minimum being the singleton, and
`identical(dcr[[9]], con_sd^2)` was `TRUE` bitwise at all three. So
`pmax(., 0)` (`:1365`) never fired in these runs and the singleton
reaches `sqrt(0) == 0`.

The guard is present and correct. One note, not a defect: `pmax()`
also silences a LARGE negative difference, which would be the signature
of a real mismatch (an `aux_car` whose `con_sd` no longer matches the
`kappa0` its precision was built from, say, after a hand-edited fit
object). Cheap to make loud; not required.

## (c) con_sd invariance, and the other CAR types

### Invariance, reproduced here

Singleton graph (`{1..4}`, `{5..8}`, `{9}`), `con_sd` in
`{1e-1, 1e-3, 1e-4}` including the default and one large. Pairwise, over
all 72 rows and all 9 levels:

    con_sd pair      max|se^2 diff|   max|se/se-1|   max|condSD^2 diff|
    1e-1 vs 1e-3        1.5e-15         3.0e-14          2.75e-15
    1e-1 vs 1e-4        3.54e-13        6.97e-12         9.05e-13
    1e-3 vs 1e-4        3.54e-13        6.96e-12         9.03e-13

For scale, the leak the fix removes is `1e-3^2 - 1e-4^2 = 9.9e-7` in the
variance, six to eight orders of magnitude above these residuals. The
residual tracks `max|dtheta|` (5.5e-13 for the 1e-4 pair, 2.4e-15 for the
1e-1 pair), so the lane's reading of it as optimizer noise rather than a
`con_sd` leak is right: it scales with how far apart the two optima land,
not with `con_sd^2`.

`con_sd = 1e-1` behaves no worse than the default, which is the useful
result: `con_sd^2 = 1e-2` is a tenth of the field variance there, so the
old code moved those standard errors by percent, not by `1e-5`.

### escar, icar, bym2 unchanged, against a base-commit core

Built a second private library from `68d6782` and fitted all four types
on the same 4x4 lattice data under both cores, comparing twelve outputs
each: `logLik`, `theta`, `beta`, `b`, in-sample `fit`/`se.fit`, newdata
`fit`/`se.fit`, `ranef`, `condSD`, `VarCorr`, and the whole
`frm_joint_cov()` matrix.

    escar : identical(whole list) = TRUE   (all 12, bitwise)
    icar  : identical(whole list) = TRUE   (all 12, bitwise)
    bym2  : identical(whole list) = TRUE   (all 12, bitwise)
    esicar: identical = FALSE, and ONLY in se, newdata se, condSD
            (max|d| 2.6e-6, 2.6e-6, 2.69e-6; everything else bitwise)

So the claim holds, and the esicar row is the positive control that the
harness can see a change at all. The change surface is exactly the three
reported quantities; point estimates, the joint covariance and `VarCorr`
are untouched, which is what "the delta method, not the model" should
mean. Relative size at the default is 1.35e-5 in the standard error,
matching the lane's number.

## (d) The rr interaction

Fitted `y ~ x + car(W, gr = loc, type = "esicar") + rr(spp + 0 | site,
d = 2)` (3x3 lattice, 10 sites, 3 species, 60 rows), so both exceptional
block types are live in one frame.

Checked `predict(se.fit = TRUE)` against `A = d eta / d(beta, b, theta)`
from `numDeriv::jacobian()` over the WHOLE linear predictor (perturbing
`fit$estimates` and re-running `predict(type = "link")`), so the
reference reuses none of `rr_jacobians()`:

    max|se - se_ref|   = 7.52e-11
    max|se/se_ref - 1| = 3.52e-10        (se in [0.202, 0.286])

That is numDeriv Richardson accuracy, not a discrepancy. The splice is
correct.

Three structural checks alongside:

- `Jb[c_idx, b_idx]` on the car block equals `P` to `max|d| = 0`; on the
  rr block it has 50 nonzeros, which is `10 sites x 5` free loadings for
  `dim = 3, rank = 2`. Neither block overwrote the other's triplets.
- `th_cols` has 5 entries, exactly the rr block's `theta_idx = 1..5`.
  The esicar block's `theta_idx = 6` is correctly absent: measured
  `max |d eta / d theta_6| = 0` EXACTLY, against 1.45 to 1.75 for the rr
  columns. Centering is tau-free, so an esicar theta column would have
  been wrong, and the lane's decision to build `th_cols` only inside the
  rr branch is right rather than lucky.
- `ranef(condVar = TRUE)` returns `NA` condSDs for the rr block and real
  ones for the esicar block in the same call, and the esicar block's
  condSDs match `sqrt(diag(P Vb P'))` rebuilt from `frm_joint_cov()` to
  7.96e-14.

## (e) newdata and conditional_effects

All three paths work. No regression: the icar control is bitwise
identical between the two cores on every one of these outputs.

1. **Regions absent from newdata.** `newdata` keeping 5 of 16 regions,
   with and without `droplevels()`: `fit` and `se.fit` match the
   corresponding in-sample rows to `max|d| = 0` (bitwise), both spellings
   agree bitwise, all se finite and positive.

2. **A region unseen in the fit.** Without `allow_new_levels` the error
   still names the level (`New levels in grouping factor 'loc': ZZZ. Use
   allow_new_levels = TRUE`). With it, prediction succeeds and all se are
   finite. The unseen rows' se (0.052 to 0.112) are NARROWER than a seen
   region's (0.19 to 0.25), which looks wrong until you read
   `lp_extra_var()` (`R/predict.R:1608`): `car` is on the deliberate skip
   list, "the levels ARE the structure there, so there is no marginal
   variance to hand an unseen one". So an unseen region gets the
   population band by design. Unchanged by this lane, and bitwise
   identical to the base core on the unseen rows; the esicar difference
   is confined to the SEEN rows (elementwise `wt - main` is
   `-2.49e-6, -2.29e-6, -2.34e-6, 0, 0, 0`).

3. **`conditional_effects(re_formula = NULL)`.** Runs, band ordered and
   finite, half-widths 0.098 to 0.345. On a car-only fit with no
   `conditions`, it equals the population band to `max|d| = 0`, which
   follows from the same skip: the grouping column goes NA, `add_b_cols()`
   is still called with the new `Jb` but every Z row for the NA level is
   zero. So that call exercises the branch without discriminating it.

   Conditioned on a real region it does discriminate, and this is the
   sharpest single check in the review:

       conditional_effects(effects = "x", re_formula = NULL,
                           conditions = data.frame(loc = "L1"))

       band se vs Z P V P' Z' rebuilt here : max|d| = 1.11e-16
       same rebuild with dc/db = I         : max|d| = 2.54e-06
       se_id^2 - se_ref^2                  = 1.000000e-06 == con_sd^2

   and `predict(newdata)` over the same grid matches the rebuild to
   `max|d| = 0`.

## (f) Tests

Run in fresh processes against the worktree core installed in a private
library, with `NOT_CRAN=true` (without it `test-predict-lp-basis.R` skips
all five tests on `skip_on_cran()`) and with the test environment
parented on the namespace (the file calls internal `car_scale_factor()`
unqualified, which errors under a plain `library()` harness).

    test-car-spde.R           23 tests  135 pass  0 fail  0 warn  0 skip
    test-predict-newdata.R     6 tests   12 pass  0 fail  0 warn  0 skip
    test-predict-lp-basis.R    5 tests   30 pass  0 fail  0 warn  0 skip
    test-message-uniqueness.R  1 test     6 pass  0 fail  0 warn  0 skip

`23/135` matches the lane's claim exactly.

### Do the new assertions discriminate?

Yes, decisively. The SAME test file run against the base-commit core
gives 11 failures (the reporter caps at 10 shown), and they fall only in
the four tests that assert new behavior:

    14  esicar constrains a disconnected graph PER COMPONENT   1 fail
    16  con_sd leaves the esicar fit AND its standard errors   6 fail
    17  the esicar delta method is Z P V P' Z', not Z V Z'     3 fail
    18  esicar handles a SINGLETON component                   1 fail

Everything else stays green, which also re-confirms (c) from the test
side. So the `1e-12` bound is NOT reachable by construction: under the
old code the same assertion fails at `9.9e-7`, and the singleton condSD
assertion fails at `2.687e-6`.

### No no-error-only assertions

Every added assertion is numeric (`expect_vector_equal`, `expect_lt`,
`expect_gt`, `expect_identical`). Nothing added is a bare
`expect_silent`/`expect_no_error`. Good.

### Singleton coverage

Present, at `test-car-spde.R:641-642`, and it is the assertion that most
directly encodes the new mathematics (`condSD == 0` exactly, and the
other eight bounded away from zero).

### Margin on the bounds (measured on the test's own configuration)

    max|se3^2 - se4^2|  = 3.779e-13   vs bound 1e-12   margin 2.6x
    max|sd3^2 - sd4^2|  = 2.866e-13   vs bound 1e-12   margin 3.5x
    max|se2^2 - se3^2|  = 2.422e-15   vs bound 1e-12   margin 413x
    max|sd2^2 - sd3^2|  = 1.318e-15   vs bound 1e-12   margin 759x
    max|se3/se4 - 1|    = 3.844e-12   vs bound 1e-9    margin 260x
    max|sd3/sd4 - 1|    = 3.648e-12   vs bound 1e-9    margin 274x

The lane's reported 3.8e-13 reproduces. The two tight ones are exactly
the pairs involving the `con_sd = 1e-4` fit, whose `kappa0` is a hundred
times stiffer and whose optimum lands 5.5e-13 away in theta. See finding
F1.

## (g) Docs

**NEWS.** Under `# frmtmb (development version)` (`NEWS.md:1`), which is
correct, and it states the behavior change in the first sentence
("standard errors no longer carry `con_sd`") rather than burying it. It
says which reported numbers move and by how much. I checked the three
magnitudes it quotes and they reproduce: 1.348e-5 relative at the 1e-3
default (I measure 1.35e-5), ~1.4e-3 at `con_sd = 0.01`, and 7.7 to 13
percent at `con_sd = 0.1` (the entry says 7.5 to 12.7). It also amends
the 0.52.0 entry in place (`NEWS.md:251`) with a pointer forward rather
than leaving a now-false claim standing. Its "`logLik()`, the estimates,
`ranef()`'s values and `VarCorr()` are bit-identical, and so is every
reported number for escar, icar and bym2" is exactly what I measured in
(c).

**Vignette.** `vignettes/frmtmb.Rmd:206-217` is accurate: the new text
puts the standard errors and conditional SDs on the invariant list and
keeps a past-tense paragraph explaining what the old numbers were. The
"7 to 13 percent at `con_sd = 0.1`" figure matches my measurement.

**The stale sentence.** `R/compat.R:1231` still ships

    con_sd changes no esicar estimate, but it does scale that type's
    prediction standard errors and ranef() conditional SDs as con_sd^2
    in the variance, so leave it at the default there.

which is now false and, unlike the 0.52.0 NEWS entry, is not tensed as
history: `frm_compat_rules()` returns it as current guidance. The lane
left it because wt-protocol owns the file this round and put the
replacement in `dev/car-jacobian-findings.md:235-238`. I checked that
replacement and it is correct as written, and I grepped `R/`, `man/`,
`vignettes/`, `inst/` and `README.md`: `R/compat.R:1231` is the ONLY
remaining stale claim. See F3.

Worth noting the neighboring rule at `R/compat.R:1240` already documents
the behavior (e2) surfaced ("Locations outside the fitted set have no
marginal variance of their own; allow_new_levels predicts them at the
population level with no block contribution"), so that is documented
design, not a gap.

## (h) Cross-lane surface

Three hunks in `R/predict.R`, all of them above base line 1506:

    @@ -629,9  +629,20   doc block of rr_jacobians()   (base 629-637)
    @@ -658,6  +669,14   the esicar branch INSIDE rr_jacobians()
    @@ -1475,19 +1494,31  doc block + head of lp_delta_A()

Functions changed, with ranges in both revisions:

    rr_jacobians()   base 637-679   worktree 648-690   (doc from 632)
    lp_delta_A()     base 1485-1562 worktree 1512-1589 (doc from 1494)

Net effect on everything below: a uniform **+31 line** offset.

**wt-spline-span's territory is untouched.** `frm_lp_basis()` is at base
3068 (worktree 3099) and `lp_basis_nl()` at base 3164 (worktree 3195);
both are far below the last hunk and neither their bodies nor their
`lp_delta_A()` call sites appear in the diff. The merge is a pure offset
shift, no textual conflict.

**No call site changed signature, and none changed at all.**
`lp_delta_A(object, lp, ed, newdata, use_re, jc, has_rr, rrj)` has the
same eight formals in both revisions. All five call sites are
byte-identical between base and worktree; they only moved:

    base 1364 -> worktree 1383
    base 1696 -> worktree 1727   (inside a per-dpar loop)
    base 2078 -> worktree 2109
    base 3100 -> worktree 3131   (inside frm_lp_basis)
    base 3198 -> worktree 3229   (inside lp_basis_nl)

What DID change is the contract, and a merging lane has to know it:
`rrj = NULL` no longer means "use the identity". `lp_delta_A()` now
derives the need for a Jacobian from `frame_needs_expand(object$frame)`
and builds one itself when the caller passed none. Any new call site
that wants the identity must not have an `rr` or `esicar` block, rather
than simply passing `NULL`.

Verified the seam the spline lane owns inherits the fix without an edit,
on an esicar fit and on an `s(x, k = 6) + car(esicar)` fit (the shape
where the two lanes actually meet):

    frm_lp_basis() sqrt(diag(A V A')) vs predict(se.fit) : max|rel| = 0
    frm_lp_basis() eta vs predict(type = "link")         : max|d|   = 0

for both, i.e. bitwise.

## Findings, ranked by severity

Nothing I found is a correctness defect in this diff. The mathematics is
right, and it is right for a stronger reason than the lane gives. What
follows is one cross-lane documentation item and three test/quality nits.

### F3 (medium, cross-lane, NOT this lane's file) `R/compat.R:1231` ships a statement this diff makes false

`frm_compat_rules()` returns, as current guidance, "con_sd ... does
scale that type's prediction standard errors and ranef() conditional SDs
as con_sd^2 in the variance, so leave it at the default there." After
this change that is false, and it directly contradicts the new NEWS
entry. Unlike the amended 0.52.0 NEWS entry, it is not tensed as
history.

Reproduction:

    grep -n "does scale that type" R/compat.R      # 1231, still present
    # then, on any esicar fit, (c) above: con_sd invariance to 1e-15

The lane handled this correctly given the ownership split: it did not
touch a file wt-protocol holds, and it wrote the replacement into
`dev/car-jacobian-findings.md:235-238`. I checked that replacement text
and it is accurate. I also confirmed by grep over `R/`, `man/`,
`vignettes/`, `inst/` and `README.md` that this is the only remaining
stale claim. Action is for the coordinator, not this lane: land the
compat.R clause with, or immediately after, this merge.

### F1 (low, test fragility) the invariance bound has a 2.6x margin, and the tight pair is the one that did not need to be there

`tests/testthat/test-car-spde.R:557` and `:565`:

    expect_lt(max(abs(se3^2 - se4^2)), 1e-12)
    expect_lt(max(abs(sd3^2 - sd4^2)), 1e-12)

Measured on the test's own data: 3.779e-13 and 2.866e-13, i.e. 2.6x and
3.5x of the bound. The residual is optimizer-path noise, not leakage:
the `con_sd = 1e-2` vs `1e-3` pair in the same test sits at 2.4e-15 and
1.3e-15 (413x and 759x margin), and I measured `1e-1` vs `1e-3` at
1.6e-15. What makes the `1e-4` pair tight is that its `kappa0` is a
hundred times stiffer, so its optimum lands 5.5e-13 away in theta and
the standard errors follow.

The assertion needs to exclude a leak of 9.9e-7. A bound of 1e-10 is
still 4 orders below that and 260x above the observed residual. As
written, a different BLAS or a different nlminb path is a plausible
flake.

Reproduction: fit `car_lattice_data(42)` at `con_sd` 1e-2/1e-3/1e-4 and
print `max(abs(se3^2 - se4^2))` and `max(abs(sd3^2 - sd4^2))`.

Suggested: loosen both to `1e-10`, or drop the `1e-4` fit from the two
absolute-variance assertions and keep it in the relative ones (which
already carry 260x and 274x margin).

### F2 (low, test fragility) `expect_identical(csd[[9]], 0)` is a bitwise assertion on a float round trip

`tests/testthat/test-car-spde.R:641`. The singleton's condSD is
`sqrt(pmax(V_99 - con_sd^2, 0))` with `V_99 = 1/(1/con_sd^2)`. It is
exactly zero on this build at all three `con_sd` I tried (`identical`
TRUE at 1e-1, 1e-3, 1e-4), because `1/(1/x) == x` holds for those `x`.
It does not hold for every double; where it fails by one ulp the
assertion sees `1.4e-11` and fails while nothing is wrong.

Note the neighboring `expect_identical(re[[9]], 0)` at `:637` is
genuinely exact (`b_i - b_i` in the AD tape) and should stay identical.
Only the condSD one is a round trip.

Suggested: `expect_lt(csd[[9]], 1e-10)`, keeping `expect_gt(min(csd[-9]),
1e-3)` beside it so the assertion still says "zero here, not there".

### F4 (nit, coverage) nothing pins esicar together with rr

The esicar branch and the rr branch are the only two consumers of
`rr_jacobians()`'s triplet accumulator, they write into the same `ii`/
`jj`/`xx`, and no test puts both blocks in one frame. I verified by hand
that the splice is correct (see (d): 3.5e-10 against a numDeriv
Jacobian, `th_cols` exactly the 5 rr thetas, the esicar theta column
exactly zero), but a future edit to either branch could break the other
with the suite still green.

Cheap to pin: one fit with `car(esicar) + rr(...)`, asserting the car
block of `Jb` equals `P`, that `length(th_cols)` equals the rr block's
`theta_idx` length, and that the esicar theta contributes no column.

### F5 (nit, performance) `rr_jacobians()` is rebuilt per linear predictor on esicar fits

The four call sites still read `rrj <- if (has_rr) rr_jacobians(object)`
(`R/predict.R:1381`, `:1712`, `:2108`, `:3130`), so on a pure-esicar fit
`rrj` arrives NULL and `lp_delta_A()` rebuilds it every call. `:1727`
sits inside a per-dpar loop, so a multivariate or distributional esicar
fit rebuilds it once per dpar.

Harmless in practice: the esicar path is sparse triplet assembly with no
finite differences, unlike the rr path. Mentioned only because the
obvious tidy (widen the call sites to `frame_needs_expand()`) would also
remove the last reason for `has_rr` to exist as a parameter.

## Not defects (checked, and they look wrong at first)

- **An unseen CAR region gets a NARROWER band than a seen one** (0.05 to
  0.11 against 0.19 to 0.25). `lp_extra_var()` (`R/predict.R:1608`) puts
  `car` on a deliberate skip list, "the levels ARE the structure there,
  so there is no marginal variance to hand an unseen one", and
  `R/compat.R:1240` documents it. Bitwise unchanged by this lane (icar
  control identical between the two cores).
- **`conditional_effects(re_formula = NULL)` equals the population band**
  on a car-only fit with no `conditions`. Same cause: the grouping column
  goes NA and every Z row for it is zero. Conditioned on a real region it
  does carry the b uncertainty, correctly (1.11e-16 against the hand
  rebuild).
- **`pmax()` in `car_center_condsd()` hides a large negative difference.**
  It cannot arise: `var(m_j) = con_sd^2` is an exact eigenvalue identity
  that holds off the mode too.

## Extra edge case: two esicar blocks in one frame

Not in the brief and not in the suite, but it is the case where the
triplet accumulator in `rr_jacobians()` has to map two DIFFERENT
projections through two different `c_idx`/`b_idx` pairs, and getting the
offsets wrong there would be silent. Fitted
`y ~ x + car(WA, gr = a, esicar) + car(WB, gr = b, esicar)` on a
5-node and a 4-node chain:

    block a (n=5)  max|Jb block - P| = 0
    block b (n=4)  max|Jb block - P| = 0
    Jb nonzeros    = 41  = 5^2 + 4^2   (exactly, no spill)
    sum(field a)   = 5.55e-17     sum(field b) = -2.78e-17
    se.fit vs a rebuild carrying BOTH projections : max|d| = 0
    the identity rebuild would be off by 2e-06 in the variance

The last line is the pleasing one: `2e-6` is `2 * con_sd^2`, one per
block, which is what the additive structure predicts.

## Counts I measured

    test-car-spde.R            23 tests   135 pass  0 fail  0 skip
    test-predict-newdata.R      6 tests    12 pass  0 fail  0 skip
    test-predict-lp-basis.R     5 tests    30 pass  0 fail  0 skip
    test-message-uniqueness.R   1 test      6 pass  0 fail  0 skip
    full suite                117 files  1161 tests  6542 pass  0 fail
                                          1 warn    91 skip

    test-car-spde.R against the BASE core: 11 fail, in tests 14/16/17/18

`23/135` and the 117-file count reproduce the lane's claims exactly, and
so does `0 fail`. My total pass count is 6542 against the lane's 6648;
the 106 difference is skips, not failures. The skipped tests are
brms/Stan-gated (`test-brms-likelihood.R` 30, `test-brms-methods.R` 8,
`test-brms-agreement.R` 2 among the 40 I listed), i.e. they need a
working Stan toolchain this reviewer environment does not have. No
`car`, `predict` or `covstruct` test skipped.

Not independently verified: the lane's `tier 33/373` and `as-cran OK`.
The as-cran surface is low risk here, since both new functions are
`@noRd` and the diff carries no `NAMESPACE` or `man/` change, which is
internally consistent.

## Verdict: MERGE

The change is correct, and correct for a stronger reason than the lane
argues. `dc/db = P` is exact, `car_center_jacobian()` reproduces `P`
bitwise on every graph I tried, a singleton's empty rows are the right
answer rather than a convenient one, and the closed form
`sqrt(diag(V) - con_sd^2)` rests on an exact eigenvector identity
(`H s_j = s_j / (con_sd^2 n_j)`) that holds at ANY parameter value and
for any family, not merely at a converged mode. The `pmax()` guard is
therefore belt-and-braces rather than load-bearing, and I measured the
sqrt argument bottoming out at exactly `0`, never below.

The change surface is exactly what it should be: `escar`, `icar` and
`bym2` are bitwise identical to a base-commit core across twelve outputs
including the whole joint covariance, and `esicar` moves in precisely
three places (`se.fit`, newdata `se.fit`, `condSD`) and nowhere else.
The rr splice is right, including the non-obvious part (no theta column
for the esicar block, verified as exactly zero). The newdata and
`conditional_effects` paths still work, and conditioned on a real region
the band matches a hand-built `Z P V P' Z'` to 1.11e-16 while the old
identity rebuild is off by exactly `con_sd^2`. The tests discriminate:
11 failures against the base core, none against the new one. The docs
are accurate and the numbers in them reproduce.

No punch list blocks the merge. Four follow-ups, in order:

1. **F3, coordinator, land with or right after this merge.**
   `R/compat.R:1231` still asserts the opposite of this change to any
   user who calls `frm_compat_rules()`. The replacement text is written
   and correct at `dev/car-jacobian-findings.md:235-238`; wt-protocol
   holds the file. This is the only remaining stale claim in the repo.
2. **F1, this lane or a follow-up.** Loosen the two `1e-12`
   absolute-variance bounds at `tests/testthat/test-car-spde.R:557` and
   `:565` to `1e-10`. Measured margin is 2.6x and 3.5x; `1e-10` keeps
   four orders of discriminating power against a `9.9e-7` leak.
3. **F2.** Replace `expect_identical(csd[[9]], 0)`
   (`tests/testthat/test-car-spde.R:641`) with `expect_lt(., 1e-10)`.
   Keep the `expect_identical(re[[9]], 0)` above it, which is genuinely
   exact.
4. **F4.** Add one `car(esicar) + rr()` test. The two branches share a
   triplet accumulator and nothing currently pins that they coexist.

F5 (rebuilding `rr_jacobians()` per dpar on esicar fits) is a nit; leave
it unless someone is already touching those call sites.

For the merge with wt-spline-span: no conflict. All three hunks are
above base line 1506, `frm_lp_basis()` and `lp_basis_nl()` are at base
3068 and 3164, `lp_delta_A()`'s signature is unchanged and no call site
was edited. Everything below shifts by +31 lines. The one thing the
other lane must know is the contract change: `rrj = NULL` no longer
means "use the identity", because `lp_delta_A()` now derives the need
from `frame_needs_expand()` and builds its own.
