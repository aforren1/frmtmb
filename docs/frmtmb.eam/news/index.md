# Changelog

## frmtmb.eam 0.5.0

[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
refuses a boundary given twice, one export for a sibling package, and
the compatibility table says refused where it used to say untested.
Requires frmtmb 0.54.0 for the exclusivity declaration.

- [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  declares `dec()` and `vint1` mutually exclusive, so a model supplying
  both is refused by name instead of fitted with the second column
  unread. Measured before the change on 120 rows: `rt | dec(u) ~ 1`,
  `rt | vint(u) ~ 1` and `rt | dec(u) + vint(1 - u) ~ 1` all gave a
  log-likelihood of -76.0486443897369, the third with the two columns
  CONTRADICTING each other. `ddm_indicator()` reads `dec` and falls back
  to `vint1`, and the declaration is that precedence written down. Needs
  frmtmb’s new `frmtmb_family(exclusive_aterms =)`.

  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)
  deliberately declares no such set and is why the rule is opt-in: there
  `dec()` and `vint1` are two data, the boundary and the condition
  index, not two spellings of one.

- Every family x addition-term cell in
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  is now decided. The refusals each family has already declared in
  `accepts_aterms` are derived through frmtmb’s
  [`compat_aterm_rules()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.html)
  rather than written out again, so `trials()`, and `vreal()` where the
  family does not read it, stop reading `untested`. Rows written by hand
  keep their own notes: the derivation defers to any pair already
  refused.

  The allow-lists themselves move to one `ddm_accepts` list that the
  five constructors and the compatibility rows both read, so the table
  cannot promise a term frame assembly refuses.

- [`wiener_lpdf()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_lpdf.md)
  is now exported: the Wiener first-passage log density with
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)’s
  parameterization and
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)’s
  tape safety.
  [`frmtmb.learn::rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.html)
  is a delta learning rule whose value difference drives the drift rate
  of this density, and until now the only route to it was
  `frmtmb.eam:::ddm_lpdf_both()`, which is a promise nobody made. This
  is that promise, made deliberately and kept to one function: the
  series truncations, the blend between them, the across-trial
  variability integrals and the CDF all stay internal. Nothing about the
  package’s own behavior changes.

## frmtmb.eam 0.4.0

[`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
gains across-trial variability with the go branch identical to
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)’s;
[`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
and
[`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
take `cens()`; every family declares the addition terms it reads.
Requires frmtmb 0.53.0 for the allow-list.

- [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  gains `variability =`, and takes Ratcliff’s `sv`, `sz` and `st` under
  the same names, links and argument \[wiener()\] takes them. The go
  branch is \[wiener()\]‘s averaged density and is not merely equivalent
  to it: evaluated at the same parameters, the two families’ log
  densities are BIT-IDENTICAL on every row, for every combination of the
  three. Making that true is why the family now has a builder that is
  handed the non-decision-time bound and the unreachable-row margin, as
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  has; at 0.3.0 it floored the decision time at a flat 1e-12 where
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  used `1e-9 * min(y)`, and no identity survives two different margins.

  The no-go branch is the part that had to be written, because it
  averages a DISTRIBUTION FUNCTION rather than a density. `sz` goes
  through Gauss-Legendre nodes, `st` shifts the deadline rather than
  cutting a range, and `sv` needs a Gauss-Hermite quadrature of its own:
  the drift enters the density only as an exponential-quadratic, which a
  normal average integrates in closed form, and it enters the
  distribution function through the eigenvalues as well, where nothing
  does.

  Verified against a 200-bit `Rmpfr` integration of the distribution
  function, against 40000-trial simulation from the process (the
  observed no-go rate is within 1.7 standard errors at every setting,
  and a ten-cell chi-square on the go response times reaches 16.0
  against a 0.999 critical value of 27.9), and against EMC2, whose
  `DDMGNG` carries the same three: the no-go probability matches
  `1 - EMC2:::pDDM()` to 3.9e-13 or better and the go density matches
  `EMC2:::dDDM()` to 1.6e-15. EMC2’s `st0` is a uniform on
  `[t0, t0 + st0]` where this family’s `st` is centred on `ndt`,
  following brms; the comparison applies that shift and
  [`?wiener_gng`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  records it.

- `wiener_gng(nogo_nodes =)` is a second node-count argument, because
  the two branches are two integrands. The density’s `st` count is 21
  because its range is cut at the response time and its integrand turns
  on sharply at the cut; the probability has no cut and saturates at 7.
  That would be a curiosity if the two cost the same, and they do not:
  the probability’s rule is a three-dimensional product, so a count
  carried over from the density is multiplied by every other count in
  the grid. With the density’s own counts a 500-row three-variability
  model exhausted memory outright.

- **`sv` is weakly identified in a go/no-go design, and
  [`?wiener_gng`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  now says so.** On 3000 simulated trials with a true `sv` of 0.6 this
  family returns 0.001 while
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
  on the same generative parameters and seeing both boundaries, returns
  0.612. It is not an optimizer failure: the fit reaches a HIGHER log
  likelihood at `sv` near zero by trading it against the drift and the
  boundary separation. Read a small fitted `sv` here as “the data did
  not pin it”.

- A start point that `sz` pushes past a boundary is now the boundary
  case it is, rather than a `NaN`.
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  documents that a wide `sz` at a biased start can leave the boundaries
  and calls the density there a barrier; that is true of the density,
  which stays finite, and it was NOT true of the no-go probability,
  which took `log1p(-w)` of a negative above one and returned a
  probability ABOVE one below zero. A `NaN` is not a barrier, it is the
  end of the tape, and one node of one row took the whole fit with it.

  The clamped value is the correct limit for THE BRANCH IT IS APPLIED
  TO. It does not repair the pair, and the documentation no longer says
  it does: the go branch is left unclamped so that the bit-identity with
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  survives, so past a boundary the two branches average different
  start-point distributions and the go mass plus the no-go probability
  falls short of one - measured, 0.970639 at `sz` = 0.5 with `bias` =
  0.85, and 0.843370 at `sz` = 0.9 with `bias` = 0.90, so up to 16
  percent of the mass. It is tolerable because the region is strictly
  downhill: profiled at `bias` = 0.85, the best point inside the
  boundaries beats the best point outside by 189 log units and the
  surface is monotone across the crossing, so the clamp is a barrier an
  optimizer walks away from rather than a corner it can be pulled into.

  Inside the boundaries, where the model is defined, the two branches
  still sum to one - to 8.9e-16 in the plain family, and otherwise to
  whatever the quadrature gives, which is 3.9e-14 at `sv` = 0.6 and
  8.1e-05 at `sv` = 2.0 with default nodes.
  [`?wiener_gng`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  carries both tables.

- [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  declares an `lccdf` and an `lcdf`, so `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) both work through core
  0.52.0’s slot. Neither needed new algebra: a race is unfinished
  exactly when every accumulator is, so the log survivor is the sum of
  the same per-accumulator survivals the density already forms for the
  losers of an observed trial, over all `n` instead of `n - 1`. Checked
  against the likelihood written out by hand at the fitted parameters:
  2.2e-15 relative on a right-censored data set, 5.0e-16 with all four
  censoring codes, and 9.4e-16 on a left-truncated fit. A censored row
  still needs a `vint()` winner, which the likelihood does not read, and
  does not read EXACTLY: moving every censored row’s winner changes the
  log likelihood by zero.

- [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  declares an `lccdf` and, deliberately, no `lcdf`. Right censoring is
  the same statement the family already makes, so it is exact rather
  than close: the same rows scored as no-go trials at the deadline and
  as trials right-censored at the deadline give log likelihoods that
  differ by no bits at all. Left censoring, interval censoring and
  [`trunc()`](https://rdrr.io/r/base/Round.html) are refused by name,
  because this likelihood is a defective density plus a point mass and a
  window normalizer on the response scale would renormalize the density
  while saying nothing about the mass.

- [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  and
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)
  declare the decision indicator through core’s any-of `required_aterms`
  instead of checking for it by hand. Both read the boundary from
  `dec()` or from `vint1`, which core 0.51.0 spells
  `list(c("dec", "vint1"))`, so frame assembly refuses a model that
  supplies neither BEFORE the frame is built rather than after. The
  refusal a user sees changes: it names the term values `dec` and
  `vint1` rather than the spellings `dec(decision)` and `vint(upper)`,
  and writes `rt | dec(<column>) ~ ...` as the example. Four pinned
  expectations changed with it. **This retires the 0.2.0 note below that
  “one hand-rolled check remains, and is not `required_aterms`’s
  fault”** - core grew the seam, and the check is gone.
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)’s
  CONDITION index is still checked by hand, and still cannot be
  declared: which slot carries it moves with the boundary’s spelling, so
  the requirement is a disjunction of conjunctions and no declaration
  says that.

- `ddm_cdf_ks`, the half-width of the no-go distribution function’s
  image sum, is 4 rather than 12. The blend gives the small-time route a
  non-zero weight only below `u = 0.197`, where `tanh` has not yet
  saturated, and at that `u` the `|j| = 2` term is already 9e-23.
  Measured over 840 rows spanning `t` in 0.05 to 15, every truncation
  from 2 to 12 gives bit-identical values AND bit-identical gradients.
  It is a pure cost change and it is worth 2.3x on a go/no-go
  variability fit, which evaluates that function on a three-dimensional
  node grid. `ddm_cdf_kl` is not reducible the same way and is
  unchanged.

- [`wiener_gng_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng_simulate.md)
  gains `sv`, `sz` and `st`, and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) on a fitted
  go/no-go model follows the variability, by drawing each trial’s
  parameters before it runs the process. The rejection is on the
  OUTCOME, so it reweights all three at once where
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  has to reweight the drift by the boundary probability it implies.

- Every family declares `frmtmb_family(accepts_aterms = )`, the core’s
  new addition-term allow-list, and the hand-written `dec()` refusal the
  two race families shared is deleted.
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)
  and
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  take `vint()` and [`weights()`](https://rdrr.io/r/stats/weights.html),
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  adds `dec()`,
  [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  takes `dec()` and `vreal()`, and
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)
  takes all four. A term outside a family’s list is refused by name at
  frame assembly.

- BEHAVIOR CHANGE:
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  refuses `vreal()`. The compatibility table recorded that pair as
  working on the ground that a model supplying one fitted and gave the
  same answer as one that did not, which is the defect rather than the
  feature: the column travelled into the fit and changed nothing. The
  row now reads refused.

- The `dec()` refusal for
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)
  and
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  is now the core’s generic sentence, which names the term and lists
  what the family does take, in place of the family-specific one. Same
  refusal, same families, one seam instead of a check per package.

## frmtmb.eam 0.3.0

- RENAMED from frmtmb.ddm. The package holds the linear ballistic
  accumulator, the racing diffusion model and the go/no-go diffusion
  beside the Wiener and generalized drift-diffusion families, so its
  name now says evidence accumulation models. Family constructors and
  the `ddm_` helper names are unchanged;
  [`library(frmtmb.eam)`](https://aforren1.github.io/frmtmb/frmtmb.eam)
  replaces [`library(frmtmb.ddm)`](https://rdrr.io/r/base/library.html),
  and the site is at frmtmb/frmtmb.eam/ with the old address
  redirecting. EMC2 leaves Suggests: no shipped code uses it.

- New `tests/testthat/test-gddm-reference.R`:
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)’s
  generalized components are now checked against PyDDM, the reference
  implementation of Shinn, Lam and Murray (2020), rather than only
  against this package. Nothing outside the solver checked them before:
  the analytic Wiener comparison covers a constant drift and fixed
  bounds and no more, the gradient tests show the derivative matches the
  value it differentiates, and recovery is circular because
  [`gddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_simulate.md)
  draws from the solver’s own density on purpose. Ten cases, covering
  leak at both signs, exponential and linear collapse, the coherence
  nonlinearity at two coherences, start-point variability and a lapse,
  are frozen into a fixture under `tests/testthat/fixtures/`; no Python
  runs at test time, and `dev/gddm-pyddm-reference.py` in the source
  repository regenerates it. Each case compares the density at both
  boundaries, the two boundary masses and the log-likelihood of a small
  fixed dataset, and the test also requires the disagreement to shrink
  as the grid is refined, so the reference is what the solver converges
  to and not merely something it lands near. One case is also checked
  against an Euler-Maruyama simulation of the equation itself, so that
  two grid solvers cannot be wrong together. The tolerances are
  measured, not guessed, and
  [`vignette("gddm")`](https://aforren1.github.io/frmtmb/frmtmb.eam/articles/gddm.md)
  records what the reference covers, what it does not, and the one case
  where the two disagree.

- [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
  now works for
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md).
  The family installs its density and its simulator in
  `family_finalize()`, so before frame assembly the family object
  carries neither, and frmtmb’s
  [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
  read the family as written and refused the model for having no
  simulator. [`simulate()`](https://rdrr.io/r/stats/simulate.html) on a
  fitted
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)
  was unaffected, because a fit carries the finalized family. The fix is
  in frmtmb rather than here: nothing in this package’s wiring was
  wrong, and
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  and
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
  which pass `sim` in the constructor, were never affected.

- New `tests/testthat/test-simulate-density.R`: the agreement tier for
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)
  and
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
  matching the one frmtmb holds its built-in families to. The other
  tests here check the density against RWiener and the exported
  [`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ddm_simulate.md)
  against its generative process; neither of those is the seam frmtmb
  uses, which is the family’s `sim` slot. The check is shaped for a
  choice-RT model: the density is defective, so the reference is that
  density renormalized on the row’s own boundary, and the two
  boundaries’ masses are asserted to sum to one.

### Two more families

- New
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md),
  the racing diffusion model of Tillman, Van Zandt and Logan (2020): `n`
  independent Wiener accumulators with their own positive drifts, a
  start point uniform on `(0, A)` and a common threshold `A + k`, so
  each finishing time is an inverse Gaussian averaged over the start
  point. It is
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)’s
  geometry with the ballistic assumption removed, and it takes
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)’s
  spelling for that reason: a model moves between the two by changing
  one word. The drifts have a LOG link rather than an identity one,
  which is the one place the two part company, because a
  racing-diffusion drift is the rate itself and an accumulator with a
  rate of zero never finishes.
  [`rdm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm_simulate.md)
  draws from the generative process.
- New
  [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md),
  the go/no-go diffusion of Gomez, Ratcliff and Perea (2007), which EMC2
  calls `DDMGNG`: a two-boundary diffusion where only the upper boundary
  produces an observable response. A go trial contributes the ordinary
  upper-boundary density; a no-go trial contributes the probability of
  no upper crossing before a deadline. The deadline goes on the family
  when every trial shares it and through `vreal()` when it does not, and
  `dec()` says which trials produced a response.
  [`wiener_gng_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng_simulate.md)
  draws from the generative process. Across-trial variability is
  deliberately not offered;
  [`?wiener_gng`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  says why.
- The Wiener defective distribution function, which this package did not
  have.
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  declares no `lcdf` and its compatibility table says so, but the
  go/no-go family’s no-go branch needs one, so it is written in
  `R/wiener-cdf.R` as two series blended in `log(u)` the way the
  density’s two are. The blend is centred at `u = 0.02` rather than the
  density’s 0.35, and that is the substantive choice: the small-time
  route reaches the no-go probability as `1 - F_upper`, which cancels
  exactly where a no-go trial is surprising, and the large-time route
  computes it directly and never subtracts. Handing over as early as the
  large-time route is accurate holds 2.2e-12 relative against a 260-bit
  reference over 1200 points, with one row worse than 1e-12 and none
  worse than 1e-9.

### Measured against EMC2

- Both families agree with EMC2 where EMC2 has digits, and the
  comparison lives in `dev/rdm-gng-emc2-reference.R` rather than in the
  suite, because every EMC2 function that computes either likelihood is
  internal. Composed exactly as `EMC2:::log_likelihood_ddmgng` does, the
  go/no-go log likelihood agrees to 6.1e-16 per go row and 1.3e-15 per
  no-go row. Composed as `EMC2:::log_likelihood_race` does, the
  racing-diffusion log likelihood agrees to 1.1e-12 over 200 rows at
  two, three and four accumulators.
- Where they disagree, the disagreement is adjudicated rather than
  asserted. EMC2 writes a race loser’s survival as `1 - pWald(...)`,
  which returns EXACTLY ZERO on 35 of 315 grid rows, the first where the
  true survival is near 1e-13;
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  writes it out directly and returns no zeros anywhere. Against
  `statmod`, which is party to neither, the single-accumulator density
  is 4.6e-13 from the truth for this package and 2.1e-04 for EMC2 on the
  rows where they differ. Similarly `WienR`, which is EMC2’s own
  distribution function, is 4.4 percent wrong at a no-go probability of
  1.19e-13 at its default precision and at every setting down to 1e-12,
  and 0.033 percent wrong at its tightest; EMC2 calls it with
  `precision = 0.005`, looser than any of those, so EMC2 sits at the 4.4
  percent end. This family’s large-time route is 3.0e-15 there.
- The suite itself uses only EXPORTED references: `statmod`’s inverse
  Gaussian averaged over the start point for
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md),
  and `WienR` and `RWiener` for the go/no-go distribution function. New
  in `Suggests`: `EMC2`, `statmod` and `WienR`.

### Fixed

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  `predict(type = "response")` and `residuals(type = "response")` now
  refuse on both new families instead of returning a drift rate. A
  family that declares no `post$mean_fn` gets frmtmb’s fallback, “the
  first primary dpar on the response scale is the mean”, and for these
  two that number is a drift: `rdm(3)` returned 3.51 for data whose
  response times average 0.36, and
  [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  returned a constant 1.05 for data whose go response times average 0.6,
  with nothing to signal it. Both now declare a mean that stops with a
  reason. For
  [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  the refusal is also the right answer on the merits: a go/no-go trial
  produces a pair, and the no-go rows have no response time to average.
- [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  refuses `dec()` by name rather than ignoring it. Measured on
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
  the sibling race family, a `dec()` term supplied alongside `vint()` is
  silently dropped and the model fits; a model ported over from
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  would otherwise fit while quietly meaning something else.

### Also

- `tests/testthat/test-simulate-density.R` covers the two new families
  on the same terms as the other three: each family’s `sim` slot against
  that family’s own log density, and the defective masses asserted to
  sum to one. For
  [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  the go branch’s support ends at the deadline rather than running to
  infinity, and a no-go row’s draw is its deadline.
- The compatibility rows for both families were RUN rather than reasoned
  about, and several first guesses were wrong: `REML` and `quadrature`
  both work where they had been written down as refused or untested,
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
  refuses
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  outright because its components need a dpar called `mu`, and neither
  family makes its addition terms mandatory on newdata for a link-scale
  prediction.
- [`vignette("ddm")`](https://aforren1.github.io/frmtmb/frmtmb.eam/articles/ddm.md)
  gains a section per family, each with a fitted example, a recovery
  check and the comparison against EMC2.

## frmtmb.eam 0.2.0

Three families where there was one: Ratcliff’s full diffusion model as
an extension of
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
the generalized drift-diffusion model
[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md),
and the linear ballistic accumulator
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
plus the queued defects. The three share one floor idiom and one
compatibility table.

### The full diffusion model

- `wiener(variability = )` adds across-trial variability to the family
  rather than forking it. Naming any of `"sv"` (drift rate), `"sz"`
  (start point) and `"st"` (non-decision time) turns that one into an
  ordinary distributional parameter, with its own link and its own
  formula.
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  with no arguments is the model it always was: on a grid of 11520
  parameter combinations its log density is byte-for-byte identical to
  the previous release everywhere the density is defined, and its tests
  are unchanged. The one difference is below the non-decision time,
  where it now returns `-Inf` instead of `NaN`, which is the deliberate
  fix below.
- The likelihood is the analytic Wiener density averaged over those
  distributions, and the three are done three different ways. The drift
  integral is Gaussian against an exponential-quadratic and is evaluated
  in CLOSED FORM: it agrees with adaptive quadrature of the same thing
  to better than 1e-13 relative and takes no nodes, so estimating `sv`
  is FREE relative to the plain density. Measured at the fit level over
  three seeds, `sv` costs 0.97 to 1.07 times plain, and at the density
  level both are 0.0002 s per call. The start-point and
  non-decision-time integrals are uniform and use fixed-node
  Gauss-Legendre quadrature, and those do cost: 0.0056 s per call for
  `sz` at 7 nodes and 0.0170 s for `st` at 21.
- Node counts are the `nodes` argument and the defaults are measured:
  `sz` reaches machine precision at 7 nodes, `st` reaches 1e-9 at 21.
  They differ because the non-decision-time range is cut by the response
  time on a fast trial, and the integrand turns on sharply at the cut.
  Node positions and counts are fixed when the family object is built,
  because an automatic-differentiation tape cannot record a branch on a
  parameter; a parameter only rescales the interval they map onto.
- Variability parameters at zero reproduce the plain Wiener density to
  better than 1e-13 in the log density, which is floating-point rounding
  on a differently associated sum rather than a quadrature error.
- Recovery on simulated data, 8 replicates of 1500 trials with a normal
  drift rate and a uniform non-decision time: every parameter’s Monte
  Carlo mean is within 2.2 of its own standard errors of the value it
  was generated from.
- Across-trial variability is NOT frmtmb’s `quadrature = TRUE`. That
  marginalizes random effects by Gauss-Kronrod, is wired to the
  random-effect coefficient vector by name, and refuses a model with no
  random-effect block. This integral shares nothing between trials, has
  no level to estimate, and exists in models with no grouping factor at
  all, so it lives inside the density.
- [`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ddm_simulate.md)
  takes `sv`, `sz` and `st`, drawing each trial’s own parameters and
  then running the ordinary process, so the simulator states the model
  independently of the density.
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) follow the
  variability rather than ignoring it. Conditioning on the boundary a
  row ended at reweights which per-trial parameters that row could have
  had, so the fitted mean is a ratio of two quadratures and the
  simulator accepts a drawn drift rate and start point with the boundary
  probability it implies. The plain closed form would not do: at an
  unbiased start point it returns the same mean for both boundaries,
  where 40000 simulated trials put them 0.06 s apart. Both are checked
  against simulated data, which knows nothing about how either is
  computed.

### The generalized drift-diffusion model

- [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm.md)
  is the generalized drift-diffusion family of Shinn, Lam and Murray
  (2020): a drift that may depend on the accumulator’s own level and on
  a covariate, boundaries that may collapse within a trial, and a
  starting distribution that may be a point or an interval. There is no
  closed-form first-passage density, so every likelihood evaluation
  solves the Fokker-Planck equation forward in time and reads the
  probability flux through each boundary. It needs no change to core
  frmtmb.
- The components are chosen by argument and are extensible.
  [`gddm_drift_constant()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-drift.md),
  [`gddm_drift_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-drift.md)
  and
  [`gddm_drift_leak()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-drift.md)
  are summed to make a drift;
  [`gddm_bound_constant()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-bound.md),
  [`gddm_bound_exponential()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-bound.md)
  and
  [`gddm_bound_linear()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-bound.md)
  give the boundary;
  [`gddm_start_point()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-start.md)
  and
  [`gddm_start_uniform()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm-start.md)
  give the start.
  [`gddm_drift_term()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_drift_term.md),
  [`gddm_bound_term()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_bound_term.md)
  and
  [`gddm_start_term()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_start_term.md)
  are the documented seams for writing more. Every free quantity is a
  dpar that takes a formula, as in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md).
- `bs` is the boundary SEPARATION and `bias` the relative start point,
  both as in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
  so estimates are directly comparable between the analytic family and
  the generalized one.
- The substitution `y = x / B(t)` pins the moving boundaries at plus and
  minus one, which keeps the grid fixed while the boundary collapses and
  is what makes the likelihood differentiable: nothing on the taped path
  branches on a parameter. With the walls stationary the scheme is
  Crank-Nicolson, which a solver that chases a moving bound cannot use.
- The likelihood is an ordinary rowwise family, not a
  [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.html).
  frmtmb calls `lpdf` once per objective evaluation with full-length
  vectors, and the condition a trial belongs to is data, so a rowwise
  density does one solve per condition, which is all a structure would
  have bought, and it keeps
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) rather than
  defaulting them to refused. Measured: 6 solves for 2400 trials over 6
  conditions.
- [`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_control.md)
  carries the grid and what is done with the answer: `dt`, `ny`,
  `t_max`, `max_ndt`, `renormalize` and `tridiagonal`. Renormalizing the
  defective density is on by default and should stay on: the discretized
  solve loses mass in a parameter-dependent way, so a likelihood that
  does not divide it out rewards fast absorption. On data simulated from
  the model, turning it off more than doubles the fitted leak and
  shrinks the boundary separation by a fifth.
- `tridiagonal` picks how the solve inside each step reaches the tape.
  `"recorded"`, the default, builds a large tape that runs in compiled
  code; `"atomic"` collapses the solve into one node with a hand-written
  adjoint, building about twelve times faster and evaluating about
  twelve times slower. Both give the same derivative to machine
  precision.
- The published coherence nonlinearity has no derivative at zero
  coherence, which a motion design normally contains. The coherence is
  data, so the zero condition is resolved once when the tape is built
  and never reaches it; the gradient in the exponent is finite, and is
  exactly zero there, because the drift is zero whatever the exponent
  is. Signed coherences are supported for stimulus coding.
- The family admits exactly two responses, because one accumulator
  between two absorbing boundaries has two walls. A decision indicator
  with more than two levels is refused at frame assembly, naming how
  many levels the data has and pointing at
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
  which fits the racing accumulators that more than two alternatives
  need, rather than being folded into one of the two.
- The boundary is read from `dec()` when it is there and from `vint()`
  otherwise, so the spelling brms uses works on this family as it does
  on
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md).
  `vint()` numbers its values positionally, so the condition index is
  the first `vint()` value alongside `dec()` and the second inside
  `vint(upper, cond)`; the family reads whichever it is. Neither can be
  declared through `required_aterms`, which names the terms a density
  needs ALL of, so both refusals are written out.
- [`gddm_conditions()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_conditions.md)
  builds the condition index.
  [`gddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_simulate.md)
  draws choices and response times from the model’s own solved density.
- Validated in the package’s own suite: against this package’s Wiener
  density with a constant drift and fixed bounds, to better than 0.01 in
  the log density at the shipped grid over decision times from 0.2 s on,
  degrading on a coarser grid and improving on a finer one; automatic
  gradients against numDeriv at parameter points including a collapsing
  boundary; parameter recovery with Monte Carlo standard errors; the
  size of the renormalization bias; and the zero-coherence gradient.
- The density is floored before it is logged, as the Wiener density in
  this package already was. Where the solved density at a trial’s own
  response time underflows, the log density is a large finite negative
  number instead of `NaN`: the optimizer gets a value it can use, and a
  mixture’s log-sum-exp is not poisoned by one component. The floored
  row is flat, so its gradient is exactly zero; `-Inf` would not do,
  because `-Inf` differentiates to `NaN`.
- `gddm_floored(fit)` is where that goes to be read. It returns the
  number of rows answered by the floor rather than by the solver at the
  fitted parameters, with the row indices in the `"rows"` attribute.
  Zero is the ordinary case and means the grid represented every
  observation. A few rows means a few trials sit within a few time steps
  of the fitted non-decision time, where a fixed grid cannot resolve a
  density climbing through orders of magnitude. Many rows means the fit
  is not to be trusted: shrink `dt`, or add a lapse component. The count
  replaces what used to surface as repeated optimizer warnings about
  `NaN` function evaluations, which a user could not act on and which
  masked real warnings in a test run.
- Known limits, measured and stated rather than hidden: the density at
  decision times of only a few time steps is far larger than the truth,
  because an implicit scheme spreads a little mass everywhere at once
  where the true density is exponentially small. A lapse component
  (`gddm(lapse = "uniform")`) floors it in the model rather than in the
  arithmetic. Cost scales with the number of conditions, so a design
  with many distinct parameter settings is where this becomes painful.
- [`vignette("gddm")`](https://aforren1.github.io/frmtmb/frmtmb.eam/articles/gddm.md)
  fits one of the paper’s models end to end and says plainly where this
  is slower than
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  and why someone would pay that.
- That vignette now carries four figures, under a `tinyplot` gate: the
  three shipped boundaries over a trial, the fit against the data it was
  fitted to as defective cumulative distributions, the solved density
  against the analytic Wiener density with the log-density error on two
  grids beside it, and the leading edge of the fitted density at three
  time steps, which is what
  [`gddm_floored()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/gddm_floored.md)
  counts. They are drawn from the vignette’s own simulated data and its
  one fit, and together they cost about two seconds of a knit that the
  fit dominates.

### The linear ballistic accumulator

- `lba(n)` is a race between `n` accumulators, each rising in a straight
  line from a uniform start point on `(0, A)` to a common threshold at a
  normally distributed rate. It is the family for choices with more than
  two alternatives: a diffusion between two absorbing boundaries admits
  exactly two responses, and no reparameterization of
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  reaches a third. The likelihood is closed form for any `n`.
- Each accumulator’s drift mean is its own distributional parameter
  (`v1`, …, `vn`, identity link) and all of them are primary, so the
  main formula reaches every drift and each gets its own coefficients.
  Give one its own formula and a covariate moves that alternative alone,
  which is the capability neither two-choice family can offer.
- The remaining parameters are `A` (start-point range, log), `k`
  (threshold above the start-point range, log) and `ndt` (non-decision
  time, a logit scaled onto `(0, max_ndt)` as in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)).
  The threshold is `A + k` rather than a parameter of its own, so
  `b > A` is structural: a threshold inside the start-point range, where
  trials would begin already finished, is not a state the optimizer can
  reach.
- The drift standard deviation is fixed, because the model is identified
  only up to a common rescaling of `A`, the threshold, the drifts and
  the drift standard deviation. It is the family argument `sd_v` rather
  than a hidden default, and is carried on the family object.
- Drift rates are truncated at zero by default, following `rtdists` and
  its `posdrift = TRUE`. Every accumulator then arrives eventually and
  the choice probabilities sum to one. `posdrift = FALSE` gives the
  untruncated convention, under which the response distribution is
  defective; the two are different models, not a rescaling of each
  other, which matters when comparing against another package.
- The single-accumulator density agrees with
  [`rtdists::dlba_norm()`](https://rdrr.io/pkg/rtdists/man/single-LBA.html)
  to better than 1e-11 relative wherever `rtdists` is itself accurate,
  and the race with
  [`rtdists::n1PDF()`](https://rdrr.io/pkg/rtdists/man/LBA-race.html)
  for two, three and four accumulators. In the fast tail the two diverge
  by design and this one is the better: `rtdists` writes the normal
  difference `Phi(g) - Phi(h)` as a subtraction of two lower tails,
  which returns exactly zero once `pnorm(h)` saturates. Written in log
  space instead, as `exp(la) * -expm1(lb - la)` on the two upper-tail
  logs, nothing saturates and no comparison is needed, which matters
  because RTMB refuses comparison on AD types. Against a 200-bit Rmpfr
  reference on 1144 points, the subtractive form returns exactly zero on
  68 of the 80 tail rows and is already 2.7e-3 wrong on 8 rows outside
  the tail; the log-space form holds 5.0e-14 in the bulk and 3.6e-4 in
  the tail. Four rows still exceed 1e-6, so it is an improvement, not a
  proof.
- That mattered beyond accuracy. The subtractive form errs in the value
  only, and a tape differentiates the function that was written, so
  value and gradient described different surfaces near the
  non-decision-time bound. It did not move the point estimate, which the
  likelihood keeps out of that region, but Hessian-based standard errors
  were up to 16.5 percent off before the change.
- The survival function is written out directly rather than as
  `1 - plba_norm()`, which returns exactly zero, and so a
  log-contribution of `-Inf`, while the true survival is still around
  1e-19; the direct form agrees with a quadrature of the density down to
  survivals of 1e-23. Under `posdrift = FALSE` the never-arriving mass
  it adds back is spelled `Phi(-v/s)` rather than `1 - Phi(v/s)`, which
  is the same saturation trap one line further on.
- The choice reaches the density through `vint()`, as the Wiener
  family’s boundary indicator does, and is declared with
  `required_aterms`, so omitting it is refused by name. A choice outside
  `1..n`, a non-positive response time and a non-decision-time bound
  above the fastest response are each refused with their own message.
- The non-decision-time bound is derived from the response through
  `family_finalize()`, so the family object carries the link it will
  actually use rather than an environment filled in later.
- [`lba_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba_simulate.md)
  draws choice and time jointly from the generative process.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) on a fitted
  object redraws times conditional on each row’s observed choice, since
  the choice is data.
- Deliberate omissions: no `lcdf`, so `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) are refused; no
  `post$mean_fn`, because the mean of the race has no closed form and a
  quadrature per row was not worth writing for
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

### Defects fixed

- **The density returns `-Inf` below the non-decision time, not `NaN`.**
  A likelihood of zero is the right answer there; `NaN` is not an answer
  and propagates through every component of a mixture. The normalized
  time is held at the smallest positive double, which is bit-for-bit
  inert everywhere the density was already right.
- **`wiener(max_ndt = )` above the smallest response time is usable in a
  mixture.** The refusal is correct for the family on its own and wrong
  inside a mixture, where the contaminant component is exactly what
  covers the trials the diffusion cannot produce, so
  `allow_unreachable = TRUE` lifts it and the refusal now names it. A
  trial below the non-decision time then gets a log density that
  exponentiates to zero AND differentiates to zero, which a true `-Inf`
  does not: `-Inf` yields a `NaN` gradient and stops the fit.

### What frmtmb 0.49.0 let this package delete

- **`dec()` is the spelling now.** frmtmb gained
  [`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.html),
  this package registers `dec` when it loads, and
  `rt | dec(response) ~ x` works and takes a factor, a character vector
  or a logical the way brms does. `vint()` carries the same thing as a
  0/1 integer and is unchanged.
- **The environment the link closures read is gone.** The bound on the
  non-decision time is a property of the response, and the family object
  is built before
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) has
  any; this package used to have `valid_y()` write the bound into an
  environment, which worked only for as long as an undocumented slot
  order held. `family_finalize()` is the documented slot for it and the
  family now derives itself from the data there.
- One hand-rolled check remains, and is not `required_aterms`’s fault:
  that argument names the terms a density needs ALL of, and this family
  needs EITHER `dec()` or `vint()`. See `dev-findings.md`.

## frmtmb.eam 0.1.0

First release. A Wiener first-passage time family for two-choice
response times, written entirely against frmtmb’s exported extension
API.

- [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  is the drift-diffusion family, with brms’s dpar names and links: `mu`
  (drift rate, identity), `bs` (boundary separation, log), `ndt`
  (non-decision time, bounded) and `bias` (relative start point, logit).
  A model written for
  [`brms::wiener()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  reads the same here, and the parameterization is pinned against brms
  in the test suite.
- The density is the Navarro and Fuss (2009) pair of series. Both are
  evaluated at a fixed truncation and their logs blended with a logistic
  weight in the normalized time, because an AD tape cannot choose
  between them on a parameter. It agrees with
  [`RWiener::dwiener()`](https://rdrr.io/pkg/RWiener/man/wienerdist.html)
  to better than 1e-11 relative over normalized times from below 1e-3 to
  50, where a fixed truncation of the small-time series alone is wrong
  by tens of percent past about 8.
- The decision indicator reaches the density through `vint()`, because
  frmtmb’s addition terms are a closed set and brms’s `dec()` is not one
  of them. Omitting it is refused with a message naming both spellings,
  rather than silently producing a log likelihood over no rows.
- `ndt` uses a logit scaled onto `(0, max_ndt)`, so the support
  constraint is structural rather than something the optimizer has to
  discover. `max_ndt` defaults to the smallest response time in the
  data, found at frame assembly.
- `post$mean_fn` gives the mean response time in closed form,
  conditional on the boundary the row ended at, with the zero-drift
  limits handled explicitly. `sim` draws conditionally by inverse
  transform through
  [`RWiener::qwiener()`](https://rdrr.io/pkg/RWiener/man/wienerdist.html),
  with a discretized forward simulation as the fallback when RWiener is
  absent.
- [`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ddm_simulate.md)
  draws response times and boundary choices jointly, for building
  example and test data sets.
- The package registers its compatibility rows with the core at load
  through
  [`frmtmb::frmtmb_register_compat()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.html),
  so `frm_compat("wiener")` states what was exercised and what was not.
- Deliberate omissions: no `lcdf`, so `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) are refused; no
  variance function and no unit deviance, so
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) answers
  `type = "response"` only.
- `dev-findings.md` in the package source records what building this
  from outside frmtmb cost, as an acceptance test of the extension API.
