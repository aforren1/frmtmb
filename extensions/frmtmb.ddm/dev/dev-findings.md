# Findings: building a family from outside frmtmb

This is the other half of the deliverable. `frmtmb.ddm` was written
against frmtmb 0.47.0's exported API only, with no edits to core and no
`:::` calls from package code, as an acceptance test of the extension
API from an outsider's position. Every place the work was forced to
reach past the documented surface, replicate private machinery, or work
around a missing hook is recorded here with the API change that would
fix it.

Severity is `blocker` (the family cannot be written correctly without
it), `friction` (the family works, but the route to it is undocumented,
silent, or hand-rolled) or `cosmetic`.

---

## 1. `dec()` is unreachable, and the aterm set is closed

**Severity: blocker (for parameterization fidelity), friction (in practice)**

brms spells the drift-diffusion decision indicator as an addition term:

```r
brms::brm(rt | dec(decision) ~ condition, family = wiener(), data = dat)
```

frmtmb refuses it:

```
Addition term `dec()` is not supported yet (currently supported:
weights(), trials(), cens(), trunc(), se(), vint(), vreal(), mi())
```

The accepted set is a local character literal inside `parse_response()`
(`R/parse.R:16-17`) and nothing consults a registry. There is no
exported hook that adds an addition term, and the three registration
seams that do exist serve other purposes:
`frmtmb_register_frame_check()` validates an assembled frame,
`frmtmb_register_compat()` writes documentation rows, and
`frmtmb_register_prior_defaults()` supplies priors.

The family is still writable, because `vint()` carries an integer data
column through to the density, and that is the route this package took.
The cost is paid by the user, not the package: the spelling every
reference on the model uses does not work, the indicator must be
hand-coded to 0/1 rather than given as a factor, and the failure mode
for the natural attempt is an error that names eight addition terms
without hinting that `vint()` is the one to reach for.

**API change that would fix it.** Either
`frmtmb_register_aterm(name, arity, coerce)`, letting a family
contribute a term the same way it contributes compat rows; or, much
cheaper, an alias table so `dec(x)` parses as `vint(as.integer(x))` and
the closed set gains the one name that a whole family class needs.
Failing both, the refusal in `parse.R` should name `vint()` as the
escape hatch for a custom family, since that is what the message is
really telling an extension author.

**CLOSED in frmtmb 0.49.0**, by the first of those and by the third as
well. `frmtmb_register_aterm(name, arity, coerce)` is exported and
documented on its own page; the refusal in `parse.R` now lists the
registered terms alongside the core ones, names `vint()` and `vreal()`
as the general route, and points at the registry.

This package registers `dec` from `.onLoad()` with a `coerce` that
follows brms's own rule -- a factor or character vector read on its
levels, second level the upper boundary -- so `rt | dec(response) ~ x`
now works and takes a factor, which was the part of the cost the user
was paying. `vint()` is unchanged and still works; the two spellings
are pinned as one model in `test-surface.R`.

What the registry does NOT reach is the shape of the seam: an entry is
keyed by its own name, so a registered term delivers its value under
that name and cannot be an ALIAS that lands in `vint1`. That is what
leaves finding 2 half open below.

---

## 2. A family cannot declare that it requires an addition term

**Severity: friction, with a silent-wrong-answer failure mode**

`vint1` is not optional for this family: it is the boundary the trial
ended at, and the density is meaningless without it. But
`frmtmb_family()` has no slot that says so. If the user omits `vint()`,
`aterms[["vint1"]]` is `NULL`, and `NULL` in the density's arithmetic
propagates to `numeric(0)`: the log-likelihood becomes a sum over no
terms, and the fit does not fail. It returns.

Confirmed directly: a custom family whose `lpdf` reads a `vint1` that
was never supplied still produced a `frmtmb_fit`.

The workaround is to hand-roll the check in `valid_y`, which this
package does, and the message it raises has to explain the whole
`dec()` situation from finding 1. That is a per-family reimplementation
of something the framework knows how to express.

**API change that would fix it.** A `required_aterms` argument to
`frmtmb_family()`, checked at frame assembly with a message built from
the family name and the missing term. It is a few lines in core and it
removes a class of silent failure from every custom family that uses
`vint()`/`vreal()`, which is every custom family that needs per-row
data.

**PARTLY CLOSED in frmtmb 0.49.0.** `frmtmb_family(required_aterms =)`
exists, is checked at `R/frame.R:1150` before every other guard is
handed the addition-term values, and builds its refusal from the family
name and the spelling that supplies the missing term. There is a second
backstop in `R/objective.R:43-59` for a family that declares nothing:
a zero-length log-likelihood over a non-empty response is refused by
LENGTH, which resolves at tape-build time and leaves no branch on the
tape. Both are exactly right, and the silent-fit failure this finding
reported is gone from core.

It is not usable *here*, and the reason is a real gap rather than an
oversight. `required_aterms` is a conjunction: it names the terms a
density needs ALL of. Since finding 1 closed, this family needs EITHER
`dec()` or `vint1`, because `dec()` is the spelling to use and
`vint()` is the spelling that already worked and must keep working.
Declaring `"dec"` refuses every model written the old way; declaring
`"vint1"` refuses the spelling the package now documents; declaring
both refuses everything. So the check stays written out in `valid_y`,
and it is now the only hand-rolled piece left in this package.

`frmtmb_register_frame_check()` can express the disjunction and was
considered. It was not used: it runs at `R/frame.R:2298`, after
`family_finalize()` and after the aterm guards, so a model missing its
indicator would take a different and worse error first. `valid_y` is
earlier and its message is better.

**API change that would close the rest.** Let `required_aterms` carry
alternatives -- a list element naming several terms means "at least one
of these" -- so `required_aterms = list("dec", c("dec", "vint1"))`
reads the way it sounds. The refusal already names spellings through
`aterm_spelling()`, so it would print "write one of: dec(<column>),
vint(<column>)" with no new message machinery. This is not exotic: any
family that gains a better-named term while keeping `vint()` for
compatibility has it, which is every family that outlives one release
of the aterm registry.

---

## 3. A family cannot see the response when it builds its links

**Severity: friction (relies on undocumented call order)**

The Wiener density is zero at and below the non-decision time, so the
likelihood has a hard edge at `ndt = min(rt)`. A log link lets the
optimizer walk over it; the right answer is a link bounded above by
`min(rt)`, which makes the constraint structural.

But links are fixed when the family object is constructed, and the
family object is constructed before `frm()` ever sees `data`. Core's own
worked example in `vignette("case-studies")` resolves this by making the
user pass the bound:

```r
frm(bf(rt | vint(upper) ~ x), family = wiener_family(max_ndt = min(dat$rt)))
```

which asks the user to compute a quantity the framework already has, and
gives a wrong answer with no error if they compute it from the wrong
column.

This package instead has `valid_y()` write the bound into an environment
the link closures read at call time, so `wiener()` needs no argument. It
works, and it works only because `valid_y()` happens to run before any
link function does. That order was verified by instrumenting a family
and reading the call sequence (`valid_y -> linkfun -> linkinv`), not by
reading a document: nothing in `?frmtmb_family` states it. A refactor in
core that moved the init pass ahead of the response check would break
this package silently, in the sense that the bound would be `NA` and the
error would come from somewhere unrelated.

There is a second, quieter trap in the same area. `init_dpars` values are
pushed through the link and then **silently dropped** if the result is
not finite (`R/fit.R:1668`, `if (is.finite(val))`). For a bounded link
an init at or above the bound produces `Inf`, the init is discarded
without a word, and the fit starts from the template default instead.
The symptom is a slow or failed optimization with no indication that the
starting values were ignored.

**API changes that would fix it.**
1. Document the slot call order in `?frmtmb_family` -- specifically that
   `valid_y` runs once at frame assembly and before any link is used. It
   costs one sentence and turns this package's mechanism from a bet into
   a contract.
2. Better: a `link_data` or `family_finalize(y, aterms)` slot that
   returns a modified family, called at assembly. The
   data-dependent-bound problem is not exotic -- every shifted family
   has it -- and an explicit hook beats an environment side effect.
3. Warn rather than silently drop a non-finite `init_dpars` value.

**CLOSED in frmtmb 0.49.0**, by both 1 and 2. `?frmtmb_family` now has
a "Slot call order" section that states the sequence as measured on an
instrumented family, and a "Deriving a family from the data" section
whose worked example is this exact bounded-link problem.
`family_finalize(fam, y, aterms)` runs at `R/frame.R:1356`, after the
response is validated and before any link is used, and the frame layer
refreshes the per-dpar link copies from whatever it returns, so a
replaced link reaches the starting values, the tape and every post-fit
method.

The environment is gone from this package. `wiener()` now builds a
family whose bounded links raise if used before the data arrives, and
`family_finalize()` rebuilds the family with the bound, the
non-decision-time range's own bound, and the unreachable-row margin,
all derived from `min(y)`. The family object no longer lies about what
it is, and nothing depends on an undocumented order.

The variability work took a second, unplanned dependency on the same
slot: the margin that holds an unreachable row's decision time off the
singularity has to be scaled to the response's units, and there is no
other place a family can learn them. A `family_finalize` that returns a
family with a fresh `lpdf` closure gets that for free, which is a use
the slot's documentation does not mention and is worth mentioning.

Point 3 was not retested and is left open.

---

## 4. The compat seam is documented in the wrong place, and defers to the source

**Severity: cosmetic**

The brief said to read the docs for `frmtmb_register_compat()`. It has
no help page of its own: its roxygen block is `@noRd` and still carries
the stale comment "Not exported yet: the two contributors are still in
this package. Moving them out is what makes it public API." It *is*
exported, via a bulk `@rawNamespace export(...)` of forty-odd internals
attached to `?"frmtmb-sampling-api"` (`R/sampling-api.R:189-201`), and
it is documented there in a section called "The compatibility registry".

Two problems. A family extension has nothing to do with the sampling
API, so an author looking for the seam has no reason to open that page.
And that page's parameter documentation is:

> `...` Arguments of the individual functions; see the sections above
> and the source, which is the reference for these.

"The source is the reference" is a fair statement about forty internal
helpers bundled for one downstream package. It is not a fair statement
about `frmtmb_register_compat()` and `compat_rule_builder()`, which the
extension guide actively tells outside packages to call.

**API change that would fix it.** Split the two compat functions out of
the sampling-api bundle into their own documented page, cross-referenced
from `?frmtmb_family` and `?"frmtmb-extension-api"`, with real `@param`
entries and the status vocabulary. Delete the stale "not exported yet"
comment.

---

## 5. Custom link objects are accepted without validation

**Severity: friction**

A family may pass a list instead of a link name, and `get_link()`
(`R/links.R:107`) returns any list untouched. Nothing checks that the
list has `name`, `linkfun`, `linkinv` and `mu_eta`. A link missing
`mu_eta` fits, summarizes and predicts happily, and then fails inside
`predict(se.fit = TRUE)` at a call site far from the family that caused
it.

**API change that would fix it.** Validate the four fields in
`get_link()` when it is handed a list, and say which field is missing
and which dpar's link is at fault. Ten lines, once, instead of a
debugging session per extension author.

---

## 6. The `cens()`/`trunc()` refusal hides the way in

**Severity: cosmetic**

Asking for truncation on a wiener model gives:

```
cens()/trunc() need a family with a CDF (currently: gaussian,
lognormal, poisson, exponential, weibull, inverse.gaussian, cox)
```

The list is accurate and the refusal is correct. But to an extension
author it reads as a closed set of blessed families, when the actual
requirement is one slot: a family that declares `lcdf` gets `cens()` and
`trunc()`, whoever wrote it. This package legitimately cannot supply one
and the refusal is the right outcome here, but the message sent me to
read `families.R` to find out whether the door was closed or merely shut.

**API change that would fix it.** Add a clause: "a family supplies one
through the `lcdf` argument of `frmtmb_family()`". One sentence.

---

## 7. Core's own worked example has a range limit a user can exceed

**Severity: friction**

`vignette("case-studies")` section 11 builds this same family, and
`tests/testthat/test-vignette-wiener.R` pins it. It is honest work: it
states that the truncation is fixed at `K = 10` because the tape cannot
branch, and it measures what that buys, pinning agreement to 1e-9 out to
a normalized time of 4 and recording that past 6 the error is
cancellation rather than truncation, so more terms do not help.

The gap is between that measurement and what response-time data does.
The normalized time is `u = (rt - ndt) / bs^2`. A two-second response
with a boundary separation of 0.5 gives `u = 8`, which is past the
documented range, and narrow boundaries are exactly what a speed-stressed
condition produces. Nothing in `frm()` warns; the likelihood is simply
wrong, and it is wrong smoothly, so the fit converges to something.

Measured here: at `u = 16` the small-time series alone is off by more
than 1% relative in the log density, and raising `K` to 40 does not
recover it. At `u ≈ 15.9` its error reaches 72%.

This package's answer is the reason it evaluates both series: the
small-time and large-time series are both computed at a fixed truncation
and their logs blended with a logistic weight in `log(u)`, which is
Navarro and Fuss's own criterion with the step smoothed so it
differentiates. Measured agreement with `RWiener::dwiener()` is below
1e-12 relative in the log density over 1125 parameter combinations
spanning `u` from 8e-4 to 50.

**Change that would fix it.** Either carry the two-series blend into the
vignette, or state the applicability bound where a reader will hit it:
"valid for `(rt - ndt) / bs^2` below about 4; check your data". The
second is cheap and would have been enough.

---

## 8. "Every function slot" is three of them

**Severity: cosmetic**

`?frmtmb_ad_overload` says:

> frmtmb applies it to every function slot of `frmtmb_family()` and
> `frmtmb_structure()`

Measured on a constructed family, by reading the stored bodies back:

| slot | carries the `ADoverload` bindings |
|------|-----------------------------------|
| `lpdf` | yes |
| `post$mean_fn` | no |
| `sim` | no |
| `valid_y` | no |

The behavior is right. `mean_fn`, `sim` and `valid_y` run on doubles,
off the tape, and wrapping them would buy nothing. Only the tape-side
slots (`lpdf`, and `lcdf` when present) are wrapped, which is the
correct rule. The sentence just does not say that rule.

It matters because the neighboring "Tape-safe scope" section of
`?frmtmb_family` is precise and correct about the thing that actually
bites -- that a helper the density CALLS inherits nothing -- so a reader
has one exact statement and one loose one on the same subject.

**API change that would fix it.** Say which slots: "frmtmb applies it to
the tape-side slots, `lpdf` and `lcdf`." One word changed, and the
contract becomes checkable.

---

## 9. A family cannot constrain two dpars jointly

**Severity: friction, with a silent-wrong-answer failure mode**

Added while building the across-trial variability family, which needs
two constraints frmtmb has no way to express.

The uniform start point runs over `bias +/- sz / 2` and the model is
only defined while that range stays inside `(0, 1)`. The uniform
non-decision time runs over `ndt +/- st / 2` and the usual
parameterization asks for `ndt - st / 2` above zero. Both are
constraints on a PAIR of distributional parameters. A link is a
property of one parameter, and it is the only constraint mechanism a
family has, so neither can be made structural the way `ndt < min(rt)`
was.

What this package does instead is pick links that make the common case
safe and document the rest: `sz` gets a logit, so the width is below 1
and the range is inside `(0, 1)` whenever `bias` is 0.5, which is what
it is in almost every fit of this model; `st` gets a logit scaled onto
`(0, 2 * max_ndt)`, which keeps the width in the units the data can
support without pinning the lower edge.

The failure mode outside those is not an error. Measured on the series
directly: at a relative start point of exactly 0 the log density is
about -45 where it is -0.34 at 0.5, so the boundary is a barrier -- but
just PAST it, at -0.02, the series come back up to -3.06. The density
is not monotone across the boundary, so a wide enough excursion has a
spurious mode on the far side and an optimizer that reaches it gets a
number rather than a refusal.

**API change that would fix it.** A `constraints` slot on
`frmtmb_family()` taking functions of the whole dpar vector, checked at
the starting values and reported at the optimum, would catch it after
the fact. Making it structural needs more: a way for a family to
declare a transformed parameter -- `sz` estimated as a fraction of the
room `bias` leaves -- which is a reparameterization hook rather than a
link, and is a larger change than this one family justifies asking for.

---

## 10. What worked without friction

Worth recording, because an acceptance test that only lists complaints
is not a measurement.

- **The tape-safe scope contract is right and is documented.**
  `?frmtmb_family`'s "Tape-safe scope" section states that `lpdf` gets
  `c()`, `[<-` and `diag<-` bound automatically and that a helper the
  density calls does not. That is exactly the trap, stated plainly, in
  the right place. This package's density needs no bindings at all, by
  construction: it is pure arithmetic with no `c()` and no `[<-`, so it
  runs on advectors and doubles alike. `frmtmb_ad_overload()` was read
  and deliberately not needed.
- **`check_custom_family()` is the right tool in the right place.** One
  exported call taped the density, compared taped against plain, and
  compared the AD gradient against central differences. An extension
  author's first question is "is my density even differentiable
  correctly", and the answer is one function.
- **`vint()` reaches everywhere it needs to.** The same values arrive in
  `lpdf`, in `post$mean_fn`, in the `sim` slot, and are required rather
  than silently dropped on `newdata`. Once the `dec()` problem is
  accepted, the data path itself is complete and consistent.
- **A dpar can be fixed on the response scale.** `bf(..., bias = 0.5)`
  puts `bias` at 0.5 and its linear predictor at 0, through the link,
  with a clear refusal when the constant is outside the link's range.
  For this family that is the common case, since bias is usually held at
  0.5.
- **Custom link objects work.** Given finding 5, they work once you know
  the four fields, and the scaled-logit link is what makes the
  non-decision time constraint structural.
- **`frm_sample()` sampled a family it had never heard of.** A short
  NUTS chain ran on a wiener fit with no special handling: an
  out-of-tree family, a custom link object, and a `vint()` payload all
  passed through `frmtmb.sample` untouched, and the brms default priors
  were applied to the dpars by name. Nothing was written to make that
  work. That is what a clean seam looks like.
- **The whole post-fit surface either worked or refused for a stated
  reason.** `summary()`, `fixef()`, `logLik()`, `AIC()`, `fitted()`,
  `predict()` on training rows and on `newdata`,
  `residuals(type = "response")`, `simulate()`, `par_template()`,
  `get_prior()`, `set_prior()`, `weights()` and a random-effect block
  all worked on a family core has never seen. The three refusals
  (`pearson`, `deviance`, `cens()`/`trunc()`) each named the missing
  slot rather than failing obscurely. No method silently produced a
  wrong number.
- **`frm_compat()` accepted contributed rows from `.onLoad()` and
  merged them into the core table.** `frm_compat("wiener")` reads
  exactly like a core family's rows.

---

## Summary

| # | Severity | Finding | Status |
|---|----------|---------|--------|
| 1 | blocker / friction | `dec()` unreachable; aterm set is closed with no registration hook | CLOSED in 0.49.0 |
| 2 | friction | No `required_aterms`; a missing `vint()` fails silently | PARTLY: exists, but cannot express "either of two terms" |
| 3 | friction | No hook for a data-dependent link bound; relies on undocumented slot order. `init_dpars` silently dropped when non-finite | CLOSED in 0.49.0 (the link hook and the documented order; the `init_dpars` half untested) |
| 4 | cosmetic | `frmtmb_register_compat()` documented inside the sampling-api bundle, deferring to the source | not retested |
| 5 | friction | Custom link objects accepted without validation | not retested |
| 6 | cosmetic | `cens()`/`trunc()` refusal does not mention the `lcdf` slot | not retested |
| 7 | friction | Core's own wiener example has a normalized-time range real data can exceed | this package's two-series blend still stands |
| 8 | cosmetic | `?frmtmb_ad_overload` says "every function slot"; measured, it is the tape-side slots only | not retested |
| 9 | friction | Nothing lets a family constrain two dpars jointly | open; see below |

Nothing in this list stopped the family from being written. Finding 1 is
the only one that changed what the user has to type, and findings 2, 3
and 9 are the ones where the failure mode is silence rather than an
error. No core file was edited and no frmtmb internal was reached with
`:::` from package code.

Findings 1 and 3 were acted on in frmtmb 0.49.0 and this package now
uses both seams: `frmtmb_register_aterm()` for `dec()`, and
`family_finalize()` for the data-derived bounds, which between them
deleted a hand-rolled coercion, an environment side effect and a bet on
an undocumented call order. Finding 2's seam exists and is the right
one; it is only the conjunction-versus-disjunction gap that keeps this
family from using it. That is a good ratio for one release, and it is
the evidence that the extension API is being maintained as an API
rather than as whatever the in-tree families happened to need.

---

## Corrections after independent review

Three measurements in this package's own write-up did not survive an
independent reviewer repeating them. They are corrected here rather
than quietly edited away, because a ledger that only records other
people's defects is not a measurement either.

**The drift-variability timing claim was wrong.** This lane reported
"0.8 s with `sv` against 3.5 s for plain", i.e. that switching on the
drift closed form made a fit *faster* than the plain family. It does
not reproduce. Measured at the fit level over three seeds, `sv`/plain
is 1.04, 0.97 and 1.07; in fresh processes, plain-first gives 2.62
then 2.86, and `sv`-first gives 2.42 then 2.41. At the density level
both are 0.0002 s per call. There is no configuration in which `sv` is
four times faster than plain.

Two errors combined. The 3.5 s figure matches the `st` path, which the
reviewer measured at 4.3 to 5.5 s, so the comparison mislabeled which
family was being timed. And the timings were taken in one process, one
after another, where an 8x warm-up swing (plain fell from 1.29 s to
0.16 s between the first and second repetition) will manufacture any
ratio asked of it.

The defensible statement, which is the interesting one anyway, is that
the drift closed form is FREE relative to the plain density: it adds a
square root and a few flops and no nodes, so the full Ratcliff drift
distribution costs nothing over a fixed drift. That is what the NEWS
entry and the vignette now say.

**How this package times things from now on.** A timing comparison
runs each arm in a FRESH process, and warms both arms before the
measured repetition. Numbers taken from consecutive runs in one warm
session are not reported as a ratio between arms. A ratio is quoted
only with the number of seeds behind it and the spread across them.

**"A separate, literally unchanged route" was the wrong wording.**
`ddm_lpdf_lower()` IS edited: `u <- t / (a * a)` became
`ur <- t / (a * a); u <- ddm_floor(ur, ddm_u_floor)`. The substantive
claim survives, and the reviewer's check is stronger than the one this
lane made: over a grid of 11520 combinations of `t`, `v`, `a`, `w` and
the boundary, the log density is byte-for-byte identical to the
previous release on the whole positive-`t` segment (`identical()` TRUE,
maximum absolute difference 0). The only difference anywhere is at
`t <= 0`, where `NaN` became `-Inf`, which is the deliberate fix. The
right wording is "byte-for-byte identical wherever the density is
defined", not "unchanged".

The one caveat the reviewer added and this lane had not stated: adding
`ddm_u_floor` is a rounding no-op only down to about 1e-300, and at
`u <= 1e-307` it perturbs. `u = t / a^2` is that small only where the
density has already underflowed, so nothing depends on it, but the
claim is "inert on the support", not "inert for every double".

**The negative-span frequency was a point estimate, not a bound.**
This lane measured 9 to 14 percent of fully-past rows rounding their
span negative and wrote it as though it characterized the phenomenon.
It characterizes one configuration. Across `ndt` and `st` the reviewer
measured 0 to 40 percent: `ndt = 0.50, st = 0.02` gives 0.0 percent,
`ndt = 0.40, st = 0.10` gives 5.0 to 11.0 percent over three seeds,
`ndt = 0.35, st = 0.14` gives 14.6, `ndt = 0.45, st = 0.06` gives 17.0,
and `ndt = 0.30, st = 0.20` gives 40.1. The defect is not rare and its
rate is a function of where the optimizer is standing, which is why the
regression test pins a specific point rather than relying on a
frequency.

**A contributed compat row naming an unknown feature is silently
dropped.** Found while adding the `gddm` and `lba` rows, and it had
already cost this package two rows without anyone noticing. A rule
whose `feature_b` is not in the registry's vocabulary does not warn and
does not error: `frm_compat()` simply returns without it. Two of the
`wiener` rows written in the previous round were being discarded that
way. `mixture()` is spelled `mixture` in the vocabulary, with no
parentheses, because it is a `structure` rather than a callable
addition term; and `dec()` was in no vocabulary at all, because
registering an addition term with `frmtmb_register_aterm()` does not
also make it a compat FEATURE. The fix on this side is one more entry
in the `features` vector, `"dec()" = "aterm"`, which is easy once you
know, and invisible until you check the output rather than the input.

**API change that would fix it.** `frmtmb_register_compat()` should
warn when a rule names a feature the vocabulary does not contain, at
registration or at first `frm_compat()` call. The registry already
knows the full vocabulary at that point, the check is a `setdiff()`,
and the current behavior means the seam's failure mode is a table that
looks complete and is not. A contributor's first sign of trouble
should not be counting rows.
