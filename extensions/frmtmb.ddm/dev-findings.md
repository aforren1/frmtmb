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

## 9. What worked without friction

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

| # | Severity | Finding |
|---|----------|---------|
| 1 | blocker / friction | `dec()` unreachable; aterm set is closed with no registration hook |
| 2 | friction | No `required_aterms`; a missing `vint()` fails silently |
| 3 | friction | No hook for a data-dependent link bound; relies on undocumented slot order. `init_dpars` silently dropped when non-finite |
| 4 | cosmetic | `frmtmb_register_compat()` documented inside the sampling-api bundle, deferring to the source |
| 5 | friction | Custom link objects accepted without validation |
| 6 | cosmetic | `cens()`/`trunc()` refusal does not mention the `lcdf` slot |
| 7 | friction | Core's own wiener example has a normalized-time range real data can exceed |
| 8 | cosmetic | `?frmtmb_ad_overload` says "every function slot"; measured, it is the tape-side slots only |

Nothing in this list stopped the family from being written. Finding 1 is
the only one that changed what the user has to type, and findings 2 and
3 are the only ones where the failure mode is silence rather than an
error. No core file was edited and no frmtmb internal was reached with
`:::` from package code.
