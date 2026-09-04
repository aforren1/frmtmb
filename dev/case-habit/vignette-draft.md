# Case study draft: replicating a response-preparation model

> **Superseded.** This draft became `vignettes/habit.Rmd`. Read that for the
> shipping version. The draft is kept because it records reasoning and
> alternatives the vignette has no room for.

**Status: draft narrative, not a vignette.** This is deliberately not an `.Rmd`
in `vignettes/`. The data it describes cannot ship until licensing is settled,
so this file holds the story in a form that can become a vignette on the day it
is. See `model-notes.md` for the licensing position and the model derivation.

The intended audience is someone who already knows `frm()` for linear models
and wants to see how far the non-linear grammar goes against a real, published,
independently implemented model.

---

## 1. The pitch

Nine years ago, in July 2017, the frmtmb maintainer opened an issue on the
analysis repository for a study they had co-authored. It is titled
"Hierarchical model", and it opens "Work in progress...". It sketches, in
pseudo-code, exactly the model this case study builds:

```r
response ~ categorical(theta)
theta = alpha * phi
alpha = f(guess_rate, upper_asymptote_a, upper_asymptote_b, habit)
phi   = g(preparation time, mu_a, mu_b, sigma_a, sigma_b)

mu_a[subject] ~ normal(mu_a_mean, mu_a_sd)   # partially pooled
mu_b[subject] ~ normal(mu_b_mean, mu_b_sd)   # partially pooled
```

The issue notes that a Stan prototype existed and "takes at least an hour to
fit". It was never closed.

That is the shape of the case study. The published analysis fit every
participant in isolation, once per model, with a MATLAB optimizer. The question
is whether frmtmb can (a) reproduce that per-subject analysis exactly enough to
be trusted, and then (b) do the thing the original procedure could not: fit one
model to a whole group, with correlated per-subject effects on the
preparation-time parameters, in one pass.

## 2. The science, briefly

Participants learn to map four symbols onto four keys. Then the mapping is
revised. In a forced-response test, a sequence of tones obliges a response at a
fixed moment, so the experimenter controls *preparation time* directly rather
than measuring it. Sweeping preparation time from near zero to about 1.2 s
traces out the whole time course of response preparation.

The signature of a habit is a bump. At short preparation times the response is
a guess. At long preparation times it is correct. In between, if the original
mapping has become habitual, there is a window where the *old* response is
prepared and emitted before the new one is ready. The model formalizes that as
a race between two processes with Gaussian preparation-time distributions, and
the whole scientific claim reduces to whether the data need the second process.

The paper's answer: after four days and after twenty days of practice, yes;
after minimal practice, no.

## 3. The model as a formula

Two independent processes are ready with probabilities \(\Phi_A(t)\) and
\(\Phi_B(t)\), giving four preparation states. A coefficient matrix maps those
states to the three response categories. Everything is in `model-notes.md`
section 1; what matters here is the shape.

The natural spelling would be a `categorical()` family with an identity link
and `nl = TRUE` bodies. **At v0.47.0 both halves of that are refused**, for
reasons worth understanding rather than working around blindly:

* `nl = TRUE` requires a family with a single `mu` location parameter
  (`R/parse.R:1225`). A categorical family has one predictor per non-reference
  category, so the combination is refused outright. The refusal is registered
  in the compatibility table at `R/compat.R:609`.
* `categorical()` accepts the logit link only (`R/families.R:4297`) and forms
  the softmax internally with `RTMB::logspace_add`. There is no identity-link
  hook that would let a body hand it probabilities.

The way through is `nlf()`, which names a distributional parameter directly and
so never passes through the single-`mu` gate, plus a change of variable. If the
model produces category probabilities \(c_1, c_2, c_3\), then feeding
\(\log(c_2/c_1)\) and \(\log(c_3/c_1)\) to the two non-reference categories
makes the softmax return the \(c_i\) unchanged. The softmax is not in the way;
it is invertible, and the log-ratio is its inverse.

That gives a formula that reads like the paper:

```r
bf(resp ~ 1) +
  nlf(PhiA ~ RTMB::pnorm((t - muA) / sgA)) +
  nlf(PhiB ~ RTMB::pnorm((t - muB) / sgB)) +
  nlf(c1 ~ qI*(1-PhiA)*(1-PhiB) + ((1-qA)/3)*PhiA*(1-PhiB) + qB*PhiB) +
  nlf(c2 ~ qI*(1-PhiA)*(1-PhiB) + qA*PhiA*(1-PhiB) + ((1-qB)/3)*PhiB) +
  nlf(c3 ~ (0.5-qI)*(1-PhiA)*(1-PhiB) + ((1-qA)/3)*PhiA*(1-PhiB) +
           ((1-qB)/3)*PhiB) +
  nlf(muhabit ~ log(c2 / c1)) +
  nlf(muother ~ log(2 * c3 / c1)) +
  lf(muA ~ 1, sgA ~ 1, qA ~ 1, muB ~ 1, sgB ~ 1, qB ~ 1, qI ~ 1)
```

Three things in that block are worth pointing out to a reader.

**`nlf()` bodies compose.** `PhiA` and `PhiB` are not distributional parameters
of any family. They are intermediate names, and frmtmb evaluates the bodies in
dependency order, so the algebra can be laid out in the order a person would
write it on paper instead of being inlined into one expression.

**The factor of two in `muother`.** The reference implementation's third
category covers two keys but is scored with the probability of one, so its
category probabilities do not sum to one. Doubling it restores a proper
distribution. This is the single most important detail for anyone comparing
numbers against the original: it shifts the log-likelihood by exactly
\(n_3 \log 2\). See `model-notes.md` section 1.

**`RTMB::pnorm`, not `pnorm`.** A bare `pnorm()` in a non-linear body resolves
lexically to **stats** and silently strips the automatic-differentiation class.
This is already recorded as a known trap in `NEWS.md`.

> **When the `nlenv` lane lands**, the two `RTMB::` prefixes above come off and
> the bodies become plain R. That is a real gain for a vignette like this one:
> the formula block is the thing a reader is meant to compare against the
> equations, and a namespace prefix on a normal CDF is noise that has nothing
> to do with the model. Nothing else in this case study changes; the
> qualification is required today only for the two CDF calls, since the
> arithmetic operators are already overloaded.

## 4. Reproducing the published fits

The case study earns trust before it claims anything, in three steps.

**Step one: a reference likelihood in plain R.** `replicate.R` contains a
hand-written translation of the MATLAB likelihood, with no frmtmb involvement.
Evaluated at the authors' own archived parameter estimates it reproduces their
stored log-likelihoods for every model, group and participant, to a maximum
absolute difference of 4.5e-13, in both the `fmincon` and the BADS archive. The
ridge penalty on the preparation-time standard deviations reproduces to the
same precision. Those checks also settle a stale comment in the original source
that gives the parameter vector in the wrong order.

**Step two: the same numbers out of `frm()`.** The ridge penalty on the
preparation-time standard deviations is a Gaussian MAP penalty, so it goes
straight into the prior vocabulary:

```r
set_prior("normal(0.07, 0.02236)", nlpar = "sgA")   # 1/(2 s^2) = 1000
```

With that and the original optimizer's box constraints in place, the
per-subject fits agree closely with the archive. The `no-habit` model matches
to below 1e-4 for every participant. The `habit` model agrees to a maximum
log-likelihood difference of 0.0027, with every preparation-time and asymptote
parameter within 1e-3.

A note the vignette should carry, because it is the sort of thing that costs an
afternoon: `logLik()` on a fit carrying a `set_prior()` penalty returns the
log-likelihood **plus the log prior density**. It is therefore not comparable
to an unpenalized published number. The comparison is made by evaluating the
hand-written reference likelihood at frmtmb's own estimates.

**On bounds.** Bounds are carried in the prior vocabulary,
`set_prior("", nlpar = ..., lb = , ub = )`, or as a transform inside the body.
The hierarchical model in section 6 does the latter; the per-subject fits here
do the former.

That was not always possible. Before v0.49, `set_prior("", nlpar = "muA",
lb = 0)` was refused, because the bounds path keyed the bound on the nonlinear
parameter's raw sub-formula column name, `"(Intercept)"`, rather than on the
qualified `"muA_(Intercept)"` that bound resolution looks for. A `normal()`
prior placed through the same `nlpar` argument worked, so the asymmetry was
specific to `lb`/`ub`. Details in `model-notes.md` section 5.

One thing is worth saying out loud even now that the gap is closed, because it
is easy to assume a prior means a transform. It does not: `lb`/`ub` is a hard
box, the same box the optimizer would get from any other spelling. The purpose
of this stage is
to reproduce a box-constrained run whose optima sit **on** their bounds: the
habitual asymptote is at its limit for most participants by construction,
because it is not separately identified from the habit-strength parameter. A
logit transform turns an attainable boundary into an asymptote, so the
replication would stop comparing like with like. Hard boxes are the right tool
for reproducing hard boxes; transforms are the right tool once random effects
need an unbounded scale to live on.

**Step three: explain the disagreements rather than hiding them.** The
`flex-habit` model disagrees by up to 2.3 log-likelihood units, in both
directions. `frm_allfit()` plus restarts from the archived estimates show these
are multimodality, not a porting error: after multi-start, frmtmb matches or
beats the archived optimum on nine of the ten worst cases, including one
participant where it finds a solution 2.3 units better than the published fit.
The remaining case is worse by 0.0013.

A note for the vignette: this is the honest shape of a replication, and it is
worth showing rather than tidying away. `frm_allfit()` exists precisely for
this question.

## 5. What the replication turned up

Two findings came out of the exercise that are not about frmtmb at all.

**The trial subsetting in the original preprocessing matches on values, not
indices.** Remapped and unchanged trials are separated with
`ismember(recodedX, revisedX)`, comparing millisecond-rounded preparation
times. Collisions are common, so a trial whose preparation time happens to
occur in the other stimulus class is admitted to both datasets. Replaying that
operation against the raw structures reproduces the archived per-participant
trial counts exactly, which confirms the mechanism. Refitting from
identity-based subsetting *strengthens* the paper's conclusion in both
directions, because the contamination was diluting the habit signal. Details
and figures in `model-notes.md` section 4.

**Only one of the three model pairs is nested.** `flex-habit` contains both
`habit` and `no-habit` at opposite ends of its habit-strength parameter, but
`habit` and `no-habit` are not nested in each other. The paper compares those
two by AIC, which is correct, and reserves the likelihood-ratio test for the
nested pair, which is also correct. Worth flagging for a vignette reader:
`anova()` on frm fits checks that the fits share a number of observations, not
that they are nested. Choosing a legitimate pair is the analyst's job.

## 6. The value-add

Everything above is replication. The new part is one model per group instead of
one model per participant, with correlated subject-level effects on the
preparation-time parameters, spelled with the `|ID|` grouping syntax:

```r
lf(lmuB ~ 1 + (1 | p | subject), ldmu ~ 1 + (1 | p | subject),
   lsgA ~ 1 + (1 | p | subject), lsgB ~ 1 + (1 | p | subject),
   lqA ~ 1, lqB ~ 1, lqI ~ 1)
```

Random effects need an unbounded scale, so the bounded parameters are
transformed inside their own bodies rather than box-constrained, which is both
the house spelling and a requirement here:

```r
nlf(sgA ~ exp(lsgA))
nlf(qB  ~ 0.5 + 0.4999 / (1 + exp(-lqB)))
```

The headline is that the paper's group-level conclusion falls out of three
fits rather than 165: the no-habit model wins outright after minimal practice,
and the habit model wins by more than 200 AIC points after both four and
twenty days. The 20-day group recovers the textbook habit geometry unprompted,
with the habitual process prepared earlier and more sharply than the
goal-directed one. Full results, including the subject-level correlation
structure that the per-subject procedure could not estimate at all, are in
`hierarchical-results.md`.

Two practical notes for the eventual vignette. Non-linear fits of this size
need their starting values placed deliberately; starting the correlated block
from a generic guess stalls the optimizer, and starting it from the pooled
solution does not. And a tempting reparameterization can be actively harmful:
writing `muA = muB - exp(ldmu)` to carry over the ordering constraint that
`fmincon` imposed as a linear inequality works for a fixed-effect fit and
breaks as soon as that offset carries a random effect, with the variance
component pinned at a round number and correlations of order 1e-12. That
failure signature is worth teaching.

## 7. Bootstrap bands

The open TODO on the cleaned-data pull request, from September 2017, reads
"Still TODO: bootstrap CIs on curves". `frm_bootstrap()` answers it in one
call, because `FUN` may return any numeric vector, including the fitted
response probabilities evaluated on a grid of preparation times. The bootstrap
distribution is then a band around the entire curve rather than around a single
parameter.

The result is the paper's central figure with intervals attached. Probability
of emitting the habitual response at its peak: 0.243 after minimal practice,
which is just the guess rate and occurs at time zero, so there is no bump at
all; 0.342 at 0.45 s after four days; and 0.516 [0.471, 0.572] at 0.40 s after
twenty days, more than double the guess rate. The bump is absent without
training and grows with it.

Two things make this cheap. `FUN` is arbitrary, so the curve is computed by the
same hand-written reference function used for validation, not by a second
implementation. And the band is taken from the complete-pooling fit, which has
no random effects to re-integrate on each of the refits, so 100 draws per group
is quick. A useful correctness check falls out for free: the three probabilities
sum to 1.000000 at every grid point, confirming the log-ratio parameterization.

## 8. What a reader should take away

* The non-linear grammar reaches a real published cognitive model, and the
  resulting formula block is close enough to the paper's equations to be read
  side by side with them.
* `nlf()` is the more capable of the two non-linear spellings. `nl = TRUE`
  cannot express a multi-category likelihood at all.
* An invertible link is not an obstacle. A hard-wired softmax can be driven
  through its own inverse.
* Replication discipline is the point: an independent reference likelihood
  first, then the package, then an explicit account of every disagreement.

## 9. Before this can ship

1. Licensing. The reference repository has no license file. Nothing from it is
   in the frmtmb tree, and this draft carries no data.
2. A decision on which data source to ship. The archived `.mat` is what the
   published fits used; the cleaned CSVs are better documented and better
   behaved. `model-notes.md` section 5 lays out where they disagree.
3. Whether to report the subsetting finding here, upstream, or both. It
   strengthens the paper rather than undermining it, but it is the original
   authors' call to make.
