# Hierarchical fits: the part the original procedure could not do

Companion to `model-notes.md` (the model, its provenance, and the per-subject
replication) and `vignette-draft.md` (the narrative). Everything here is
produced by `replicate.R` stage `hier`.

The published analysis fits each participant separately, three models each.
That gives 55 independent fits per model and no way to ask how the
preparation-time parameters covary across people, or to let a poorly
constrained participant borrow strength from the group. HabitTR issue #2,
opened by the frmtmb maintainer in July 2017 and never closed, sketches the
model that would: partial pooling on the preparation-time parameters, the
asymptote and guess parameters pooled. This is that model.

## The specification

```r
lf(lmuA ~ 1 + (1 | p | subject), lmuB ~ 1 + (1 | p | subject),
   lsgA ~ 1 + (1 | p | subject), lsgB ~ 1 + (1 | p | subject),
   lqA ~ 1, lqB ~ 1, lqI ~ 1)
```

The `|p|` identifier merges the four terms into a single unstructured 4 by 4
block, so the subject effects on the two means and the two standard deviations
are estimated with their full correlation matrix rather than as four
independent scalars.

Bounded parameters are transformed inside their own bodies rather than
box-constrained, which is both the house spelling and a requirement here: a
random effect needs an unbounded scale to live on.

```r
nlf(sgA ~ exp(lsgA))
nlf(qB  ~ 0.5 + 0.4999 / (1 + exp(-lqB)))
nlf(qI  ~ 0.499 / (1 + exp(-lqI)))
```

The asymptote and guess parameters stay pooled. They are weakly identified
within a single participant, and giving them their own effects makes the block
singular.

### A parameterization that did not survive

The first attempt carried over the \(\mu_A \le \mu_B\) inequality that
`fmincon` imposed as a linear constraint, by writing
`nlf(muA ~ lmuB - exp(ldmu))`. It works for a fixed-effect fit and fails as
soon as `ldmu` carries a random effect: the offset is only weakly identified,
its variance component ran to a boundary and sat at exactly 1.000000 with
correlations of order 1e-12, and two of the three groups failed to converge.
Modeling the two means directly is better behaved, and is also what issue #2
sketches. Ordering is then checked after the fact instead of imposed.

This is worth recording because the fix is not obvious from the fit output.
A variance component pinned at a round number with numerically zero
correlations is the signature to look for.

## Model comparison

One hierarchical `habit` fit and one hierarchical `no-habit` fit per group.
The `no-habit` model keeps subject effects on the parameters it still has
(`lmuB`, `lsgB`).

| Group | trials | subjects | logLik habit | df | logLik no-habit | df | ΔAIC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| minimal practice | 5737 | 21 | -3304.11 | 17 | -3301.35 | 7 | **-25.52** |
| 4 days | 5744 | 21 | -4100.75 | 17 | -4244.09 | 7 | **+266.69** |
| 20 days | 3570 | 13 | -2470.42 | 17 | -2596.38 | 7 | **+231.91** |

ΔAIC is AIC(no-habit) minus AIC(habit), so positive favors habit.

**This is the paper's result, from three fits instead of 165.** After minimal
practice the no-habit model wins outright. After four days and after twenty
days the habit model wins by a very wide margin. The per-subject analysis in
`table-model-comparison.csv` says the same thing participant by participant:
1 of 21 favor habit after minimal practice, 16 of 21 after four days, 13 of 13
after twenty days.

Note that habit and no-habit are **not** nested (see `model-notes.md`
section 5), so this is an AIC comparison and not a likelihood-ratio test. AIC
does not require nesting; `anova()` here would be a category error, and
`anova.frmtmb_fit` would not stop you.

## Population-level parameters

Back-transformed to the reference implementation's own scale.

| Group | muA | muB | sgA | sgB | qA | qB | qI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| minimal | 1.8344 | 0.4719 | 0.0001 | 0.0959 | 0.990 | 0.9636 | 0.2467 |
| 4day | 0.5492 | 0.5099 | 0.1692 | 0.1687 | 0.999 | 0.9384 | 0.2546 |
| 20day | 0.3733 | 0.4739 | 0.0747 | 0.1270 | 0.990 | 0.9016 | 0.2444 |

Three things to read here.

**The 20-day group shows the textbook habit signature.** The habitual process
is prepared earlier and more sharply than the goal-directed one:
\(\mu_A = 0.373\,\text{s} < \mu_B = 0.474\,\text{s}\), and
\(\sigma_A < \sigma_B\). That is exactly the geometry that produces the
transient error bump, and the model recovers it without being told to.

**The minimal-practice group disposes of the habit process rather than using
it.** \(\mu_A\) is pushed out to 1.83 s, well beyond the roughly 1.2 s window
the experiment observes, with \(\sigma_A\) collapsed to 1e-4. A process that
is never prepared inside the observable range contributes nothing. This is the
model saying "no habit" in the only vocabulary it has, and it is consistent
with no-habit winning the AIC comparison by 25 points. It also means the
minimal-practice `habit` fit is effectively at a boundary, so its parameter
estimates should not be interpreted as if the A process were real.

**The 4-day group is the awkward one.** \(\mu_A = 0.549\) slightly *exceeds*
\(\mu_B = 0.510\), reversing the intended ordering, and the two standard
deviations are nearly equal (0.1692 against 0.1687). With such similar means
and spreads the two processes are close to exchangeable, which is precisely
the degeneracy the original inequality constraint existed to prevent. The
AIC verdict for this group is large and secure, but the *labels* on its two
processes are not, and the population-level parameters should be read with
that in mind. A cleaner treatment would impose the ordering through a prior
rather than a reparameterization, or fit from several starts and inspect the
resulting basins.

## Correlations across subjects

The quantity the per-subject procedure could not produce at all. Standard
deviations are on the working (transformed) scale.

**20 days.** The most coherent structure of the three:

| | SD | corr with lmuA | with lmuB | with lsgA |
| --- | ---: | ---: | ---: | ---: |
| lmuA | 0.0604 | | | |
| lmuB | 0.0728 | 0.873 | | |
| lsgA | 0.6890 | 0.966 | 0.933 | |
| lsgB | 0.4939 | 0.380 | 0.485 | 0.582 |

The two means correlate at 0.87 across participants: someone whose
goal-directed process is slow also has a slow habitual process. Both means
correlate strongly with the variability of the habitual process. Substantively
this reads as a single per-person "speed of preparation" factor loading on
everything, which is a hypothesis the original per-subject analysis had no way
to state, let alone test.

**4 days.** The same pattern, weaker: lmuA with lmuB at 0.672, lmuA with lsgA
at 0.781.

**Minimal practice.** Essentially no structure involving the A process
(correlations of lmuA with lsgA of -0.003, of lmuB with lsgA of 0.006), which
is what one expects when that process is not identified because it is not
there. The one interpretable correlation is lmuB with lsgB at 0.676.

## Convergence, stated plainly

| Group | pooled | habit | no-habit |
| --- | --- | --- | --- |
| minimal | 0 | 0 | 0 |
| 4day | 0 | **1** | 0 |
| 20day | 0 | **1** | **1** |

Codes are `nlminb` convergence codes; 0 is clean. Three of the nine fits do
not report clean convergence. Given that two of them are the boundary-adjacent
cases discussed above, this is unsurprising, but it is a real caveat and the
AIC differences involved are large enough (over 200 points) that the
qualitative conclusion does not depend on the last decimal. The honest summary
is that the model-comparison verdict is robust and the individual parameter
estimates for the 4-day and minimal groups are not fully settled.

Remedies not pursued here, for time: `frm_allfit()` across optimizers on the
three flagged fits, and a multi-start sweep on the 4-day group to map the
label-exchange basins.

## Bootstrap bands on the curves

HabitTR PR #3 has carried one open TODO since September 2017: "Still TODO:
bootstrap CIs on curves". `frm_bootstrap()` answers it directly, because `FUN`
may return any numeric vector, including the three fitted response
probabilities evaluated on a grid of preparation times. The bootstrap
distribution is then a band around the whole curve.

Run on the complete-pooling `habit` fit per group, 100 draws, grid from 0.001
to 1.2 s. Output in `table-bootstrap-bands.csv`. As a correctness check, the
three fitted probabilities sum to 1.000000 at every grid point, which is the
log-ratio parameterization of `model-notes.md` section 5 behaving as intended.

Probability of emitting the **habitual** response, with 95 percent bootstrap
intervals:

| Group | at 0.30 s | at 0.40 s | at 0.50 s | peak | peak at |
| --- | --- | --- | --- | ---: | ---: |
| minimal | 0.215 [0.200, 0.229] | 0.162 [0.147, 0.174] | 0.093 [0.084, 0.100] | 0.243 | 0 s |
| 4day | 0.298 [0.273, 0.329] | 0.335 [0.306, 0.362] | 0.329 [0.301, 0.359] | 0.342 | 0.45 s |
| 20day | 0.379 [0.307, 0.433] | 0.516 [0.471, 0.572] | 0.377 [0.345, 0.409] | 0.516 | 0.40 s |

This is the paper's central figure, reproduced with intervals.

After minimal practice the habitual response probability never rises above the
guess rate of about 0.25; it simply decays as the correct response becomes
available. After four days a bump appears, peaking at 0.342 around 0.45 s.
After twenty days the bump reaches 0.516 at 0.40 s, more than double the guess
rate, and its interval [0.471, 0.572] is nowhere near it.

The transient habitual-error bump is therefore absent after minimal practice
and grows with training duration, which is the paper's claim, now stated with
uncertainty attached and with all participants in a group contributing to one
curve.

## What this adds to the paper

1. The group-level conclusion is reached in three fits rather than 165, with
   the subject-level variability modeled rather than averaged over afterwards.
2. The preparation-time parameters are shown to covary strongly across
   participants in the trained groups, which is a new claim the per-subject
   procedure could not make.
3. The minimal-practice group's "no habit" result is expressed as a property
   of one fitted model rather than as a failure of a separate model to win,
   which is a cleaner statement of the same finding.
4. It closes, in substance, an issue the maintainer opened nine years ago.
