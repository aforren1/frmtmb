# The response-preparation ("habit") model of Hardwick et al. (2019)

Working notes for a frmtmb replication of the model in:

> Hardwick RM, Forrence AD, Krakauer JW, Haith AM (2019). Time-dependent
> competition between goal-directed and habitual response preparation.
> *Nature Human Behaviour* 3, 1252-1262.

Reference implementation: <https://github.com/adrianhaith/HabitTR>

## Licensing status: NOT REDISTRIBUTABLE PENDING PERMISSION

The reference repository has **no license file**. Neither its MATLAB code nor
its data may enter the frmtmb tree. This directory therefore contains only:

* original R code written for this replication, and
* numeric result tables derived from running that code.

No trial-level data and no MATLAB source is committed here. `replicate.R`
clones the reference repository into a scratch directory at run time and
derives everything from that clone.

Two sources exist, both under the same unlicensed upstream repository:

| Source | Clone | Commit |
| --- | --- | --- |
| Upstream MATLAB and `HabitData.mat` | `adrianhaith/HabitTR` (master) | `acaaa2b` |
| Cleaned per-trial CSVs | `aforren1/HabitTR` (branch `clean_data`, PR #3) | `06bd3ab` |

The cleaned CSVs are the frmtmb maintainer's own preprocessing of a study the
maintainer co-authored. That makes the permission question near-trivial for
the maintainer to settle, but it is not settled yet, so the stance above holds
and the vignette question stays deferred.

## 1. The model

### Response categories

Each remapped-stimulus trial yields one of three categories
(`habit_lik.m:5-7`):

| Code | Meaning |
| --- | --- |
| 1 | correct, that is the revised (new) mapping, "mapping B" |
| 2 | habitual, that is the original (old) mapping, "mapping A" |
| 3 | any other key |

There are four response keys in total, so category 3 covers two keys.

### Preparation-time processes

Two independent processes race to be *prepared* by the time of the forced
response. Each has a Gaussian preparation-time distribution, so the
probability that it is ready at prep time \(t\) is a normal CDF
(`getResponseProbs.m:8-9`):

\[
\Phi_A(t) = \Phi\!\left(\frac{t - \mu_A}{\sigma_A}\right), \qquad
\Phi_B(t) = \Phi\!\left(\frac{t - \mu_B}{\sigma_B}\right)
\]

Process A is the habitual (original) mapping, process B the goal-directed
(revised) mapping. Independence gives four mutually exclusive preparation
states, whose probabilities form the columns of

\[
\Phi(t) = \bigl[(1-\Phi_A)(1-\Phi_B),\;\; \Phi_A(1-\Phi_B),\;\;
(1-\Phi_A)\Phi_B,\;\; \Phi_A\Phi_B\bigr]^{\!\top}
\]

(`getResponseProbs.m:50`). Read the states as: neither prepared, A only,
B only, both.

### Category probabilities

A state-by-category coefficient matrix \(\alpha\) maps preparation states to
response probabilities, and the model is the plain matrix product
(`getResponseProbs.m:51`):

\[
p(t) = \alpha \, \Phi(t)
\]

The free quantities in \(\alpha\) are the asymptotic accuracy of each process,
\(q_A\) and \(q_B\), and the guess rate \(q_I\) (`initAE` in their code).
\(q_I\) is the probability of hitting the correct key, and separately of
hitting the habitual key, when nothing is prepared; the two remaining keys
split the rest, hence the \(0.5 - q_I\) entries.

**`habit`** (`getResponseProbs.m:22-27`):

| category | neither | A only | B only | both |
| --- | --- | --- | --- | --- |
| 1 correct | \(q_I\) | \((1-q_A)/3\) | \(q_B\) | \(q_B\) |
| 2 habit | \(q_I\) | \(q_A\) | \((1-q_B)/3\) | \((1-q_B)/3\) |
| 3 other | \(0.5-q_I\) | \((1-q_A)/3\) | \((1-q_B)/3\) | \((1-q_B)/3\) |

Once B is prepared it wins outright, whether or not A is also prepared: the
"B only" and "both" columns are identical. That is the whole habit story. If A
alone is prepared the participant emits the habitual response with probability
\(q_A\), which is what produces the transient error bump at intermediate prep
times.

**`no-habit`** (`getResponseProbs.m:29-34`) is the same table with the "A
only" column replaced by a copy of the "neither" column. Preparing A alone
then does nothing, \(\Phi_A\) drops out algebraically, and the model collapses
to a single process: baseline guessing until B is ready.

**`flex-habit`** (`getResponseProbs.m:36-42`) mixes the two "A only" columns
with weight \(\rho = \) `params(8)`, the habit-strength parameter:

\[
\alpha^{\text{flex}}_{\cdot,2} = \rho\,\alpha^{\text{habit}}_{\cdot,2}
  + (1-\rho)\,\alpha^{\text{no-habit}}_{\cdot,2}
\]

So \(\rho = 1\) recovers `habit`, \(\rho = 0\) recovers `no-habit`, and
`no-habit` is nested in `flex-habit` at a parameter-space boundary. That
boundary matters for the likelihood-ratio test, see below.

Rows 4 and 5 of \(\alpha\) describe unchanged (non-remapped) stimuli. They are
never reached by the likelihood, which indexes rows 1 to 3 only
(`habit_lik.m:16-18`), and exist for plotting.

### Parameter vector

`params` is length 8, ordered

\[
[\;\mu_A,\; \sigma_A,\; q_A,\; \mu_B,\; \sigma_B,\; q_B,\; q_I,\; \rho\;]
\]

The header comment at `habit_lik.m:10` claims `[sigmaA muA qA ...]`. That
comment is **stale and wrong**; the code at `getResponseProbs.m:8` passes
`paramsA(1)` as the mean of `normcdf`. The ordering above is confirmed
empirically in section 3.

### Objective

The fitted objective is not the log-likelihood. It carries a ridge penalty
pulling both preparation-time standard deviations toward 0.07 s
(`habit_lik.m:45-47`):

\[
\text{nLL} = -\sum_i \log p_{r_i}(t_i)
  \; + \; 1000\,(\sigma_A - 0.07)^2 \; + \; 1000\,(\sigma_B - 0.07)^2
\]

`LLopt` in the stored fits is this penalized quantity; `LLactual`
(`fit_habit_model.m:102`) is the clean, unpenalized log-likelihood, and is the
one used for AIC.

### An unnormalized likelihood

Category 3 covers two keys, but its row of \(\alpha\) gives the probability of
**one** such key. Column sums of rows 1 to 3 are therefore not 1; they become
1 only when row 3 is counted twice. Their own sliding-window code halves the
observed "other" rate to match (`preprocess_data.m:72`,
`habit_lik.m:70`), so the convention is deliberate and internally consistent
for plotting.

The likelihood, however, indexes row 3 directly (`habit_lik.m:16-20`). Each
category-3 trial therefore contributes \(\log(p_3/2)\) rather than
\(\log p_3\), and the reported log-likelihood sits below the properly
normalized one by exactly \(n_3 \log 2\), where \(n_3\) is that subject's
count of category-3 trials. Because \(n_3\) is fixed for a given dataset, the
offset is identical across the models being compared, so every AIC difference
and likelihood-ratio statistic they report is unaffected. It matters here only
because a normalized categorical family in frmtmb will not reproduce their
absolute log-likelihood without adding this offset back. See section 3.

## 2. Design, data, and fitting

### Conditions

Three groups, indexed `c = 1,2,3` in their code:

| `c` | Name in code | Paper | Source |
| --- | --- | --- | --- |
| 1 | `minimal` | minimal practice | experiment 1, untrained |
| 2 | `4day` | 4 days of practice | experiment 1, trained |
| 3 | `4week` | 20 days of practice | experiment 2 |

Conditions 1 and 2 are the **same people** measured before and after training
(`preprocess_data.m:19-27`, and PR #3: subjects in the 1-day and 5-day
experiments share an ID). Condition 3 is a separate cohort, given IDs offset
by 100 in the cleaned CSVs. Their per-subject fitting procedure ignores this
pairing; a hierarchical model does not have to.

### Session and file structure

From issue #1 ("File IDs & notes", aforren1, 2017-06-07):

*Experiment 1, untrained (1-day)*

| Session | Content |
| --- | --- |
| 001 | hand practice (1-200 free RT, 201-400 forced) |
| 002 | criterion assessment, original mapping; one block free RT |
| 003 | criterion assessment, revised mapping; one block free RT |
| 004 | timed-response assessment (100 trials per block) |

*Experiment 1, trained (5-day)*

| Session | Content |
| --- | --- |
| 001 | 200 hand practice, 200 hand forced, 401+ symbol practice |
| 002-004 | symbol practice |
| 005-008 | replicate the 1-day protocol 001-004 after training |

*Experiment 2 (4 week)*

| Session | Content |
| --- | --- |
| 001 | hand practice |
| 002 | hand forced |
| 003 | criterion, original |
| 004 | symbol forced |
| 005-008, 010-013, 015-018, 020-023 | symbol practice |
| 009, 014, 019, 024 | symbol forced |
| 025 | criterion revised |
| 026 | 500 trials forced to revised, rest practice |

The timed-response assessment that the model is fit to is therefore session
004 (1-day), 008 (5-day), and 026 (20-day).

### Known data corrections

Issue #1 records that for subject 4 of experiment 2 (ID 104 in the cleaned
CSVs) the screen was minimized on three trials, which must be set to `NA`:

| Session | Trial |
| --- | --- |
| 010 | 1073 |
| 015 | 790 |
| 026 | 758 |

`replicate.R` applies these. Note that only the session 015 and session 026
rows are present in the cleaned CSVs; the session 010 row is already absent.

### Raw column layout

Issue #1 gives the 13 raw columns: `Subject`, `Session`, `Stamp`, `Trial`,
`DesOnset`, `RecOnset`, `Stimulus`, `TargetButton`, `ResponseNum`,
`ResponseButton`, `RT`, `Correct`, `Tally`. The cleaned CSVs of PR #3 drop
`Stamp`, `DesOnset`, `RecOnset` and `RT`, and add `IsFreeResp`, `IsHand`,
`IsToCriterion`, `PrepTime`, `IsCatchTrial`, `WillRemap`, `IsRemapped`,
`IsRemappedSession`, `OldButton` and `NewButton`.

### Cleaning rules for the cleaned CSVs

From PR #3 and the maintainer's own `adrian_doc.Rmd:23-24`, filter out
mishits (`ResponseButton > 4`), user errors (`PrepTime > 5`) and coding
errors (`PrepTime` far below 0, present as a `-999.98` sentinel). Logicals are
exported as 0/1 integers. Only the initial response of each trial is
retained. Recoding of a remapped trial follows `adrian_doc_2.Rmd:27-32`:
`ResponseButton == NewButton` is correct, `== OldButton` is habitual,
otherwise other; `!WillRemap` marks the unchanged control stimuli.

Note that `clean.r` and `misc.r` are `source()`d by those documents but were
never committed to the branch, so the loader in `replicate.R` is written from
scratch against the documented rules.

### Subject exclusions

Two filters, applied in sequence:

1. `preprocess_data.m:81-86` blanks a subject's remapped `RT` vector when
   accuracy for prep times above 0.8 s falls below 0.7. The unchanged vector
   survives, which is why some subjects have unchanged data but no remapped
   data.
2. `HabitTR_analysis2.m:5-26` documents a further pass on asymptotic correct
   below 0.7 or asymptotic habit above 0.2, evaluated above a 0.9 s threshold.
   This one only marks points on a figure; it does not remove anyone from the
   stored fits.

Subjects surviving filter 1 in `HabitData.mat`, per condition: 21, 21, 13.

### Fitting procedure

Per subject, per condition, per model, by maximum penalized likelihood
(`fit_habit_model.m:39-74`). Both a `fmincon` and a `bads` (Bayesian Adaptive
Direct Search) run are stored, in `HabitModelFits.mat` and
`HabitModelFits_bads.mat` respectively.

Start values `[.4 .05 .99 .5 .05 .95 .25 .95]` (`fit_habit_model.m:24`) and
box bounds (`fit_habit_model.m:13-21`):

| Parameter | LB | UB |
| --- | --- | --- |
| \(\mu_A\) | 0 | 0.75 |
| \(\sigma_A\) | 0.01 | 100 |
| \(q_A\) | 0.99 | 0.999 |
| \(\mu_B\) | 0 | 10 |
| \(\sigma_B\) | 0.01 | 100 |
| \(q_B\) | 0.5 | 0.9999 |
| \(q_I\) | 0.0001 | 0.499 |
| \(\rho\) | 0.0001 | 0.9999 |

\(q_A\) is deliberately pinned into a narrow band just below 1
(`fit_habit_model.m:20-21`), because it trades off against \(\rho\) and is not
separately identified. There is also an inequality constraint
\(\mu_A \le \mu_B\) (`fit_habit_model.m:26-27`) and, in some code paths, an
equality constraint \(q_A = 0.99\) (`fit_habit_model.m:50-51`).

The stored \(q_A\) values span 0.99 to 0.999 in both fit files, so the
equality constraint was **not** live for the archived runs; the box bound was.
`fit_model_constrained.m` is a near-duplicate driver that comments out the
\(q_A\) box lock (`fit_model_constrained.m:25-27`) and writes to the same two
filenames, so the two drivers cannot both have produced the archive. The
\(q_A\) spread identifies `fit_habit_model.m` as the source.

Model 4 is a separate `no-habit` fit to the unchanged stimuli
(`fit_habit_model.m:83-86`), used as a skill reference rather than in the
habit comparison.

### Model comparison as published

Parameter counts (`fit_habit_model.m:98-100`): 4 for `no-habit`
(\(\mu_B, \sigma_B, q_B, q_I\)), 6 for `habit` (adding \(\mu_A, \sigma_A\)),
7 for `flex-habit` (adding \(\rho\)). AIC uses the unpenalized `LLactual`
(`fit_habit_model.m:102-103`). The headline claim is that
\(\Delta\text{AIC}\) favors `habit` over `no-habit` in the 4-day and 20-day
groups but not after minimal practice.

## 3. Validation of the R reference

`replicate.R` contains a hand-written R translation of `habit_lik.m` and
`getResponseProbs.m`, with no frmtmb involvement, following the house pattern
of an independent reference likelihood in the comparison script.

Evaluated at the authors' own stored `paramsOpt`, for every model, condition
and included subject, it reproduces their stored `LLactual` to machine
precision:

| Fit file | comparisons | max abs difference |
| --- | --- | --- |
| `HabitModelFits.mat` (`fmincon`) | 165 | 4.3e-13 |
| `HabitModelFits_bads.mat` (BADS) | 165 | 4.5e-13 |

The ridge penalty is confirmed separately: for every stored fit,
`LLopt + LLactual` equals \(1000[(\sigma_A - 0.07)^2 + (\sigma_B - 0.07)^2]\).

Both checks would fail under any other permutation of the parameter vector, so
they also settle the ordering question raised by the stale comment at
`habit_lik.m:10`.

## 4. Finding: value-matched trial subsetting inflates the fitted datasets

`preprocess_data.m:59-66` splits the recoded trials into remapped and
unchanged sets like this:

```matlab
revised_trials = ismember(recodedX,revisedX);
recodedXr = recodedX(revised_trials)';
recodedYr = recodedY(revised_trials)';

non_revised_trials = ismember(recodedX,unchangedX);
recodedXnr = recodedX(non_revised_trials)';
recodedYnr = recodedY(non_revised_trials)';
```

`recodedX` is a vector of prep times and `revisedX` is the prep times of the
remapped trials. `ismember` therefore selects **by value, not by trial
index**. Prep times are rounded to the millisecond and several hundred of them
fall in a range of about 1.3 s, so collisions between the two sets are common.
A trial whose prep time happens to occur anywhere in the other stimulus class
is admitted to both datasets.

Replaying that `ismember` against the raw `tmp.mat` structures reproduces the
length of `data(subject,c).RT` **exactly, for all 58 subject-conditions**,
which confirms the mechanism rather than merely being consistent with it. The
full audit is in `table-subsetting-audit.csv`:

| Group | subjects | trials run | true remapped | remapped as fit | double-counted | share of fitted set |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| minimal | 22 | 493 to 500 | 250 | 265 to 281 | 1052 | 17.45% |
| 4day | 22 | 492 to 500 | 250 | 263 to 285 | 1052 | 17.44% |
| 20day | 14 | 498 to 500 | 250 | 269 to 281 | 693 | 18.00% |

That is 2797 double-counted trials out of 28940 run. No trial is lost: every
trial lands in at least one dataset, and about a sixth land in both.

Because roughly half of the wrongly admitted trials are genuinely unchanged
trials, which carry no habitual response option at all, the contamination
dilutes the habit signal. The direction of the bias is therefore **against**
the paper's own conclusion, which makes the published result conservative
rather than inflated.

Refitting from identity-based subsetting off the cleaned CSVs confirms the
predicted direction. Mean per-subject ΔAIC favoring habit:

| Group | value-matched, as published | index-based, corrected | change |
| --- | ---: | ---: | ---: |
| minimal | -4.80 | -5.87 | -1.06 |
| 4day | +8.21 | +8.49 | +0.27 |
| 20day | +19.94 | +20.58 | +0.64 |

The groups with habit show more of it, and the group without shows less. The
defect was working against the authors.

This is a preprocessing defect, not a modeling one. The likelihood itself is
sound.

## 5. Mapping the model onto frmtmb v0.47.0

### Two refusals shape the spelling

The shape the brief anticipated, a categorical or Bernoulli family with an
identity link and `nl = TRUE` bodies, is not available at v0.47.0. Both halves
are refused, for separate reasons:

| Wanted | Refused by | Why |
| --- | --- | --- |
| `bf(..., nl = TRUE) + categorical()` | `R/parse.R:1225-1229` | `nl = TRUE` requires a family whose `primary_dpars` is exactly `"mu"`. `categorical()` has one predictor per non-reference category (`R/families.R:3457`), so the gate rejects it. Registered at `R/compat.R:609`. |
| `categorical(link = "identity")` | `R/families.R:4297` | The constructor accepts `"logit"` only. The softmax is built internally with `RTMB::logspace_add` (`R/families.R:3359`), so there is no hook for supplying probabilities. |

`bernoulli("identity")` does exist and does take a probability directly
(`R/families.R:1200-1204`), but the response here has three categories, so it
does not apply.

### The route that works

`nlf()` names a distributional parameter and therefore bypasses the single-`mu`
gate (`R/parse.R:1211-1232`). Combined with a change of variable it reaches the
model exactly. If the model yields \(c_1, c_2, c_3\), then setting

\[
\texttt{muhabit} = \log(c_2/c_1), \qquad
\texttt{muother} = \log(2 c_3 / c_1)
\]

makes the internal softmax return \((c_1, c_2, 2c_3)\), which is the properly
normalized version of the reference model's categories. The softmax is
invertible, so hard-wiring it costs nothing as long as its inverse can be
written in a body.

The factor of two restores the second "other" key that the reference likelihood
drops. Consequently

\[
\ell_{\text{frmtmb}} = \ell_{\text{reference}} + n_3 \log 2
\]

exactly, where \(n_3\) counts category-3 trials. Verified at the authors' own
archived estimates for condition 1 subject 1: the difference from that
prediction is 0.

To my knowledge this is the first use of `nlf()` on a categorical dpar in this
codebase; no test, vignette or man page exercises it. It works, and this
directory is evidence, but it is an unpinned code path.

### Constraints and the penalty

| Reference construct | frmtmb equivalent |
| --- | --- |
| `fmincon` box bounds (`fit_habit_model.m:13-21`) | `set_prior("", nlpar = ..., lb =, ub =)`, one specification per parameter |
| ridge `1000*(sigma - 0.07)^2` (`habit_lik.m:47`) | `set_prior("normal(0.07, s)", nlpar = ...)` with \(1/(2s^2) = 1000\), so \(s = 1/\sqrt{2000}\) |
| linear inequality \(\mu_A \le \mu_B\) (`fit_habit_model.m:26-27`) | no box equivalent; attempted as \(\mu_A = \mu_B - e^{\delta}\), which does not survive a random effect (see `hierarchical-results.md`) |

### Finding: prior-carried bounds cannot address a nonlinear parameter

**Resolved in v0.49.** The bound is now keyed by the template position's name
rather than by the design-matrix column, so a bound addressed by `nlpar`,
`dpar` or `resp` lands where the distribution beside it lands. `replicate.R`
sets `BOUNDS_VIA_PRIOR <- TRUE`, and the record below is kept as the
diagnosis. `frm(lower =, upper =)` was then removed outright in the same
release, once `thetaac_*` and `thetar_*` gained `set_prior()` classes of their
own (`"ar"`, `"ma"`, `"cosy"`, `"cortime"`, `"rescor"`, and `class = "theta"`
with `coef = "thetaac_1"` for one raw internal parameter).

House style is to spell bounds through the prior vocabulary rather than
`frm(lower =, upper =)`. For a nonlinear parameter that was not possible
before v0.49.

```r
set_prior("", nlpar = "guess", lb = 0, ub = 1)
#> Error: Unknown parameter(s) in bounds: (Intercept).
#>   Available: eta_ability, guess_(Intercept)
```

Reproduced on the package's own documented nonlinear example (the brms
guessing-floor model from `dev/brms-vignettes/brms_nonlinear.R`), so it is not
specific to this case study. The four relevant spellings:

| Spelling | Result |
| --- | --- |
| `frm(lower = c(guess = 0), upper = c(guess = 1))` | worked then; the arguments were removed in 0.49 |
| `set_prior("", nlpar = "guess", lb = 0, ub = 1)` | **refused** |
| `set_prior("normal(0.33, 0.05)", nlpar = "guess")` | works |
| `set_prior("", class = "b", coef = "x", lb = 2.5)` on a linear model | works |

Cause: `R/priors.R:861-862` key the bound as `lower[tg$name]`, where `tg$name`
is `cn[k]` (`R/priors.R:831`), the design-matrix column name of the nonlinear
parameter's own sub-formula. For `guess ~ 1` that is `"(Intercept)"`.
`resolve_bounds()` (`R/priors.R:1297-1322`) matches names against
`outer_par_names()`, which spells the same parameter `"guess_(Intercept)"`, so
the lookup fails. The distribution path is unaffected because it keys on the
index-based `nm_of(tg$comp, tg$idx)`, which is why a `normal()` prior placed
through `nlpar` works while a bound placed the same way does not. A model with
several nonlinear parameters has the further problem that every bound collides
on the single key `"(Intercept)"`.

Consequence for the plan to retire `frm(lower =, upper =)`: those arguments
were the **only** way to bound a nonlinear parameter, so they could not be
removed until `tg$name` was qualified with the nonlinear parameter's name on
the bounds path. That is the v0.49 fix, and it is why this stage now writes
its bounds as priors. The arguments themselves went in the same release.

Separately, a transform is still not a substitute in this stage, and `lb`/`ub`
is not one either: it is a hard box, the same box the argument built. The
point of the per-subject fits is to reproduce
a box-constrained `fmincon` run whose optima sit **on** the bounds: \(q_A\) is
at its limit for the large majority of participants, by construction, because
it is not separately identified from \(\rho\). A logit transform converts an
attainable boundary into an asymptote. Hard boxes are the right tool for
reproducing hard boxes, and `replicate.R` uses them here, written in the prior
vocabulary, with the reasoning recorded at the point of use. The hierarchical
model, where nothing sits on a bound and the random effects need an unbounded
scale, uses transforms.

A caution that cost some time: `logLik()` on a fit carrying a prior returns
the log-likelihood **plus the log prior density**, normalizing constants
included. It is therefore not the reference implementation's `LLactual`. The
comparable number is obtained by evaluating the hand-written reference at
frmtmb's own estimates, which is what `frm_LL()` in `replicate.R` does.

### Per-subject agreement

Every one of the 165 per-subject fits converged. Agreement with the archived
`fmincon` estimates, measured as the reference log-likelihood at each side's
own estimates:

| Model | n | max abs difference | within 1e-4 |
| --- | --- | --- | --- |
| `no-habit` | 55 | 0 | 55 |
| `habit` | 55 | 0.0027 | 16 |
| `flex-habit` | 55 | 2.28 | 10 |

For `habit`, every parameter except \(q_A\) agrees to better than 1e-3 for
every participant. \(q_A\) differs by up to 0.008, which is most of the width
of its own permitted band, and is a flat direction: it is confined to
\([0.99, 0.999]\) precisely because it is not separately identified from
\(\rho\).

The `flex-habit` disagreements are real but are optimizer basins, not a porting
error. Restarting frmtmb from the archived estimates and running
`frm_allfit()` recovers an optimum at least as good as the archive for nine of
the ten worst cases; the tenth is worse by 0.0013. On one participant frmtmb
finds a solution 2.28 log-likelihood units better than the published fit.

### Nesting, and what `anova()` may legitimately be asked

* `flex-habit` contains `habit` at \(\rho = 1\), a boundary.
* `flex-habit` contains `no-habit` at \(\rho = 0\), also a boundary, and
  \(\mu_A\) and \(\sigma_A\) become unidentified there, so the
  likelihood-ratio statistic has no standard null distribution.
* `habit` and `no-habit` are **not** nested. Collapsing the "A only" column of
  the habit model onto the "neither" column needs \(q_A = q_I = 0.25\), and
  \(q_A\) is bounded at or above 0.99.

The paper handles this correctly: AIC for habit against no-habit, and a
likelihood-ratio test only for habit against flex-habit (`LR_test.m:6-13`,
which also clamps the statistic at zero in acknowledgement of the boundary).

Worth recording as a package observation: `anova.frmtmb_fit` requires only that
the fits share `n_obs` (`R/confint.R:1231-1237`). It does not and cannot verify
nesting, so a user can obtain a confident-looking p-value for a non-nested
pair. That is standard for the method, and the documentation already carries a
boundary caveat at `R/confint.R:1141-1148`, but a case study is a good place to
say it out loud.

## 6. Cross-check of the two data sources

The cleaned CSVs and `HabitData.mat` describe the same measurements but do not
agree trial for trial:

* `modelCoded.x` in the upstream `tmp.mat` has the same length as the cleaned
  CSV's session block (499 for condition 1 subject 1), so the two **raw**
  sources agree.
* `HabitData.mat` disagrees with both, for the reason given in section 4.
* The cleaned loader drops `PrepTime <= 0`, whereas the upstream pipeline
  retains small negative prep times (down to about -0.9 s in condition 1).
  These are responses issued before stimulus onset. Retaining them is
  defensible, since the model's normal CDFs are defined there and such trials
  are informative about the guess rate, but the two sources differ on it.

For the replication proper, `HabitData.mat` is the required input: it is what
the archived fits were computed from, and any per-subject comparison against
those fits must use the same rows. The cleaned CSVs are used for the corrected
refit of section 4 and are the better base for a future shipped vignette.

## 7. Deliverables and what does not ship

Files in this directory:

| File | Contents |
| --- | --- |
| `model-notes.md` | this file: the model, licensing, reference map, findings |
| `replicate.R` | clone to comparison, runnable end to end |
| `hierarchical-results.md` | the hierarchical fits and the model-comparison story |
| `vignette-draft.md` | case-study narrative, for the day licensing is settled |
| `table-reference-validation.csv` | R reference against archived MATLAB likelihoods |
| `table-persubject.csv` | per subject, per model: frmtmb against MATLAB |
| `table-pooled-summary.csv` | the same, pooled by model and group |
| `table-model-comparison.csv` | habit against no-habit and flex-habit, per subject |
| `table-corrected-subsetting.csv` | refit from identity-based trial subsetting |
| `table-subsetting-audit.csv` | the `ismember` audit, all 58 subject-conditions |
| `table-bootstrap-bands.csv` | bootstrap bands on the fitted curves |

**Draft NEWS bullet: none.** Nothing ships this round. This is a development
case study that produced findings, not a package change. The two package-facing
findings, the nonlinear-parameter bounds addressing gap in section 5 and the
absence of a nesting check in `anova()`, are reported for the maintainer to act
on rather than fixed here, because this lane may not edit `R/`.

**No data is committed.** The brief permitted a tidy per-trial `.rds` in this
directory provided the licensing note flagged it. The note does flag it, but
the standing instruction is that nothing from the unlicensed upstream
repository enters the frmtmb tree pending permission, so the per-trial data
stays in the scratch directory and only the code that regenerates it from a
clone is committed here.
