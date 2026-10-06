# Lane ordmix: ordinal mixtures; hurdle_cumulative() thres(gr = ) and cs()

Base: 9e902909 (frmtmb 0.67.0), branch `wt-ordmix`, uncommitted. Base
build for "before" measurements: `rellib-r5`. Lane library:
`C:/Users/adf44/source/r/wt-ordmix-lib`. Reference: brms 2.23.0, run
(not recalled) through `stancode()`, `standata()`, `default_prior()`,
its compiled programs and its R-side densities. Scripts and logs are
`dev/ordmix-*`; logs are local, and every number below is pasted from
the log named beside it.

## Punch round 2b (the final check's B5 and c1)

Logs `dev/ordmix-p2b-*`; driver `dev/ordmix-p2b-final.sh`, run on the
final build after the last R edit and install.

**B5. Weights.** `mixture_ord_degeneracy()` summed the case weights as
given, so the cut on the sharpening (an absolute 0.1) moved with their
scale. The weights now enter divided by their mean. Tests ("constant
weights do not change the degenerate verdict", `test-ordinal-mixture.R`):
- a sound fit with weights 1/n (the review's `w_sum1`, seed 3) does
  not warn;
- weights 1 and 10 give the same verdict on a degenerate fit (seed 19,
  n = 300) and a sound one (seed 3);
- at one fit's parameters, the sharpening with its weights scaled by
  10 and by 1/300 equals the sharpening at 1. The two optimizer paths
  for weights 1 and 10 end within 1e-5 of each other, not at the same
  doubles, so that comparison is made at one fit.

Seen to fail on a copy of the final source without the normalization
(`dev/ordmix-suite-p2b-mut/`: 5 failures, lines 707, 731 and 732), and
on the round-2 build (`dev/ordmix-suite-p2b-before/`).

The review's driver `dev/ordmix-rev3-b3.R` (`dev/ordmix-p2b-log/`;
before: the review's `dev/ordmix-rev3-b3-log/`):

    design    fits  sound  warn before  warn after  sound fits warning after
    w_sum1     20    11        20           1               0
    w_tenth    20    12         0           0               0

The one `w_sum1` fit that warns (seed 8) is not sound: its standard
errors are NaN, and it warns on a collapsed gap (its sharpening is
-24.7 and -7.1). The sharpenings for weights 1/n and 0.1 now agree
to within 2.2% on every seed; they differ only where the two fits end a
little apart (seed 19: -20.9 and -10.1, against -20.9 and -9.89).

The review's driver `dev/ordmix-rev3-wdegen.R`
(`dev/ordmix-p2b-log/wdegen.txt`), the review's degenerate fits with
weights 1, 10 and 1/n:

    seed 19 : w=1 warn=TRUE sharpen=8.12e-07,-33 ll=-389.0058 | w=10 warn=TRUE sharpen=6.64e-07,-33 ll=-389.0058 | w=0.003333 warn=TRUE sharpen=3.17e-06,-33 ll=-389.0058
    seed 33 : w=1 warn=TRUE sharpen=3.39e-06,-32.3 ll=-389.5565 | w=10 warn=TRUE sharpen=2.37e-06,-32.3 ll=-389.5565 | w=0.003333 warn=TRUE sharpen=2.52e-06,-32.3 ll=-389.5565
    seed 35 : w=1 warn=TRUE sharpen=-39.5,2.19e-05 ll=-384.5550 | w=10 warn=TRUE sharpen=-39.5,9.06e-06 ll=-383.8555 | w=0.003333 warn=TRUE sharpen=-40.7,8.29e-06 ll=-384.9821
    seed 38 : w=1 warn=TRUE sharpen=-30.1,-0.185 ll=-398.7081 | w=10 warn=TRUE sharpen=-30.1,-0.185 ll=-398.7081 | w=0.003333 warn=TRUE sharpen=-31.8,2.99e-06 ll=-398.3846
    seed 40 : w=1 warn=TRUE sharpen=1.17e-05,-26.7 ll=-386.3642 | w=10 warn=TRUE sharpen=2e-05,-26.7 ll=-386.3642 | w=0.003333 warn=FALSE sharpen=-0.392,-26.4 ll=-390.1929

Every fit with weights 1 and 10 warns, and the sharpening no longer
scales with the weights. Seed 40 with weights 1/n is a different fit
(log-likelihood per unit weight -390.19 against -386.36): the
optimizer, not the check, ends elsewhere at that scale. Its sharpening
is -0.39 there, and it does not warn. Seeds 35 and 38 also end
elsewhere at 1/n, and they still warn.

**c1. The identified count of a hurdle mixture with `thres(gr = )`**:
one zero proportion for all groups, not one per group, since `hu<k>`
and `theta<k>` are shared (`mixture_ord_ident_check()`, `n_id <-
sum(nth) + isTRUE(extra_cat)`). The review's case is now 6, its
Hessian rank. Added to "with no predictor, priors are counted against
the proportions": priors on mu1's 5 thresholds and on hu1 leave 7
free, so the fit warns. Round 2 counted 7 and stayed silent: seen to
fail on the round-2 build (`dev/ordmix-suite-p2b-before/`, line 680).

**Nothing else moved** (`dev/ordmix-p2-crit.sh` rerun, summary
`dev/ordmix-p2b-log/crit-sum.txt`): the installed build's warning
agrees with the rule on all 537 fits. Margin: 0 of 167 sound and 0 of
43 poor fits warn, 7 of 8 degenerate fits warn, and 0 of 19 probit
saturation fits warn. Review's 160: 20 of 22 and 0 of 138. Identified
140: 3 of 8 and 0 of 132. These are round 2's numbers. Unweighted
fits have weights 1 with mean 1, so their sharpening is unchanged.

**Tests.** `test-ordinal-mixture.R` pass 100, `test-hurdle-cum-thres-cs.R`
pass 47, both with fail 0, error 0 and warning 0 (`dev/ordmix-suite-p2ba/`).
Core suite, one file per process, every gate on
(`dev/ordmix-suite-p2bcore/`): 206 files. 205 have fail 0, error 0 and
warning 0. The other, `test-ordinal-mixture.R`, failed twice there on
the first form of the constant-weights test, which compared two fits
to the default relative tolerance. That test now compares at one fit,
and the file's rerun above is clean. The core total is pass 16670
(16666 - 96 + 100), fail 0. No documentation changed (the edits are
`@noRd` comments; `man/` is as in round 2), so R CMD check was not
rerun.

Functions touched: R/families.R `mixture_ord_degeneracy()` (the
weights), `mixture_ord_ident_check()` (the count); tests in
`test-ordinal-mixture.R`.

## Punch round 2 (the review's "Re-check after punch round 1")

Scripts and logs are `dev/ordmix-p2-*`; the review's scripts
`dev/ordmix-rev2-*` were rerun unchanged, into `dev/ordmix-p2-*-log/`.

### B3. The degenerate warning on sound strong-predictor fits (fixed)

The round-1 cut, a latent distance above 50, is replaced
(`mixture_ord_degenerate_check()`, `mixture_ord_degeneracy()` in
R/families.R). A component warns when:

- **it is a step function**: doubling all its latent distances (its
  discrimination times 2) costs the mixture's log-likelihood less than
  0.1 (`sharpen > -0.1`). With `cs()`, each threshold's distances are
  also doubled alone, through the cs slot (`tau_j - cs'_ij - eta_i =
  2 (tau_j - cs_ij - eta_i)`), because one threshold's own slope can be
  the step function while the rest are placed;
- **a threshold has run off**: no row lies within a latent distance of
  50 of it (`near > 50`), thresholds above the data left out; or
- **a gap is collapsed**: two ordered thresholds are the same double.

The sharpening test needs no link-specific scale and no cut on the
size of the predictors: a sound component loses likelihood when made
sharper, and one that has run off does not. The reach is still
reported, and no longer decides. The doc comment's "0 or 1 for every
link" is gone (the cauchit's F(50) is 0.9936).

Candidates measured (`dev/ordmix-p2-crit.R`, driver
`dev/ordmix-p2-crit.sh`, summary `dev/ordmix-p2-crit-sum.R`, log
`dev/ordmix-p2-log-crit-sum.txt`; the round-1 split of the margin
degenerate fits from `dev/ordmix-p2-r1split.R`): the fraction of rows
a component predicts as a point mass, the sharpening, and the nearest-row
distance. The point-mass fraction misses every cauchit step function
(its tails never reach 1 - 1e-6). The rule taken, by population, the
review's label as the truth (margin designs: sound = all standard
errors finite and both slopes within 50%; poor = finite but a slope
off; the 300: `|estimate| > 30` or a non-finite standard error):

    population                      round 1 (reach > 50)   round 2
    margin, 12 designs, 237 fits
      sound 167                       43 warn               0 warn
      poor 43                         23 warn               0 warn
      degenerate 8 (not probit)        7 of 8               7 of 8
      probit saturation 19             1 warn               0 warn
    review's 160                   20 of 22, 0 of 138    20 of 22, 0 of 138
    identified designs, 140         3 of 8, 0 of 132      3 of 8, 0 of 132
    the 6 suite fires (all degenerate)   6                5

The margin designs are the review's ten (`dev/ordmix-rev2-margin.R`)
plus probit and cloglog at c = 10; three probit c = 10 fits stop at
"NA/NaN gradient evaluation". Among all 480 sound and poor fits of the
three populations the largest `sharpen` is -0.31 and the largest
`near` 20.5. "Probit saturation" is 19 probit fits
whose standard errors are NaN at a latent distance of about 38, where
the probit's log-odds form underflows (dev/test-backlog.md): no
component ran off (`near` at most 0.15), and the sharpened likelihood
is itself NaN there, so the check stays silent. Of the 6 suite fires
of round 1 (`dev/ordmix-p2-log-suitefires.txt`), the `thres(x = )`
fit stays silent as in round 1's final build (its run-off threshold
bounds categories no row takes), and the other 5 warn.

The review's margin script rerun on the final build
(`dev/ordmix-p2-margin-log/`; round 1: the review's
`dev/ordmix-rev2-margin-log/`), sound fits that warn:

    design            sound   round 1   round 2
    logit 1 5000        20        0         0
    logit 5 500         19        0         0
    probit 5 500        17        0         0
    cloglog 1 500       16        0         0
    cloglog 5 500       18        0         0
    logit 10 500        16       13         0
    logit 10 2000       20       20         0
    cauchit 1 500        6        0         0
    cauchit 5 500       12        0         0
    cauchit 10 500      11        3         0

Misses: the 8 of the 300 and the margin set are components that
stopped using an interior category without running off (`near` 9.4 to
33.8, `sharpen` -2.3 to -13.5), the class the review accepted as
recorded.

`cs()` designs, new here because the per-threshold probe is new (20
seeds each; `cs_sr_acat`: `y ~ cs(x)` with `mixture(sratio(),
acat())`, n = 300; `cs_mu2`: `mu2 ~ cs(x)` with `mixture(cumulative(),
sratio())`, n = 400). Under the review's label the warning fires on 1
of 11 and 1 of 7 sound fits; both have a standard error that means
nothing (3.2e4, and 42 on an estimate of -9.6), which the label does
not read. With a label that also calls an SE above 30 degenerate: 0 of
16 sound fits warn, and 22 of 24 degenerate ones do. The two misses
are `cs_sr_acat` seed 1 (an estimate of 31.6 with an SE of 25.7, at
the label's edge) and seed 9 (two thresholds of 31.5 and -31.4 moving
together, SE 3e5: a ridge, not a step; doubling either alone costs 301
and 954, `dev/ordmix-p2-log-csdetail.txt`).

Tests (`test-ordinal-mixture.R`): "strong predictors are not a
degenerate boundary" (logit c = 10 seed 1 silent; a cauchit step
function warns). Seen to fail on the round-1 build
(`dev/ordmix-suite-p2-before/`).

### B4. A threshold prior silenced the flat-direction warnings (fixed)

`mixture_ord_prior_holds()` is replaced by `mixture_ord_prior_held()`
(the template keys a prior holds) and `mixture_ord_n_held()`:

- the hurdle branch counts only the free `hu<k>` and `theta<k>` a
  prior holds, and warns while `n_free - n_held > K`;
- the no-predictor branch warns while the parameters no prior holds
  outnumber the free proportions the data identify. That count is one
  per threshold of each threshold group, plus one per group for a
  hurdle's zero: `ncat - 1` without `thres(gr = )`, as the review
  says, and the sum over groups with it (`ordinal_ncat()` gives the
  largest group's, which undercounts there).

The review's script rerun on the final build (`dev/ordmix-p2-b1-log/`;
before: the review's `dev/ordmix-rev2-b1-log/`). The script's own
`flatwarn` column is FALSE on every hurdle fit on both builds: its
pattern holds "P(Y = 0)" as a regular expression, whose parentheses
are a group, so it never matches the hurdle text. Counted by the
warning text in the logs instead, fits that warn of 10:

    case                                   before   after
    thres_only (prior on mu1's thresholds)    0       10
    nopred_mu1 (y ~ 1, prior on mu1's)        0       10
    none (hurdle, no prior)                  10       10
    hu1_only                                  0        0
    beta11 (beta(1, 1) on hu1 and hu2)        0        0
    disc_differs                              0        0
    disc_same (control)                      10       10

Tests (`test-ordinal-mixture.R`): "a threshold prior does not place
the hu and theta direction" (and its absent case, a prior on hu1) and
"with no predictor, priors are counted against the proportions" (and
its absent case, priors on both components' thresholds). Both seen to
fail on the round-1 build (`dev/ordmix-suite-p2-before/`: fail 3 in
all, one per new test).

The round-1 false-alarm table rerun on the final build
(`dev/ordmix-p1-fa.sh`, `dev/ordmix-p1-fa-log/`; round 1's copy
`dev/ordmix-p1-fa-log-round1/`): the flat-direction counts are round
1's (0 everywhere but the control `mu_same_theta`, 20 of 20). The
degenerate warning fires on 7 fits: `gr_cum2` 4 and `hurdle_prior` 2,
as in round 1, each with a non-finite standard error, and newly
`mu_cum_sr_theta` seed 8, which stopped at singular convergence with a
largest standard error of 1.1e5 and gains 2.3 in log-likelihood when
component 2 is sharpened.

### For the consolidator (record only)

- Lane fixes may add a post-fit Hessian null-direction warning. It
  needs the same guard as `nl_flat_message()`: silent when the
  Hessian, or the gradient, is not finite. Otherwise, on a fit whose
  gradient is not finite, it repeats what `check_convergence()`'s
  non-finite-gradient warning says about the same standard errors.
- The review's rerun list for the merged tree, one file per process:
  - `test-compat.R`: fixes asserts the mixture cells "refused" and
    that `mixture(cumulative(), cumulative())` errors; both must change.
  - `test-ordinal-mixture.R`, `test-hurdle-cum-thres-cs.R`,
    `test-nonfinite-gradient.R`, gated `test-ordinal-mixture-brms.R`.
  - `test-nl.R`, `test-brms-priors.R`, `test-brms-suite-priors.R`,
    `test-brms-formula-priors.R`, `test-prior-compat.R`: the
    per-threshold rows under `dpar = mu<k>`.
  - `test-emmeans.R`, `test-case-studies.R`.
  - frmtmb.sample `test-ordinal-mixture-draws.R`,
    `test-ordinal-draws-names.R`, `test-draws-spellings.R`: the
    `incl_thres` refusal for mixtures.
  - Then all eight suites with `dev/ordmix-rev2-profile.R` as
    `R_PROFILE_USER`, to recount the non-finite-gradient firings and
    to see any fit where both warnings fire.
  - Then the false-alarm drivers (`dev/ordmix-p1-fa.sh`,
    `dev/ordmix-rev2-margin.R`, `dev/ordmix-rev2-b1.R`).

Functions this round touched: R/families.R `mixture_ord_ident_check()`,
`mixture_ord_prior_held()` and `mixture_ord_n_held()` (new, replacing
`mixture_ord_prior_holds()`), `mixture_ord_degenerate_check()`,
`mixture_ord_degeneracy()`, `mixture()` (docs); NEWS.md.

### Tests and checks after the round

All on the final lane build. The warning drivers above ran after the
last R edit and install (`dev/ordmix-p2-final.sh`).

The core and frmtmb.sample suites, one file per process, every gate on
(`dev/ordmix-run-par.sh lane p2final dev/ordmix-p2-jobs-suites.txt
20`, log `dev/ordmix-suite-p2final/p2final.log`): 252 of 252 files
ran. frmtmb 206 files, pass 16662, fail 0, error 0, skip 0, warning 0;
frmtmb.sample 46 files, pass 2549, fail 0, error 0, skip 1
(`test-scale.R`, the scale tier), warning 0. The affected files in it:
`test-ordinal-mixture.R` 92, `test-hurdle-cum-thres-cs.R` 47,
`test-nonfinite-gradient.R` 5, `test-ordinal-mixture-brms.R` (gated)
22, `test-xbeta-zibb-hurdle-cum.R` 606; frmtmb.sample
`test-ordinal-mixture-draws.R` 30, `test-ordinal-draws-names.R` 31.

R CMD check --as-cran of frmtmb (`dev/ordmix-check/p2-frmtmb/`):
Status 1 NOTE, V8 unavailable for the HTML manual's math, as in rounds 0
and 1.

## Punch round 1 (review `dev/reviews/2026-10-05-ordmix.md`)

Every item of the review, in its order. Scripts and logs are
`dev/ordmix-p1-*`. "Seen to fail" means the new test was run against
the build before the fix and failed there.

### B1. The flat-direction warnings on identified models (fixed)

`mixture_ord_ident_check()` (R/families.R):

- The shared-threshold branch now needs identical components: one
  family, one link, one threshold structure, and every `disc<k>` held
  at the same value.
- The no-predictor branch and the hurdle branch are silent when a
  prior entry holds a `theta<k>`, an `hu<k>` or a threshold of the
  mixture (new `mixture_ord_prior_holds()`, the `pinned` pattern of
  `ord_fit_check()`).
- Test "the flat-direction warnings stay off identified mixtures"
  (`test-ordinal-mixture.R`): different families, different links,
  priors on `hu<k>`, plus the truly flat control, which still warns.
  Seen to fail: 4 failures on the build before the fix
  (`dev/ordmix-p1-log-t-mixture-before.txt`).

The review's own script, every model it defines (its driver ran 11 of
the 13; `dev/ordmix-p1-fa.sh`, logs `dev/ordmix-p1-fa-log/`, table by
`dev/ordmix-p1-fa-sum.R`). Fits that raised each kind of warning, of
20 per model, before (the review's logs) and after:

    model                 flat before  flat after  degenerate  other
    cum2                       0           0           0         0
    probit2                    0           0           0         0
    cum_sr_clog                0           0           0         0
    mu2                        0           0           0         0
    hurdle_huz                 0           0           0         0
    sr_acat_th                 0           0           0         0
    gr_cum2                    0           0           4         0
    mu_cum_sr_theta           20           0           0         1
    mu_cum_links_theta        20           0           0         1
    mu_same_theta (control)   20          20           0         0
    hurdle_prior              20           0           2         0
    plain_cum                  0           0           0         0
    plain_probit               0           0           0         0

"other" is nlminb's "singular convergence (7)", one fit each, as
before. Each of the 6 degenerate warnings is on a fit with a
non-finite standard error (`gr_cum2` seeds 5, 8, 14, 20; `hurdle_prior`
seeds 10, 13); `gr_cum2` seed 10 also has one and does not warn (see
B2, misses). No errors in 260 fits.

### B2. Degenerate fits, the infinite gradient, the record

**The degenerate-component warning** (new
`mixture_ord_degenerate_check()` and `mixture_ord_degeneracy()`, run
from `mixture_ord_fit_check()` at fit end). Per component: `reach`,
the largest `|disc (tau - cs - eta)|` over rows and threshold groups,
and `collapsed`, two ordered thresholds that are the same double. It
warns when `reach > 50` (the CDF is 0 or 1 to double precision for
every link there) or a gap is collapsed, names each such component and
what it found, says that the likelihood's supremum is often there, and
names the remedies (several starts, priors as brms has them, fewer
components). Measured on 300 fits (`dev/ordmix-p1-degen.R`, logs
`dev/ordmix-p1-degen-log/`), the review's label (`|estimate| > 30` or
a non-finite standard error) as the truth:

    fits                                 degenerate  warns on them  false alarms
    review's 160 (4 x 40 seeds)              22          20         0 of 138
    review's identified designs, 7 x 20       8           3         0 of 132

The largest `reach` of a fit the label calls sound is 41.1 (the
threshold 50 is above it). The 7 misses are fits whose component
stopped using a category without running off (`reach` 18.7 to 41.1, no
collapsed gap): `dev/ordmix-p1-seed10.R` shows one, a log-increment of
-17.7 and a threshold at -15.2, with two Hessian eigenvalues of 1e-5
and 8e-9. The statistic `cover` (the largest probability a component
gives its least used category) does not separate them: a sound fit
reaches 1.7e-7 and a missed one 0.13. Five of the seven are seed 10 of
the identified designs, one data set.

Tests: "a component run off to a degenerate boundary warns" (seed 19
of the review's process warns, seed 3 does not), seen to fail
(`dev/ordmix-p1-log-t-mixture-before.txt`); "the degenerate check says
only what it found" (three components on two-class data): the first
version of the check also warned "no non-missing arguments to max"
when a component gave `NaN` for a category on every row, seen as 4
escaped warnings on that build (`dev/ordmix-suite-p1f-before/`).

**The non-finite gradient** (`check_convergence()`, R/fit.R). A
gradient that is not finite at the reported optimum skipped the
gradient test, because that test needs a finite number. It now warns,
naming the parameters. Test file `test-nonfinite-gradient.R`: a
gaussian fit with its gradient made `Inf`, then `NaN`, and the
review's three-component mixture (`r_three_mixlink`). Seen to fail on
rellib-r5 (2 failures, 1 error: ordinal mixtures do not exist there;
`dev/ordmix-suite-p1e-base/`).

False alarms over every suite of the eight packages, by an
instrumented copy of the lane source (`dev/ordmix-p1-instr.R`, a log
line per fit to `dev/ordmix-p1-check/suite-log/`; suite logs
`dev/ordmix-suite-p1instr/`): 343 files, pass 23448, fail 0, error 0,
skip 15 (the scale tier), 0 escaped warnings. 9155 fits reached
`check_convergence()`. The non-finite-gradient warning fired on 3,
all in `test-nonfinite-gradient.R`, which builds them on purpose; on
no other fit of the eight suites. The degenerate warning ran on all 47
ordinal-mixture fits and fired on 6. Refitted with the review's label
(`dev/ordmix-p1-suitefires.R`, `dev/ordmix-p1-log-suitefires.txt`),
all 6 are degenerate (`|estimate|` up to 2.2e7, or a non-finite
standard error). One of them, the `thres(x = )` test's fit, ran off
only in a threshold above the data, which `thres_fit_check()` already
reports and which runs off whatever the component does; the check now
leaves such thresholds out (`th$unident`), and that fit no longer
fires (its `reach` is 5.8 and 18.8 without that threshold). That test
had allowed "bound", which also muffled "boundary"; it now allows its
threshold warning by name, seen as an escaped warning on the build
before (`dev/ordmix-suite-p1g-before/`). The same build's check also
warned "no non-missing arguments to max" on a fit with an all-`NaN`
category (above); the instrumented run was of that build, and both
fixes come after it, so the final suite run below is of the fixed
build.

**The record.** The identifiability section below is corrected with
the review's numbers and this rerun. `mixture()`'s documentation
("Ordinal components") now says that maximum likelihood often has its
supremum at a degenerate boundary (22 of 160, 16 of 22 above the best
other start), to compare several starts, the more so with three or
more components, and to hold the components with priors.

**Multi-start: filed, not built** (`dev/test-backlog.md`, "Filed by
lane ordmix after punch round 1"). In 16 of the review's 22 degenerate
fits the degenerate point has the higher likelihood, so a "best of n
starts" option would select it; a useful option ranks starts by
something else and reports them all, which is not small.

**Probit underflow past `|eta| = 38.2`: filed, not fixed.** The
candidate `pnorm(eta, log.p = TRUE) - pnorm(-eta, log.p = TRUE)` is not
the current form bit for bit where that is finite: 483 of 2001 values
on `[-38, 38]` differ (at most 4.3e-12 relative) and 1903 of 1977
derivatives (at most 1.4e-13), through RTMB's tape
(`dev/ordmix-p1-probit.R`, `dev/ordmix-p1-log-probit.txt`). It would
move every plain probit fit, so the 484 plain outputs could not stay
identical. Filed with the repro.

**The gap floor binds only at equal doubles**, in the docs and in the
warning: a collapsed gap is one of the two things the degenerate
warning names ("two of its thresholds are the same number, where the
likelihood holds the category between them at about exp(-692)").

### Minors

- m1. NEWS: `order = "mu"` with ordinal components moved from Breaking
  changes into the ordinal-mixture feature bullet.
- m2. `test-ordinal-mixture.R`: the equidistant check is now
  `abs(diff(diff(tk))) <= 4 eps max|tk|` (measured 0, and 0.33 of
  that unit on the other component); the stuck slopes are
  `expect_identical()` (they are the same double;
  `dev/ordmix-p1-m2.R`, `dev/ordmix-p1-log-m2.txt`).
- m3. `test-ordinal-mixture-draws.R`: both `suppressWarnings()` are
  `allow_warnings(..., omd_stan)`, Stan's sampler diagnostics by name
  (divergent transitions, maximum treedepth, the pairs() note, R-hat,
  ESS), measured per fit in `dev/ordmix-p1-m3.R`
  (`dev/ordmix-p1-log-m3.txt`). The "none" fit's degenerate warning is
  allowed by name at the fit.
- m4. The reason for keeping the floor out of plain fits, measured
  (`dev/ordmix-p1-m4.R`, `dev/ordmix-p1-log-m4.txt`): on 100000 random
  pairs the floored value is the unfloored one bit for bit, as the
  review found, but the taped gradient is not: 8675 of 200000 partial
  derivatives differ (median 5.6e-16 relative, up to 4e-3 near `b =
  -40`). The comment on `ord_log_interior()` now says that.
- m5. `cs()` outside the latent predictor gives brms's sentence first
  on every family (`sigma ~ cs(w)` on a gaussian too); the mixture's
  `mu1, mu2, ...` are named only on an ordinal mixture (R/frame.R; test
  in "cs() outside mu is refused").
- m6. NEWS and the `thres() x mixture` compat row qualify the 2.5 ulp:
  the 20 shapes of `dev/ordmix-lpcheck.R`, and a probit component past
  a latent distance of about 38 is `NaN` where brms is finite.
- m7. The Links comment in `test-ordinal-mixture.R` says brms leaves
  out `theta1` without a predictor and listing it is frmtmb's
  convention.
- m8. `draws_to_natural()` (frmtmb.sample R/draws-brms.R) used
  `x[-seq_len(R)]` in both directions; with `R = 0` that selects
  nothing. Now `x[seq_along(x) > R]`. Test "a threshold block with no
  free parameter keeps its column" (`test-ordinal-draws-names.R`), seen
  to fail: 1 failure and 1 error (`dev/ordmix-suite-p1b-before/`).
- m9. A multivariate model, `mixture(cumulative(), sratio())` with
  `theta1 ~ z` beside a gaussian response, against brms's compiled
  program at matched parameters (`dev/ordmix-p1-mv.R`,
  `dev/ordmix-p1-log-mv.txt`): 0 ulp at the optimum, 3.3, 0.8 and 1.6
  ulp at three perturbed points; brms's gradient at frmtmb's optimum is
  8.99e-5, frmtmb's own.

### User decisions, not changed

`cs()` on `cumulative()` stays refused, and the compat row's reason is
now the crossing-threshold hazard rather than "not identified". The
hurdle-beside-non-hurdle refusal stays.

### For the consolidator (merge with lane fixes; not merged here)

The three `R/priors.R` conflicts the review lists:

1. `prior_table()`: this lane adds `dpar = tdp` to the `Intercept`
   rows, fixes adds the per-threshold `coef = "1"..` rows. brms lists
   the per-threshold rows under `dpar = mu<k>` for `order = "none"` and
   with no `dpar` for `order = "mu"`, so the merged code takes the
   count from `ob$comp` (not `extra_tpl_name(..., "tau_raw")`, empty in
   a mixture) and passes `dpar = tdp`.
2. `ordinal_threshold_entry()`: keep fixes' form of the early return
   (no `nzchar(s$coef)`), or a per-threshold row never reaches a
   mixture.
3. The returned entry carries both this lane's `offset` (none under
   `order = "mu"`, `dpar` otherwise) and fixes' `coef_k`.

Clean textual merges that interact:

- fixes' `test-compat.R` asserts the `disc`, `equidistant` and
  `sum_to_zero` x `mixture` cells "refused" and that
  `mixture(cumulative(), cumulative())` errors; after the merge these
  are "works" (equidistant refused under `order = "mu"`, as in brms)
  and the fit succeeds.
- fixes' `ord_thres_linpred()` receives an ordinal mixture (type
  "ordinal"); add brms's "'incl_thres' is not supported for mixture
  models." branch first.
- fixes' `ord_internal_labels()` filters on `extra_tpl_name(...,
  "tau_raw")` and `dpar == "mu"`, so a mixture's `confint()` keeps
  `tau_raw1_1`; a port through `ord_lp_blocks()` must label a shared
  vector once.
- `ordinal_center_offset()` gained a `dpar` argument defaulting to
  `NULL`; fixes' two-argument calls keep their meaning.

Functions this round touched: R/families.R `mixture_ord_fit_check()`,
`mixture_ord_ident_check()`, `mixture_ord_prior_holds()` (new),
`mixture_ord_degenerate_check()` (new), `mixture_ord_degeneracy()`
(new), `ord_log_interior()` (comment), `mixture()` (docs); R/fit.R
`check_convergence()`; R/frame.R the `cs()` refusal; R/compat.R two
row texts; frmtmb.sample R/draws-brms.R `draws_to_natural()`.

### Tests and checks after the round

Affected files, one per process, on the final lane build:
`test-ordinal-mixture.R` pass 86, `test-nonfinite-gradient.R` 5,
`test-hurdle-cum-thres-cs.R` 47, `test-ordinal-mixture-brms.R` (gated)
22, `test-brms-parity-defects.R` 170, `test-xbeta-zibb-hurdle-cum.R`
606, `test-compat.R` 712, `test-compat-register.R` 103,
`test-prior-compat.R` 197; frmtmb.sample `test-ordinal-draws-names.R`
31, `test-ordinal-mixture-draws.R` 30; every one fail 0, error 0,
warning 0 (`dev/ordmix-suite-p1c/`, `-p1d/`, `-p1g/`).

Every suite of the eight packages on the final lane build, one file
per process, every gate on (`dev/ordmix-run-par.sh lane p1final
dev/ordmix-p1-jobs-all.txt 20`, log
`dev/ordmix-suite-p1final/p1final.log`): 343 of 343 files ran, pass
23449, fail 0, error 0, skip 15 (`test-scale.R`, the scale tier, as in
round 0), warning 0. Per package: frmtmb 206 files 16656, coupling 11
542, eam 29 1743, latent 10 359, learn 15 500, ode 11 547, sample 46
2549, spline 15 553.

R CMD check --as-cran (`dev/ordmix-check.ps1`, `dev/ordmix-check/p1-*`):
frmtmb 1 NOTE (V8 unavailable for the HTML manual's math, as in round
0); frmtmb.sample Status OK.

## What brms does (measured)

`dev/ordmix-brms-code.R` and `dev/ordmix-brms-code2.R` (outputs
`dev/ordmix-log-brms-code.txt`, `-code2.txt`), `dev/ordmix-brms-links.R`,
`dev/ordmix-brms-csdisc.R`:

- **Ordinal mixtures.** `mixture()` of ordinal families is allowed; a
  mix of ordinal and non-ordinal is refused ("Cannot mix ordinal and
  non-ordinal families", for poisson and categorical; gaussian is
  refused earlier as real against integer support). brms prints
  "Setting order = 'none' for mixtures of ordinal families": the
  default is `order = "none"`.
- **`order = "none"`**: each component has its own thresholds,
  `Intercept_mu<k>` (`ordered` for cumulative and hurdle, `vector`
  otherwise), its own centering of `X_mu<k>`, and draws
  `b_mu<k>_Intercept`. `disc<k>` is fixed at 1 in transformed
  parameters unless modeled. Prior rows: `Intercept` and `b` per
  `dpar = mu<k>`, `delta` per `mu<k>` under equidistant, one row per
  `thres(gr = )` level per component, `dirichlet(1) theta`.
- **`order = "mu"`** (or `TRUE`): brms's `fix_intercepts()`. One
  parameter `fixed_Intercept` (per level under `thres(gr = )`), ordered
  when any component is cumulative or hurdle (`vector` for sratio +
  acat), assigned to every `Intercept_mu<k>`; a sum-to-zero component
  centers it (`Intercept_mu2_stz = Intercept_mu2 - mean(...)`), and no
  design is centered. The prior row is class `Intercept` with no
  `dpar`. With an equidistant component brms refuses: "Cannot use
  equidistant and fixed thresholds at the same time."
- **Equidistant under `order = "none"`** (lane ordinal's defect 4,
  confirmed): the parameters block declares `real<lower=0> delta;`
  once per component, and the transformed parameters read `delta_mu1`,
  `delta_mu2` (`delta_mu1_1` under `thres(gr = )`, beside a declared
  `delta_1`). stanc refuses the program.
- **cs()**: per component, `bcs_mu<k>`, with the thresholds read as
  `Intercept_mu<k> - transpose(mucs_mu<k>[n])`; brms accepts it on a
  cumulative component too, with its "experimental" warning.
  `disc ~ cs(z)` is refused: "Category specific effects are only
  supported for the main parameter 'mu'."
- **Hurdle components** are allowed, also beside a component with no
  hurdle (`mixture(cumulative, hurdle_cumulative)` builds).
- **Links line**: `brms:::summarise_links()` gives
  `mu1 = probit; mu2 = logit`, and lists only predicted dpars
  (`disc1 = log; theta1 = identity` when modeled).
- **`hurdle_cumulative()` with `thres(gr = )`**: the hurdle stays in
  each group's density, `hurdle_cumulative_*_merged_lpmf`, with
  `nthres` counted per group over the categories above the hurdle;
  prior rows one `Intercept` per level, `delta` per level under
  equidistant. With `cs()`: `hurdle_cumulative_logit_lpmf(Y | mu, hu,
  disc, Intercept - transpose(mucs[n]))`. `cs()` with `thres(gr = )` is
  refused ("Cannot use category specific effects in models with
  multiple thresholds"), as for every ordinal family.

## What changed

Core, by file and function (for the merge with lane fixes):

- `R/families.R`
  - `mixture()`: ordinal components go to `mixture_ordinal()`;
    `mixture_check_order(order, ordinal)` returns `"none"`/`"mu"` and
    accepts `"mu"` for ordinal components only.
  - New: `mixture_ordinal()` (the refusals), `mixture_build_ordinal()`
    (the family), `mixture_logspace_add()`, `mixture_ord_spread()`,
    `mixture_ord_start_spread()`, `mixture_ord_shared_start()`,
    `mixture_ord_fit_check()`, `mixture_ord_ident_check()`,
    `cs_target_family()`, `ord_cumulative_logpmf_cs()`,
    `ord_log_interior()`, `ord_gap_floor()`, `hurdle_cum_draw()`,
    `hurdle_above_aterms()`.
  - Changed: `ord_family_tail()` (grouped density and simulator
    factories `ord_glpdf_make`, `ord_gsim_make`), `fam_cumulative()`
    (reads `.gap_floor`), `ord_cumulative_logpmf()` (`cs`,
    `gap_floor`), `fam_hurdle_cumulative()` (cs, grouped factories),
    `hurdle_cum_probs()` (cs), `hurdle_thres_finalizer()` (no refusal
    of `gr =`), `ord_fit_check()` (component `k`),
    `family_link_str()` (a mixture's `mu<k>` reports its cdf); the
    roxygen of `mixture()` (new section "Ordinal components") and of
    the hurdle family.
- `R/thres.R`: new `ord_lp_block()`, `ord_lp_blocks()`,
  `ord_prior_block()`, `thres_center_slices()`; changed
  `thres_finalizer()` (uses the family's grouped factories),
  `thres_lpdf()` and `thres_lpdf_cumulative()` (gap floor),
  `thres_pin_of_fit()` (a mixture's count), `thres_fit_check()`
  (`fam`, `comp` arguments).
- `R/compat.R`: rows `thres() x hurdle_cumulative` (works),
  `thres() x mixture` (works), new `cs_pred() x hurdle_cumulative`,
  `cs_pred() x mixture`.
- `R/confint.R`: `ord_extra_comps()`, `hyp_put_ordinal()`,
  `ord_delta_info()` iterate the threshold blocks.
- `R/brms-shapes.R`: `brms_extra_fixef()` iterates the blocks.
- `R/priors.R`: the default-row branch of the ordinal thresholds,
  `ordinal_threshold_entry()`, `ordinal_delta_entry()` (by `dpar`),
  `ordinal_center_offset()` (`dpar` argument). The per-threshold
  `Intercept` coef rows are lane fixes' and are not added here.
- `R/predict.R`: `ordinal_ncat()`, `ord_probs()` (a mixture goes to the
  new `ord_mix_probs()`), `cs_offsets_add()`, `with_cs_offsets()`, new
  `cs_slot()`.
- `R/objective.R`, `R/simulate-newdata.R`: the `cs()` offsets ride in
  `cs_slot(dpar)`, `.cs` for `mu` and `.cs_mu<k>` for a component.
- `R/frame.R`: the `cs()` check reads the component family and
  refuses `cs()` outside the latent predictor.
- `R/fit.R`: `make_start()` calls a family's `post$start_spread`.
- `R/methods-fit.R`: `coef()` finds a component's thresholds.
- `R/interop.R`: `emm_target_one()`: `mu<k>` of an ordinal mixture is
  latent.
- `R/brms-names.R`: `brms_par_labels()` (bug fix below).

frmtmb.sample: `R/draws-brms.R` (`draws_ordinal_cols()` joins the
blocks that read one component, new `draws_join_blocks()`),
`R/sample.R` (`default_priors_for()`, `default_prior_notes()`, new
`ord_disc_dpar()`). Nothing in the `posterior_linpred(incl_thres = )`
path, which is lane fixes'.

### The model

`mixture_build_ordinal()` builds the family twice: once in `mixture()`
and again in its `family_finalize()`, always from the components as
`mixture()` received them, after each component's own finalizer has
resolved `thres()` and its structure. Under `order = "none"` component
`k` reads its own block `tau_raw<k>` through its own (finalized)
density; under `order = "mu"` every component reads `tau_raw`, mapped
once (ordered when any component is) and centered by a sum-to-zero
component, through a density built on the thresholds themselves. The
mixture density is the log-sum-exp of the components' log densities
plus the log mixing weights. `fitted()` and `simulate()` read the
category probabilities out of that density one category at a time,
which is the theta-weighted sum. Each component's view of the dpars
carries its `.cs_mu<k>` as `.cs`, its `.eta_<dp><k>` as `.eta_<dp>`,
and `.gap_floor`.

## Measurements

### Log density against brms's compiled program

`dev/ordmix-lpcheck.R` (data seed `20261005 + case number`, n = 400;
perturbations seed `99 + case`), driver `dev/ordmix-lpcheck.ps1`,
logs `dev/ordmix-lpcheck-log/`, summary `dev/ordmix-lpcheck-sum.R` ->
`dev/ordmix-log-lpcheck-sum-final.txt`. Each case fits frmtmb, builds
brms's program with every prior flat, and compares `-obj$fn(p)` with
`rstan::log_prob(adjust_transform = FALSE)` at the optimum and at three
perturbed vectors (`opt + N(0, 0.3)`), plus brms's gradient at frmtmb's
optimum. A simplex `theta` keeps brms's dirichlet(1), whose density is
the constant Gamma(K); it is added (`log(2)` for `cum3`). The
generated summary:

    cases run: 30; with all four points: 29
    mixtures  cases 20 of 20, largest ulp 2.5, largest relative 5.55e-16
    hurdle    cases 7 of 8, largest ulp 1.6, largest relative 3.59e-16
    single    cases 2 of 2, largest ulp 37.2, largest relative 8.26e-15
    points where both densities are NaN: 3
    and frmtmb NaN at each of them too: TRUE

Per case (max ulp over the four points; brms's and frmtmb's largest
absolute gradient at frmtmb's optimum):

    acat_probit      2.5  NaN       0.000279
    cs_sratio_acat   0.9  3.31e-05  3.31e-05
    cum2             1.0  0.000169  0.000169
    cum3             0.7  0.000377  0.000377
    disc_theta       1.0  5.72e-05  5.88e-05
    equi_flex        1.0  2.23e-05  2.23e-05
    gr_cratio_acat   1.0  2.44e-05  2.44e-05
    gr_mu            0.9  0.00023   0.00023
    hurdle2          0.7  0.000219  0.000199
    hurdle_mu        0.8  3.09e-05  3.09e-05
    mix_hu_cs        0.0  0.000374  0.000374   (2 points NaN in both)
    mix_hu_gr        0.8  0.000257  0.000239
    mu2z             0.9  5.24e-05  5.24e-05
    mu_cum_sratio    0.9  0.000167  0.000167
    mu_disc          0.9  0.000208  0.000208
    mu_flex_stz      1.0  5.33e-05  5.33e-05
    mu_stz           0.9  7.4e-05   7.4e-05
    probit_sratio    1.0  0.00011   0.000104
    stz_flex         0.9  7.57e-05  5.35e-05
    thres5           0.9  0.000189  0.000189
    hu_cs_disc       0.8  0.000291  0.000291
    hu_cs_probit     0.0  0.000212  0.000212   (1 point NaN in both)
    hu_gr_disc       1.6  0.000368  0.000368
    hu_gr_equi       0.8  0.00082   0.00082
    hu_gr_logit      0.8  0.000319  0.000319
    hu_gr_probit     1.5  0.000758  0.000764
    hu_gr_stz        0.8  5.93e-05  5.55e-05
    hu_cs_logit      brms stops: "index 4 out of range" (brms-1)
    acat_probit_1    0.9  0.00137   0.00137    (plain acat probit)
    cum_cloglog_1   37.2  0.000163  0.000163   (plain cumulative cloglog)

The 20 mixture shapes cover both orders, cumulative/sratio/cratio/acat
and hurdle components, the logit, probit and cloglog, `disc1 ~ 0 + z`,
`disc2 ~ 0 + z` under shared thresholds, `theta1 ~ z`, `mu2 ~ z`,
three components, `thres(x = 4)`, `thres(gr = )` under both orders,
`cs()` per component, sum-to-zero beside flexible under both orders,
and equidistant components (`equi_flex`, against brms's program with
its defect 4 repaired the way its body reads it: the k-th `delta`
declaration renamed `delta_mu<k>`). brms's gradient equals frmtmb's
wherever both are finite (the fit stops at `grad_tol = 1e-8` on its
relative-function test, which leaves gradients of 1e-5 to 1e-3), so the
translation is right and the optimum shared. The NaN points are rows
whose `cs()` offsets cross two thresholds, where both densities are
undefined. The two single-family cases are no code of this lane: they
localize brms's NaN gradient in `acat_probit` (below); the 37.2 ulp of
the plain cumulative cloglog is 8.3e-15 relative, the precision of
frmtmb's log-odds form of the cloglog against brms's.

The gated test file `test-ordinal-mixture-brms.R` runs the same
translation on 11 shapes at two points each: pass=22 fail=0
(`dev/ordmix-log-g-mixture-brms.txt`).

### Against brms's R-side densities (ungated tests)

`test-ordinal-mixture.R` evaluates `-obj$fn(p)` at the estimates and
at a fixed perturbation against the theta-weighted sum of
`brms:::dcumulative()`, `dsratio()`, `dcratio()`, `dacat()` with the
thresholds rebuilt by hand from the internal vector, for 4 structure
mixes, `cs()`, `thres(gr = )` and two hurdle components; `fitted()` in
sample and on newdata to `1e3 * .Machine$double.eps`.
`test-hurdle-cum-thres-cs.R` does the same against
`brms:::posterior_epred_hurdle_cumulative()`'s formula
(`cbind(hu, (1 - hu) * dcumulative())`), per group and per row of
`cs()`. frmtmb.sample's `test-ordinal-mixture-draws.R` checks
`posterior_epred()` and `log_lik()` at the first and last draw against
the same sums at the stored columns.

### The hurdle's post-fit outputs the refusals guarded

Tested: `fitted()` in sample and on newdata, `simulate()` (each group's
own categories), `conditional_effects()` per level of `g`, the draws'
names and `posterior_epred()`. `dev/ordmix-hurdle-postfit.R` (data
seed 20261060, sampler seed 3, `dev/ordmix-log-hurdle-postfit.txt`)
adds `predict()` and frmtmb.sample's `log_lik()` and
`posterior_predict()` for `thres(gr = g)` and for `cs(x)`, each with
`hu ~ z`:

    gr: predict() 300 x 5, max |z| of the proportions against fitted() = 3.32 over 1500 cells
    gr draw 1: max |log_lik - log epred[y]| = 4.44e-16
    gr draw 150: max |log_lik - log epred[y]| = 4.44e-16
    cs: predict() 300 x 5, max |z| of the proportions against fitted() = 3.33 over 1500 cells
    cs draw 1: max |log_lik - log epred[y]| = 4.44e-16
    cs draw 150: max |log_lik - log epred[y]| = 4.44e-16

(`propagate_error = FALSE`, 4000 draws; the largest of 1500 standard
normals is about 3.4.) `disc` under `thres(gr = )`: held at 1 it is
hidden (`Links: cdf = logit; hu = logit`, no `disc_(Intercept)` in
`fixef(flatten = TRUE)`), and `disc ~ 0 + z` is fitted and shown
(`disc_z` 0.265, se 0.065; `dev/ordmix-hudisc.R`), its density checked
against brms's program (`hu_gr_disc`, 1.6 ulp).

### Fits against brms with seeds

`dev/ordmix-brmsfit.R` (data seed `20261040 + case`, n = 800; brms's
default priors, one chain of 1500, seed 1, started at frmtmb's
estimates so that both describe one labeling), logs
`dev/ordmix-brmsfit-log/`. `fixef()` row names and order are brms's in
all five cases (`ROWS identical names and order: TRUE`).

    case  model                                rows max|z|  se ratio      rhat
    hgr   hurdle, thres(gr = g), hu ~ z          9  0.094  [0.983, 1.064] 1.003
    hcs   hurdle probit, cs(x)                   6  0.126  [0.999, 1.046] 1.007
    mu    cumulative x2, order = "mu"            8  0.278  [0.855, 1.083] 1.004
    none  cumulative + sratio                    8  0.647  [0.259, 1.149] 1.159
    none2 cumulative x2                          8  1.231  [0.691, 1.171] 1.014

`z` is (frmtmb estimate - brms posterior mean) / brms posterior sd,
and the se ratio is frmtmb's standard error over brms's posterior sd.
The two `order = "none"` cases are weakly identified: in `none` the one
chain has not mixed (rhat 1.159); in `none2` the maximum likelihood
estimate puts `mu2`'s first two thresholds together (-1.8954 both: the
component gives category 2 no probability), where brms's student_t
prior keeps them 0.77 apart, and that is the row at |z| = 1.23.

### Identifiability and label switching

`dev/ordmix-ident.R` (replicate seeds 1..20, n = 500, starts seed
`1000 + s`), logs `dev/ordmix-log-ident-none.txt`, `-mu.txt`:

    SUMMARY mode=none R=20 default_at_best(gap<1e-4)=17 max_gap=0.7589 comp1_is_class1=0 max_swap_diff=1.14e-13
    SUMMARY mode=mu R=20 default_at_best(gap<1e-4)=10 max_gap=1.294 comp1_is_class1=20 max_swap_diff=1.14e-13

- Label switching is an identity: the components swapped (and
  `theta1` negated) give the same objective to 1.1e-13 in every
  replicate, under both orders. Under `order = "mu"` the two
  components share thresholds and differ only in `mu`, which is the
  label-switching case the brief names; brms's `order = "mu"` fixes
  the thresholds and orders nothing, so it does not identify the
  labels either.
- What the start values do. A continuous mixture here starts its `mu`
  intercepts at spread quantiles, which orders the components by
  location; an ordinal `mu` has no intercept, so the lane shifts each
  component's starting thresholds by `-F^-1(k / (K + 1))`
  (`mixture_ord_spread()`), and under `order = "mu"` spreads the
  starting slopes instead (`mixture_ord_start_spread()`). The labels
  follow: component 1 took the lower-location class in 20 of 20
  replicates under `"none"`, and the positive-slope class in 20 of 20
  under `"mu"`.
- Without the slope spread, shared thresholds start the components
  identical, a stationary point the optimizer does not leave: on
  `dev/ordmix-explore1.R` (seed 20261005) both slopes stopped at
  0.2552104 with `theta1` 0.5 and no standard error.
  `test-ordinal-mixture.R` pins both halves.
- Multimodality is real. The default start reached the best of 12
  starts in 17 of 20 replicates under `"none"` (largest gap 0.76
  log-likelihood units) and in 10 of 20 under `"mu"` (largest 1.29).
  The `"mu"` gaps are distinct local maxima, not early stops: on seed
  11 the default fit continued by BFGS stays at -663.8317 with a
  gradient of 2.3e-7 and a positive definite Hessian (smallest
  eigenvalue 7.7), while another start ends at -662.5381
  (`dev/ordmix-dbg6.R`, `dev/ordmix-dbg7.R`).
- **Correction (punch round 1): degenerate fits are common, and this
  record missed them.** `dev/ordmix-ident.R` ran seeds 1..20 at n =
  500, which hold none of the degenerate seeds, and its `quiet()`
  wrapper suppressed every warning. The review (`dev/ordmix-rev-degen.R`,
  the same data-generating process, seeds 1..40, n = 300 and 500,
  `cum + cum` and `cum + sratio`) found 22 of 160 default fits
  degenerate (an `|estimate| > 30` or a non-finite standard error), 15
  of the 22 with no warning at fit time, and in 16 of the 22 the
  degenerate point had a higher log-likelihood than the best of 6
  other starts: the supremum of the likelihood is at the boundary,
  where a component is a step function of `x`. The rerun of the same
  160 fits and of 140 fits of the review's identified designs
  (`dev/ordmix-p1-degen.R`, logs `dev/ordmix-p1-degen-log/`) gives the
  same 22 of 160, and 8 of 140. See "Punch round 1", B2, for the
  warning this now raises.

### Defects of the base build found here (fixed)

1. **`cs()` outside `mu` moved the thresholds.** The objective put the
   offsets of every predictor in the one `.cs` slot the ordinal
   densities read, so `bf(y ~ x, disc ~ cs(z))` with `sratio()` fitted
   `y ~ x + cs(z)`: logLik -549.1128459074 against -549.1128459085,
   with the coefficients reported as `disc_z[1..3]`
   (`dev/ordmix-base-behavior.R`, `dev/ordmix-log-base-behavior-base.txt`).
   Refused now in brms's words; the slots are per predictor
   (`cs_slot()`).
2. **The draws of a model with no location column were named one
   column off.** `brms_par_labels()` took `lab[-seq_len(n_beta)]`,
   which selects nothing when `n_beta` is 0, so every betad label was
   dropped: `bf(y ~ 1, disc ~ 0 + z)` with `cumulative()` stored its
   draws as `b_Intercept[1] b_Intercept[2] b_Intercept[3] tau_raw[3]`,
   the first column holding disc's slope (`dev/ordmix-emptybeta.R`,
   `dev/ordmix-log-emptybeta-base.txt`). The lane's `y ~ cs(x)` on the
   hurdle met it first, losing `hu`. The lane log shows the same until
   the fix was installed; the test pins the fix.

### The interior-category gap floor

An ordinal mixture's component can stop using a category; its
threshold increment then runs to 0 until the two thresholds are the
same double, and `log(F(a) - F(b))` is `log(0)` with an infinite
derivative, which reaches the mixture's gradient as NaN and stops the
optimizer. Measured: `mixture(cumulative(), sratio())` with
`thres(gr = g)` on `test-ordinal-mixture.R`'s data seed 20261011
stopped at "NA/NaN gradient evaluation" with thresholds 25.76 and
25.76 + exp(-32.85) (`dev/ordmix-dbg4.R`). `ord_log_interior()` holds
the gap at 1e-300 or above with the exact `min()` of
`dev/rtmb-pitfalls.md` item 11; it is on only inside a mixture
(`.gap_floor`), so a plain fit's tape is unchanged and its arithmetic
order is kept. With it the case fits (the simulate test). A plain
cumulative fit with an empty middle category stops on its own at an
increment of exp(-21.58) on both builds (`dev/ordmix-emptycat.R`), so
the floor is not needed there.

Three cumulative components on two-class data (`dev/ordmix-cum3two.R`)
now finish instead of stopping at NaN, at a degenerate optimum (a
slope of -4468): the third component captures a few rows exactly. That
is the likelihood's, not the floor's. Since punch round 1 the fit
warns there, naming component 1 (latent distance 21200) and component
2 (two thresholds the same number) (`dev/ordmix-p1-log-cum3two.txt`).

## Tests

New: `tests/testthat/test-ordinal-mixture.R`,
`tests/testthat/test-hurdle-cum-thres-cs.R`,
`tests/testthat/test-ordinal-mixture-brms.R` (gated),
`extensions/frmtmb.sample/tests/testthat/test-ordinal-mixture-draws.R`.
Changed: `tests/testthat/test-xbeta-zibb-hurdle-cum.R` (the refusals of
`cs()` and `thres(gr = )` on the hurdle are gone; their pair stays
refused).

On the lane build, one process per file:

    test-ordinal-mixture.R          pass=75 fail=0 err=0 skip=0 warn=0
    test-hurdle-cum-thres-cs.R      pass=47 fail=0 err=0 skip=0 warn=0
    test-ordinal-mixture-brms.R     pass=22 fail=0 err=0 skip=0 warn=0 (gated)
    test-ordinal-mixture-draws.R    pass=30 fail=0 err=0 skip=0 warn=0
    test-xbeta-zibb-hurdle-cum.R    pass=606 fail=0 err=0 skip=0 warn=0

The first two gained their `conditional_effects()` blocks after the
suite below ran (69 and 34 there); they were run again alone, with the
counts above. After the suite the R sources changed only in the text
of one compat row (`thres() x mixture`: its ulp count and number of
shapes, from the final lpcheck).

Seen to fail on rellib-r5 (`dev/ordmix-log-base-*.txt`):

    test-ordinal-mixture.R          pass=0 fail=1 err=12 skip=0 warn=1
    test-hurdle-cum-thres-cs.R      pass=1 fail=2 err=5 skip=0 warn=0
    test-ordinal-mixture-draws.R    pass=1 fail=1 err=5 skip=0 warn=0
    test-xbeta-zibb-hurdle-cum.R    pass=603 fail=0 err=1 skip=0 warn=0

The errors are the base's refusals (the weak form). The behavioral
failures on the base: `disc ~ cs(z)` fits instead of being refused
(defect 1 above); the draws of `y ~ 1, disc ~ 0 + z` are named
`b_Intercept[1..3], tau_raw[3]` (defect 2); and the hurdle density
handed `.cs` with crossing thresholds returns a finite value and its
simulator a draw, because the base ignores the offsets.

Guards and their absent cases, each in a test: equidistant beside
`order = "mu"` refused, equidistant under `"none"` fits;
`hurdle_cumulative()` beside `cumulative()` refused, two hurdles fit;
`cs()` on a `cumulative()` component refused, on an `sratio()`
component fits; `cs()` with `thres(gr = )` refused, each alone fits;
the flat-direction warnings fire on their models and not on the same
models with a predictor (`hu1 ~ z`, `hu1 = 0.1`, `y ~ x` under both
orders); class `Intercept` with no `dpar` refused under `"none"`, with
`dpar = "mu1"` under `"mu"`; `order = "mu"` for gaussian components
still refused.

### Plain fits are unchanged to the bit

`dev/ordmix-bitwise.R` (data seed 20261050), base against lane, on
nine models: cumulative, cumulative probit with `thres(gr = )`,
cumulative with `disc ~ 0 + z`, sratio with `cs()`, acat with
`thres(gr = )`, cratio sum-to-zero, the hurdle with `hu ~ z`, the
equidistant hurdle and a gaussian mixture; `opt$par`, `logLik`,
`fitted()`, `simulate(seed = 1)`, `variables()`, `fixef()` and
`print()` (`dev/ordmix-log-bitwise.txt`):

    BITWISE 63 of 63 outputs identical

## Every suite

`dev/ordmix-run-par.sh` over `dev/ordmix-jobs-all.txt` (all 342 test
files of frmtmb and its seven extensions, one process each, the lane
library first and rellib-r5 behind it for the six extensions this lane
did not change; `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`,
`FRMTMB_DRMTMB_FIT_TESTS` and `FRMTMB_FUZZ` set). Logs
`dev/ordmix-suite-all/`, totals from `dev/ordmix-suite-sum.R`
(`dev/ordmix-log-suite-sum.txt`):

    RAN 342 of 342
                 pkg files  pass fail err skip warn
              frmtmb   205 16621    0   0    0    0
     frmtmb.coupling    11   542    0   0    5    0
          frmtmb.eam    29  1743    0   0    3    0
       frmtmb.latent    10   359    0   0    2    0
        frmtmb.learn    15   500    0   0    2    0
          frmtmb.ode    11   547    0   0    1    0
       frmtmb.sample    46  2545    0   0    1    0
       frmtmb.spline    15   553    0   0    1    0
    all: 342 files, pass 23410 fail 0 err 0 skip 15 warn 0

    logs: 342  loading frmtmb from wt-ordmix-lib: 342

Every skip is a `test-scale.R` (the scale tier is not set).

## R CMD check --as-cran

`dev/ordmix-check.ps1`, built with vignettes and the manual inside
`dev/ordmix-check/<pkg>/`, `R_LIBS` = lane library, rellib-r5, user
library, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, from the final
sources:

- frmtmb: `Status: 1 NOTE`, the environmental "Skipping checking math
  rendering: package 'V8' unavailable"; tests 480 s, OK.
- frmtmb.sample: `Status: OK`.

## Ported brms-suite rows

None flips. The gated brms-suite files ran in the suite above with no
stale "now HOLDS" failure, and no ledger row of
`dev/brmsport-ledger.tsv` reaches an ordinal mixture or the hurdle's
`thres(gr = )` or `cs()`: brms's own tests of them are `stancode()`
string matches (`tests.stancode.R`, the `fixed_Intercept` and
`*_merged_lpmf` lines), which the port does not carry.

## Decided, and why

- **Equidistant mixtures fit brms's intended model** (tiebreaker:
  brms is clearly wrong; stanc refuses its own program). The names are
  the ones its transformed parameters use, `delta_mu<k>`
  (`delta_mu<k>_<g>` under `thres(gr = )`), with prior rows
  `delta, dpar = mu<k>` as `default_prior()` lists them.
- **A hurdle component beside one without a hurdle is refused.** brms
  builds it, but the two read the response differently (0..K against
  1..K): its cumulative component cannot evaluate a 0 response, so the
  model only runs on data with no 0, where every hurdle probability is
  pushed to 0.
- **`groups = ` is refused with ordinal components.** brms has no such
  argument; a latent-class ordinal mixture would need the extras passed
  through `mixture_structure()` and its simulator, which nothing tests.
- **`cs()` on `cumulative()` stays refused**, as is (and on a
  cumulative component of a mixture), although brms fits it and the
  density path `ord_cumulative_logpmf_cs()` now supports it for the
  hurdle. Lifting it is a scope decision for the user: the refusal's
  stated reason ("not identified") is not accurate (brms and
  `ordinal::clm(nominal = )` fit it); the real hazard is the negative
  probability of crossing thresholds, which the hurdle now documents.
- **Crossing `cs()` thresholds give NaN**: the density, as brms's; the
  row's `fitted()`, where brms's `posterior_epred()` returns the
  negative difference; and `NA` from `simulate()`, which has no
  distribution to draw from.
- **The gap floor is on inside a mixture only**, to keep every plain
  ordinal fit's tape bit for bit; the floored form is not bitwise
  equal to `logspace_sub(-b, -a)`.
- **`residuals(type = "osa")` is refused for ordinal mixtures** by
  name: a one-step residual needs the mixture's CDF, and the components'
  OSA paths step through their own.
- **Flat-direction warnings.** brms fits these models and lets its
  priors place the flat direction; here the standard errors would mean
  nothing, so the fit says so (measured: hurdle `hu1`, `hu2` standard
  errors 18.2 and 10.4 on the logit scale, `dev/ordmix-log-smoke.txt`).

## Not done, and why

- **No multi-start for `order = "mu"`.** The default reaches the best
  of 12 starts in 10 of 20 replicates; a start from a one-component fit
  spread by its slopes may do better, and is not measured. Filed in
  `dev/test-backlog.md` in punch round 1 (B2).
- **The per-threshold `Intercept` prior rows** (`coef = "1"`, ...) and
  the `theta` dirichlet row are not listed: the first is lane fixes',
  and no mixture lists the second today.
- **Sum-to-zero components have no class `Intercept` row**, as for one
  sum-to-zero family; brms lists one.
- **`insight` parameter scale** of an `order = "mu"` mixture: the block
  list holds `tau_raw` twice, so `insight_param_scale()` keeps the
  internal scale there.
- ~~A degenerate mixture optimum does not warn~~: it warns now (punch
  round 1, B2).
- ~~`draws_to_natural()` with an empty block~~: fixed in punch round 1
  (m8).
- **The scale tier** (`FRMTMB_SCALE_TESTS`) was not run; it times the
  extensions' own objectives.

## brms 2.23.0 defects met

- Defect 4 of lane ordinal, confirmed for mixtures: one bare `delta`
  per component, read as `delta_mu<k>` (and `delta_1` declared,
  `delta_mu1_1` read, under `thres(gr = )`).
- brms-1 of `dev/upstream-bugs.md` on the `cs()` path: `yh ~ cs(x)`
  with `hurdle_cumulative()` (logit) stops at "index 4 out of range;
  expecting index to be between 1 and 3" (`hu_cs_logit`).
- New: brms's gradient is NaN at frmtmb's optimum of
  `mixture(acat("probit"), cumulative("cloglog"))` while its log
  density equals frmtmb's to 2.5 ulp and frmtmb's gradient is 2.8e-4.
  Each component alone has a finite brms gradient (0.00137 and
  0.000163), so it comes from the mixture's `log_sum_exp` over a
  component term that underflows on some row, the pattern of lane
  ordinal's defect 1. Not filed.

## Version

Minor bump for frmtmb (new features, two behavior changes marked
BREAKING). frmtmb.sample: minor bump, floor the frmtmb development
version (ordinal mixtures, `ord_lp_blocks()` through
`brms_fixef_rows()`, the label fix).
