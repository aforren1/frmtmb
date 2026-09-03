# rOpenSci statistical software submission: draft

Working text for the submission issue at
<https://github.com/ropensci/software-review>, template
`F-submit-statistical-software-for-review.md` (fetched 2026-09-02).

This is a working file. It is not part of the package.

## Header block

```
Submitting Author Name: Alex Forrence
Submitting Author Github Handle: <!--author1-->@aforren1<!--end-author1-->
Other Package Authors Github handles: (delete if none)
Repository:  <!--repourl-->https://github.com/aforren1/frmtmb<!--end-repourl-->
Version submitted: 0.39.0
Submission type: <!--submission-type-->Stats<!--end-submission-type-->
Badge grade: <!--statsgrade-->bronze<!--end-statsgrade-->
Editor: <!--editor--> TBD <!--end-editor-->
Reviewers: <!--reviewers-list--> TBD <!--end-reviewers-list-->
Archive: TBD
Version accepted: TBD
Language: <!--language-->en<!--end-language-->
```

Paste the full `DESCRIPTION` into the code block the template gives.

## Scope

Check one box:

- [x] Regression and Supervised Learning

The template offers no "General" box. The General standards apply to
every statistical submission, so they are not selected; they are
claimed in `R/srr-stats-standards.R` and shown by `srr_report()`.

The seven other statistical categories are declined, each with a
reason, in the `NA_standards` block and the header of
`R/srr-stats-standards.R`. The Bayesian and Monte Carlo entry is the
one to raise with the editor; see "Open questions" below.

## Pre-submission inquiry

Not yet opened. The stats guide asks for one, and the Bayesian
category question in "Open questions" is a good reason to open one
before the full submission.

## Target audience and scientific applications

Applied statisticians and quantitative researchers who already write
brms formulas and who need a likelihood-based answer. Four situations
recur:

1. A maximum-likelihood workflow is wanted or required, with
   confidence intervals and likelihood-ratio tests.
2. A model must be refit many times, inside a simulation, a
   bootstrap, a multi-start sweep, or a leave-one-out loop. frmtmb
   re-tapes in milliseconds and compiles nothing.
3. A model has to be screened before a Bayesian fit. The formula
   moves to brms without change.
4. No toolchain can be installed, so Stan is unavailable. Teaching
   settings and locked-down analysis machines are the usual cases.

Fields where this comes up: psychology and psychophysics
(drift-diffusion and psychometric-function fits), ecology
(zero-inflated and spatial counts), pharmacometrics (population PK
through `frm_ode()`), and meta-analysis (`se()` addition terms).

## Statement of need

brms gave R a formula grammar that covers a very wide model space. It
reaches that space only through Bayesian estimation in Stan. A user
who wants the same grammar under maximum likelihood must divide one
model across several packages, accept a smaller model, or wait for a
sampler. frmtmb closes that gap: it fits the brms grammar by maximum
likelihood with the Laplace approximation, and compiles nothing at
run time.

The existing frequentist packages each stop at a different place.
lme4 fits mixed models in a small set of GLM families; only the mean
gets a formula, and each grouping factor gets one unstructured
covariance. glmmTMB adds families, structured random-effect
covariances, and separate formulas for dispersion and zero inflation,
but its predictor grammar is not brms's: no nonlinear formulas, no
multivariate responses with residual correlation, no monotonic or
measurement-error terms, no correlation of effects across formulas.
gamlss, and its successor gamlss2, give every distributional
parameter its own additive predictor, which is the part lme4 and
glmmTMB lack, but random effects there are additive terms rather than
an lme4 grammar with structured covariances.

frmtmb supplies the combination, under ML or REML: distributional
regression on every parameter of the family, nonlinear formulas,
multivariate responses with `rescor`, effects correlated across
formulas with `|ID|`, custom families written as plain R
log-densities, and the lme4 random-effect grammar with structured and
spatial covariances. The second need is migration: brms code ports by
removing the priors and changing `brm()` to `frm()`, and a measured
audit of the brms vignettes puts about 7 of 10 of their model calls
through that transform unchanged.

The same text is the "Statement of need" section of `README.md`.

## General Standard G1.1

frmtmb is *an improvement on other implementations of similar
algorithms in R*. It introduces no new estimator. Every method it
implements has a primary reference, listed in the package-level
`@references` and summarized in the "Algorithm provenance" section of
`?frmtmb-package`:

- Laplace approximation and AD: [TMB / RTMB](https://cran.r-project.org/package=RTMB),
  Kristensen et al. (2016); Tierney and Kadane (1986).
- Formula grammar: [brms](https://paulbuerkner.com/brms/),
  Buerkner (2017).
- Random-effect syntax and structured covariances:
  [lme4](https://cran.r-project.org/package=lme4), Bates et al. (2015);
  [glmmTMB](https://glmmtmb.github.io/glmmTMB/), Brooks et al. (2017).
- Smooth terms: [mgcv](https://cran.r-project.org/package=mgcv),
  Wood (2011, 2017).
- Approximate Gaussian processes: Riutort-Mayol et al. (2023).
- Spatial fields: Lindgren et al. (2011); Riebler et al. (2016).
- One-step-ahead residuals: Thygesen et al. (2017).
- Multiple-imputation pooling: Rubin (1987).

The contribution is architectural. Existing frequentist packages put
a formula front end on one fixed likelihood; frmtmb generates the
objective from the formula as an R closure and differentiates it on
the RTMB tape. Feature combinations that are structural dead ends
elsewhere become ordinary code paths, and they compose.

Other relevant software: [gamlss](https://cran.r-project.org/package=gamlss),
[GLMMadaptive](https://cran.r-project.org/package=GLMMadaptive),
[BayesRTMB](https://github.com/norimune/BayesRTMB).

## Ethics, data privacy and human subjects research

Not applicable. The package fits models to data the user already
holds. It connects to no service, collects nothing, and writes no
files.

## Badging

Aiming for bronze at submission. Silver needs at least one of the
four aspects in the Guide for Authors; the two nearest are the
documentation tier (seven vignettes, a queryable feature
compatibility map, and a full API reference) and the extended testing
tier (per-file suite of more than 5000 tests, exact comparison
against nine reference packages, a pairwise grammar fuzzer with
metamorphic invariants, and an opt-in tier that checks estimates
against the mode of the Stan programs brms generates). Raise silver
with the editor after the first review round.

## The RTMBode dependency

`RTMBode` is in `Suggests`, and `DESCRIPTION` carries

```
Additional_repositories: https://kaskr.r-universe.dev
```

`RTMBode` is not on CRAN. It supplies the ODE solvers that
`frm_ode()` uses, and nothing else in the package needs it. Every
code path that reaches it is guarded by `requireNamespace()`, every
example and vignette chunk that uses it is conditional, and the tests
skip without it. The package installs, checks and runs with it
absent.

This is the pattern CRAN already accepts for an r-universe suggested
dependency, and the precedent is directly comparable: brms
(2.23.0), bayesplot (1.16.0) and bridgesampling (1.2-1) are all on
CRAN today with `cmdstanr` in `Suggests` and
`Additional_repositories: https://stan-dev.r-universe.dev/` in
`DESCRIPTION`. Verified against those three packages' installed
`DESCRIPTION` files on 2026-09-02.

`README.md` tells the user how to install `RTMBode` from
`https://kaskr.r-universe.dev`.

## Technical checks

- [ ] rOpenSci packaging guide read.
- [ ] Author guide read; maintenance for at least 2 years intended.
- [ ] Statistical Software Peer Review Guide for Authors read.
- [ ] `autotest` run with no failures. (Owned by the separate
      input-validation lane; fill in from that lane's report.)
- [ ] `srr_stats_pre_submit()` confirms the package may be submitted.
      (Fill in from the last run; the report renders to
      `dev/srr-report.html`.)
- [ ] `pkgcheck()` confirms the package may be submitted.
      (Fill in from `dev/srr-audit.md`, "pkgcheck" section.)
- [x] Violates no Terms of Service. The package contacts no service.
- [x] CRAN and OSI accepted license: GPL (>= 2).
- [x] README carries installation instructions for the development
      version.

## srr statement

Standards are tagged in place with `@srrstats` at the code or prose
that implements them, in `R/`, in `tests/testthat/`, and in
`vignettes/inputs.Rmd`, which carries the input, terminology,
attribute, scaling and benchmark-reproduction standards that a
document rather than code satisfies. Standards that do not apply are
collected in the `NA_standards` block of `R/srr-stats-standards.R`,
each with a reason. There are no `@srrstatsTODO` tags.

A tag audit before submission read every tag against the code it
sits on and recorded the verdict; see `dev/srr-audit.md`.

## Use of Generative AI

- [x] Generative AI tools were used to produce some of the material
      in this submission.

Fill in a description before posting: which tools, on which parts of
the package and of the submission text, and what human review each
received. rOpenSci asks for transparency here, not for a minimal
answer. Link the relevant commits or the `dev/` records.

## Publication options

- [x] CRAN intended.
- [ ] Bioconductor: no.

## Open questions for the editors

1. **The Bayesian and Monte Carlo category.** The package declines
   it, and the reason is written out in the header of
   `R/srr-stats-standards.R`. In short: the documented inferential
   surface is maximum likelihood with the Laplace approximation, and
   every estimate, standard error, interval and test in the
   vignettes is a maximum-likelihood quantity. `frm_sample()` is an
   opt-in bridge to NUTS through tmbstan. On a fitted model it
   explores the likelihood under flat priors, which is what makes
   `check_laplace()` meaningful. From a formula it does sample a
   posterior, under weakly informative default priors that match
   brms, with an LKJ default on correlations and a non-centered
   parameterization, and it reports `n_eff` and `Rhat`.

   That second route is the honest edge of the position, and it has
   grown since the position was first written. A reviewer reading
   `NEWS.md` will see brms-matched default priors, the LKJ default,
   `loo()` and `waic()`, a large draws method surface, and
   non-centered sampling arrive one after another. The question for
   the editor: does a documented, opt-in bridge whose purpose is
   brms-migration parity and Laplace diagnosis cross into the
   Bayesian category, and if so, is the right answer to claim the
   category (a tagging project across the BS standards, roughly the
   size of the existing RE tagging) or to narrow the bridge?

2. **brms in `Suggests`.** brms is suggested for three separate
   purposes: reference validation of the model-building layer
   against `brms::make_standata()`, the opt-in tier that compares
   estimates against the mode of a compiled brms Stan program, and
   coexistence handling so that a session with both packages
   attached behaves predictably (`bf()` masking detection). Is a
   Suggests dependency on the package a submission is positioned
   against acceptable, or would the editors prefer the reference
   comparisons move behind an environment-variable gate that the
   default check never reaches? They are already gated with
   `skip_on_cran()` and `requireNamespace()`.

3. **The r-universe `Suggests` dependency.** See "The RTMBode
   dependency". Confirming that the brms/bayesplot/bridgesampling
   precedent is acceptable to rOpenSci would settle it early.

4. **Badge grade.** Bronze at submission, with a case for silver
   sketched above. Would the editors prefer the silver case be made
   at submission instead?
