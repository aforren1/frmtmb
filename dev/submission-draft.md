# rOpenSci statistical software submission: draft

Working text for the submission issue at
<https://github.com/ropensci/software-review>, template
`F-submit-statistical-software-for-review.md`. The template was read
again on 2026-09-04, together with the "Package Development" chapter
of the statistical software peer review guide
(<https://stats-devguide.ropensci.org/pkgdev.html>).

This is a working file. It is not part of the package.

## Header block

```
Submitting Author Name: Alex Forrence
Submitting Author Github Handle: <!--author1-->@aforren1<!--end-author1-->
Other Package Authors Github handles: (delete if none)
Repository:  <!--repourl-->https://github.com/aforren1/frmtmb<!--end-repourl-->
Version submitted: 0.49.1
Submission type: <!--submission-type-->Stats<!--end-submission-type-->
Badge grade: <!--statsgrade-->bronze<!--end-statsgrade-->
Editor: <!--editor--> TBD <!--end-editor-->
Reviewers: <!--reviewers-list--> TBD <!--end-reviewers-list-->
Archive: TBD
Version accepted: TBD
Language: <!--language-->en<!--end-language-->
```

0.49.1 is the version in `DESCRIPTION` as this draft is written. The
package is in active development, so the number at submission will
probably be later. Set it from `DESCRIPTION` on the day you post, and
check that nothing below has gone stale with it.

Paste the full `DESCRIPTION` into the code block the template gives.

## What is submitted, and what is beside it

The repository holds five R packages. Review covers one of them.

**The package submitted for review is `frmtmb`, the core package at
the repository root.** It fits models and reports on them. Its source
is `R/`, `tests/`, `vignettes/`, `man/` and `DESCRIPTION` at the root.
Its `.Rbuildignore` excludes `extensions/`, so the built tarball a
reviewer installs contains the core alone.

Four companion packages live under `extensions/`, one directory each.
Each is a separate R package with its own `DESCRIPTION`, its own
tests, its own vignette and its own R CMD check workflow. Each depends
on the core and reaches it only through the core's public extension
API. No companion package's R sources use `:::` to reach a core
internal; a few of their test files do, which is white-box testing and
not a runtime dependency. None of the four is submitted for review.

| package | adds |
|---|---|
| `frmtmb.sample` | NUTS sampling through tmbstan, and the posterior method surface |
| `frmtmb.latent` | discrete latent states: `hmm()` and `lca()` |
| `frmtmb.ode` | ordinary differential equation dynamics: `frm_ode()` |
| `frmtmb.ddm` | drift-diffusion response times: `wiener()` |

A user installs each one separately:

```r
remotes::install_github("aforren1/frmtmb")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.sample")
```

### Why the split was made

The design record is `dev/structured-family-protocol.md`, which holds
the packaging decision, and `NEWS.md` for 0.45.0 through 0.48.0, which
records each step. The reasons:

1. **The core's dependencies resolve entirely from CRAN.** The ODE
   surface needed `RTMBode`, which is not on CRAN. It left with
   `frm_ode()`, and the `Additional_repositories` line left the core
   `DESCRIPTION` with it. `Additional_repositories` now sits in
   `extensions/frmtmb.ode/DESCRIPTION`, where it belongs.
2. **The core's promise of no compilation stays exact.** Sampling
   goes through tmbstan, and so through Stan, which compiles. That
   promise now holds for the whole of the submitted package.
3. **One package uses one engine.** The core maximizes a likelihood.
   Each companion package brings a different engine or a different
   literature.
4. **One package is reviewable.** rOpenSci reviews one package. The
   split makes the reviewed thing the part that matters, and gives the
   reviewer a boundary that a test enforces rather than a convention.

A test in the core suite polices that boundary. No core file may name
a structured family or an ODE symbol anywhere. It asserts zero, with no
exempt file, and its positive control lives in the extension sources
(`NEWS.md`, 0.48.0).

### Why the split is a monorepo, not five repositories

One git repository, one directory per package. The precedent is
kaskr's own layout, where RTMB and TMB are each a subdirectory of one
repository. Development installs are one documented line each.
r-universe builds every package in a monorepo natively. CI partitions
by path filter. See `dev/structured-family-protocol.md`, "Packaging".

## The extension API is part of the reviewed surface

The seam between the core and the companion packages is public, and it
is in the submitted package. A reviewer who reviews the core reviews
this API.

What it consists of:

- `frmtmb_family()`, which builds a family from a plain R log-density.
- `frmtmb_structure()`, which declares a likelihood that does not
  factorize over the rows of the data.
- The registration functions: `frmtmb_register_aterm()`,
  `frmtmb_register_compat()` and `frmtmb_register_frame_check()`.
- The read-only accessors an out-of-tree family may use:
  `single_response()`, `eval_dpars()`, `fit_extras()`,
  `dpar_linpred()`, `response_mean()`, `as_frmtmb_family()`,
  `frame_block_of()`, `structure_supports_all()`,
  `mixture_posterior()`, `mixture_multimodal_refusals()`,
  `latent_probs()` and `frmtmb_ad_overload()`.
- Two reference pages that document the contract as one thing:
  `?"frmtmb-extension-api"` and `?"frmtmb-sampling-api"`. Both are on
  the package's pkgdown reference index.

Every name above is exported in the core `NAMESPACE`, and both
reference pages are built from the core `R/` sources.

**Its consumers are named, and they are outside the reviewed tree.**
The four companion packages are written against this API. That gives
the API something most extension seams do not have, which is a set of
independent users that a reviewer can read. `frmtmb.ddm` is the
strongest case: it was built against the exported API alone, with no
core change and no reach into internals, as an acceptance test of the
interface from an outsider's position (`NEWS.md`, 0.48.0).

The API carries a stability promise, stated in the README under "Life
cycle": the extension interface is public and settled, and internal
fields of the fitted object that no exported method reaches are not
part of the API.

## Scope

Check one box:

- [x] Regression and Supervised Learning

The template offers no "General" box. The General standards apply to
every statistical submission, so they are not selected. They are
claimed in `R/srr-stats-standards.R` and shown by `srr_report()`.

The seven other statistical categories are declined, each with a
reason, in the `NA_standards` block and the header of
`R/srr-stats-standards.R`. The Bayesian and Monte Carlo entry is the
one to raise with the editor. See "Open questions" below.

## Pre-submission inquiry

Not yet opened. The stats guide asks for one. Two things make it worth
opening first: the Bayesian category question, and the monorepo
question of what exactly is under review.

## Target audience and scientific applications

Applied statisticians and quantitative researchers who already write
brms formulas, and who need a likelihood-based answer. Four situations
recur:

1. A maximum-likelihood workflow is wanted or required, with
   confidence intervals and likelihood-ratio tests.
2. A model must be refit many times, inside a simulation, a bootstrap,
   a multi-start sweep, or a leave-one-out loop. frmtmb re-tapes in
   milliseconds and compiles nothing.
3. A model has to be screened before a Bayesian fit. The formula moves
   to brms without change.
4. No toolchain can be installed, so Stan is unavailable. Teaching
   settings and locked-down analysis machines are the usual cases.

Fields where this comes up: psychology and psychophysics, ecology
(zero-inflated and spatial counts), meta-analysis (`se()` addition
terms), and pharmacometrics. The last two of those reach further with
a companion package installed: population PK through
`frmtmb.ode::frm_ode()`, and drift-diffusion response times through
`frmtmb.ddm::wiener()`.

## Statement of need

brms gave R a formula grammar that covers a very wide model space. It
reaches that space only through Bayesian estimation in Stan. A user
who wants the same grammar under maximum likelihood must divide one
model across several packages, accept a smaller model, or wait for a
sampler. frmtmb closes that gap. It fits the brms grammar by maximum
likelihood with the Laplace approximation, and it compiles nothing at
run time.

The existing frequentist packages each stop at a different place.
lme4 fits mixed models in a small set of GLM families; only the mean
gets a formula, and each grouping factor gets one unstructured
covariance. glmmTMB adds many families, structured random-effect
covariances, and separate formulas for dispersion and for zero
inflation. Its predictor grammar is not brms's: there are no nonlinear
formulas, no multivariate responses with residual correlation, no
monotonic or measurement-error terms, and no correlation of effects
across formulas. gamlss, and its successor gamlss2, give every
distributional parameter its own additive predictor, which is the part
lme4 and glmmTMB lack. Random effects are additive terms there rather
than an lme4 grammar with structured covariances, and the spelling is
not brms's.

frmtmb supplies the combination, under ML or REML: distributional
regression on every parameter of the family, nonlinear formulas,
multivariate responses with `rescor`, effects correlated across
formulas with `|ID|`, custom families written as plain R
log-densities, and the lme4 random-effect grammar with structured and
spatial covariances.

The second need is migration. brms code ports by changing `brm()` to
`frm()`. The priors can stay where they are: `frm(prior = )` takes
brms's own spelling, and a prior object brms itself built is
translated. The fit is then MAP, so a prior is a penalty on the
likelihood rather than a posterior. A measured audit of the brms
vignettes, recorded in `dev/brms-vignette-audit.md`, puts about 7 of
10 of their model calls through that transform unchanged.

The same argument is the "Statement of need" section of `README.md`.

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

The contribution is architectural. Existing frequentist packages put a
formula front end on one fixed likelihood. frmtmb generates the
objective from the formula as an R closure and differentiates it on
the RTMB tape, so it does not select from a fixed set of likelihoods.
Feature combinations that are structural dead ends elsewhere become
ordinary code paths, and they compose with each other.

The extension API is the same argument carried one step further. A
family written outside the package reaches the core the same way a
built-in family does. The four companion packages are the evidence
that this works, because each one adds a model class to the grammar
without a core change.

Other relevant software: [gamlss](https://cran.r-project.org/package=gamlss),
[GLMMadaptive](https://cran.r-project.org/package=GLMMadaptive),
[BayesRTMB](https://github.com/norimune/BayesRTMB),
[qbrms](https://github.com/Tony-Myers/qbrms).

## Validation

Validation has three layers. The README states them under "Status".

1. **Comparison against an exact external reference.** Each model
   class for which an existing implementation is available is compared
   with a package that implements the same likelihood: glmmTMB, lme4,
   mgcv, nlme, MASS, survival, nnet, GLMMadaptive, quantreg or mice.
   For the classes with no such package, the comparison is against a
   closed-form marginal or a hand-written maximum likelihood built
   with `RTMB::MakeADFun()`. References are called live in the test
   run. They are never transcribed.
2. **Comparison against brms itself.** Design matrices, random-effect
   structures and special-term data are compared with
   `brms::make_standata()` to near machine precision. An opt-in tier,
   gated on `FRMTMB_BRMS_FIT_TESTS`, compiles the Stan program brms
   generates and compares our estimates against its posterior mode.
3. **A pairwise grammar fuzzer.** It sweeps feature combinations
   against metamorphic invariants. The resulting compatibility map is
   queryable with `frm_compat()` and is published as the "feature
   compatibility" article. The companion packages contribute their own
   rows to that map through `frmtmb_register_compat()`.

The heavy tiers are gated by environment variable, and `CONTRIBUTING.md`
documents each gate: `NOT_CRAN` for the reference validation,
`FRMTMB_FUZZ` for the fuzzer, and `FRMTMB_BRMS_FIT_TESTS` for the brms
fit tier. Nothing is downloaded by any of them.

**Continuous integration.** `.github/workflows/` holds R CMD check for
the core on four platforms, test coverage, and pkgcheck. Each
companion package has its own R CMD check workflow, filtered on its
own path: `check-frmtmb-sample.yaml`, `check-frmtmb-latent.yaml`,
`check-frmtmb-ode.yaml` and `check-frmtmb-ddm.yaml`. A change to one
extension checks that extension.

## Ethics, data privacy and human subjects research

Not applicable. The package fits models to data the user already
holds. It connects to no service, collects nothing, and writes no
files. The `habit_prep` data set that `vignette("habit")` uses is
redistributed under the MIT license of the original authors' OSF
deposit; see `inst/COPYRIGHTS`.

## Technical checks

- [ ] rOpenSci packaging guide read.
- [ ] Author guide read; maintenance for at least 2 years intended.
- [ ] Statistical Software Peer Review Guide for Authors read.
- [ ] `autotest` run with no failures. (Owned by the separate
      input-validation lane. Fill in from that lane's report, and see
      `dev/autotest-triage.md`.)
- [ ] `srr_stats_pre_submit()` confirms the package may be submitted.
      (Fill in from the last run. The report renders to
      `dev/srr-report.html`.)
- [ ] `pkgcheck()` confirms the package may be submitted. (Fill in
      from `dev/pkgcheck-docker.md`.)
- [x] Violates no Terms of Service. The package contacts no service.
- [x] CRAN and OSI accepted license: GPL (>= 2).
- [x] README carries installation instructions for the development
      version, for the core and for each companion package.

Run each of the three tool checks against the core package at the
repository root. `.Rbuildignore` excludes `extensions/`, so none of
them reaches a companion package. Run them again inside each extension
directory only if the editor asks for it.

## srr statement

Standards are tagged in place with `@srrstats` at the code or prose
that implements them, in `R/`, in `tests/testthat/`, and in
`vignettes/inputs.Rmd`, which carries the input, terminology,
attribute, scaling and benchmark-reproduction standards that a
document rather than code satisfies. Standards that do not apply are
collected in the `NA_standards` block of `R/srr-stats-standards.R`,
each with a reason. There are no `@srrstatsTODO` tags.

A tag audit before submission read every tag against the code it sits
on and recorded the verdict. It is `dev/srr-audit.md`. The audit also
logs the follow-up work it found and did not do, which is code rather
than prose.

All tags are in the submitted package. The companion packages carry no
srr tags, and claim none.

## Badging

Aiming for bronze at submission.

Silver needs at least one of the four aspects in the Guide for
Authors. Two are near:

- **Documentation.** A vignette for each part of the grammar and the
  inference surface, a queryable feature compatibility map, a full API
  reference, and a case study that replicates a published model.
- **Extended testing.** Every model class compared against an exact
  external reference, a pairwise grammar fuzzer with metamorphic
  invariants, and an opt-in tier that checks estimates against the
  mode of the Stan programs brms generates.

Raise silver with the editor after the first review round.

## In flight

Two changes are being written as this draft is prepared. Neither is
merged into the branch this draft describes. Both are stated here as
pending, not as facts about the submitted version. **Update this
section, and any claim it touches, when they merge.**

- **[PENDING] `get_prior()` will stop depending on which packages are
  attached.** Today it reports `(flat)` with frmtmb alone and the
  sampling defaults once `frmtmb.sample` registers them (`NEWS.md`,
  0.47.0); the change adds a `route` argument that selects the fit
  defaults or the sample defaults explicitly, so one call gives one
  answer whatever is on the search path.
- **[PENDING] An importance-sampling correction of the Laplace
  approximation.** `frm(importance = N)` will correct the Laplace
  approximation of the marginal likelihood by importance sampling,
  which gives the user a further remedy beside `quadrature = TRUE` and
  the measurement in `frmtmb.sample::check_laplace()`.

If either has merged by the time you post, move its sentence into the
body, cite the `NEWS.md` entry, and delete the bracket.

## Use of Generative AI

- [x] Generative AI tools were used to produce some of the material in
      this submission.

**[USER TO APPROVE]** The paragraph below is drafted from what the
repository shows, and from nothing else. The maintainer decides the
final wording and what to claim. Do not extend it with process claims
that the tree does not evidence.

> Generative AI was used throughout the development of this package
> and in preparing this submission text. The tool was Claude Code, in
> an agent configuration driven by the maintainer. The commit history
> records this. Most commits carry a `Co-Authored-By` trailer naming
> the Claude model that worked on them, and the trailers name more
> than one model, because the work spans several model releases. The
> command
> `git log --format="%(trailers:key=Co-Authored-By,valueonly)" | sort -u`
> lists them. The `dev/` directory holds the working records of that
> process: design records approved before implementation
> (`dev/structured-family-protocol.md`), feasibility memos written
> before a feature was attempted, and audits run against the finished
> code rather than against its description (`dev/srr-audit.md`,
> `dev/brms-vignette-audit.md`, `dev/bracket-sweep.md`,
> `dev/autotest-triage.md`). The maintainer is the sole author in
> `DESCRIPTION` and is responsible for the package. The validation
> described under "Validation" above is the check on the output: every
> model class is compared against an external reference implementation
> or a hand-written likelihood, and the model-building layer is
> compared against brms itself. Those comparisons do not care who
> wrote the code.

Points for the maintainer to settle before posting:

- Whether to name the specific models, or only the tool.
- Whether to say anything about human review of AI-written code and
  text. The repository does not record a review procedure, so this
  draft claims none.
- Whether to link the `dev/` records. They are in the repository, but
  `.Rbuildignore` excludes them from the built package.

## Publication options

- [x] CRAN intended, for the core and for the companion packages.
- [ ] Bioconductor: no.

## Code of Conduct

- [ ] I agree to abide by rOpenSci's Code of Conduct during the review
      process and in maintaining my package should it be accepted.

The repository carries `CODE_OF_CONDUCT.md` and `CONTRIBUTING.md` at
its root.

## Open questions for the editors

1. **What is under review in a monorepo.** The submitted package is
   the core at the repository root, and its built tarball excludes
   `extensions/`. The public extension API is in the core, so it is
   under review. Its four consumers are in the same repository, but
   outside the reviewed tree. Should a reviewer read those consumers
   as evidence about the API? If so, how should that be framed, so
   that it does not become a review of four unsubmitted packages?

2. **The Bayesian and Monte Carlo category.** The package declines it.
   The reason is written out in the header of
   `R/srr-stats-standards.R`. In short: the documented inferential
   surface of the submitted package is maximum likelihood with the
   Laplace approximation, and every estimate, standard error, interval
   and test in its vignettes is a maximum-likelihood quantity.

   The split makes this position cleaner than it was. The sampler is
   no longer in the reviewed package at all. `frm_sample()`,
   `check_laplace()`, `as_tmbstan()` and the whole draws method
   surface moved to `frmtmb.sample` at 0.47.0. What remains in the
   core is `frm()`, which optimizes, and `set_prior()`, which
   penalizes.

   The question for the editor is therefore narrower than it was. Does
   a package that declares a public sampling API
   (`?"frmtmb-sampling-api"`), and whose companion package uses that
   API to sample a posterior, take on the Bayesian standards? Our
   reading is that it does not. The standards describe a package whose
   front door is the sampler, and the front door here is `frm()`.

   **Before posting:** the header of `R/srr-stats-standards.R` still
   argues this case as though `frm_sample()` and `check_laplace()`
   were in the core package. They are not. That header needs rewriting
   to the post-split position. It is code, not prose, and is out of
   this draft's scope.

3. **brms in `Suggests`.** brms is suggested for three separate
   purposes: reference validation of the model-building layer against
   `brms::make_standata()`, the opt-in tier that compares estimates
   against the mode of a compiled brms Stan program, and coexistence
   handling, so that a session with both packages attached behaves
   predictably (`bf()` masking detection). Is a Suggests dependency on
   the package a submission is positioned against acceptable? Or would
   the editors prefer the reference comparisons move behind an
   environment-variable gate that the default check never reaches?
   They are already gated with `skip_on_cran()` and
   `requireNamespace()`.

4. **Badge grade.** Bronze at submission, with a case for silver
   sketched above. Would the editors prefer the silver case be made at
   submission instead?

**Withdrawn.** An earlier draft asked about an r-universe `Suggests`
dependency on `RTMBode`, which is not on CRAN, and about the
`Additional_repositories` line it needed. The question no longer
applies to the submitted package. Both left the core with `frm_ode()`
at 0.47.0, and they now sit in `extensions/frmtmb.ode/DESCRIPTION`.
Every dependency of the submitted package resolves from CRAN.
