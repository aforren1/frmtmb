# Contributing to frmtmb

Thank you for your interest in frmtmb. This document tells you how to
report problems and how to send changes.

Everyone who takes part in this project must follow the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Report a bug

Open an issue at <https://github.com/aforren1/frmtmb/issues>.

Include a minimal reproducible example. A good example:

- uses a small simulated data set, or a data set from a package that
  the issue already needs (for example `lme4::sleepstudy`);
- runs in a fresh R session;
- shows the output you got and the output you expected.

Add the output of `sessionInfo()`.

If frmtmb disagrees with a reference package (brms, glmmTMB, lme4,
mgcv, and others), show both calls and both results. Disagreements
with a reference are treated as defects until proven otherwise.

## Suggest a feature

Open an issue first. Describe the model you want to fit and write the
formula you expect to work. If brms can fit the model, give the
`brm()` call. The formula grammar follows brms, so a brms call is the
clearest specification.

## Send a change

1. Fork the repository and make a branch from `main`.
2. Make your change.
3. Add tests. See "Tests" below.
4. Run the checks. See "Checks" below.
5. Open a pull request. Describe what changed and why.

Keep each pull request to one topic. Small pull requests are reviewed
faster.

Do not edit these files by hand:

- `NAMESPACE` and `man/*.Rd`. These are generated. Run
  `roxygen2::roxygenise()`.
- `codemeta.json`. Run `codemetar::write_codemeta()`.

## Style

Follow the style of the code around your change. The package uses
base R and a small set of imports. Do not add a dependency without
discussing it in an issue first.

Write comments that say why, not how.

## Tests

Every change to model code needs a test. The test suite is
`testthat` edition 3, in `tests/testthat/`.

Validation is layered. Put your test in the layer that fits:

- **Reference agreement.** Fit the same model with an exact reference
  (glmmTMB, drmTMB, lme4, mgcv, MASS, survival, nnet, GLMMadaptive,
  quantreg, or a closed-form result) and compare estimates within an
  explicit tolerance. This is the strongest test and is preferred.
  The drmTMB tier is gated on `FRMTMB_DRMTMB_FIT_TESTS`, because every
  test in it fits about 50 models; see
  `tests/testthat/test-drmtmb-agreement.R`.
- **Parameter recovery.** Simulate data from known parameters with a
  fixed seed, fit, and check the estimates are inside tolerance.
- **Invariants.** Check a property that must hold, such as agreement
  between two code paths that must give the same answer.
- **Edge conditions.** Check the error or warning for bad input.

Set a seed with `set.seed()` in any test that uses random numbers.

### Where test data comes from

A validation test must generate its reference independently of
frmtmb. Write the generative process by hand with `rnorm()`,
`rpois()` and the other base samplers, or take the data from another
package. Do not use `frm_simulate()` there. A simulator and a density
that share a convention error still agree with each other, so a
validation test that draws its data from the package under test cannot
see that error.

Every other test uses `frm_simulate()`. Tests of frame assembly,
refusals, method surfaces, printing and ergonomics only need data of
the correct shape, and so do the examples and the vignettes. A
`frm_simulate()` call states the model once, in the same vocabulary as
the `frm()` call below it. Hand-written `rnorm()` lines state the
model a second time, and the second statement can drift away from the
formula beside it.

Give the draw a seed of its own, distinct from the one that made the
covariates. `frm_simulate(seed = )` restarts the stream from that
seed, so reusing the covariates' seed hands the covariate draws back
as the residuals and the test fits a noiseless line. The convention in
this repository is the covariate seed plus 1000.

`tests/testthat/test-simulate-density.R` is what makes the second rule
safe. It is the one file that closes the loop on purpose. For every
family with a simulator it draws through the family's `sim` slot and
tests the draws against that same family's `lpdf`, evaluated
numerically: a goodness-of-fit test on cells whose probabilities are
summed or integrated from the density, and the first two moments taken
from the same measure. A family that declares no simulator must be
refused by name at both entry points. So a convention error shared by
a simulator and a density fails in that one file, instead of passing
quietly in every test that used the simulator to make its data.

Add a family, and add it to that file.

Tests that need a suggested package must call
`skip_if_not_installed()`. Slow tests must call `skip_on_cran()`.

Run the suite:

```r
devtools::test()
```

The default run skips the slow tests. Three environment variables
switch on the extended tiers. They all use the same testthat
framework as the regular tests:

| Variable | What it adds |
| --- | --- |
| `NOT_CRAN=true` | The heavy reference-validation files, which compare every model class against glmmTMB, lme4, mgcv, MASS, survival, and the rest. |
| `FRMTMB_FUZZ=true` | The pairwise grammar fuzzer. Set `FRMTMB_FUZZ_N` to cap the plan size for a quick smoke run. |
| `FRMTMB_BRMS_FIT_TESTS=true` | The tier that compiles brms Stan programs and checks that our log-likelihood equals the Stan program's log density at the estimate. |

To run everything:

```r
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(FRMTMB_FUZZ = "true")
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
devtools::test()
```

Requirements for the extended tiers:

- The packages in `Suggests`. The brms fit tier also needs a working
  `rstan` toolchain, because it compiles Stan models.
- Time. The full run takes much longer than the default run, and the
  brms fit tier takes the longest because of Stan compilation.
- No special hardware, and no large memory. Nothing is downloaded:
  every test simulates its own data from a seed or uses a data set
  from an installed package.
- No artifact needs manual inspection. A failure prints the
  reproduction case; the fuzzer also writes a runnable reproduction
  for any new finding.

## Checks

Before you open a pull request, run:

```r
devtools::document()
devtools::test()
devtools::check()
```

`R CMD check` must give no errors, no warnings, and no notes.

The package carries rOpenSci statistical software standards as
`@srrstats` roxygen tags. If your change affects a tagged behavior,
keep the tag with the code. To see the current state:

```r
srr::srr_report()
```

## Documentation

The website at <https://aforren1.github.io/frmtmb/> is one GitHub Pages
site that holds eight pkgdown sites: the core package at the root and
each extension in a subdirectory, which is how the Extensions menu in
the navbar works.

GitHub Actions builds and deploys it. You do not commit a built site.
The workflow is `.github/workflows/pkgdown.yaml`, and it runs on a push
to `main` that touches the documentation, on a release, and on demand.

To build the same site on your own machine, install all eight packages
and run:

```
Rscript dev/release/build-docs.R --dest=/tmp/frmtmb-site
```

On Windows, `dev/release/run-docs.ps1` is a wrapper that sets the paths
of the development box and calls the same script.

The core site is built first, then the seven extensions, and then the
root search index and sitemap are rebuilt. The last step is not
tidiness: pkgdown builds those two by walking the destination, so the
search box at the root finds extension pages only if the extension sites
are already there when the index is written. The script asserts that the
finished index holds a path under every extension.

The script empties the destination before it starts, which is how a
renamed package's directory stops being served, and it exits non-zero
if:

- fewer than eight sites build, or a site has no `reference/index.html`;
- a `href` in the Extensions menu of the root `_pkgdown.yml` names a
  subdirectory that did not build;
- a `destination` in an extension `_pkgdown.yml` does not end in the
  package name;
- an installed package version differs from its `DESCRIPTION`. pkgdown
  renders the *installed* package, so a stale install publishes stale
  output.

To rebuild one site while you work, name it and lower the count:

```
Rscript dev/release/build-docs.R --dest=/tmp/one \
  --pkgs=frmtmb.eam --min-sites=1
```

A partial tree has dead navbar links, so the script says so and skips
the whole-tree checks. Do not deploy one.

### Articles that do not run their code

Each article gates its fitted output on `requireNamespace()`, so a
machine without a `Suggests` package renders the page with the numbers
dropped and reports success. Pass `--require-articles` to turn a
missing one into a failure before the build starts; CI always passes it.
The script reads the list of gate packages out of the article sources,
so adding a gate needs no change here.

One of those packages is not on CRAN: `RTMBode`, which the whole
`frmtmb.ode` article depends on, comes from
<https://kaskr.r-universe.dev>.

**No article fits a model with brms or Stan.** `vignettes/brms-migration.Rmd`
and the two brms articles in `frmtmb.sample` set `eval = FALSE` for the
whole document, and `sampling.Rmd` runs only the chunks that need no
sampler and says so on the page. CI therefore installs no Stan
toolchain, and the site is not degraded by that. The brms and Stan
agreement evidence is in `.github/workflows/brms-likelihood.yaml`. If
you make an article fit with brms, reuse that workflow's Stan setup in
the pkgdown workflow; do not let the article render short.

## Questions

Open an issue with the question. There is no separate mailing list.
