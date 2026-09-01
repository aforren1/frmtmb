# Contributing to frmtmb

Thank you for your interest in frmtmb. This document tells you how to
report problems and how to send changes.

Everyone who takes part in this project must follow the [Code of
Conduct](https://aforren1.github.io/frmtmb/CODE_OF_CONDUCT.md).

## Report a bug

Open an issue at <https://github.com/aforren1/frmtmb/issues>.

Include a minimal reproducible example. A good example:

- uses a small simulated data set, or a data set from a package that the
  issue already needs (for example
  [`lme4::sleepstudy`](https://rdrr.io/pkg/lme4/man/sleepstudy.html));
- runs in a fresh R session;
- shows the output you got and the output you expected.

Add the output of
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

If frmtmb disagrees with a reference package (brms, glmmTMB, lme4, mgcv,
and others), show both calls and both results. Disagreements with a
reference are treated as defects until proven otherwise.

## Suggest a feature

Open an issue first. Describe the model you want to fit and write the
formula you expect to work. If brms can fit the model, give the `brm()`
call. The formula grammar follows brms, so a brms call is the clearest
specification.

## Send a change

1.  Fork the repository and make a branch from `main`.
2.  Make your change.
3.  Add tests. See “Tests” below.
4.  Run the checks. See “Checks” below.
5.  Open a pull request. Describe what changed and why.

Keep each pull request to one topic. Small pull requests are reviewed
faster.

Do not edit these files by hand:

- `NAMESPACE` and `man/*.Rd`. These are generated. Run
  [`roxygen2::roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html).
- `codemeta.json`. Run
  [`codemetar::write_codemeta()`](https://docs.ropensci.org/codemetar/reference/write_codemeta.html).

## Style

Follow the style of the code around your change. The package uses base R
and a small set of imports. Do not add a dependency without discussing
it in an issue first.

Write comments that say why, not how.

## Tests

Every change to model code needs a test. The test suite is `testthat`
edition 3, in `tests/testthat/`.

Validation is layered. Put your test in the layer that fits:

- **Reference agreement.** Fit the same model with an exact reference
  (glmmTMB, lme4, mgcv, MASS, survival, nnet, GLMMadaptive, quantreg, or
  a closed-form result) and compare estimates within an explicit
  tolerance. This is the strongest test and is preferred.
- **Parameter recovery.** Simulate data from known parameters with a
  fixed seed, fit, and check the estimates are inside tolerance.
- **Invariants.** Check a property that must hold, such as agreement
  between two code paths that must give the same answer.
- **Edge conditions.** Check the error or warning for bad input.

Set a seed with [`set.seed()`](https://rdrr.io/r/base/Random.html) in
any test that uses random numbers.

Tests that need a suggested package must call `skip_if_not_installed()`.
Slow tests must call `skip_on_cran()`.

Run the suite:

``` r

devtools::test()
```

The default run skips the slow tests. Three environment variables switch
on the extended tiers. They all use the same testthat framework as the
regular tests:

| Variable | What it adds |
|----|----|
| `NOT_CRAN=true` | The heavy reference-validation files, which compare every model class against glmmTMB, lme4, mgcv, MASS, survival, and the rest. |
| `FRMTMB_FUZZ=true` | The pairwise grammar fuzzer. Set `FRMTMB_FUZZ_N` to cap the plan size for a quick smoke run. |
| `FRMTMB_BRMS_FIT_TESTS=true` | The tier that compiles brms Stan programs and compares our estimates against their posterior mode. |

To run everything:

``` r

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
- No special hardware, and no large memory. Nothing is downloaded: every
  test simulates its own data from a seed or uses a data set from an
  installed package.
- No artifact needs manual inspection. A failure prints the reproduction
  case; the fuzzer also writes a runnable reproduction for any new
  finding.

## Checks

Before you open a pull request, run:

``` r

devtools::document()
devtools::test()
devtools::check()
```

`R CMD check` must give no errors, no warnings, and no notes.

The package carries rOpenSci statistical software standards as
`@srrstats` roxygen tags. If your change affects a tagged behavior, keep
the tag with the code. To see the current state:

``` r

srr::srr_report()
```

## Questions

Open an issue with the question. There is no separate mailing list.
