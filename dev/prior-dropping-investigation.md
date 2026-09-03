# Why the pkgcheck container drops priors

Investigated 2026-09-03 against `0afcd27`, in the container the other lane
committed as `frmtmb-pkgcheck:deps`.

## Summary

The priors are not dropped. Nothing reaches the sampler at all.

`tmbstan` 1.2.0 built against `StanHeaders` 2.39.1 samples
`y ~ std_normal()` instead of the `TMB` objective. The likelihood, the
priors and every bound the objective encodes are absent from the sampled
density. The chains are a standard normal in the unconstrained space.

This is an upstream bug in `tmbstan`, not in `frmtmb`. It is triggered by
`StanHeaders` 2.39.1, which CRAN published on 2026-09-02, one day before
the container run. It is silent: the build prints one warning that no
installer shows, and the run produces plausible-looking draws.

`dev/pkgcheck-docker.md` concludes that "`StanHeaders` alone does not
explain the drift" because the macOS runner also holds 2.39.1. That
reasoning does not hold: what matters is the `StanHeaders` version present
when `tmbstan` is COMPILED, and the macOS runner installs a `tmbstan`
binary compiled earlier. See "Who is affected" below. That section of
`dev/pkgcheck-docker.md` needs a correction.

## The symptom, restated

Reported by the previous lane: in the container, sampling with
`priors = "flat"` and sampling with the default priors give
`max |flat - default| = 0`, and `prior_summary()` reports the prior in
both. Everything else in that table follows from the same cause. Read the
container numbers again with "the sampler ignores the model" in mind:

| Quantity | Container | Reading |
| --- | --- | --- |
| posterior mean of `x`, ML gives 0.454 | -0.037 | the mean of a standard normal |
| `sd_ratio` for `x` | 10.9 | posterior sd near 1 over a Wald se near 0.09 |
| mixture separation | 0.047 | the components never see the data |
| `sd(flat) / sd(default)` | 1.000 | both chains are the same standard normal |

## Step 1. The augmented objective is correct

`dev/pdi-repro.R` builds the objective `frm_sample()` hands to `tmbstan`
(`frmtmb:::prior_augmented_obj()`) and evaluates `obj$fn` directly, which
separates taping from sampling. It probes the priored coordinate at four
values, so a prior taped to a constant would show as a constant
difference.

Model B, `y ~ x + (1 | g)`, `prior_normal(0, 0.01)` on `theta_1`,
`aug - plain`:

| `theta_1` | Container | Host |
| --- | --- | --- |
| -1.0 | 4996.313768347216 | 4996.313768347216 |
| 0.0 | -3.686231652783 | -3.686231652783 |
| 0.5 | 1246.313768347217 | 1246.313768347217 |
| 2.0 | 19996.313768347220 | 19996.313768347220 |

The two platforms agree to the last printed digit (the final digit at
theta = 2.0 is one ulp at magnitude 2e4 and moves with the frmtmb build;
review re-measured the host at ...220), the difference varies
with the parameter, and it equals the untaped R-side
`neg_log_prior_fn()` at the same points. The gradient gains exactly
`5000` in the priored coordinate and nothing elsewhere. `RTMB`, the
prior resolution and the tape are all correct in the container.

So the fault is between `MakeADFun()` and the draws.

## Step 2. `tmbstan` ignores the objective

`dev/pdi-tmbstan.R` removes `frmtmb`. It builds two one-parameter `RTMB`
objectives that differ only by a tight `N(0, 0.05)` prior term and samples
both with the same seed. Analytic posteriors: flat is `N(3.4432, 0.2236)`,
prior is `N(0.1640, 0.0488)`.

| | flat mean | flat sd | prior mean | prior sd |
| --- | --- | --- | --- | --- |
| Analytic | 3.44324636708 | 0.22360679775 | 0.163964112718 | 0.0487950036474 |
| Host | 3.44372398232 | 0.216727529914 | 0.15819851925 | 0.0453407835105 |
| Container | 0.0362052934462 | 1.02771600421 | 0.0362052934462 | 1.02771600421 |

The container returns the same chain for both objectives, with mean 0 and
sd 1. It is sampling a standard normal. `tmbstan(debug = TRUE)`, which
keeps the R closures instead of the `TMB` tape pointer, gives the same
wrong answer, so the fault is below both plugin routes.

## Step 3. The root cause

`tmbstan` has no Stan model of its own. `inst/model.stan` is a
placeholder:

```stan
parameters { vector<lower=...,upper=...>[N] y; }
model { y ~ std_normal(); }
```

`tools/autogen.R`, run by `configure` at install time, translates it with
`rstan::stanc()` and then rewrites the placeholder line to call the `TMB`
objective:

```r
pattern <- "lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));"
replace <- "lp_accum__.add(custom_func::custom_func(y));"
searchReplace(pattern, replace)
```

`searchReplace()` takes the FIRST match:

```r
i <- grep(pattern, mod, fixed = TRUE)
stopifnot(length(i) >= 1)
if (length(i) > 1) warning("More than one match; Using first")
i <- i[1]
```

`stanc` 2.32.2 emits ONE `log_prob_impl`, so the first match is the only
match. `stanc` 2.39.0 emits TWO, one per scalar type:

| Overload | Guard | Used for | Patched |
| --- | --- | --- | --- |
| double | `require_not_st_var<VecR>` | value only | yes, line 226 |
| reverse-mode | `require_st_var<VecR>` | HMC value AND gradient | NO, line 273 |

`autogen.R` patches the double overload and leaves the reverse-mode one on
`std_normal_lpdf`. NUTS evaluates the log density through the reverse-mode
overload, so it samples the placeholder.

Running the stock `tools/autogen.R` in the container prints the warning
that says so, then the build continues:

```
Warning messages:
2: In searchReplace(pattern, replace) : More than one match; Using first
```

Generated `model.hpp` after that run: one `custom_func` call and one
surviving `std_normal_lpdf`. On this machine, `stanc` 2.32.2: one
`custom_func` call and zero survivors.

`pak` and `R CMD INSTALL` do not surface a `configure` warning, so the
install log is clean and the package works, for the wrong density.

### Which `stanc`

`rstan` does not carry its own `stanc`. Its `.onLoad()` reads
`stanc.js` out of `StanHeaders` AT RUN TIME:

```r
if (packageVersion("StanHeaders") == "2.26.28") {
  stanc_js <- system.file("exec", "stanc.js", package = "rstan", ...)
} else stanc_js <- system.file("stanc.js", package = "StanHeaders", ...)
```

So the `stanc` that generates `tmbstan`'s `model.hpp` is whatever
`StanHeaders` was installed when `tmbstan` was compiled. The `rstan`
version is not involved: both machines run `rstan` 2.32.7.

| | `StanHeaders` | its `stanc.js` | `rstan::stan_version()` |
| --- | --- | --- | --- |
| Container | 2.39.1 | 2.39.0 | 2.39.0 |
| This machine | 2.32.10 | 2.32.2 | 2.32.2 |

CRAN has no `StanHeaders` between 2.32.10 (2024-07-15) and 2.39.1
(2026-09-02), so the boundary is that single release.

## Step 4. The proof

Change `tools/autogen.R` so the log-density patch applies to every
remaining match, then rebuild `tmbstan` in the SAME container, against the
same `StanHeaders` 2.39.1, the same `rstan` and the same `stanc` 2.39.0:

```r
mod <- sub("lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));",
           "lp_accum__.add(custom_func::custom_func(y));", mod, fixed = TRUE)
```

`dev/pdi-tmbstan.R` against that build:

| | flat mean | flat sd | prior mean | prior sd |
| --- | --- | --- | --- | --- |
| Container, stock | 0.0362052934462 | 1.02771600421 | 0.0362052934462 | 1.02771600421 |
| Container, patched | 3.44372398232 | 0.216727529914 | 0.15819851925 | 0.0453407835105 |
| This machine | 3.44372398232 | 0.216727529914 | 0.15819851925 | 0.0453407835105 |

The patched container reproduces this machine to every printed digit. One
line of `autogen.R` is the whole difference.

The same flip at the `frmtmb` level, `dev/pdi-frm-sample.R`, default
priors against `priors = "flat"` on `y ~ x + (1 | g)`:

| Parameter | Stock sd (default) | Stock sd (flat) | Patched sd (default) | Patched sd (flat) |
| --- | --- | --- | --- | --- |
| `Intercept` | 0.97963927 | 0.97963927 | 0.20562114 | 0.20096960 |
| `x` | 0.92835674 | 0.92835674 | 0.09968708 | 0.08962403 |
| `sigma_Intercept` | 0.95997532 | 0.95997532 | 0.14189301 | 0.13151154 |
| `theta_1` | 1.01031905 | 1.01031905 | 0.65236765 | 0.39683889 |

Stock: every posterior sd is 1, and `identical(default, flat)` is `TRUE`.
Patched: the sds are on the scale of the model, and
`max |default - flat|` is 32.59.

## The version diff

471 packages are installed in both places. 27 versions differ, 21 of them
newer in the container:

```
class 7.3-23/7.3-24      cluster 2.1.8.2/2.1.8.3   frmtmb 0.28.0/0.40.1
KernSmooth 2.23-26/-27   lattice 0.22-9/0.23-1     MASS 7.3-65/7.3-66
Matrix 1.7-5/1.7-6       nlme 3.1-169/3.1-170      nnet 7.3-20/7.3-21
pkgcheck 0.2.0.003/.012  spatial 7.3-18/7.3-19     survival 3.8-6/3.8-11
StanHeaders 2.32.10/2.39.1
```

`RTMB` 1.9, `TMB` 1.9.25, `rstan` 2.32.7 and `tmbstan` 1.2.0 are equal on
both, as the earlier report claimed. `StanHeaders` is the only numerically
relevant delta, and Step 4 proves it is the one. `Matrix` 1.7-6 is
excluded by the same proof: it did not change when the bug flipped.

`OpenBLAS` against reference BLAS is not the cause either. It cannot
produce a chain that agrees with a standard normal to three digits.

## Who is affected

The trigger is the `StanHeaders` version present when `tmbstan` is
COMPILED, so the split is between source installs and binary installs.

| Install | `tmbstan` `model.hpp` from | Affected |
| --- | --- | --- |
| Linux, any CRAN install since 2026-09-02 | local `stanc` 2.39.0 | YES |
| macOS or Windows, `type = "source"`, since 2026-09-02 | local `stanc` 2.39.0 | YES |
| macOS or Windows CRAN binary built before 2026-09-02 | CRAN's `stanc` 2.32.2 | no |
| Any machine pinned to an RSPM snapshot before 2026-09-02 | `stanc` 2.32.2 | no |

CRAN publishes no Linux binaries, so every Linux user who installs
`tmbstan` from now on gets the broken build. macOS and Windows binaries
protect their users only until CRAN rebuilds them against 2.39.1.

That explains the whole platform table the previous lane assembled:

- The pkgcheck container: `pak` takes newest CRAN, Linux, source. Broken.
- This machine: `StanHeaders` 2.32.10, `tmbstan` binary. Correct.
- The Ubuntu `R-CMD-check` runner: RSPM snapshot, `StanHeaders` 2.32.10.
  Correct.
- The macOS `R-CMD-check` runner: `StanHeaders` 2.39.1 but a `tmbstan`
  BINARY whose `model.hpp` was generated earlier. Correct, and only by
  timing.

The blast radius is not `frmtmb`. It is every `tmbstan` user: plain `TMB`
models are affected identically, and the failure is silent.

## Whose bug it is, and the fix

Upstream, `kaskr/tmbstan`, `tools/autogen.R`. `frmtmb` does everything
right. There is no open issue for this; the newest open issue in that
repository is #31, from 2026-06-04.

The one-line fix, verified in Step 4, is to make the log-density
replacement apply to every match rather than the first.

### What `frmtmb` should do, but not in this lane

`frmtmb` cannot fix `tmbstan`, but it must not return standard-normal
draws as a posterior. Two guards were tested in `dev/pdi-guard.R`, and
both work. Do NOT implement either here.

1. Static, and preferred. `tmbstan` installs its generated `model.hpp`
   beside itself, so the miscompilation is readable without sampling:

   ```r
   src <- readLines(system.file("model.hpp", package = "tmbstan"))
   sum(grepl("std_normal_lpdf<propto__>(y)", src, fixed = TRUE)) > 0
   ```

   Container: 1, so miscompiled. This machine: 0. The check costs one
   file read, has no side effects, and names the exact defect, so
   `frm_sample()` can refuse with a message that tells the user to
   reinstall `tmbstan` after downgrading `StanHeaders`, or to wait for the
   upstream fix.

2. Numeric, for the test suite. Compare Stan's own log density against the
   objective on the returned `stanfit`. `rstan::log_prob()` and
   `rstan::grad_log_prob()` both route through the reverse-mode overload,
   so both discriminate:

   | `mu` | `stan_lp` | `-obj$fn` | `stan_gr` | `-obj$gr` |
   | --- | --- | --- | --- | --- |
   | 0 | 0.00000000 | -151.49847477 | 0.00000000 | 68.86492734 |
   | 1 | -0.50000000 | -92.63354742 | -1.00000000 | 48.86492734 |
   | 2 | -2.00000000 | -53.76862008 | -2.00000000 | 28.86492734 |

   The container's `stan_lp` is `-mu^2 / 2` exactly. On this machine the
   columns agree to the last digit.

`tests/testthat/helper-reference.R` also needs its explanation corrected:
the chain-agreement gates are not hiding seeded-chain irreproducibility,
and `test-reparam.R:627` was a true positive.

## Reproducing this

```
# the discriminating measurement, host and container
Rscript dev/pdi-repro.R
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/adf44/source/r/frmtmb-wt-dockergate:/work:ro" -w /work \
  -e R_PROFILE_USER= --entrypoint Rscript frmtmb-pkgcheck:deps \
  /work/dev/pdi-repro.R

# tmbstan alone, no frmtmb
Rscript dev/pdi-tmbstan.R          # and the same docker run

# the end-to-end symptom
Rscript dev/pdi-frm-sample.R

# the guards
Rscript dev/pdi-guard.R
```

`-e R_PROFILE_USER=` is required. `.github/pkgcheck.Rprofile` ends with
`Sys.setenv(FRMTMB_SAMPLER_GATES = "false")`, which runs after `-e` and
wins. See `dev/pkgcheck-docker.md`.

The rebuild in Step 4 needs the `tmbstan` source, which the warm
`frmtmb-pkgcheck-cache` volume already holds at
`/root/.cache/R/pkgcache/pkg/src/contrib/tmbstan_1.2.0.tar.gz`. Mount that
volume, patch `tools/autogen.R`, then `R CMD INSTALL --library=`
somewhere off the default path.

## Issue skeleton for kaskr/tmbstan

> Title: `tmbstan` silently samples `std_normal()` instead of the `TMB`
> objective when built against `StanHeaders` 2.39.1
>
> `tools/autogen.R` rewrites the `std_normal_lpdf` placeholder in the
> `stanc` output into a call to `custom_func::custom_func()`.
> `searchReplace()` replaces only the first match and warns
> "More than one match; Using first" for the rest.
>
> `stanc` 2.32.2 emitted one `log_prob_impl`. `stanc` 2.39.0, which
> `StanHeaders` 2.39.1 ships as `stanc.js` and `rstan` sources at run
> time, emits two: a `require_not_st_var<VecR>` overload for doubles and a
> `require_st_var<VecR>` overload for reverse-mode `var`. `autogen.R`
> patches the first and leaves the second on the placeholder. HMC uses the
> reverse-mode overload for both the value and the gradient, so the
> sampler explores `y ~ std_normal()` and never touches the `TMB` model.
>
> There is no error. `configure` succeeds, the package loads, chains run
> and `rstan` reports the usual diagnostics. The draws are a standard
> normal in the unconstrained space.
>
> Reproduce, with `StanHeaders` 2.39.1 and a source install of `tmbstan`
> 1.2.0:
>
> ```r
> library(RTMB); library(tmbstan)
> set.seed(7); y <- rnorm(20, 3, 1)
> obj <- MakeADFun(function(p) -sum(dnorm(y, p$mu, 1, log = TRUE)),
>                  list(mu = 0), silent = TRUE)
> obj$fn(obj$par)
> s <- tmbstan(obj, chains = 1, iter = 1000, seed = 42, init = "0")
> mean(extract(s)$mu); sd(extract(s)$mu)
> # broken : 0.036, 1.028   (should be 3.444, 0.217)
> rstan::log_prob(s, 2)     # -2, that is -mu^2/2, not -53.77
> ```
>
> Static check on an installed copy:
>
> ```r
> sum(grepl("std_normal_lpdf<propto__>(y)", fixed = TRUE,
>           readLines(system.file("model.hpp", package = "tmbstan"))))
> # 1 when miscompiled, 0 when correct
> ```
>
> Fix, verified: make the log-density replacement apply to every match.
> Adding this before `## Write` in `tools/autogen.R` restores the correct
> posterior under `StanHeaders` 2.39.1, byte for byte against a 2.32.10
> build:
>
> ```r
> mod <- sub("lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));",
>            "lp_accum__.add(custom_func::custom_func(y));", mod,
>            fixed = TRUE)
> ```
>
> A cheap regression test would be to fail the build when any
> `std_normal_lpdf<propto__>(y)` survives in the generated `model.hpp`,
> and to turn the "More than one match" warning into an error for this
> pattern.
>
> Affects every Linux install from CRAN since `StanHeaders` 2.39.1
> (2026-09-02), and any source install on macOS or Windows. Binary
> installs still carry a `model.hpp` generated by `stanc` 2.32.2.
