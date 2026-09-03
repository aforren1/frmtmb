# Run pkgcheck in its container before you push

`.github/workflows/pkgcheck.yaml` runs rOpenSci's `pkgcheck` on every push to
`main`, and the run takes 35 to 50 minutes. Run the same container on this
machine first. The container is the pre-push gate: it returns the same check
lines as the workflow, so a check that fails here will fail there.

This is a how-to guide. For why the chain-agreement gates are off in that
workflow, read `tests/testthat/helper-reference.R` and the "Sampler gates"
section below.

## Before you start

- Docker Desktop must be running. Start it and wait for the daemon:
  `until docker info >/dev/null 2>&1; do sleep 5; done`.
- Pull the image once. It is 21.3 GB, and a full run then commits a
  `frmtmb-pkgcheck:deps` image on top of it. Keep 30 GB free. The pull took
  about 50 minutes here:

  ```
  docker pull ghcr.io/ropensci-review-tools/pkgcheck-action:latest
  ```

- `gh auth status` must succeed. `pkgcheck` asks GitHub whether the package
  name is free and whether continuous integration passes.

## Run the gate

```
dev/run-pkgcheck-docker.sh
```

The script clones the committed `HEAD` into a temporary directory and mounts
that clone, so commit your work first. It writes a timestamped log and a copy
of the check summary to `dev/`, which `dev/.gitignore` keeps out of git.

Measured on 2026-09-03 at `0afcd27`: 36 minutes in total, of which 21 minutes
went to `pak`, 5 minutes to `covr` and 9 minutes to `rcmdcheck`. A warm
`frmtmb-pkgcheck-cache` volume removes most of the `pak` time.

Read the `Check Results` block at the end of the log, or the copied
`dev/pkgcheck-docker-<stamp>-summary.md`. Every line that starts with a cross
is a failure that the workflow turns into a red run.

### It agrees with continuous integration

The container run at `0afcd27` returned the same 17 lines as the GitHub run on
the same commit, in the same order, from the same `pkgcheck` 0.2.0.12:

```
✔ Package name is available          ✔ Package has continuous integration checks.
✔ has a 'contributing' file.         ✖ Package coverage failed
✔ uses 'roxygen2'.                   ✔ complies with all applicable standards
✔ 'DESCRIPTION' has a URL field.     ✖ Not all authors have ORCID IDs
✔ 'DESCRIPTION' has a BugReports     ✖ R CMD check found 1 error.
✔ at least one HTML vignette         ✔ R CMD check found no warnings.
✔ All functions have examples.       ℹ Some goodpractice linters failed.
✔ Repository has a website           ℹ Function names duplicated in other packages
                                     ℹ Examples should not use `\dontrun`
```

The `R CMD check` error is the same failing test file in both. The test counts
differ by one: the container gave `FAIL 9 | SKIP 14 | PASS 5559`, and GitHub
gave `FAIL 8 | SKIP 14 | PASS 5560`. The extra local failure is
`test-reparam.R:627`, "a prior on an sd applies identically under both routes",
which is not a gated test. See the next section for why that one matters.

Treat a one-test disagreement as normal and a changed check line as a real
difference.

## The invocation

The script builds this command. Use it directly if you must:

```
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/path/to/staged/frmtmb:/github/workspace" \
  -v "frmtmb-pkgcheck-cache:/root/.cache" \
  -w /github/workspace \
  -e R_PROFILE_USER=.github/pkgcheck.Rprofile \
  -e FRMTMB_SAMPLER_GATES=false \
  -e NOT_CRAN=true -e CI=true -e ROPENSCI=true \
  -e GITHUB_TOKEN="$(gh auth token)" \
  ghcr.io/ropensci-review-tools/pkgcheck-action:latest
```

### Windows path rules

Git Bash rewrites arguments that look like absolute POSIX paths. Without
`MSYS_NO_PATHCONV=1`, `-w /github/workspace` becomes
`C:/Program Files/Git/github/workspace` and the daemon refuses it. Two
spellings work:

| Form | Works |
| --- | --- |
| `-v /c/src/frmtmb:/github/workspace -w /github/workspace` | No. MSYS mangles the container path. |
| `MSYS_NO_PATHCONV=1 ... -v "//c/src/frmtmb:/github/workspace" -w //github/workspace` | Yes |
| `MSYS_NO_PATHCONV=1 ... -v "C:/src/frmtmb:/github/workspace" -w /github/workspace` | Yes. Use this one. |

With `MSYS_NO_PATHCONV=1` set, write the host side as a native Windows path
and the container side with one leading slash.

### Which image

Use `pkgcheck-action:latest`, not `pkgcheck:latest`. The action image adds
three things that `pkgcheck:latest` does not have:

- the `pkgcheck` R package itself, installed by the action's `install.R`;
- `/check.R`, which is the `ENTRYPOINT`, run as `Rscript /check.R`;
- `octolog`, which formats the results.

`pkgcheck:latest` is only the base layer. It has the system libraries, R,
`rstan`, `cmdstan` and Julia, but no `pkgcheck`.

### What the entrypoint does

`/check.R` runs in the working directory, so the working directory must be
the package root. In order:

1. `pak::lockfile_create("local::.", dependencies = "all")`, then
   `lockfile_install`. This installs every dependency and `frmtmb` itself.
2. `usethis::use_build_ignore(".pkgcheck")`, which edits `.Rbuildignore`.
3. `pkgstats::ctags_install(sudo = TRUE)`.
4. `pkgcheck()`, which runs `covr`, `cyclocomp`, `lintr`, `rcmdcheck`,
   `roxygen2`, `srr` and the rest.
5. `summary.md` and `full.md`, written to the workspace root and copied to
   `.pkgcheck/`.

Steps 2 and 5 write into the mounted directory. This is why the script mounts
a throwaway clone and not the real tree.

### Environment variables

| Variable | Set by | Why |
| --- | --- | --- |
| `NOT_CRAN`, `CI`, `ROPENSCI` | the image, as `ENV` | Already `true`. The script repeats them so the invocation is self-describing. |
| `GITHUB_TOKEN` | the action, from `github.token` | The name-available and continuous-integration checks call the GitHub API. |
| `R_PROFILE_USER` | `.github/workflows/pkgcheck.yaml` | Points at `.github/pkgcheck.Rprofile`. See below. |
| `FRMTMB_SAMPLER_GATES` | the same workflow | `false` mirrors continuous integration. |
| `CMDSTAN_PATH` | the image | `/root/.cmdstan`. `frmtmb` does not use it. |

### How RTMBode resolves

`RTMBode` is a `Suggests` that only exists on r-universe, and `pak` does not
read `Additional_repositories` from `DESCRIPTION`. `.github/pkgcheck.Rprofile`
adds `https://kaskr.r-universe.dev` to `options(repos)`, and `pak` forwards the
session's `repos`. The container then gets `RTMBode` as a source package and
compiles it, which takes a few seconds. Verified in the 2026-09-03 run:
`Got RTMBode 1.0 (source) (162.98 kB)` then `Installed RTMBode 1.0`, from
`https://kaskr.r-universe.dev`. The trick is still needed; nothing in the image
resolves `RTMBode` on its own.

`R_PROFILE_USER` is a path relative to the working directory, and the working
directory is the package root, so the relative path in the workflow works
unchanged in a local run.

Confirm the resolution in the log:

```
grep -a RTMBode dev/pkgcheck-docker-*.log
```

You want `Got RTMBode 1.0 (source)` and `Installed RTMBode 1.0`. If instead
you see `Can't find package called RTMBode`, the profile did not load: check
that `-e R_PROFILE_USER=.github/pkgcheck.Rprofile` is present and that the
working directory is the package root.

## Sampler gates

The chain-agreement assertions compare NUTS output against a reference
quantity. `sampler_gates_on()` in `tests/testthat/helper-reference.R` turns
them off when `FRMTMB_SAMPLER_GATES=false`, which the `pkgcheck` workflow
sets. Run the tests with the gates on to study chain behavior:

```
dev/run-pkgcheck-docker.sh --tests 'v07|v15|sample-direct' --gates on --repeat 2
```

This needs a full run first, because it reuses the `frmtmb-pkgcheck:deps`
image that the full run commits. It does not run `pkgcheck` or
`R CMD check`; it runs `testthat::test_dir` against the installed package.

Two traps make a gates-on run silently become a gates-off run. Both are
handled by the script, and both will catch you if you write the `docker run`
by hand:

- `.github/pkgcheck.Rprofile` ends with
  `Sys.setenv(FRMTMB_SAMPLER_GATES = "false")`, which runs after `-e` and wins.
  A gates-on run must pass `-e R_PROFILE_USER=` to suppress the profile.
- `docker commit` copies the committed run's `-e` values into the new image's
  config. A `frmtmb-pkgcheck:deps` image built without care therefore carries
  `R_PROFILE_USER`, the gate switch, and your `GITHUB_TOKEN`. The script clears
  all three with `--change` at commit time. Check any image you did not build
  this way:

  ```
  docker inspect frmtmb-pkgcheck:deps --format '{{range .Config.Env}}{{println .}}{{end}}' | grep -iE 'TOKEN|PROFILE|FRMTMB'
  ```

### The container does reproduce the drift, and it is not chain luck

Run on 2026-09-03 at `0afcd27`, `test-v07`, `test-v15` and `test-sample-direct`
with the gates on. Result: 6 failures, and the same 6 with the same numbers on
a second run in the same container.

| Quantity | Gate | Container | This machine |
| --- | --- | --- | --- |
| `sd_ratio` for `x`, `test-v07.R:117` | \|r - 1\| < 0.5 | 10.904044 | 0.992817 |
| `z_shift` for `x` | \|z\| < 0.75 | -0.504778 | -0.070278 |
| `ess_bulk` for `x` | >= 100 to judge | 743.14 | 332.71 |
| posterior mean of `x` (ML gives 0.454158) | | -0.037404 | 0.447927 |
| mixture separation, `test-sample-direct.R:477` | > 2 | 0.047361 | 3.863931 |
| flat `theta_1` 2.5 percent quantile, `:440` | < -3.5 | -1.935587 | -5.135134 |
| sd(flat) / sd(default) `theta_1`, `:443` | > 1.3 | 1.000000 | 1.361089 |

Every container number repeats exactly across runs, so the container is not
drawing a random chain each time. It draws one particular chain, always the
same one, and that chain is not the one every other platform draws.

The last two rows say what is really happening. Fit the same model twice, once
with `priors = "flat"` and once with the default half-t on the standard
deviation, and compare the `theta_1` draws:

```
container:  max |flat - default| = 0          identical(flat, default) = TRUE
this machine: max |flat - default| = 6.109209 identical(flat, default) = FALSE
```

`prior_summary()` reports `student_t(3, 0, 2.5) class=sd` for the default fit
in both places, so the prior is declared in both. The follow-up investigation
(`dev/prior-dropping-investigation.md`) sharpened this measurement's reading:
it is not that the prior is dropped, it is that NOTHING reaches the sampler.
A tmbstan built under StanHeaders 2.39 samples a standard normal instead of
the model, likelihood and priors alike, because its install-time generator
patches only the first of two generated log-density overloads. The flat chain
and the half-t chain are one chain because both are that standard normal.

That explains the rest of the table: a standard normal's spread is 10.9 times
this model's Wald standard error, and the mixture "components" are two reads
of the same N(0, 1) coordinate.

So the gates were not hiding a seeded-chain reproducibility problem. They were
partly masking a real upstream defect, which frm_sample() now refuses
statically. The non-gated `test-reparam.R:627` fails in the container for the
same reason, and no switch hides that one.

Do not read a green `pkgcheck` run as evidence that the priors work. Use
`--gates on` when you touch prior or reparameterization code.

### The numeric stack, container against this machine

Measured on 2026-09-03 with `pkgcheck-action` image `411580451793`.

| | Container | This machine | ubuntu R-CMD-check | macOS R-CMD-check |
| --- | --- | --- | --- | --- |
| R | 4.6.1, `x86_64-pc-linux-gnu` | 4.6.1, `x86_64-w64-mingw32` | 4.6.1 | 4.6.1 |
| BLAS | OpenBLAS pthread 0.3.26 | R reference BLAS | OpenBLAS | Accelerate |
| LAPACK | 3.12.0, from OpenBLAS | 3.12.1, R internal | | |
| `rstan` | 2.32.7 | 2.32.7 | 2.32.7 | 2.32.7 |
| `StanHeaders` | 2.39.1 | 2.32.10 | 2.32.10 | 2.39.1 |
| `tmbstan` | 1.2.0, built in the run | 1.2.0 | 1.2.0 | 1.2.0 |
| `RTMB` | 1.9 | 1.9 | 1.9 | 1.9 |
| Cores | 8 | 16 | 4 | 3 |

Two differences stand out. The container gets a threaded OpenBLAS, and this
machine gets R's single-threaded reference BLAS. The container also builds
`tmbstan` against `StanHeaders` 2.39.1, which is a different Stan Math
release from the 2.32.10 that this machine and the Ubuntu runners hold.

`StanHeaders` at tmbstan COMPILE time is the whole explanation, which the
follow-up investigation proved by a one-line patch flip: the macOS job also
holds 2.39.1 but installs a prebuilt tmbstan binary generated under stanc
2.32.2, so its runtime StanHeaders version is irrelevant. `pak` always takes
the newest CRAN build and compiles from source on Linux, so the container is
exactly the affected configuration. Re-read this table after any `pak`
install; it will change.

## Container against action: the deltas

- The action pins `docker://ghcr.io/ropensci-review-tools/pkgcheck-action:latest`,
  so both track the same moving tag. A local image can be older than the one
  the runner pulls. Re-pull before you trust a clean local run.
- GitHub mounts the workspace at `/github/workspace` and sets `HOME` to
  `/github/home`. Locally `HOME` stays `/root`, so `pkgcheck`'s cache lands in
  `/root/.cache/R/pkgcheck` instead. The script mounts a named volume there so
  the cache survives between runs.
- The runner gives the job 4 cores and about 16 GB. Docker Desktop on this
  machine gives 8 cores and 15.6 GB, and `pak` uses all of them: `docker stats`
  showed 657 percent CPU during the dependency build. More cores make the local
  run faster, and they also make the BLAS thread pool wider, which is a
  numerical difference and not only a speed difference.
- The action's last step is `exit ${{ steps.pkgcheck.outputs.status }}`, and
  `status` counts only the crossed lines in the summary. Read the log, not the
  exit code.
- `octolog` reports a failed check as an Actions annotation and continues only
  when `GITHUB_ACTIONS` is set. Without it, `/check.R` aborts at the first
  failed check and never prints `Check Results`, although it has already
  written `summary.md` and `full.md`. The script sets `GITHUB_ACTIONS=true` for
  this reason.
- The action's `Fail if pkgcheck found problems` step is the only place the run
  turns red. Nothing in the container does that, so a local run that ends in
  `Execution halted` is not by itself worse than a green GitHub run.

## Clean up

```
docker image rm frmtmb-pkgcheck:deps
docker volume rm frmtmb-pkgcheck-cache
```
