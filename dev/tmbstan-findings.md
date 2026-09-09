# Lane tmbstan: the gated sampler tier on a broken tmbstan

Worktree `frmtmb-wt-tmbstan`, branch `wt-tmbstan`, base 92e9330.
Nothing was committed.

The red CI run of 2026-09-09 was `check_tmbstan_build()` doing its job.
The guard is correct and was not changed. What this lane changed is
what a TEST does when it meets such a build, and what CI installs.

## Summary of the numbers

| Measurement | Value |
| --- | --- |
| `test_that` blocks in the frmtmb.sample suite | 192 (was 188) |
| blocks that reach `frm_sample()` or `as_tmbstan()` | 111 |
| of those, gated by the shared `skip_sampler()` BEFORE | 30 |
| gated by a file-LOCAL shadowing copy of it | 47 |
| gated only by `skip_if_not_installed("tmbstan")` | 34 |
| gated by the shared `skip_sampler()` AFTER | 111 |
| blocks that skip on a broken build, after | 114 |
| blocks that still FAIL loudly on a broken build | 1 |

## The green-while-broken claim is TRUE, with artifact evidence

Stated in the brief as something to verify rather than repeat. It
holds, and the evidence is the artifact rather than an inference.

**The last green run of this job sampled a standard normal.** Run
`34258865676`, 2026-09-08 17:44 UTC, frmtmb.sample at 0.55.0:

    [ FAIL 0 | WARN 0 | SKIP 10 | PASS 890 ]     Status: OK

None of its 10 skips concerns the build; they are Linux-only, gates
off, brms-fit and scale-tier gates. It resolved
`tmbstan_1.2.0.tar.gz` from
`https://packagemanager.posit.co/cran/__linux__/noble/latest`.

**That artifact is broken, and it is a BINARY.** Downloaded and read
directly (headers `X-Package-Type: binary`,
`X-Package-Binary-Tag: 4.6-noble`):

| RSPM path | `Built:` | unpatched `std_normal_lpdf<propto__>(y)` | `log_prob_impl` |
| --- | --- | --- | --- |
| `__linux__/noble/latest` | 2026-09-03 08:40:31 UTC | **1** | 4 |
| `__linux__/noble/2026-09-02` | 2026-08-28 08:22:33 UTC | 0 | 3 |
| `__linux__/noble/2026-09-01` | 2026-08-28 08:22:33 UTC | 0 | 3 |
| `__linux__/noble/2026-08-20` | 2026-07-31 07:43:09 UTC | 0 | 3 |

RSPM REBUILT the noble R-4.6 binary on 2026-09-03, the day after
StanHeaders 2.39.1 reached CRAN, and the rebuild carries the defect:
two `log_prob_impl` overloads at lines 193 and 242, the first patched
to `custom_func::custom_func(y)` at line 226, the second left with
`stan::math::std_normal_lpdf<propto__>(y)` at line 273.

So the 2026-09-08 green run and the 2026-09-09 red run downloaded the
same broken binary. The only difference between them is that the guard
had shipped.

**And garbage satisfies the structural assertions.** Measured, not
argued. `dev/tmbstan-probe-06-garbage.R`, noise SEED 909, model and
seeds from `ce_case()` in `test-conditional-effects-draws.R`. A correct
`frm_sample()` run has its draws matrix replaced column by column with
standard normal noise, `exp()` applied to the columns whose correct
draws are strictly positive, which is what the defect produces: a
standard normal in the unconstrained space put through Stan's
constraint transforms.

| block | on garbage draws |
| --- | --- |
| the draws frame has the fit frame's columns and grid | PASS |
| the default effect list and its grids match the fit method | PASS |
| the bands exist, are finite and are ordered | PASS |
| the draws accessor surface answers | PASS |
| **the drawn curve tracks the ML curve** | **FAIL** |

Posterior means, correct against garbage, same fit:

    Intercept      x       z   sigma_Intercept
    0.871      0.565   0.147   0.078            correct
    -0.055     1.733  -0.028  -0.023            garbage

Four of four structural blocks pass on noise. Only the correctness
assertion discriminates, and this workflow sets
`FRMTMB_SAMPLER_GATES: "false"`, which is what turns the
chain-agreement assertions off. That is the whole mechanism of the
green run.

## Half one, the local half

### What changed

`skip_sampler()` in
`extensions/frmtmb.sample/tests/testthat/helper-sampling.R` now asks
whether the sampler WORKS, not only whether it is installed. It calls
the package's own detector, `frmtmb.sample:::tmbstan_build_broken()`,
rather than reimplementing the file read, so the skip and the refusal
cannot drift apart.

Then every block that reaches a sampler was put behind it. Three
routes existed and all three are gone:

- **A file-LOCAL `skip_sampler()`** shadowed the shared one in
  `test-draws-methods.R`, `test-loo.R` and `test-reparam.R`. Three
  identical two-line copies, covering 47 sampling blocks, that a grep
  for `skip_sampler` makes look covered. Removed; the files now
  inherit the helper.
- **An inline `skip_if_not_installed("tmbstan")` pair**, 34 sampling
  blocks in `test-sample-direct.R` (11), `test-sampling-ported.R`
  (21), `test-parallel-chains.R` (1) and `test-scale.R` (1).
  Converted.
- `test-conditional-effects-draws.R`, the file named in the brief,
  already used the shared helper. Its 12 blocks were fixed by the
  helper change alone.

Four sites in `test-sampling-ported.R` were not simple adjacent pairs
and were converted by hand: the `as_tmbstan` block whose `rstan` skip
sat eight lines lower, the two `frmtmb.latent` blocks, and the
brms-agreement block with no `rstan` skip at all.

Six blocks that use `skip_sampler()` do NOT sample: they assert a
REFUSAL, and `sample_preflight()` runs before `check_tmbstan_build()`
in both doors, so they would still pass on a broken build. They now
skip with everything else. That is a deliberate loss of six refusal
assertions in an environment that is already red, and it buys one
gate instead of two. The refusal blocks that were NOT already on
`skip_sampler()` were left on the plain installed check:
`test-parallel-chains.R:40`, `test-sample-direct.R:568`,
`test-sampling-ported.R:687` and `:869`.

### The classification was static, not guessed

`dev/tmbstan-probe-02-classify.R` parses the suite with srcrefs, walks
each block's call tree, takes a fixpoint over the top-level function
definitions so a block calling a cached case closure counts as
sampling, honors a file-local definition shadowing a helper one, and
does not count a sampler mention inside `expect_error()`. Before: 111
sampling blocks, 30 on the helper, 47 on a local copy, 34 inline, 0
ungated. After: 111 on the helper, 0 elsewhere.

### Constructing the case where the guarded thing is ABSENT

Four seams were tried. Two work and two do not, and the two that fail
are worth recording because the obvious one is among them.

1. `testthat::local_mocked_bindings(system.file =, .package =)`.
   FAILS with "Can't find binding for `system.file`". testthat 3.3.2
   will not create a binding that does not already exist in the
   namespace, and a base function imported implicitly has none.
2. Prepend a synthesized tmbstan to `.libPaths()`. Works ONLY while
   tmbstan is not loaded. Once it is, which it is by the time any of
   these tests run, `system.file()` resolves from the loaded namespace
   and ignores `.libPaths()`. Measured both ways in
   `dev/tmbstan-probe-05-seam2.R`, which prints the fake path in the
   first case and the real one in the second.
3. Bind `system.file` in the detector's own `local()` environment.
   WORKS. Used by `with_tmbstan_hpp()`.
4. Poke `cached` in the same environment. WORKS. Used to run whole
   files in the broken arm.

`tmbstan_build_broken()` is built by `local({...})`, so its enclosing
environment is not the namespace and is not locked. A `system.file`
bound there shadows base's for the one call the detector makes, and
the memo lives in the same place. `with_tmbstan_hpp()` in
`helper-sampling.R` writes a synthesized `model.hpp` with the real
two-overload shape, points the detector at it, and restores both
bindings on exit. Nothing is installed, `.libPaths()` is untouched,
and the real tmbstan is not disturbed.

`test-tmbstan-build-guard.R` (new) uses it for five blocks:

1. this installation's tmbstan is healthy. A CANARY, not a skip: it
   FAILS on a broken build, which is what keeps the tier from going
   quiet.
2. the broken `model.hpp` is detected and both callers refuse.
3. the patched `model.hpp` is NOT refused. The arm that matters more:
   a check that fires on a healthy build refuses every correct
   installation.
4. `skip_sampler()` skips on the broken build AND REACHES THE CODE
   AFTER THE SKIP on the clean one. A gate that skipped
   unconditionally would look identical on this machine without this
   second arm.
5. the detector reads `model.hpp`, not the runtime StanHeaders
   version.

It replaces a block that lived in `test-sampling-ported.R` and
asserted the refusal PATTERN by grepping a marker out of a string it
had just written. That is true of the literal and says nothing about
the guard; its own comment said the `stop()` branch "is not executed
here". It is executed now.

### Seen it fail

The BEFORE arm was reconstructed by reverting `skip_sampler()` to its
old two lines and running the file that broke CI, with the detector's
memo poisoned. `test-conditional-effects-draws.R`, one file per
process, script `dev/tmbstan-run.R`, seeds the file's own
(`set.seed(9)` and `set.seed(46)` for the data, sampler seeds 1 and 2):

| arm | blocks | passed | failed | error | skipped |
| --- | --- | --- | --- | --- | --- |
| BEFORE (old skip) x broken build | 12 | 0 | **0** | **12** | 0 |
| AFTER x broken build | 12 | 0 | 0 | 0 | **12** |
| AFTER x clean build | 12 | **57** | 0 | 0 | 0 |

The BEFORE row is the CI failure reproduced, and it is also the
counting trap from the lane rules in the flesh: `failed` is 0 for all
twelve. A runner that sums `failed` and not `error` calls that file
clean.

### The whole tier, both arms

Every touched file, one file per R process, `dev/tmbstan-sweep.sh`.
Thirteen files, the same 178 blocks in both arms.

| arm | blocks | passed | failed | error | skipped | elapsed |
| --- | --- | --- | --- | --- | --- | --- |
| clean | 178 | 1063 | 0 | 0 | 3 | 155 s |
| broken | 178 | 344 | 1 | 1 | 117 | 58 s |

The clean arm's 3 skips are pre-existing environment gates
(`FRMTMB_BRMS_FIT_TESTS` twice, `FRMTMB_SCALE_TESTS` once) and are
present in both arms, so 114 blocks skip because of the build. The
broken arm's single failing and erroring block is the canary, block 1
of `test-tmbstan-build-guard.R`: `expect_false()` fails and
`expect_silent()` errors, on purpose, in one place instead of 108.

The static count and the dynamic count agree exactly, which is the
check that the parser missed nothing. The classifier says 117 blocks
reach `skip_sampler()`: 111 that sample plus the 6 refusal blocks that
do not. The broken arm skips 117.

### R CMD check

One pass, on the built tarball, `--as-cran --no-manual`, with
`NOT_CRAN=true`, `FRMTMB_SAMPLER_GATES=false` and
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH:

    Status: OK
    [ FAIL 0 | WARN 0 | SKIP 11 | PASS 1026 ]

No ERROR, no WARNING, no NOTE. All 11 skips are pre-existing platform
and environment gates: chain-agreement gates off (4), brms-fit tier
(2), scale tier (1), RTMBode absent (1), an unhealthy-chain skip (1),
installed-package sources (1), and one empty test. None of them is the
build skip, because this machine's tmbstan is clean, which is the
point.

Checking the source DIRECTORY does not work and is worth recording:
it stops at `checking for file 'frmtmb.sample/DESCRIPTION' ... ERROR,
Required fields missing or empty: 'Author' 'Maintainer'`, because the
package declares `Authors@R` and those two fields are generated by
`R CMD build`. The tarball must be built first.

`frmtmb.eam`'s `test-sampling.R:10` is the same wall in another
package: it samples through `frmtmb.sample::frm_sample()` behind an
installed-only check. Guarded the same way. Clean arm 3 blocks, 97
passing, 0 skipped; broken arm 1 skipped, 94 passing, 0 errors.

### An instrument fault, corrected

The first sweep reported one error each in `test-loo.R` and
`test-sampling-ported.R` in BOTH arms. They were the runner's fault:
`test_file()` without `package =` parents the test environment on
`globalenv()`, so unqualified calls to the internals `stan_cores()`
and `loo_matrix()` cannot resolve, where `test_check()` does pass it.
With `package = "frmtmb.sample"` both files are clean. Every count in
this document is from the corrected runner.

## Half two, the CI half

### What the ubuntu job actually installs

The workflow's own header claimed its compiled dependencies "come as
binaries from the RSPM build". Half right, and the wrong half is the
one that matters. It IS a binary. The binary is broken, because RSPM
rebuilt it on 2026-09-03 against StanHeaders 2.39.1. Being a binary
protects nobody once the binary itself is regenerated.

### What was written

A step, `Pin a tmbstan whose generated model.hpp is patched`, between
the dependency install and `check-r-package`. It:

1. derives a dated snapshot from the repo the runner already has,
   `sub("/latest/?$", "/2026-09-01", getOption("repos")[["CRAN"]])`,
   so the distro codename is never hardcoded, and stops if the
   substitution is a no-op;
2. installs StanHeaders, rstan and tmbstan from it. The whole trio,
   not tmbstan alone, so the three stay the set RSPM built against
   each other. That snapshot has StanHeaders 2.32.10 (2.39.1 returns
   404 there), rstan 2.32.7 and the clean 2026-08-28 tmbstan binary;
3. reads the installed `model.hpp` and FAILS the job if it is still
   unpatched, with a message that says not to fix it by deleting the
   pin;
4. writes the pin, the three versions and the patched verdict to
   `$GITHUB_STEP_SUMMARY` on every run;
5. probes whether the pin is still needed, by downloading the
   UNPINNED artifact and reading its `model.hpp`. If upstream has been
   fixed it emits a `::warning::` annotation saying to delete the
   step. Non-fatal and wrapped in `try()`, so a network hiccup cannot
   fail a good run, but loud on every run so the pin cannot quietly
   outlive its reason. The probe is against the ARTIFACT because the
   version number never changed: tmbstan stayed 1.2.0 through the
   whole episode. Only the build date moved.

The probe's first draft was wrong in a way worth recording, because
it is the failure mode the lane rules name. It built the download URL
from the version it had just INSTALLED, which is the pinned one. The
day tmbstan is bumped, that URL 404s at the unpinned repository, the
`try()` swallows it, and the obsolescence warning goes quiet exactly
when upstream is most likely to have fixed the bug. It now reads the
version from `available.packages(repos = repo)` instead. Verified from
here: that call returns 24861 rows against
`__linux__/noble/latest` and yields `tmbstan 1.2.0`, giving the same
URL that was successfully downloaded and read above.

No source build of tmbstan is proposed, so its cost was not measured
and is not claimed. It is not needed: a clean binary exists at a dated
snapshot and the pin costs one `install.packages()` of three binaries.
For scale, the whole green job of 2026-09-08 took 4 m 07 s including
installing every dependency as a binary.

Why not the second option in the brief, accepting that the tier cannot
run. Because it can. The measurement above shows a clean binary is one
URL away, and a tier that runs beats a tier that announces it did not.
The visibility the second option asked for is kept anyway, by the
canary: if the pin is ever removed, or a snapshot stops serving a
clean build, block 1 of `test-tmbstan-build-guard.R` fails and the
other 114 skip, so the job is red for one legible reason instead of
green or drowned in 108 errors.

That sentence was too strong as first written. It held for the TESTS
and not for the EXAMPLES, which the review caught: 19 sampler calls in
15 `.Rd` files ran behind an installed-only guard and would have added
their own errors under `--as-cran`. Fixed in punch round 1 below, and
the sentence is true now.

### Other jobs that hit the same wall

- `R-CMD-check.yaml`, ubuntu legs. The root DESCRIPTION suggests
  tmbstan and `tests/testthat/test-lkj.R:376` samples through
  `tmbstan::tmbstan()` directly. It already carries its own detector,
  `tmbstan_samples_the_model()` at `test-lkj.R:157`, so it SKIPS. It
  skips SILENTLY, which is the state this lane is against, but the
  fix is a pin in a three-OS matrix and this lane's scope is the
  sample workflow. Not changed. Recorded here as the next thing.
- `test-coverage.yaml` runs the core suite on ubuntu, so the same one
  block is silently uncovered.
- `pkgcheck.yaml` is the container the defect was found in.
- `check-frmtmb-eam.yaml` does NOT install tmbstan: frmtmb.eam's
  DESCRIPTION lists neither tmbstan nor frmtmb.sample, so
  `test-sampling.R:10` skips on the installed check. No CI wall, but a
  developer machine has both, which is why the guard was added there.
- `brms-likelihood.yaml` deliberately omits tmbstan. No exposure.
- The other five extension workflows do not install tmbstan; no
  package but frmtmb and frmtmb.sample names it.

## What is NOT verified, and what would settle it

The workflow change cannot be verified without a push, and pushing is
not this lane's to do.

VERIFIED from this machine:

- The snapshot serves a CLEAN tmbstan binary. Downloaded, extracted
  and read; 0 unpatched overloads, `Built: 2026-08-28`.
- The pin URL derivation is right for this runner. Checked against the
  repository URL that appears in the CI logs of both runs,
  `https://packagemanager.posit.co/cran/__linux__/noble/latest`, and
  unit-checked on the trailing-slash and non-RSPM cases.
- The step's R parses. Extracted from the YAML and parsed.
- The YAML is well formed and the step is in the right place: step 7
  of 9, after the dependency install and before the check.
- The obsolescence probe's `available.packages()` resolves the live
  version from the same repository URL.

NOT verified, and only a run settles it:

- That `install.packages()` from the snapshot SUCCEEDS on the runner.
- That the snapshot's StanHeaders 2.32.10 and rstan 2.32.7 satisfy
  every other dependency resolved from `latest`. brms and bayesplot
  depend on rstan, and these are the versions CRAN had a week earlier,
  so a conflict is unlikely. Unlikely is not measured.
- That the tier goes green on CI with the pin. This is the whole point
  of a run.
- That the obsolescence probe's `download.file()` receives the BINARY
  on the runner. It works from here only with an explicit Linux user
  agent, and the assumption is that R on the runner sends one by
  itself. If it were handed the SOURCE tarball instead, the extracted
  tree would have no `tmbstan/model.hpp` at that path,
  `file.exists()` would be FALSE, and the probe would stay silent
  rather than warn wrongly. It fails safe, but it fails.

`R CMD check` was run on frmtmb.sample only. frmtmb.eam's one changed
test file was run in both arms, one file per process, but that
package's own check was not run: the change there is four lines inside
an existing `test_that` block, behind two skips that already existed.

The evidence that would settle all of them is one run of
`check-frmtmb-sample.yaml` on this branch, and its step summary. If
the pinned install fails, the fallback is to pin the whole job's
repository to the dated snapshot in `setup-r`, which trades freshness
of every other package for a working sampler.

## Defects found and not fixed

1. **`tmbstan_build_broken()` fails OPEN when `model.hpp` is absent.**
   `nzchar(hpp) && ...` means a tmbstan that ships no `model.hpp` is
   called healthy. Measured: `dev/tmbstan-probe-05-seam2.R` prints
   `no model.hpp at all -> detector = FALSE`. That is defensible, since
   an unreadable installation cannot be diagnosed and refusing every
   such build would be worse, but it IS a fail-open path and it is not
   stated in the `@noRd` block. Not changed: the brief put the guard
   out of scope.
2. **`test-lkj.R:157` duplicates the detector** in the core suite,
   with the same marker string written out a second time. If the
   upstream marker ever changes, one of the two will be updated. Not
   changed: different package, and the core suite is not this lane's.
3. **The silent skip in the core ubuntu jobs**, above.

## Files touched

    extensions/frmtmb.sample/tests/testthat/helper-sampling.R
    extensions/frmtmb.sample/tests/testthat/test-tmbstan-build-guard.R   NEW
    extensions/frmtmb.sample/tests/testthat/test-draws-methods.R
    extensions/frmtmb.sample/tests/testthat/test-loo.R
    extensions/frmtmb.sample/tests/testthat/test-parallel-chains.R
    extensions/frmtmb.sample/tests/testthat/test-reparam.R
    extensions/frmtmb.sample/tests/testthat/test-sample-direct.R
    extensions/frmtmb.sample/tests/testthat/test-sampling-ported.R
    extensions/frmtmb.sample/tests/testthat/test-scale.R
    extensions/frmtmb.eam/tests/testthat/test-sampling.R
    .github/workflows/check-frmtmb-sample.yaml

No `R/` source was changed, so nothing needed roxygenising. No version
number, DESCRIPTION version or NEWS section was touched.

Scripts, all prefixed with the lane name:

    dev/tmbstan-probe-01-state.R      the installed tmbstan and StanHeaders
    dev/tmbstan-probe-02-classify.R   which blocks sample, and their skip route
    dev/tmbstan-probe-03-sites.R      the exact skip lines to rewrite
    dev/tmbstan-probe-04-seam.R       the two seams that do not work
    dev/tmbstan-probe-05-seam2.R      the two that do
    dev/tmbstan-probe-06-garbage.R    does noise satisfy the structural tests
    dev/tmbstan-run.R                 one test file, one process, every count
    dev/tmbstan-run-eam.R             the same for frmtmb.eam
    dev/tmbstan-sweep.sh              both arms across the tier

# Punch round 1

Against `dev/reviews/2026-09-09-tmbstan.md`. The reviewer's own
instruments are `dev/rev-tmbstan-*.R` and were not touched. Numbers
already measured there are cited rather than re-run.

## B1, the pin read the wrong repository. FIXED

`r-lib/actions/setup-r@v2` with `use-public-rspm: true` writes TWO
entries, and the CRAN one is not RSPM. The runner prints the pair in
`setup-r-dependencies` under "Repo status":

    1  RSPM  https://packagemanager.posit.co/cran/__linux__/noble/latest
    2  CRAN  https://cran.rstudio.com

So `sub("/latest/?$", ...)` was applied to `https://cran.rstudio.com`,
was a no-op, and the step's own no-op guard aborted the job before
`check-r-package` ran. The `sub()` was right; its input was wrong.

The step now reads `getOption("repos")[["RSPM"]]`, falls back to
`Sys.getenv("RSPM")`, and stops with an explicit message when neither
exists.

**A second defect, found by testing the fix rather than by reading
it.** `[[` on a name a vector does not carry is an ERROR, not `NULL`.
The first version of the fix therefore died on `subscript out of
bounds` on any runner without an RSPM entry, and the `Sys.getenv()`
fallback was unreachable. The lookup is now guarded with
`"RSPM" %in% names(rp)`.

Five repository configurations, `dev/tmbstan-probe-09-pinstep.R`,
which extracts the step's source from the YAML rather than retyping
it:

| configuration | result |
| --- | --- |
| RSPM and CRAN, as the runner has them | pin `.../noble/2026-09-01` |
| CRAN only, no RSPM anywhere | STOPS, naming `use-public-rspm` |
| CRAN in options, RSPM in the environment | pin `.../2026-09-01` |
| RSPM already pinned to a date | STOPS, and says to check the date |
| RSPM with a trailing slash | pin `.../2026-09-01` |

The absent-RSPM case is a STOP, never a silent skip. The reasoning is
in the message: an unpinned tier would run against whatever the
resolver produced, which since 2026-09-03 is a build that samples a
standard normal, so "carry on without the pin" is the one outcome that
must not be available.

## B2, a verdict printed from a file that was not read. FIXED

The step borrowed the package's `nzchar(path) &&`, so an empty
`system.file()` gave `bad = FALSE` and it wrote
`model.hpp patched: TRUE` into the job summary and passed. Fail-open
is defensible in a skip decision and is not defensible in an
assertion.

`unpatched()` is replaced by `verdict()`, which returns one of four
states. Only `patched` continues:

| model.hpp | verdict | job |
| --- | --- | --- |
| path does not exist | `absent` | STOPS |
| `system.file()` returned `""` | `absent` | STOPS |
| present but zero lines | `empty` | STOPS |
| patched | `patched` | continues |
| carries the marker | `unpatched` | STOPS |

The summary line is now `- model.hpp: **<state>**` rather than a
boolean, so the report cannot say "patched" about a file nobody read.
The two stops are worded apart: one names the defect, the other says
the check could not tell, and that tmbstan's own `configure` aborts
rather than ship without the file, so reaching it means a damaged
installation.

## The obsolescence probe's trap, written down

The probe also took the wrong `repo` and is fixed with B1. The second
problem is the reviewer's, found by downloading the CRAN source
tarball: tmbstan ships `inst/model.hpp` and `src/include/model.hpp`
byte-identical at 516 lines with 0 markers, and `configure`
regenerates both at install time, so the shipped pair says nothing
about what any build produced. A maintainer who "fixed" a silent probe
by pointing it at `inst/model.hpp` would get a warning that is
permanently TRUE.

That is now a comment beside the probe, in those terms, ending with
the instruction that matters: if the path stops matching, check that
`repo` is still RSPM rather than hunting for a `model.hpp` deeper in
the tree.

## N1 and N2, the examples. FIXED

The review is right and the record was wrong: coverage was complete
for the tests and not for the examples. 15 `.Rd` files carried 19
sampler calls in `\donttest`, behind
`requireNamespace("tmbstan") && requireNamespace("rstan")` and nothing
else, which is the exact installed-only guard this lane removed from
the tests. `man/sample-as_draws.Rd` had no tmbstan guard at all (N2).

**What an example should do on a broken build.** It cannot skip, so it
must not run. The condition is the same one the tests ask, so the
examples now ask it: `!frmtmb.sample:::tmbstan_build_broken()` is a
further conjunct of the guard that was already there. The example then
does not enter its body, which is what it already did when tmbstan was
absent. The alternative, letting it run and error, is loud but adds up
to 20 errors on top of the canary and buys nothing: the canary has
already said the one thing there is to say.

The marker string is NOT written a third time. The examples call the
package's own detector, so there is one literal in the package and the
duplication the review criticized in core's `test-lkj.R` is not
repeated here.

Edited in the roxygen sources, not in `man/`:

    R/loo.R            4 guards extended
    R/methods-draws.R  7 guards extended, plus sample-as_draws (N2)
    R/sample.R         3 guards extended

After `roxygenise()`: 15 `.Rd` files carry the build guard, and the
count of `.Rd` files that hold a sampler call and do NOT carry it is
**0**. `NAMESPACE` is unchanged.

**The control, which is the arm a guard usually gets wrong.** A guard
that suppressed the examples outright would look identical on a green
run. It does not: `--run-donttest` under `--as-cran` on this clean
machine executes the same 15 example topics and prints the same 16
`frm_sample()` banners after the change as before it, and reports OK
both times. The only new lines are the 15 executions of the detector
itself, one per guarded topic.

    checking examples with --run-donttest ... [15s] OK   before
    checking examples with --run-donttest ... [21s] OK   after

The 6 s is machine load, not behavior: the two runs execute the same
examples.

**One thing recommended and NOT done.** The workflow passes
`args: 'c("--no-manual")'`, which drops the `--as-cran` that
`r-lib/actions/check-r-package@v2` includes by default, and that is
the only reason CI never met this. With the examples guarded it would
now be safe to restore `--as-cran` there, and it would put the
examples under CI for the first time. It is not done in this round:
`--as-cran` also turns on the CRAN-incoming checks, which on a package
that is not on CRAN produce findings unrelated to this lane, and I
cannot see what they would be without a push. It belongs in the same
follow-up as the `FRMTMB_SAMPLER_GATES` decision, where a red run
would mean one thing.

## The fail-open, documented. N5 DONE

`tmbstan_build_broken()`'s `@noRd` now states that `nzchar(hpp) &&`
makes an unreadable installation count as healthy, lists all three
measured shapes that return FALSE (no file, an EMPTY file, and the
defect under a renamed placeholder), and says why fail-open is right
there and only there: tmbstan's `configure` aborts on a failed
autogen, so each reachable absence breaks tmbstan itself and the user
meets a loud error rather than plausible draws. It ends with the rule
the workflow now follows, that a caller turning the answer into a
printed claim must not inherit the fail-open. The detector's behavior
is unchanged.

The EMPTY-file shape is the reviewer's measurement, not mine; I had
measured only the absent one.

## The grad_log_prob assertion. IT WORKS, and it is now a test

The reviewer's mechanism: the patched overload is selected by
`stan::require_not_st_var` and the unpatched one by
`stan::require_st_var`, so the entry points HMC uses read a standard
normal on an affected build. `rstan::grad_log_prob()` is one of them.

I could not install a broken tmbstan, so the broken arm is the
defect's CLOSED FORM rather than a chain: its density is the standard
normal kernel `-sum(u^2)/2`, so its gradient is `-u`. That is not
invented. `dev/prior-dropping-investigation.md` recorded exactly those
numbers from the affected container, `stan_lp` coming back as
`-mu^2/2` and `stan_gr` as `-mu`.

`dev/tmbstan-probe-08-grad-wide.R`, SEEDS 4021, 4023, 4024 and 4025
for the data, sampler seed 11, five unconstrained points per model
(origin, mode+0.3, mode-0.7, a random draw at seed 4030, and three
times the mode). The comparison is `rstan::grad_log_prob(sf, u)`
against `-fit$obj$gr(u)`, the sign because `obj$fn` is the negative
log posterior. Gaps are relative to `max(abs(g_obj))`, which the run
measures.

| model | pars | clean gap | defect gap | max ulp | bitwise | as_tmbstan s |
| --- | --- | --- | --- | --- | --- | --- |
| gaussian | 3 | 0 | 0.9888 | 0 | TRUE | 7.62 |
| poisson | 2 | 0 | 0.9856 | 0 | TRUE | 1.17 |
| bernoulli | 2 | 0 | 0.9631 | 0 | TRUE | 0.33 |
| gaussian, sigma ~ x | 4 | 0 | 0.9945 | 0 | TRUE | 0.31 |

Twenty comparisons. The correct gradient agrees BITWISE every time,
0 ulp, and the defect's would sit at 0.963 to 0.994 of the gradient's
own size. False-alarm rate 0 in 20. The assertion itself costs 0.02 s
across all four models; the cost is the `as_tmbstan()` call.

Written up as one block in `test-tmbstan-build-guard.R`, on the
gaussian model at four points. It asserts a ratio against a ratio,
both measured in the run, with a separation factor of a million, which
is far inside a measured separation of 0 against 0.96. No absolute
tolerance appears.

Measured, `dev/tmbstan-run-gradblock.R` run twice with `desc =`
so the block runs
alone (the fixture blocks in the same file reset the detector's memo
on exit and would otherwise un-poison it):

| arm | passed | failed | error | skipped |
| --- | --- | --- | --- | --- |
| clean | 10 | 0 | 0 | 0 |
| broken | 0 | 0 | 0 | 1 |

**What it is worth, stated honestly.** It runs only on builds the
static check passes, because `as_tmbstan()` refuses the others before
a chain starts. So it does not double-cover the marker case. Its value
is the case the marker cannot see: the renamed placeholder, row 3 of
the reviewer's absent-construction table, where the detector says
FALSE and the sampler is still wrong. It is also chain-free, seed-free
and platform-variance-free, so unlike the chain-agreement assertions
it takes no position on sampler luck and has no reason to be gated by
`FRMTMB_SAMPLER_GATES`.

One bookkeeping note, so the reviewer does not have to find it. The
final `R CMD check` ran on a tarball built just before a comment-only
edit to this block's explanation. That the edit is comment-only is
measured, not asserted: `parse(keep.source = FALSE)` of the tarball's
copy and of the worktree's give 7 expressions each and
`identical()` TRUE, deparse included. No code in the checked tarball
differs from the tree.

**What I could not do.** Run it against a genuinely miscompiled
tmbstan. The broken arm is a closed form with a recorded provenance,
not an observation. What would settle it is one run on an affected
machine, where it should fail at every one of the four points.

## FRMTMB_SAMPLER_GATES, not flipped

Left alone, per the review's own first recommendation: with B1
unfixed the job never ran, so flipping the switch would have shipped
an untested change, and ubuntu chain luck cannot be measured from
here. The reviewer's measurement of what the switch does belongs in
this file whatever is decided:

- 25 `sampler_gates_on()` sites in the sample suite;
- 4 are `skip_if_not()` on a whole block and leave a skip line;
- **21 are inline `if (sampler_gates_on())` inside blocks that
  otherwise pass, and leave no trace at all**;
- across the six files that hold all 25, one file per process, clean
  arm: gates ON 809 passing expectations and 1 block skipped; gates
  OFF 716 and 7.

So the switch removes 93 passing expectations and shows 6 of them.
The other 87 vanish silently. That is the arithmetic behind
`PASS 890, SKIP 10, Status OK` on a standard normal, and it is a
better argument for the grad_log_prob block above than for flipping
the switch: that block cannot be gated away, because it has no reason
to be gated.

The order to do it in, from the review and agreed: land the pin, see
one green run, then flip the switch in a change of its own where a red
run means one thing.

## Also corrected from the review

- The green history is FOUR consecutive green runs on the broken
  binary, not one: 34046054790, 34178185201, 34220763660 and
  34258865676, all with StanHeaders 2.39.1 and PASS 889 or 890. The
  RSPM build stamp of 2026-09-03 08:40 UTC has not moved since, so the
  whole 09-03 to 09-09 window served it.
- N3, the fixture's patched arm, is relabeled. Its comment claimed the
  two-patched-overload shape is stanc 2.32 output. It is not: 2.32
  emits ONE overload. It is what a FIXED autogen would emit under
  2.39, which is the harder case and the one worth pinning. The
  comment now says so and points at the block that covers the real
  2.32 shape.
- Core's `test-lkj.R` is left alone, as instructed. It is a core lane's
  to fix and the review says what it should do.

## R CMD check, punch round

One pass, built tarball, `--as-cran --no-manual`, same environment as
the round-1 pass (`NOT_CRAN=true`, `FRMTMB_SAMPLER_GATES=false`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH):

    Status: OK
    * checking examples with --run-donttest ... [21s] OK
    [ FAIL 0 | WARN 0 | SKIP 11 | PASS 1036 ]

No ERROR, no WARNING, no NOTE. PASS went 1026 to 1036, which is
exactly the ten expectations the grad_log_prob block adds. The 11
skips are the same pre-existing platform and environment gates as
before: chain-agreement gates off (4), brms-fit tier (2), scale tier
(1), RTMBode absent (1), an unhealthy-chain skip (1),
installed-package sources (1) and one empty test. None is the build
skip, because this machine's tmbstan is clean.

## Files touched in this round

    .github/workflows/check-frmtmb-sample.yaml
    extensions/frmtmb.sample/R/sample.R          (@noRd only)
    extensions/frmtmb.sample/R/loo.R             (@examples only)
    extensions/frmtmb.sample/R/methods-draws.R   (@examples only)
    extensions/frmtmb.sample/tests/testthat/test-tmbstan-build-guard.R
    extensions/frmtmb.sample/man/*.Rd            (15, roxygenised)

New scripts:

    dev/tmbstan-probe-07-grad.R       the gradient identity, one model
    dev/tmbstan-probe-08-grad-wide.R  four shapes, false alarms, cost
    dev/tmbstan-probe-09-pinstep.R    the pin step's own arithmetic
    dev/tmbstan-run-gradblock.R       that one block alone, both arms

No behavior in `R/` changed: the detector, the refusal and every
exported function are byte-identical to the base commit apart from the
`@noRd` prose and the `@examples` guards. `NAMESPACE` is unchanged.
