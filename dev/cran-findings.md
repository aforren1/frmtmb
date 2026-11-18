# Lane cran: what a CRAN submission of core frmtmb would take

Round of 2026-11-16, lane cran. Base b255350e (frmtmb 0.83.0 plus CI
files). This lane MEASURES and PROPOSES. It changed no package code,
no test and no vignette. Every file it wrote is under `dev/cran-*`.

## Summary

- As it is, core fails CRAN's time budget by a factor of 4. The whole
  `R CMD check --as-cran` with `NOT_CRAN` unset took 2,502 s (41.7 min)
  on this machine. The tests took 1,824 s of it (20,838 passes; 546
  blocks skipped by `skip_on_cran()`). The vignettes took 323 s.
- The check status is clean apart from one local NOTE (no V8 for the
  HTML manual's math). Installed size is 9.1 MB, reported as INFO.
- The proposed tier keeps every `test_that()` block that runs in under
  0.25 s, plus 6 blocks that are the cheapest to call 6 exports nothing
  else in the tier names. Measured in one session as `R CMD check` runs
  it: 240 s, 10,171 passes, 0 failures, 0 errors, 0 warnings.
- With that tier, three vignettes precomputed and the 8 `\donttest{}`
  examples made plain examples, the projected check is about 570 s
  (9.5 min) on this machine. The margin under 10 minutes is thin; the
  section "Buying margin" lists four more cuts with their measured
  savings.
- The largest obstacle is not time. It is the release cadence: 18 minor
  releases between 2026-09-29 and 2026-11-14, with 133 BREAKING bullets
  in 0.75.0 to 0.83.0. CRAN asks for no more than one update every one
  to two months.

## Method and environment

Scripts (all in `dev/`, all read-only towards the package):

| Script | What it does |
|---|---|
| `cran-check.sh` | `R CMD build`, then `R CMD check --as-cran --timings`, every line stamped with the wall clock |
| `cran-runall.sh` | runs `cran-vignette.R` per vignette and `cran-testfile.R` per test file, one R process each, 1 at a time while the check ran and 2 after |
| `cran-testfile.R` | one test file, with CPU and elapsed time per file and per `test_that()` block |
| `cran-vignette.R` | one vignette knit in a scratch copy, with every chunk timed |
| `cran-tier.R` | builds the tier from the block timings and maps exports to blocks and examples by name |
| `cran-tier-run.R` | runs the tier (or the full suite) in ONE session, moved blocks skipped, with an optional trace of every export |
| `cran-tier-cov.R` | joins the trace counts with the names and the Rd topics |
| `cran-suggests-scan.sh` | Suggests packages named in a file with no guard in that file |
| `cran-static.R`, `cran-urls.R` | title case, NEWS parse, `cat()` outside print methods, urlchecker |

Environment of every run: `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS` and
`FRMTMB_STAN_CACHE` unset; `OMP_THREAD_LIMIT=2`;
`_R_CHECK_LIMIT_CORES_=TRUE` (what `--as-cran` sets);
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; `_R_CHECK_TIMINGS_=0`, so that R
prints the elapsed time of every step (it prints only steps over 10 s
under `--as-cran` when the variable is unset, and on Windows it prints
elapsed time only, never CPU); the CPU-to-elapsed thresholds of CRAN's
incoming check set to 2.5 for examples and vignettes.

Library: `C:/Users/adf44/source/r/wt-cran-lib` (frmtmb installed from
this worktree; drmTMB 0.7.0, a CRAN package, copied from rellib-r22
because the user library lacks it), then the user library. rellib-r22
was NOT on the path, on purpose: it holds the eight extensions, and
CRAN, like core's CI, has none installed. This departs from the lane
environment's `R_LIBS` order for that reason.

Load: other lanes ran during the whole measurement. Every test and
vignette process timed the same fixed arithmetic control first. Over
the 349 test processes it took 0.78 s minimum, 0.80 s median, 0.86 s at
the 90th percentile and 0.92 s maximum; over the 11 vignette processes
0.78 to 0.81 s. So the load varied each time by at most a factor of 1.18.
Calibration against CI: this machine's vignette step (323 s) and
examples step (57 s) sit beside the ubuntu-24.04 runner's (291 s and
66 s on b255350e, run 44009382186), so a number here is close to what
a GitHub ubuntu runner would take. CRAN's own machines were not
measured; nothing here can say how fast they are.

## 1. The check as CRAN runs it

Log: `dev/cran-log/check.log` (stamped), `00check.log`, `testthat.Rout`,
`frmtmb-Ex.timings`, `installed-du.txt`, `build.log`.

Status: `1 NOTE` (the HTML manual: "Skipping checking math rendering:
package 'V8' unavailable", local to this machine).

`R CMD build` (with vignettes; not part of CRAN's check): 396 s.
Tarball `frmtmb_0.83.0.tar.gz`: 4,285,926 bytes (CRAN's policy:
"Source package tarballs should if possible not exceed 10MB").
Uncompressed: R 5.40 MB, tests/testthat 5.14 MB, inst/doc 1.82 MB,
man 0.93 MB, NEWS.md 0.81 MB, vignettes 0.33 MB, data 0.09 MB.

Steps, elapsed (R's own `[Ns]` and the line stamps):

| Step | s |
|---|---:|
| install | 41 |
| R code for possible problems | 83 |
| examples | 57 |
| examples with `--run-donttest` | 63 |
| tests (`Running 'testthat.R' [30m]`) | 1,824 |
| re-building of vignette outputs | 323 |
| PDF manual | 7 |
| HTML manual | 30 |
| every other step | 74 |
| **total check** | **2,502** |

CPU: R on Windows does not print CPU time per step. The parts that carry
it: the tests' own `proc.time()` in `testthat.Rout` reads user 1,715.23
s, system 100.00 s, elapsed 1,823.06 s (CPU/elapsed 1.00); the examples'
user plus system time is 24.14 s of 24.38 s elapsed (sum over
`frmtmb-Ex.timings`); the 11 vignettes rendered in separate processes
took 316.4 s CPU and 330.9 s elapsed. No step comes near CRAN's
CPU-to-elapsed ratio of 2.5, so core uses one core.

Tests in the check: `[ FAIL 0 | WARN 0 | SKIP 562 | PASS 20838 ]`. The
skips: 546 "On CRAN", 10 "package sources are not available", 2 "not
run from the source tree", 1 "not a monorepo checkout", 1 each for
frmtmb.eam and frmtmb.survival absent, 1 `FRMTMB_FUZZ`.

Examples: 129 timed. None takes 5 s CPU or elapsed. The slowest:
`frm_hazard_reads` 3.19 s, `ps` 2.45 s, `profile.frmtmb_fit` 1.80 s,
`dharma_residuals` 1.35 s, `confint.frmtmb_fit` 1.34 s; all others
under 1.0 s. The examples' own time is 24.4 s, but each examples step
takes about 60 s: the rest is the runner. The second step exists only
because 8 Rd files hold `\donttest{}` (R's check runs it when the first
pass counted a `\donttest`); each of those 8 examples runs in under
1.4 s with the `\donttest` part included (`conditional_effects` 0.56,
`confint.frmtmb_fit` 1.34, `frmtmb-emmeans` 0.17,
`frmtmb-marginaleffects` 0.44, `influence.frmtmb_fit` 0.21, `mixture`
0.87, `plot.frmtmb_influence` 0.25, `residuals.frmtmb_fit` 0.55).

Installed size (`du -sk` of the installed package, KiB): total 9,349;
R 5,128 (of it `frmtmb.rdb` 5,084); doc 1,856 (largest:
case-studies.html 420, compatibility.html 224, reinforcement-learning
.html 116); help 1,192; NEWS.md 796; data 130; html 72; bcm 64; Meta 42;
rl 20. R 4.6.1 prints it as `checking installed package size ... INFO`
("installed size is 9.1Mb; R 5.0Mb, doc 1.8Mb, help 1.2Mb"), not as a
NOTE, and the status line does not count it. The R directory alone is
over the 5 MB threshold (`_R_CHECK_PKG_SIZES_THRESHOLD_`), so no cut
short of removing code brings the size under it; leaving NEWS.md
history out of the tarball saves 0.8 MB.

## 2. Test time per file, NOT_CRAN unset

Logs: `dev/cran-log/tests/<file>.txt` (one per file),
`files-sorted.tsv`, `blocks.tsv` (one row per `test_that()` block).

All 349 files ran: pass 20,913, fail 0, err 0, warn 0; 513 of 3,429
blocks skipped. The sum of the files' own elapsed time is 2,055 s
(2,038 s CPU). A file that skips at once still costs 0.76 s, which is
the helpers being sourced; in one session (`R CMD check`) that cost is
paid once, which is why the check's 1,824 s is less than the sum.

The time is concentrated: the 10 slowest files take 698 s (34 %), the
25 slowest 1,085 s (53 %), the 50 slowest 1,348 s (66 %). 72 files take
under 1 s and 170 under 2 s. 65 files skip every block on CRAN (58 s,
all of it helper sourcing).

The 20 slowest (elapsed s | file | passes | slowest blocks):

| s | File | Passes | Slowest blocks |
|---:|---|---:|---|
| 157.7 | test-gp-start.R | 243 | body start reads every coefficient (32 s); marginal's gate (29 s) |
| 104.5 | test-latres-newrow.R | 39 | a category's Est.Error beside a gp() and a latent (85 s) |
| 94.4 | test-predict-re-uncertainty.R | 45 | propagate_error = FALSE (61 s); joint law of rows (27 s) |
| 58.9 | test-predict-draw-mean.R | 134 | draws that move the mean (27 s); a draw on a period (17 s) |
| 57.3 | test-me-newdata.R | 180 | predict()/simulate() per replicate (10 s) |
| 55.2 | test-gp-noncentered.R | 76 | replicated positions and sdgp (24 s) |
| 46.1 | test-ce-levels.R | 101 | bootstrap band with smooths (9 s) |
| 46.1 | test-nl-family.R | 333 | 87 small blocks, none over 5 s |
| 41.1 | test-marginaleffects-re.R | 59 | new level and unseen gp() position (27 s) |
| 37.0 | test-compat.R | 868 | addition term on ordinal families (11 s) |
| 36.8 | test-quadrature-nested.R | 69 | memory guard (8 s); re-anchoring loop (7 s) |
| 35.5 | test-predfix.R | 92 | cs() in sample and at newdata (27 s) |
| 33.5 | test-latres.R | 305 | fitted() on new rows (5 s) |
| 31.3 | test-ordinal-mixture.R | 88 | simulate() category distribution (5 s) |
| 30.5 | test-grad-verdict.R | 74 | importance fit under grad_tol (22 s) |
| 25.4 | test-pp-check-types.R | 216 | every bayesplot ppc type (22 s) |
| 24.0 | test-cdfs.R | 229 | censored fit per family (4 s) |
| 23.4 | test-aterm-expr.R | 64 | weights()/rate()/se()/trunc() expressions (12 s) |
| 22.9 | test-cens-windows.R | 28 | interval in the upper tail (12 s) |
| 22.6 | test-mi-rows-read.R | 69 | reader that skips every missing row (15 s) |

What they test: the exact gp() (its start, its non-centered form, its
prediction draws), new-row prediction with latent residuals, me() terms
and unseen groups (draws and their means), the draw check, nested
quadrature, and parity sweeps over every family (compat, cdfs,
nl-family). Most slow blocks are a few large fits or a loop of draws;
they are regression tests for defects the releases found.

## 3. The proposed CRAN tier

### Budget

Fixed costs measured above that the tier does not touch: install 41,
code check 83, one examples step 57, manuals 37, the rest 74: 292 s.
With the vignette and example changes of section 4 the vignette step is
about 39 s (see there) and the second examples step goes. 600 s minus
331 s leaves about 270 s for the tests on this machine.

### Rule

Keep every block that runs on CRAN and takes under 0.25 s; move every
slower block behind `NOT_CRAN`; then, for each export that no kept
block, no kept file's top-level code and no example names, add back the
cheapest moved block that names it (greedy, by time per newly named
export). Output: `dev/cran-log/tier-0.25.txt`, `tier-kept.tsv`,
`tier-moved.tsv`, `tier-perfile.tsv`.

Block thresholds measured (block time is the sum of the kept blocks'
own elapsed time):

| Threshold | Kept blocks | Kept time | Expectations kept | Moved |
|---:|---:|---:|---:|---:|
| 0.20 s | 1,861 | 143 s | 8,879 | 12,017 |
| 0.25 s | 2,003 | 174 s | 9,929 | 10,967 |
| 0.30 s | 2,153 | 214 s | 10,855 | 10,041 |
| 0.40 s | 2,328 | 273 s | 12,193 | 8,703 |
| 0.50 s | 2,423 | 315 s | 13,095 | 7,801 |

A file-level rule does worse: keeping files under 4 s keeps 244 files
and 7,614 passes for 234 s; the block rule at 0.30 s keeps 10,855 for
214 s.

The 6 greedy additions (8.5 s):

| File | Block | s | Export it adds |
|---|---|---:|---|
| test-portability.R | FN-11: frm_multiple refuses what it cannot pool | 0.36 | as_draws_array, nchains |
| test-families.R | compois fits and recovers the mean model | 0.27 | compois |
| test-hurdle-nb-zoib.R | the registry says simulate() works where a simulator exists | 0.58 | sim_can |
| test-ordinal-thres-names.R | incl_thres on an ordinal mixture is refused in brms's words | 0.61 | ord_thres_linpred |
| test-covstruct-seam.R | a registered bar structure fits as its density says | 3.10 | log_sd_theta_index |
| test-gp-start.R | the warning at the nugget names it, and only there | 3.56 | gp_nugget_message |

### The tier, measured

`dev/cran-tier-run.R` ran the whole suite in one session, with the 907
moved blocks skipped by file and description (log
`dev/cran-log/tierrun-time.txt`):

    COUNTS  pass=10171 fail=0 err=0 skipblocks=1414 warn=0 fileerrors=0
    SESSION mode=time files=349 helpers_real=0.78 total_real=240.28 total_cpu=239.21

So the tier is 2,009 blocks in 273 files, 240 s, 10,171 passes (49 % of
the check's 20,838). Moved behind `NOT_CRAN`: 907 blocks, 10,931
expectations, 1,557 s of block time, in 219 files; 11 files move whole:
test-quadrature-defects.R (21 s), test-gp-predict-draw.R (15 s),
test-message-uniqueness.R (9 s), test-fd-chain.R (7 s),
test-reml-shapes.R (5 s), test-gp-field.R, test-influence-plot.R (2 s
each), test-multiple-pooling.R, test-functional.R,
test-bcm-psychophysics.R, test-stationary-passes.R (1 s each).

The 240 s is 57 s more than the kept blocks' 183 s: top-level code
outside the blocks, which runs even when every block of the file is
skipped. The largest: test-autocor-sigma.R 5.0 s, test-ordinal-mixture.R
4.8 s, test-fcor.R 2.5 s, test-cdfs.R 2.1 s, test-rl-example.R 1.6 s,
test-xbeta-zibb-hurdle-cum.R 1.4 s, test-compat.R 1.3 s.

The files that lose the most time (kept/moved blocks, kept s/moved s):
test-gp-start.R 3/22, 3.8/153.0; test-latres-newrow.R 3/7, 0.5/103.2;
test-predict-re-uncertainty.R 2/6, 0.3/93.2; test-predict-draw-mean.R
7/12, 0.6/57.4; test-me-newdata.R 10/21, 1.1/55.3; test-gp-noncentered.R
6/17, 1.0/53.2; test-ce-levels.R 2/20, 0.3/44.6; test-nl-family.R 27/60,
2.6/42.3; test-marginaleffects-re.R 2/12, 0.3/39.9;
test-quadrature-nested.R 3/11, 0.2/35.7. The full list, block by block
with times, is `dev/cran-log/tier-moved.tsv`.

### What the tier exercises, export by export

NAMESPACE exports 334 names. Two measures:

1. By name (`cran-tier.R`, `all.names()` of each block's parsed code):
   206 exports are named by some core test block; 205 by a block that
   runs on CRAN; 195 by the 0.25 s tier before the additions; 135 by an
   example. After the additions, every export that a CRAN-run block
   names is named by a tier block or by an example.
2. By execution (`cran-tier-run.R trace`, log `tierrun-trace.txt`,
   counts `tierrun-trace-calls.tsv`, joined in `tier-cov.txt`): a
   tracer on each of the 288 exports that have a plain function binding
   (the active-binding generics such as `fixef`, `loo` and `VarCorr` and
   the re-exports cannot be traced). The tier reached 261 of the 288.

Split by kind (extension API = the aliases of `frmtmb-sampling-api`,
`frmtmb-extension-api`, `frmtmb-shared-generics` and the
`\keyword{internal}` topics):

| Kind | Reached | Not reached | Untraced |
|---|---:|---:|---:|
| user-facing | 132 | 5 | 46 |
| extension API | 129 | 22 | 0 |

The same traced run of the FULL CRAN-mode suite, nothing moved (log
`fullrun-trace.txt`): 1,859 s, pass 20,913, fail 0, err 0, warn 0, and
267 of the 288 traced exports reached. So the tier loses 6 exports by
execution (calls in the full suite in brackets; `trace-lost.tsv`):

- user-facing: `frmtmb_autocor_pars` (2) and `update_adterms` (15).
  Both have examples, which run on CRAN.
- extension API: `ce_plan_has_new` (70), `predict_category_props` (14),
  `resp_refuse_single` (5), `hyp_tail_p` (2).

Adding back, for each, the cheapest moved block that reaches it would
close the gap; the name-based greedy step missed these because the
blocks reach them through other functions, not by name. The 21 exports
the full CRAN-mode suite does not reach either are not a tier loss.
Tracing cost 3.6 % (249 s traced, 240 s untraced).

User-facing exports the tier does not reach by execution: `fcor` and
`unstr` (formula terms; the parser reads them as symbols and never calls
them; 30 and 12 tier blocks name them), `frmtmb_autocor_pars` (3 tier
blocks name it in test-autocor-seam.R, the example runs),
`update_adterms` (the example runs; no tier block), `prior_logistic`
(no core test names it; no example).

Untraced user-facing exports (generics owned by brms, loo, posterior
or nlme), by name in the tier: `autocor` 24 blocks, `conditional_effects`
57, `fixef` 142, `hypothesis` 34, `ranef` 38, `VarCorr` 25, `variables`
43, `pp_check` 16, `ngrps` 7, `prior_summary` 6, `refit` 5,
`conditional_smooths` 3, `loo` 3, `waic` 2, `add_criterion`,
`as_draws_array`, `log_lik`, `nchains`, `parnames`, `stancode`,
`standata` 1 each. Named by no tier block but run by an example:
`as_draws`, `as_draws_df`, `bayes_R2`, `expose_functions`, `kfold`,
`loo_predict`, `loo_subsample`, `ndraws`, `posterior_summary`,
`pp_mixture`. Named by no tier block and no example: `as_draws_list`,
`as_draws_matrix`, `as_draws_rvars`, `LOO`, `WAIC`, `loo_compare`,
`loo_epred`, `loo_linpred`, `loo_model_weights`, `loo_moment_match`,
`loo_predictive_interval`, `loo_R2`, `niterations`, `nvariables`,
`reloo`. These are the methods for draws objects that frmtmb.sample
makes: no core test names `as_draws_list`, `as_draws_matrix`,
`as_draws_rvars`, `LOO`, `WAIC`, `loo_compare`, `loo_epred`,
`loo_linpred`, `loo_predictive_interval`, `niterations` or `nvariables`
at all (frmtmb.sample's tests do), so the tier loses nothing there.

### How to ship it (proposal, user decision)

The measurement skipped blocks through a `test_that()` placed ahead of
testthat's. The package should not ship that. Two options:

1. Generate the edit: a dev script inserts `skip_on_cran()` as the
   first line of each of the 907 blocks in `tier-moved.tsv`, and a guard
   test asserts that the tier's time stays under a budget measured
   against a control. Plain, visible, and what `dev/lane-rules.md`
   already asks of a test over 10 s; but the 0.25 s line is about 40
   times stricter, so every new test needs the same judgment.
2. One list, `tests/testthat/cran-tier.txt`, read by a helper
   `skip_unless_cran_tier()` called by the slow blocks. Same effect, one
   place to regenerate.

Both rest on per-block times that move with the machine. Regenerate
the list from `cran-runall.sh` before each CRAN submission, not at
every release.

### Buying margin

Measured savings, on this machine:

- Move top-level fits into their blocks (or skip at file level when no
  kept block needs them): up to 57 s, of which the seven files above
  hold 18.7 s.
- Threshold 0.20 s instead of 0.25 s: 31 s less block time, 1,050
  fewer expectations.
- Precompute diagnostics.Rmd too: 11.5 s.
- Precompute every vignette: 34 s (the 8 live vignettes' render time).

## 4. Vignettes

Per-vignette render time, NOT_CRAN unset, one process each (logs
`dev/cran-log/vig/*.txt`; render CPU / elapsed s):

| Vignette | Chunks | CPU | Elapsed | Heaviest chunks (elapsed s) |
|---|---:|---:|---:|---|
| habit.Rmd | 17 | 197.8 | 201.7 | boot 100.8, persubject 73.0, hier 24.1 |
| case-studies.Rmd | 90 | 66.0 | 67.6 | smocc-fit 24.6, animal-boot 15.9, animal-mv 6.4, animal-mv-sens 5.6, fosr-fs-check 4.4 |
| reinforcement-learning.Rmd | 17 | 26.4 | 27.5 | importance 24.0 |
| diagnostics.Rmd | 14 | 10.7 | 11.5 | imp-fit 4.7 |
| spectral.Rmd | 18 | 4.2 | 5.1 | |
| compatibility.Rmd | 32 | 2.9 | 4.5 | |
| curve-inference.Rmd | 17 | 3.2 | 4.0 | |
| frmtmb.Rmd | 15 | 2.6 | 3.4 | |
| inputs.Rmd | 7 | 1.8 | 2.9 | |
| bayesian-cognitive-modeling.Rmd | 3 | 0.7 | 1.6 | |
| brms-migration.Rmd | 5 | 0.2 | 1.0 | |
| **sum** | | **316.4** | **330.9** | check step: 323 s |

Proposal:

- Precompute habit.Rmd whole (knit `habit.Rmd.orig` locally to
  `habit.Rmd`, the rOpenSci pattern): saves about 200 s. Its three heavy
  chunks feed the figures after them, so `eval = FALSE` on them alone
  would break the later chunks.
- Precompute case-studies.Rmd: saves about 66 s. Partial alternative:
  `eval = FALSE` on `smocc-fit` and `animal-boot` with their printed
  output in the text saves 40 s, but `smocc-fit` feeds six chunks after
  it.
- reinforcement-learning.Rmd: `eval = FALSE` on `importance` (the chunk
  is already `error = TRUE` and nothing reads `fit_i` after it), with its
  output pasted as text: saves 24 s.
- Keep the other 8 live: 34 s of render. With 3 precomputed vignettes
  at about 1.5 s each, the step is about 39 s instead of 323 s.

What is lost: the precomputed code no longer runs on CRAN. The tier
keeps blocks for every function those chunks call (frm_bootstrap,
importance, nl formulas); the vignettes stay live in CI and pkgdown.

## 5. Other items a submission would hit

### Examples

- `\dontrun{}` in 4 Rd files. CRAN's reviewers ask that `\dontrun{}` be
  used only when the example really cannot run. `frm.Rd` wraps an lme4
  sleepstudy fit: use `if (requireNamespace("lme4", quietly = TRUE))`
  (or `@examplesIf`) instead. `frmtmb_register_aterm.Rd`,
  `frmtmb_register_covstruct.Rd` and `frmtmb_register_frame_check.Rd`
  wrap a registration that lasts for the session: run it and undo it in
  the example, or keep `\dontrun{}` and say why in the submission
  comments.
- `\donttest{}` in 8 Rd files, all under 1.4 s: make them plain
  examples (guarded by `requireNamespace()` where they need a Suggests
  package). That removes the second examples step: 63 s.
- No example writes a file, calls `setwd()`, `options()` or
  `Sys.setenv()`; `plot.frmtmb_fit.Rd` resets `par()` with `par(op)`.
- One Rd file with `\usage` has no `\value`: `frmtmb-gev-xi-chart.Rd`,
  which documents 4 EXPORTED functions (`dpar_has_predictor`,
  `gev_xi_chart_terms`, `gev_xi_chart_frame`, `gev_xi_chart_start`)
  under `\keyword{internal}`. Reviewers ask for `\value` on exported
  functions.

### Cores

Core starts no parallel work: no `detectCores()`, `mclapply()`,
`makeCluster()` or `future` in `R/`; `predict(cores = )` refuses;
three test files set `mc.cores = 1`. The check ran with
`OMP_THREAD_LIMIT=2` and `_R_CHECK_LIMIT_CORES_=TRUE` and the tests'
CPU/elapsed ratio was 1.00 (policy: "it must never use more than two
simultaneously").

### Suggests used unconditionally

Policy: "A package listed in 'Suggests' or 'Enhances' should be used
conditionally in examples or tests", and Writing R Extensions
recommends always. `dev/cran-suggests-scan.sh` (per file, so it can
only under-report; log `suggests-scan.txt`) found 69 file/package
pairs. After reading each: 33 are `withr` (installed whenever testthat
is, since testthat imports it), most others are comments, and these
run on CRAN with no guard:

- lme4 data or calls: test-aliased-grouping.R (`data(sleepstudy, package
  = "lme4")` in a block), test-confint-anova.R (at TOP LEVEL, so the
  whole file errors without lme4), test-numerical-robustness.R
  (`lme4::glmer()`), test-open-issues.R, test-verbose.R.
- emmeans: test-latres-brms.R, test-resp-read.R,
  test-extension-home-mixed.R (to be confirmed per block).
- marginaleffects: test-resp-read.R.
- MASS: test-predict-draw-mean.R (`MASS::mvrnorm()`).

Each needs `skip_if_not_installed()`. The CI skip gate already fails a
skip on a missing Suggests package, so the guards cost nothing on CI.
Vignettes: every Suggests use is behind `requireNamespace()` or
`eval = FALSE`; the hits in inputs.Rmd are prose. Examples: the only
hit is the lme4 example inside `\dontrun{}` above.

### Package code

- `.onLoad()` sets `options(marginaleffects_model_classes = )` on every
  load, whether or not marginaleffects is installed. It adds to a
  user-visible option at load. marginaleffects documents this as the way
  to whitelist a class, but a reviewer may ask; setting it only when
  marginaleffects is installed narrows it.
- Two functions set a fixed seed (`mixture_starts_draws()` with
  20261111L, `draw_check_z()` with `draw_check_seed`). Both restore the
  caller's `.Random.seed` on exit (the first also the RNG kinds).
  Reviewers ask "do not set a seed to a specific number within a
  function"; the answer is the restore, and it should be in the
  submission comments.
- `cat()` outside print, format and summary methods in 6 functions:
  `cat_family_links`, `cat_formula_lines`, `diagnose`, `print_pform`,
  `print_summary_block` (`print()`), `prior_summary.frmtmb_fit`.
  `diagnose()` is user-facing; reviewers ask for `message()` for
  information that cannot be suppressed. It takes `quiet = `, which is
  the answer to give.
- Two vignettes call `frmtmb:::build_objective()`: case-studies.Rmd runs
  it, brms-migration.Rmd shows it (eval = FALSE). Reviewers object to
  `:::` in user documentation; export the function or drop the call.
- No `:::` to other packages, no `installed.packages()`, no writes
  outside `tempdir()`, no `Sys.setenv()` in `R/`.
- Active bindings (`frm_install_generics()`, `R/ad-env.R`) live in
  frmtmb's own namespace and environments, made at load; they do not
  touch another namespace.

### URLs

urlchecker 2.0.0 in the user library: `urls with a problem: 0` (log
`dev/cran-log/urls.txt`, pandoc on PATH).

### DESCRIPTION

- Title "Formula-Based Regression Models via 'RTMB'":
  `tools::toTitleCase()` returns it unchanged.
- Description quotes 'brms', 'RTMB' and 'lme4'. Missing: references for
  the methods. Reviewers routinely ask for them in the form
  `Kristensen et al. (2016) <doi:10.18637/jss.v070.i05>` (TMB) and
  Bürkner (2017) <doi:10.18637/jss.v080.i01> (brms).
- Authors@R: one person, roles aut and cre. Add "cph" (R Packages, 2e:
  "Check that Authors@R: includes a copyright holder (role 'cph')").
  inst/COPYRIGHTS names other holders whose material ships in the
  tarball: Hardwick and Haith (MIT, data/habit_prep.rda) and the
  stan-dev example-models authors, translation by Martin Smira (BSD-3,
  Stan programs in tests/testthat/test-bcm-*.R). Policy: copyright
  "should preferably be indicated via 'cph' roles in the 'Authors@R'
  field". Either add them as cph, or add `Copyright: See file
  inst/COPYRIGHTS` (Writing R Extensions 1.1.1) and say so.
- inst/COPYRIGHTS says "see the LICENSE file", and the package has no
  LICENSE file (none is needed for GPL (>= 2)). Fix the sentence.
- `Roxygen: list(... "srr::srr_stats_roclet")` names a package that is
  not in Suggests. R CMD check ignores the field.

### NEWS

NEWS.md is 813,369 bytes, 14,182 lines, 82 versions; R's NEWS.md parser
reads it with 0 bad entries (291 entries). It installs (796 KiB). CRAN
sets no limit, but a first submission usually starts the file at the
CRAN version. Proposal: keep the history on the website and in git;
start NEWS.md at the first CRAN version.

### Release cadence

CRAN's policy: once a package is established, updates "no more than
every 1–2 months". The project released 18 minor versions from 0.66.0
(2026-09-29) to 0.83.0 (2026-11-14), one every 2.6 days; NEWS.md holds
190 BREAKING bullets, 133 of them in 0.75.0 to 0.83.0 (2026-10-25 to
2026-11-14). Policy for API changes: "If an update will change the
package's API and hence affect packages depending on it", the reverse
dependencies must be checked and their maintainers told. Once the
extensions are on CRAN, each BREAKING core release means up to 8
coordinated resubmissions.

Proposal: CRAN gets a snapshot at most every two months; development
continues on GitHub and an r-universe (the ODE extension already points
users at kaskr.r-universe.dev for RTMBode). Keep the extensions off
CRAN until the extension API (137 aliases in `frmtmb-sampling-api`
alone) stops moving.

## 6. The extensions

CRAN availability of every dependency was read from
`https://cloud.r-project.org/src/contrib/PACKAGES` on 2026-11-17. Every
Imports and Suggests of every extension is on CRAN except RTMBode
(frmtmb.ode's Suggests; DESCRIPTION already has
`Additional_repositories: https://kaskr.r-universe.dev`) and the
sibling extensions themselves. Times are the `check-r-package` step of
the CI run on c6383527 (2026-11-15, ubuntu, NOT_CRAN=true, so an upper
bound for CRAN), with the tests step from the log stamps; logs
`dev/cran-log/ci/<ext>.log`, `ci-jobs.tsv`, `ci-steps.tsv`.

| Extension | Depends/Imports on extensions | Non-CRAN Suggests | Test files (skip_on_cran) | expect_ calls | CI check step | CI tests | CI vignettes |
|---|---|---|---:|---:|---:|---:|---:|
| frmtmb.survival | none | none | 16 (5) | 378 | 75 s | 26 s | 4 s |
| frmtmb.latent | none | none | 13 (7) | 447 | 104 s | 57 s | 5 s |
| frmtmb.spatial | none | none | 8 (3) | 338 | 93 s | 43 s | 3 s |
| frmtmb.coupling | none | none | 11 (5) | 260 | 88 s | 40 s | 6 s |
| frmtmb.ode | none | RTMBode | 11 (5) | 359 | 154 s | 81 s | 23 s |
| frmtmb.learn | Imports frmtmb.eam | none | 16 (6) | 335 | 99 s | 51 s | 5 s |
| frmtmb.eam | none | frmtmb.sample | 33 (11) | 1,149 | 645 s | 521 s | 38 s |
| frmtmb.sample | none | frmtmb.eam, .latent, .ode, .survival | 92 (63) | 3,097 | 566 s | 492 s | 2 s |

Order: core, then survival, latent, spatial, coupling, ode, eam, then
learn (it imports eam), then sample. The Suggests cycle between eam and
sample is allowed, but each must use the other conditionally.

- Six extensions (survival, latent, spatial, coupling, ode, learn) check
  in under 3 minutes on CI even with NOT_CRAN=true. They need no tier.
- frmtmb.ode: `frm_ode()`, its main function, needs RTMBode, which is
  not on CRAN. On CRAN its tests that need RTMBode skip
  (`skip_if_not_installed("RTMBode")`) and its example is behind
  `requireNamespace("RTMBode")`, so CRAN would check little of it. CRAN
  accepts a Suggests from `Additional_repositories`, but the reviewer
  may ask why the package is on CRAN when its main function cannot run
  from CRAN alone.
- frmtmb.eam and frmtmb.sample take 9 to 11 minutes on CI with
  NOT_CRAN=true. Their CRAN-mode time was not measured (the brief says
  not to run them all); each would need its own tier, measured the way
  section 3 measures core's.
- Each extension's `Depends: frmtmb (>= 0.8x.0)` tracks core's minor
  version, so the cadence of section 5 applies to all nine packages at
  once.

## What this lane did not do

- Changed no package code, test, vignette or DESCRIPTION; the tier is a
  measurement. The mechanism to ship is a user decision (section 3).
- Did not run the check on a CRAN machine or win-builder (posting
  nothing externally), so CRAN's own speed is unknown; this machine is
  about as fast as GitHub's ubuntu runner on the steps both measured.
- Did not run the extension suites.
- Did not run a no-Suggests check (`_R_CHECK_DEPENDS_ONLY_=true`); the
  Suggests list above is from a static scan.
- `R CMD check --as-cran` ran once, on the unchanged package; the lane
  changed only `dev/`.
