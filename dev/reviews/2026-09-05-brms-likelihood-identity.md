# Review: brms likelihood-identity lane (branch `wt-brmslp`)

Reviewer notes. Base `a233c3c`; main at `3a856e6`. Everything below is
reproduced locally unless explicitly marked otherwise.

Private review library: `<scratchpad>/rb-lib`, with the worktree core
installed there (frmtmb 0.49.1, `R CMD INSTALL --library=$RB_LIB`).
Three separate Stan caches under the scratchpad, never the lane's.

Sections are numbered by the VERIFY list and appear in the order I
finished them, not in list order: 1, mo(), 6, 5, 7, 8, 9, 3, 4, 2,
then the verdict.

---

## 1. Toolchain diagnosis: CONFIRMED, and the mechanism is exact

The lane's account is right in every particular, and the evidence is
stronger than the lane's own hand-compile.

Installed pair, verified: rstan 2.32.7, StanHeaders 2.39.1, R 4.6.1,
brms 2.23.0, lme4 2.0-6, RTMB 1.9.

Clean state first: `R_MAKEVARS_USER` is unset and `~/.R/Makevars` does
not exist, so the failing run is not the product of leftover config.

**Without the workaround.** `rstan::stan_model()` on
`parameters { real mu; } model { mu ~ normal(0, 1); }` fails in 26.8 s
with `Error: invalid connection`, exactly the reported symptom. With
`verbose = TRUE` rstan prints the make command it would use, and it
contains, in this order on one compiler invocation:

    g++ -std=gnu++17 -I"...include" ... -std=c++1y ... -O2 -Wall ... -c file.cpp

rstan sets `CXX='$(CXX17) $(CXX17STD)'`, which is where the leading
`-std=gnu++17` comes from. The `-std=c++1y` after it is the Rcpp
`cpp14` plugin's contribution through `PKG_CXXFLAGS`, and
`ALL_CXXFLAGS = $(R_XTRA_CXXFLAGS) $(PKG_CXXFLAGS) ... $(CXXFLAGS)`
(Makeconf line 284) is what puts it there. Last `-std` wins, so the
translation unit compiles as C++14 against StanHeaders 2.39.1, which
needs C++17.

**With the workaround.** `R_MAKEVARS_USER` pointing at a file holding

    CXX17FLAGS = -O2 -Wall -mfpmath=sse -msse2 -mstackrealign -std=gnu++17

compiles the identical model successfully in 55.6 s. This is the
controlled experiment: one flag changed, nothing else, fail to pass.
The lever is correct for a non-obvious reason worth recording, because
a maintainer might otherwise reach for `CXX17STD`: rstan passes
`CXXFLAGS='$(CXX17FLAGS)'`, and in the emitted command line
`$(CXXFLAGS)` lands AFTER `PKG_CXXFLAGS`. `CXX17STD` would land before
it and lose. `CXX17FLAGS` is the only one of the two that wins.

One correction to the lane's prose, cosmetic but worth fixing so a
future reader can follow it: the plan says Makeconf places the plugin
flag "after its own `-std=gnu++17`". R 4.6.1's Makeconf default is
`CXX = ... g++ -std=gnu++20` (line 118); the `gnu++17` in the actual
command line comes from rstan's own `CXX17STD` override, not from a
Makeconf default. The conclusion is unchanged.

Measured cold compile of the trivial model after the repair: 55.6 s
(the lane recorded 97 s; same order, machine variance).

**This is an environment fault, not a frmtmb fault.** Nothing under
`R/` or `tests/` works around it, which is the right call. See item 6
for whether CI is exposed.

---

## The mo() divergence: CONFIRMED, quantified, and it has prior art
   the lane missed

(Not a numbered VERIFY item; it is the lane's most consequential claim.)

### The structural claim, reproduced

`brms::make_standata` / `make_stancode` on the same 300-row frame:

| formula | Imo | Jmo | simplexes brms declares | frmtmb zetas | frmtmb npar |
|---|---|---|---|---|---|
| `y ~ mo(inc) + z` | 1 | 3 | `simo_1` | `zeta1` | 6 |
| `y ~ mo(inc):z`   | 1 | 3 | `simo_1` | `zeta1` | 5 |
| `y ~ mo(inc) * z` | 2 | 3, 3 | `simo_1`, `simo_2` | `zeta1` | 7 |

Exactly as the lane reported, down to `Imo = 2` and `Jmo = c(3, 3)`.
The source claim checks out too: `R/frame.R:1865` is
`mo_zetas <- list()   # simplexes are shared per mo() variable`, and
the key at line 1888 is `vkey <- deparse1(mexpr)`, the variable
expression, so `mo(inc)` and `mo(inc):z` collide on one `zeta1`.

### The size of the gap, measured

The shared-simplex model is nested in the per-term model, so the gap
can be measured entirely inside frmtmb by giving the interaction a
second copy of the variable (`inc2`), which deparses differently and
therefore gets its own simplex. Same data, same optimizer, both fits
converged:

    shared simplex  (frmtmb's mo(inc)*z)  logLik -419.908893  npar 7
    free simplexes  (brms's  mo(inc)*z)   logLik -419.073982  npar 9
c)*z)   logLik -419.073982  npar 9

    log-likelihood gap  0.834911 nats
    free-parameter gap  2
    LRT 2*gap = 1.6698 on 2 df, p = 0.4339

The gap is small on THIS frame only because the data were simulated
with a single monotonic shape, which is the null the restriction
imposes. The restriction is nonetheless real: frmtmb cannot express a
model brms can, and the deficiency is structural, not numeric.

### Is the sharing deliberate? Yes, but the stated reason is false

This is the part the lane did not run down, and it changes the
recommendation.

The sharing is a deliberate v0.18 design choice, recorded in two
places besides the code comment:

- `dev/feature-gaps.md:421` — "mo()/mi() interactions DONE in v0.18
  (two-way `:`/`*` with numeric terms; shared simplex per mo variable)."
- `R/parse.R:878` — "`*` also emits the main effects; mo()
  interactions share their variable's simplex."

But `NEWS.md:2364`, in the v0.18 entry, states the rationale:

> `mo()` interactions share their variable's simplex (**brms
> convention**).

That parenthetical is factually wrong. It is not the brms convention.
brms declares one simplex per special term, as the table above shows
and as brms's monotonic vignette describes. So the design was chosen
in order to match brms, and it does not match brms. That reframes the
finding: deliberate in mechanism, mistaken in premise. A maintainer
reading only `NEWS.md` would conclude the behavior is correct and
brms-compatible, and would be wrong.

### Prior art the lane presents as new

`dev/brms-vignette-audit.md:513-517` already recorded this, under the
`mo()` prior finding ranked 5:

> A second modeling difference sits next to it: `mo(income) * age`
> makes the interaction SHARE the main effect's simplex. brms fits
> two. The shape of the monotonic effect therefore cannot vary with
> age, which is what that vignette section is about.

and lists "the shared `mo()` simplex across an interaction" at line
680 among items new to that audit. The lane's plan document presents
the divergence as its own discovery and cites neither document. The
conclusion is unaffected and the lane's evidence is better (it has the
Stan parameter block, which the audit did not), but the write-up
should cross-reference the audit so the maintainer sees this is the
second independent report of the same defect, not the first.

### Recommendation to the maintainer

Not the lane's call to make, and not mine, but the evidence points one
way:

1. **Fix `NEWS.md:2364` regardless of what else is decided.** The
   "(brms convention)" claim is false and is the reason the divergence
   has survived two audits. This is a one-line documentation defect
   and it is independent of the modeling question.
2. **Prefer keying the simplex on the term, not the variable.** brms's
   behavior is the defensible one: a user writing `mo(inc) * z` after
   reading brms's monotonic vignette is asking whether the monotonic
   shape varies with `z`, and the shared simplex answers a different
   question silently. The change is local: key `mo_zetas` on the term
   rather than on `deparse1(mexpr)` at `R/frame.R:1888`.
3. **If sharing is kept**, it must be opt-in or at least announced:
   documented at `mo()`'s help page and in the monotonic vignette, and
   ideally a message when one mo() variable appears in more than one
   term. Silently fitting a stricter model than the formula implies is
   the failure mode both audits object to.
4. Either way the lane's structural assertion in the test file is the
   right holding pattern: it fails loudly if either package moves.

---

## 6. The workflow

### Shared steps against check-frmtmb-sample.yaml

Compared line by line. The five steps the two have in common are
byte-identical in substance:

| | check-frmtmb-sample | brms-likelihood |
|---|---|---|
| `runs-on` | ubuntu-latest | ubuntu-latest |
| `actions/checkout` | v4 | v4 |
| `setup-pandoc` | v2 | v2 |
| `setup-r` | v2, `r-version: release`, `use-public-rspm: true` | identical |
| `setup-r-dependencies` | v2, `needs: check` | v2, `needs: check` |
| env | GITHUB_PAT, R_KEEP_PKG_SOURCE, NOT_CRAN | same three, plus two |

Differences, all of them justified by what the job does:

- `timeout-minutes: 90` added. The reference workflow sets no timeout.
- `schedule` added. The reference has none.
- `extra-packages` is a five-package list instead of `any::rcmdcheck`.
- The reference does `R CMD INSTALL .`; this one uses
  `pkgload::load_all(".")`. Sound here: frmtmb has no `src/`, so there
  is nothing for `load_all` to fail to compile, and the tests read
  `ord_tau_from_raw()`, an internal, which `load_all` makes reachable.
- No `check-r-package` step, by design: this tier is not an R CMD check.

### Triggers and path filter: all correct

Parsed with an R YAML reader (note for anyone repeating this: a naive
`read_yaml` turns the `on:` key into boolean `TRUE` under YAML 1.1, so
address it positionally; GitHub Actions itself is unaffected).

    TRIGGERS: push, pull_request, workflow_dispatch, schedule
    push branches: main master
    paths: R/** | tests/** | DESCRIPTION | NAMESPACE |
           .github/workflows/brms-likelihood.yaml
    pull_request paths identical to push paths: TRUE
    cron: 17 4 * * 1   (Mondays 04:17 UTC)
    runs-on ubuntu-latest, timeout 90, 7 steps

Four triggers as claimed, path filter as claimed.

### Inline Rscript blocks

Both extracted and passed to `parse()`. Both parse clean. The gate
assertion in the run step is the right shape: it `stopifnot()`s rstan,
brms, lme4 and the two env vars BEFORE the run, and it fails the job on
`any(df$skipped)` afterwards, so a fully skipped suite cannot report
green. That is the strongest form of this guard and it is correct.

### FINDING: RSPM serves the exact pair that fails, and the workflow
### has no mitigation. This is the one blocking item.

The lane recorded "CI is not known to be affected, because RSPM
resolves its own pair on Linux, and this could not be verified from
here." I verified it. It is affected.

`available.packages()` against the repositories `setup-r` with
`use-public-rspm: true` actually uses:

| repository | rstan | StanHeaders | brms |
|---|---|---|---|
| RSPM `__linux__/noble` (ubuntu-24.04 = today's ubuntu-latest) | 2.32.7 | 2.39.1 | 2.23.0 |
| RSPM `__linux__/jammy` (ubuntu-22.04) | 2.32.7 | 2.39.1 | 2.23.0 |
| RSPM source | 2.32.7 | 2.39.1 | 2.23.0 |
| CRAN | 2.32.7 | 2.39.1 | 2.23.0 |

That is the identical mismatched pair that fails on this machine. Not
a similar one: the same two version numbers.

The three ingredients of the fault are all platform-generic, and I
checked each rather than assuming:

1. **The plugin emission is not platform-guarded.** `rstan:::rstanplugin()`
   branches on `.Platform$OS.type == "windows"`, and BOTH arms return
   `includes = "// [[Rcpp::plugins(cpp14)]]\n"`. The branch differs
   only in library paths. Linux gets `cpp14` too.
2. **The C++17 toolchain selection is not platform-guarded.**
   `rstan:::cxxfunctionplus()` sets `USE_CXX17 = 1` unconditionally,
   which is what makes R use `CXX='$(CXX17) $(CXX17STD)'` and
   `CXXFLAGS='$(CXX17FLAGS)'`.
3. **The flag ordering is generic Makeconf**, not a Windows quirk:
   `ALL_CXXFLAGS = $(R_XTRA_CXXFLAGS) $(PKG_CXXFLAGS) ... $(CXXFLAGS)`
   with the rule `$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) -c`. The
   plugin's `-std=c++1y` rides in `PKG_CXXFLAGS` and therefore lands
   after `CXX17STD`'s `-std=gnu++17` on every platform.

Same versions, same unguarded plugin, same ordering. I cannot execute
a Linux runner from here, so I mark this a high-confidence prediction
rather than an executed result; but every input to it is verified, and
the failure mode if I am right is a 90-minute job that goes red on its
first run, on a workflow whose entire purpose is to be trusted.

The mitigation is cheap enough that it should go in regardless of how
confident anyone is, because it is a no-op when unnecessary. Add a step
before the run step:

```yaml
      # rstan 2.32.7 emits [[Rcpp::plugins(cpp14)]], whose -std=c++1y
      # lands after CXX17STD's -std=gnu++17 and wins, while StanHeaders
      # 2.39.1 needs C++17 (std::void_t in eigen_plugins.h). CXX17FLAGS
      # is the only lever that lands last. Harmless if the pair is ever
      # fixed upstream.
      - name: Force C++17 for rstan's generated models
        run: |
          mkdir -p ~/.R
          echo 'CXX17FLAGS = -O2 -Wall -std=gnu++17' >> ~/.R/Makevars
        shell: bash
```

Pinning StanHeaders to 2.32.x is the alternative, but it is worse: it
fights RSPM's resolution and it will drift.

---

## 5. The cache: CONFIRMED, with one loose sentence in the plan

### Key sensitivity

`brms_stan_cache_key()` is md5 over the Stan code with the rstan
version appended, so by construction it tracks the code. Measured:

| variation | Stan code identical | key identical |
|---|---|---|
| same formula, different data (n and values) | TRUE | TRUE |
| different family (gaussian -> poisson) | FALSE | FALSE |
| rstan version string bumped | n/a | FALSE |

So the key is invariant to the data and sensitive to the program and
to the rstan version, which is the required behavior and is what makes
one compile serve every data variant of a row.

**Nit, in the plan not the code.** The plan says "the code is a
function of the formula and the family, never of the data". The first
half is not quite true and the test proves it: `y ~ x` and
`y ~ x + I(x^2)` generate *byte-identical* Stan code, because brms
declares `vector[K] b` and passes `K` as data. Two different formulas
therefore share one cache entry. This is correct and in fact desirable
(more hits), and it is safe precisely because the cache is addressed on
the code rather than on the formula. But the sentence as written would
mislead someone reasoning about cache collisions. Suggested rewording:
"the code is a function of the model's structure and never of the data,
so one compile serves every data variant of a row, and two formulas
that generate the same program correctly share an entry."

### Corruption fall-through

Tested at both corruption shapes, and in a fresh process, which is the
case that matters (CI restores a cache into a new session):

    fresh session, cache file overwritten with non-RDS text
      -> 70.2 s, real recompile, returns a stanmodel
      -> rewrites a loadable RDS in place

    same session, corrupted then truncated to zero bytes
      -> returns a valid stanmodel both times, no error

A corrupt cache never fails the run. Confirmed as designed.

(Aside for anyone re-running this: inside one session the "recompile"
after corruption takes 0.2 s, because rstan reuses the DSO it already
built. Only a fresh process shows the true cost. Do not mistake the
0.2 s for the cache still being hit.)

### Ignore rules

    git check-ignore -v dev/stan-cache
      .gitignore:8:dev/stan-cache/   dev/stan-cache

A probe file inside it does not appear in `git status`. The
`.gitignore` diff adds `dev/stan-cache/` and also `dev/brmslp/`; the
second is a lane scratch directory, harmless, though unmentioned in the
lane's write-up.

`.Rbuildignore` needs no new entry: it already carries `^dev$`, which
excludes the whole directory from the tarball. Confirmed below.

---

## 7. The retirement and the prose: CONFIRMED

**The deleted test was the only brms MCMC run in the suite.** At the
base commit, `git grep "brms::brm(" a233c3c -- tests/` returns exactly
one hit, `test-brms-agreement.R:589`, inside the deleted
`"distributional sleepstudy agrees with brms posterior means"` test.
After the change the only `brms::brm(` in `tests/` is
`helper-brms.R:190`, which passes `empty = TRUE` and so compiles and
samples nothing. Every other Stan call in the suite is
`rstan::sampling(..., chains = 0)`, which does not sample. The one
remaining true sampler is `tmbstan` in `test-lkj.R`, which is frmtmb's
own sampler and unrelated.

**Nothing depended on it.** The deleted block created `bform` and
`bfit` inside its own `test_that()`, and `grep -rn "bfit" tests/`
returns no hits anywhere. No fixture, no helper, no snapshot.

**The prose reads correctly.** README's validation bullet now says the
tier "verifies that our log-likelihood equals the Stan program's log
density at the estimate", which is what the tests actually assert (the
old wording, "our estimates equal the mode", described the estimand the
lane retired). CONTRIBUTING's gate row matches. The NEWS bullet is
accurate on every claim I checked: the 1e-6 tolerance, the zero
gradient, the joint-density treatment of random effects, the reason the
posterior-mean comparison is gone, the file name and the gate variable.

---

## 8. The plan document: CONFIRMED a strict superset

`diff -u` between the maintainer's untracked original on main
(`C:\Users\adf44\source\r\frmtmb\dev\brms-likelihood-tests.md`, 137
lines) and the lane's version in the worktree (491 lines):

- Removed lines whose first character is `-`: exactly **one**.
- That one line is `Status: plan, 2026-09-04. Not implemented.`,
  replaced by a Status paragraph recording what landed.
- Everything else is pure addition: the "Implementation log" section
  and everything under it, appended after the maintainer's last line.

Nothing of the maintainer's was dropped, reworded, resequenced or
silently corrected. The claim holds exactly.

---

## 9. Scope: CLEAN

Worktree, `git diff --name-only a233c3c` plus untracked:

    M .gitignore
    M CONTRIBUTING.md
    M NEWS.md
    M README.md
    M tests/testthat/helper-brms.R
    M tests/testthat/test-brms-agreement.R
    ?? .github/workflows/brms-likelihood.yaml
    ?? dev/brms-likelihood-tests.md
    ?? tests/testthat/test-brms-likelihood.R

**`R/` is untouched.** No file under `R/` appears in the diff or in the
untracked list. The tests reach `ord_tau_from_raw()` as a package
internal under `load_all`, which is why no export was needed. Confirmed.

Diffstat: 489 insertions, 42 deletions across the six tracked files;
the 42 deletions are the retired sleepstudy test (37 lines) and the
five-line comment rewrite above it.

**Main is clean and unmoved.** `git status --short` in
`C:\Users\adf44\source\r\frmtmb` reports only
`?? dev/brms-likelihood-tests.md`, the maintainer's own untracked plan.
HEAD is `3a856e6`, as given. I made no edit there and ran nothing that
writes to it, except `R CMD build` which I pointed at a scratch output
directory (see the note in my edit log: it needed a probe file under
`dev/stan-cache/` in the WORKTREE, which I created and removed).

---

## 3. Independent hand-check of three rows: ALL PASS

`helper-brms.R` was never sourced for this. The translations were
written from the plan document and from my own reading of frmtmb's
sources. Compiled programs were reused from the cache, but each one was
verified by comparing its stored `model_code` with the freshly
generated code (whitespace-normalized) before use, so nothing about the
lane's keying was taken on trust.

### Row 1, gaussian `y ~ x + z` with `sigma ~ x`

My translation: `b = fixef$mu[c("x","z")]`;
`Intercept = b0 + sum(colMeans(X)[c("x","z")] * b)`; same for sigma.

    lp   = -239.538977525
    ours = -239.538977525   (logLik)
    measured constant = -5.68e-14   tol 2.40e-04   CHECK A PASS
    max|grad|         =  2.638e-04  vs 1e-3        CHECK B PASS

The gradient matches the lane's recorded 0.000264 to three figures.

### Row 12d, `acat(y ~ x)` — the family the threshold bug hit

I confirmed the storage convention from the source myself before
translating: `R/families.R:2547` gives acat
`ord_tau_init(y, ordered = FALSE)`, and `ord_tau_from_raw()` at 2287-90
returns `raw` unchanged when `ordered` is FALSE. So `tau_raw` IS the
threshold vector for acat, with no exp/cumsum. (2436 gives sratio TRUE,
2480 cratio FALSE, and cumulative TRUE — matching the helper's
`fam %in% c("cumulative", "sratio")` exactly.)

I also confirmed the ordinal centering shape from standata rather than
from the helper: `Kc == ncol(X) == 1`, so there is no intercept column,
and the thresholds take the centering with a minus sign.

    lp   = -307.792386991
    ours = -307.792386991
    measured constant = 0.000000e+00  tol 3.08e-04  CHECK A PASS
    max|grad|         = 2.781e-05     vs 1e-3       CHECK B PASS

A machine-exact zero. (The lane records 4.36e-06 for this gradient
against my 2.78e-05; both are four decades below the threshold and the
difference is optimizer noise, not a mapping residual.)

Control, in the same script: substituting the WRONG (ordered) transform
moves `log_prob` by **-17.2373 nats**. Hold that number for item 4.

### Check C, sleepstudy `(Days | Subject)` with `sigma ~ Days`

I took brms's within-group order from standata directly rather than
from any table, and verified it: `Z_1_1` is all ones and `Z_1_2` equals
`Days`, so the order is Intercept then Days. My own map:
`sd = sqrt(diag(V))`, `L = t(chol(V / tcrossprod(sd)))`,
`z = solve(diag(sd) %*% L, t(r))` with `r = ranef()[[1]][levels, coefs]`,
and `logJ = n_levels * (sum(log(sd)) + sum(log(diag(L))))`.

    n_levels = 18   logJ = 90.491828
    z occupies 36 unconstrained slots (18 x 2), found by perturbation
    lp   = -863.177038673
    ours = -863.177038673   (-f(last.par.best) + logJ)
    measured constant = 2.27e-13   tol 8.63e-04   CHECK A PASS
    max|grad| on z    = 1.044e-14  vs 1e-3        CHECK B PASS
    sd/L reproduce VarCorr to 8.88e-16

The gradient matches the lane's recorded 1.04e-14 exactly.

**Verdict on item 3: the harness computes what it claims to compute.**
Three independent reimplementations land on the same numbers, and two
of the three gradients match the lane's table to the digit.

---

## 4. The two translator bugs

Method: revert the one fix in a scratch copy of `helper-brms.R`, run
the affected row through a measuring (not asserting) variant of
`brms_lp_check`, and run the current helper as a control.

### Bug 1, group coefficient order: CONFIRMED TO THE DIGIT

Revert = `brms_group_info()` orders its rows by `coef` (alphabetical,
as `get_prior()` returns them) instead of by `cn`.

| | check A miss | max abs gradient on z |
|---|---|---|
| reverted (alphabetical) | **-1148.7330 nats** | **294.8430** |
| current helper | 0.0000 | 0.0000 |

The lane recorded 1149 nats and 295 on z. Both reproduce exactly. This
is a clean confirmation, and it also demonstrates that check B does the
job it is there for: the value error and the gradient error appear
together, and the gradient localizes the fault to the z block.

### Bug 2, ordinal thresholds: MECHANISM CONFIRMED, MAGNITUDES DO NOT
### REPRODUCE

Revert = `brms_ord_thresholds()` passes `ordered = TRUE` unconditionally.

| family | reverted, check A miss | reverted, max grad | current helper |
|---|---|---|---|
| cumulative | 0.0000 | 0.0002 | 0.0000 |
| sratio | 0.0000 | 0.0001 | 0.0000 |
| cratio | **-0.0234** | **1.5295** | 0.0000 |
| acat | **-17.2373** | **44.8362** | 0.0000 |

What is confirmed, and it is the substance of the finding:

- Exactly the two families that store raw thresholds break, and the two
  that store `(tau_1, log increments)` are untouched. That is the
  conditional the fix introduces, and it is right.
- Both break in value AND in gradient, so check B catches it.
- The current helper is machine-exact on all four.

What does NOT reproduce: the plan document says "Rows 12c and 12d
failed check A by 0.83 and 9.6 nats, with gradients of 9 where zero
belongs." I measure **-0.0234 and -17.2373 nats**, with gradients of
1.53 and 44.84. Not the stated amounts, and not close to them.

I am confident in my numbers rather than the lane's, because the acat
figure was produced twice by independent routes: -17.2373 from the
reversion above, and -17.2373 from the control inside my hand-written
spot-check in item 3, which shares no code with the helper. The data,
seed and formula are exactly the test file's.

So the finding is real and the fix is correct, but the three numbers
recorded next to it in `dev/brms-likelihood-tests.md` are stale or
mistranscribed. That matters because those numbers are the evidence a
future reader would use to decide the conditional is load-bearing. They
should be corrected to the measured values. (Note also that the "0.83"
is numerically the same as the mo() shared-simplex log-likelihood gap
of 0.834911 measured in item 2, which raises the possibility of a
transcription slip between sections; I cannot prove that, only flag it.)

---

## 2. Running the suite: GREEN WARM, RED COLD. This is the blocking
##    finding.

Both runs: `load_all()` in one process, `FRMTMB_BRMS_FIT_TESTS=true`,
`NOT_CRAN=true`, the C++17 repair in place.

| | cold (empty cache) | warm (24 programs cached) |
|---|---|---|
| wall clock | **1477.7 s (24.6 min)** | **23.8 s** |
| tests | 19 | 19 |
| assertions | **93** | **97** |
| failed | 0 | 0 |
| errors | **2** | **0** |
| skipped | 0 | 0 |
| programs in cache at end | 24 | 24 |

The warm run reproduces the lane's headline exactly: 19 tests, 97
assertions, all green, nothing skipped. **The cold run does not.**

### The two cold errors

    3  row 1: distributional gaussian, y ~ x + z with sigma ~ x   ERROR  0.43 s
    4  row 2: monotonic effects, y ~ mo(inc) + z                  ERROR  0.29 s

Both die in `rstan::sampling` at `helper-brms.R:476`, down through
`object@mk_cppmodule` to `Rcpp::Module(module, mustStart = TRUE)`.

These are not arbitrary. Test 1 (the translator round-trip) uses the
*same Stan program* as test 3 (row 1), and test 2 uses the same program
as test 4 (row 2). They are the only two programs the file uses twice.

### Minimal reproduction, seven lines

    m1 <- brms_stan_model(code)        # compiles, saveRDS
    rstan::sampling(m1, chains = 0)    # OK
    m2 <- brms_stan_model(code)        # same session, now a readRDS
    rstan::sampling(m2, chains = 0)
      -> Failed to initialize module pointer:
         Error in FUN(X[[i]], ...): NULL value passed for DllInfo

**The rule:** within one R session, the first `brms_stan_model(code)`
compiles and works; every later call for the same code re-reads the RDS
and returns an object whose DSO cannot be initialized. A model read
from cache in a session that has *not* compiled it is fine, which is
why the fully warm run is green and why my hand-checks in item 3
worked.

### Why this is blocking, not cosmetic

The workflow's own gate turns this into a hard failure:

    if (sum(df$failed) > 0 || any(df$error)) stop(...)

- **The first ever run of this workflow has an empty cache.** It is a
  cold run. It errors on two tests and the job goes red. There is no
  configuration in which the first run is green.
- Any run whose `restore-keys` fallback misses those two programs
  (a changed formula in either test file) compiles them fresh and hits
  the same wall.
- Whether it self-heals on run 2 depends on whether `actions/cache`
  writes its post-job cache after a failed job. I did not verify that
  semantic and will not assert it; but even in the best case the
  workflow's debut is red, on a job whose entire value is being
  trusted.

It also means the lane's results table was produced against a warm
cache and the file has never been run green from scratch.

### The fix, and it is small

Memoize per session so a program compiled in this process is reused
rather than re-read. `tests/testthat/helper-brms.R:86-98`:

```r
.brms_stan_models <- new.env(parent = emptyenv())

brms_stan_model <- function(code) {
  key <- brms_stan_cache_key(code)
  hit <- .brms_stan_models[[key]]
  if (!is.null(hit)) {
    return(hit)                       # compiled or read earlier HERE
  }
  path <- file.path(brms_stan_cache_dir(), paste0(key, ".rds"))
  if (file.exists(path)) {
    mod <- try(readRDS(path), silent = TRUE)
    if (!inherits(mod, "try-error")) {
      .brms_stan_models[[key]] <- mod
      return(mod)
    }
  }
  mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(mod, path)
  .brms_stan_models[[key]] <- mod
  mod
}
```

This removes the failure and is faster: the warm run stops paying for
24 `readRDS` calls of DSO-carrying objects.

I did NOT apply this. It is more than the "minor listed fix" my brief
permits me to edit, and it is the maintainer's call whether a test
helper gets a session cache. It is the first punch item.

### The other two files

    test-message-uniqueness.R   1 test,  6 assertions, 0 failed, 0 skipped,  2.5 s
    test-bracket-access.R       3 tests, 7 assertions, 0 failed, 0 skipped,  0.9 s

Both clean. `test-bracket-access.R` scans `R/` and `extensions/` only,
so the new helper is out of its scope; the lane nonetheless followed the
convention, writing `fit$estimates[["tau_raw"]]` rather than
`estimates$tau_raw`, which is correct since `estimates` is on that
file's hazard list.

`test-brms-agreement.R` is reported separately below.

### Program count: the workflow's sizing comment is wrong

The cold run left **24** compiled programs, not the "handful" the
workflow comment claims, and not 17. At the 55-70 s per program I
measured, that is 22-28 minutes of compiling on this box, and the cold
run took 24.6 minutes end to end. 90 minutes is still an adequate
timeout on Linux, so the number is safe, but the rationale comment at
`.github/workflows/brms-likelihood.yaml:40-43` ("about four minutes per
Stan program", "a handful of programs") is wrong in both factors and
should be corrected to the measured 24 programs at about a minute each.


### The measured constants, row by row

testthat prints pass/fail, not the numbers, so I called `brms_lp_check()`
directly against the warm cache to recover what the plan's results table
records. "const" is `log_prob` minus frmtmb's log density.

| row | measured constant | admitted | max abs gradient |
|---|---|---|---|
| 1 gaussian, sigma ~ x | -5.68e-14 | 0 | 2.64e-04 |
| 2 mo(inc) + z | **0.693147181** | lgamma(3) | 2.00e-05 |
| 3b mo(inc):z | **0.693147181** | lgamma(3) | 5.82e-07 |
| 5 nonlinear | -1.42e-14 | 0 | 1.55e-06 |
| 12 cumulative | 0 | 0 | 1.63e-04 |
| 12 sratio | 0 | 0 | 5.42e-05 |
| 12 cratio | 0 | 0 | 1.75e-04 |
| 12 acat | 0 | 0 | 2.78e-05 |
| 12e sratio + cs(z) | 0 | 0 | 2.56e-04 |
| 13 categorical | -5.68e-14 | 0 | 3.42e-05 |
| 14a cens | -2.84e-14 | 0 | 3.75e-05 |
| 14c se | -1.14e-13 | 0 | 4.90e-05 |
| 15 binomial trials(n) | -2.84e-13 | 0 | 2.85e-09 |
| 16 zero-inflated poisson | 0 | 0 | 2.00e-05 |
| 20 weights(w) | 0 | 0 | 1.84e-04 |
| 21a poisson | -1.71e-13 | 0 | 7.96e-05 |
| 21b Gamma(log) | 5.68e-14 | 0 | 7.17e-06 |
| 21c negbinomial | -5.68e-13 | 0 | 2.77e-04 |
| 21d bernoulli | 0 | 0 | 6.23e-09 |

Every constant is machine zero or exactly `lgamma(3) = 0.693147181`,
the admitted flat Dirichlet on the mo() simplex. Nothing is a
suspicious near-zero and nothing needed a widened tolerance. Every
gradient is at least four decades below the 1e-3 threshold. This is a
clean reproduction of the plan's results table and it is the strongest
single piece of evidence that the identity holds.

Two caveats on this table, both mine rather than the lane's. The
individual gradient magnitudes do not all match the lane's recorded
values: about half agree to three figures (rows 1, 3b, 5, 13, 15, 16,
21a) and the rest differ by up to two decades. That is expected and is
not a finding. These numbers are frmtmb's own optimizer residuals, and
my script builds some data frames in a different RNG order than the
test file, so several rows are fitted to slightly different data. The
constants, which are the actual claim, are invariant to all of that.
The mixture row is absent because my standalone script hit brms
masking frmtmb's `mixture()`; it passes in the suite itself.

### test-brms-agreement.R

    19 tests, 187 assertions, 0 failed, 0 errors, 0 skipped, 106.3 s

Clean. It compiles its two Stan programs directly with
`rstan::stan_model()` and does not use the cache helper, so it is
unaffected by the defect above.

### Why the lane never saw the cold failure

Circumstantial but consistent. The worktree's `dev/stan-cache` holds 27
programs written between 23:28 and 00:16, and the lane's plan document
was saved at 00:21. That cache accumulated across many iterative runs,
so by the time the full file was run end to end, tests 1 and 3 always
found their shared program already on disk from an earlier session and
took the readRDS-first path, which works. The file was never run from
an empty cache. (My cold run used a separate scratch directory and
started at 00:29, so it did not contribute to that cache; the three
extra programs there are the lane's row 8 and `mo(inc) * z` attempts.)

---

# VERDICT: GO WITH FIXES

The science is sound and I verified it independently. The translator
computes what it claims: three rows reimplemented by hand from the plan
and from frmtmb's sources land on the same log densities to 1e-13 or
better, and two of three gradients match the lane's table to the digit.
Both translator bugs are real, and the group-ordering one reproduces to
four significant figures. The mo() divergence is real, correctly
diagnosed, and correctly adjudicated in brms's favor. Scope is clean,
`R/` is untouched, the plan document is a strict superset, the prose is
accurate, and the retirement removed the suite's only brms MCMC run
with nothing depending on it.

What stops this being a plain GO is that the test file has never been
run from an empty cache, and it fails when it is. That plus the CI
toolchain exposure means the workflow as written cannot go green on its
first run. Both are defects in the scaffolding, not in the claim under
test, and both fixes are small.

## Punch list

**Blocking**

1. `tests/testthat/helper-brms.R:86-98`, `brms_stan_model()`. Within one
   session the first call compiles and works; every later call for the
   same code re-reads the RDS and returns a model whose DSO cannot
   initialize (`NULL value passed for DllInfo`). The file uses two
   programs twice, so a cold run errors on "row 1" and "row 2": 93
   assertions and 2 errors instead of 97 and 0. Fix: memoize the model
   per session in an environment keyed on the cache key, before the
   `file.exists()` branch. Patch given in item 2 above.

2. `.github/workflows/brms-likelihood.yaml`, before the run step. RSPM
   serves rstan 2.32.7 with StanHeaders 2.39.1 on both noble and jammy,
   the exact pair that cannot compile a Stan model here, and the fault's
   three ingredients are all platform-generic (verified in rstan's
   source, not assumed). Add the `~/.R/Makevars` step with
   `CXX17FLAGS = -O2 -Wall -std=gnu++17`. It is a no-op if I am wrong.

**Should fix**

3. `NEWS.md:2364`. "mo() interactions share their variable's simplex
   (brms convention)" is factually wrong: brms builds one simplex per
   term. This false rationale is why the divergence survived two
   audits. Fix it regardless of what is decided about the behavior.

4. `dev/brms-likelihood-tests.md`, the ordinal-threshold finding. The
   recorded misses of 0.83 and 9.6 nats with gradients of 9 do not
   reproduce. Measured, by two independent routes: cratio -0.0234 nats
   with gradient 1.53, acat -17.2373 with gradient 44.84. Correct the
   numbers; the finding itself is right.

5. `dev/brms-likelihood-tests.md`, the divergence section. Cross-
   reference `dev/brms-vignette-audit.md:513-517`, which already
   recorded the shared-simplex divergence and that "brms fits two".
   This is the second independent report, not the first, and saying so
   strengthens the case for acting on it.

**Nits**

6. `.github/workflows/brms-likelihood.yaml:40-43`. The timeout rationale
   says "about four minutes per Stan program" and "a handful of
   programs". Measured: 24 programs at 55-70 s each. The 90-minute
   timeout is still right; the reasoning behind it is not.

7. `dev/brms-likelihood-tests.md`, cache section. "the code is a
   function of the formula and the family, never of the data" is not
   quite true: `y ~ x` and `y ~ x + I(x^2)` generate byte-identical
   programs. Harmless and in fact beneficial, but reword so nobody
   reasons about collisions from it.

8. `dev/brms-likelihood-tests.md`, toolchain section. The `-std=gnu++17`
   that the plugin flag overrides comes from rstan's own `CXX17STD`
   override, not from a Makeconf default (R 4.6.1's default is
   `gnu++20`). Conclusion unchanged.

## Edits I made to the worktree

One file, created, not modifying anything the lane wrote:

- **`C:\Users\adf44\source\r\frmtmb-wt-brmslp\dev\review-brms-likelihood.md`**
  — this document, which my brief required.

Also, transiently, and both removed before I finished:

- `dev/stan-cache/probe.rds`, created twice (once to test
  `git check-ignore`, once to test the `R CMD build` exclusion) and
  deleted both times. `dev/stan-cache` itself is the lane's, not mine.

I committed nothing, and I did not touch the main checkout. Everything
else I ran wrote only under my scratchpad, including my own private
library and three separate Stan caches.

## What the maintainer must know about the local rstan toolchain

**The installed pair cannot compile any Stan model, and the error names
nothing useful.** rstan 2.32.7 with StanHeaders 2.39.1 on R 4.6.1. What
you see is `Error in sink(type = "output") : invalid connection`, which
is rstan's output-sink unwind failing on the way out of a compile error
it never shows you. The whole failing run is 19 lines and the compiler
diagnostics are not among them.

**The cause.** rstan emits `// [[Rcpp::plugins(cpp14)]]`, which becomes
`-std=c++1y` in `PKG_CXXFLAGS`. rstan also asks for the C++17 toolchain
(`CXX='$(CXX17) $(CXX17STD)'`), so the command line opens with
`-std=gnu++17` and then, further along, contradicts it. Last `-std`
wins, the unit compiles as C++14, and StanHeaders 2.39.1 needs C++17.
The first of the 142 errors is `'void_t' is not a member of 'std'` at
`StanHeaders/include/stan/math/prim/eigen_plugins.h:13`.

**The repair.** Point `R_MAKEVARS_USER` at a file containing

    CXX17FLAGS = -O2 -Wall -mfpmath=sse -msse2 -mstackrealign -std=gnu++17

`CXX17FLAGS` specifically. rstan passes it as `CXXFLAGS`, which lands
after `PKG_CXXFLAGS` on the command line and therefore wins. `CXX17STD`
is the intuitive lever and it does not work, because it lands before.

**This is not confined to this lane.** StanHeaders 2.39 has already
cost this repo once: commit 049fcd5 skips an LKJ marginal test because
a tmbstan built against 2.39 silently samples a standard normal, since
tmbstan's generator patches only the first of the two log-density
overloads stanc now emits. Two independent breakages from the same
upgrade suggests the durable fix is to pin StanHeaders to the 2.32.x
that rstan 2.32.7 was built against, rather than to keep adding
per-symptom workarounds.

**And CI is exposed to the same thing**, which the lane could not check
and I could: RSPM serves that identical pair on Linux. See punch item 2.
