# Lane `generics`: Phase 2.5, frmtmb must not break brms

Worktree `frmtmb-wt-generics`, branch `wt-generics`, based on the Phase
2 release. Items 2.5a, 2.5b, 2.5c and 2.5d. Nothing is committed.

Round 2, after `dev/reviews/20260915-generics.md`. Four blockers and
twelve nits; what changed in response is marked **R2** below. The
round-1 mechanism, a binding SWAP plus one load hook per owner, is GONE
and replaced by an active binding, and three numbers on the round-1
page did not survive and are corrected with their constructions.

Punch round 2, after the "Round 2 re-check" in the same review file.
One BLOCKER, E, which round 1's own `hypothesis()` fix caused in
`frmtmb.sample`, and five nits; what changed is marked **R3** below.

Private library `C:/Users/adf44/source/r/generics-lib`; the base commit
read-only from `C:/Users/adf44/source/r/rellib-r3`; the review's
library `C:/Users/adf44/source/r/genrev-lib`, which still holds the
round-1 SWAP build and is used READ-ONLY as a third arm so the two
designs are timed and tested side by side rather than across runs.
R 4.6.1, brms 2.23.0, lme4 2.0.6, posterior 1.7.0, loo 2.10.1,
rstantools 2.7.1, bayesplot 1.16.0, nlme 3.1.171, generics 0.1.4.

Every count below is emitted by a summarizer into a marked block and
pasted verbatim; none is typed.

## The defect, restated from the measurement

`UseMethod()` consults the method table of the namespace where the
generic **it reached** was defined. frmtmb defined its own generic for
28 borrowed names and exported all of them, so a frmtmb generic sitting
above brms on the search path sent a `brmsfit` into frmtmb's table,
which has no entry for that class. brms's method was registered the
whole time, in brms's table, unreachable.

The measurement that matters most is about the fix that was already
there. **The base commit already carried 38
`S3method(pkg::generic, class)` directives over 36 distinct generics**
(`conditional_effects` and `hypothesis` are each registered for two
classes), and the defect reproduces there. Delayed registration puts
frmtmb's method in the owner's table. It does nothing about frmtmb's
rival generic. That is why the defect survived, and why the fix is not
more delayed registration.

## What I changed

| file | what |
|---|---|
| `R/generic-owners.R` | NEW. The ownership table, the ACTIVE BINDING that resolves each borrowed name to its owner, and the `importFrom`/`export` pair for the four static names |
| `R/scales.R` | NEW. `?frmtmb-scales`, the scale contract (item 2.5d) |
| `R/zzz.R` | `.onLoad()` calls `frm_install_generics()` |
| `R/methods-fit.R` | `fixef`, `ranef`, `VarCorr` generics DELETED, imported from nlme; `VarCorr.frmtmb_fit()` gains nlme's `sigma`; the doc blocks move onto the methods |
| `R/sugar.R` | `refit` generic DELETED, imported from generics; `refit.frmtmb_fit()` also registered on `lme4::refit` |
| `R/draws-generics.R` | `nvariables()` gains `...`; `posterior_summary()`'s first formal becomes `x`; **R2** ten new `frmtmb_fit` refusals |
| `R/loo.R` | **R2** `loo_compare.default`'s comment described the defect this change removes |
| `R/predict.R` | `@seealso [frmtmb-scales]` on `predict()`, `fitted()`, `residuals()` |
| `DESCRIPTION` | `nlme` moves Suggests to Imports; `generics` added to Imports |
| `_pkgdown.yml` | **R2** `frmtmb-scales` indexed |
| `NEWS.md` | **R2** the breaking-change bullet, no version number |
| `NAMESPACE`, `man/*` | roxygenised |
| `tests/testthat/test-generic-collision.R` | NEW, items 2.5a to 2.5c |
| `tests/testthat/test-scale-contract.R` | NEW, item 2.5d |
| `extensions/frmtmb.sample/R/methods-draws.R` | three method signatures follow core's new generic formals, and one `@param` |

## BLOCKER A. The active binding, which does dominate

The round-1 mechanism swapped frmtmb's binding for the owner's inside
`.onLoad`, and set a load hook per owner for an owner that turned up
later. Review built a third design on three throwaway packages, and it
is better on every axis I can measure, so it is what ships.

`.onLoad` removes frmtmb's own binding for each borrowed name and
installs an ACTIVE BINDING in its place, whose function returns the
owner's generic when the owner is loaded and frmtmb's own when it is
not. This can only be done from `.onLoad`: `loadNamespace()` runs
`registerS3methods()` at line 345 and `.onLoad` at line 367, so every
`S3method()` directive has already resolved against frmtmb's own
generic and landed in frmtmb's own table, which is exactly where the
fallback needs them, and the namespace is not sealed yet.

What it buys, on the SAME axes, each measured rather than argued. The
round-1 SWAP build is still installed in the review's library, so every
row below is one run with three arms, not a number quoted from another
session:

| axis | BASE | SWAP (round 1) | ACTIVE (now) |
|---|---|---|---|
| `library(brms); library(frmtmb)`, class method lost | 27 of 27 | 0 of 27 | 0 of 27 |
| `library(frmtmb); library(brms)` | 0 of 27 | 0 of 27 | 0 of 27 |
| `library(frmtmb)`, brms only LOADED | 27 of 27 | 0 of 27 | 0 of 27 |
| `library(lme4); library(frmtmb)` | 5 of 5 | 0 of 5 | 0 of 5 |
| every owner absent, frmtmb's own methods not frmtmb's **(R3: the refusals' win, not the binding's)** | 10 of 26 | 10 of 26 | 0 of 26 |
| `getS3method()` lost, brms first | 27 of 27 | 0 of 27 | 0 of 27 |
| owner unloaded and reloaded | n/a | STALE, `loo.brmsfit` unreachable | correct |
| owner exports a NON-generic under the name | n/a | adopts it | keeps frmtmb's |
| `frmtmb.sample` re-exports, brms merely LOADED | 26 of 26 | 23 of 26 | 0 of 26 |
| `hypothesis()` shadowing note, brms loaded | 1 | **0** | 1 |
| `unlockBinding()` in the namespace | no | yes | no |
| `R CMD check` NOTEs from this lane | none | one, `unlockBinding` | see "What was run" |
| load hooks | none | six | none |

Sources: `dev/generics-out2/blocks/collision.txt` for the first six
rows, `dev/generics-out2/blocks/propagate.txt` for the re-export row,
`dev/generics-out2/blocks/notice.txt` for the `hypothesis()` row, and
the collision test below for the rest.

**R3. Two rows in that table are not the binding's.** The
owner-absent row, 10 of 26 to 0 of 26, is BLOCKER B's ten refusal
methods: review registered those same ten on the SWAP build in-session
and got 0 of 26 there too (`dev/genrev-r2-attrib.R`). It holds on
either mechanism and says nothing about which one ships. The rows that
ARE the binding's are the stale owner, the non-generic owner, the
`frmtmb.sample` re-exports with brms only loaded, and the NOTE.

**On 23 against the review's 22. R3: the round-2 reason was wrong
and is withdrawn here.** Round 2 said the review's probe did not count
`posterior_summary` because a `.default` answered it. It did count it:
`posterior_summary` is on the review's own lost list
(`dev/genrev-out/p-sampleprop.txt`). The one-name difference is in the
NAME LISTS, not the method. The review's list omits `loo_compare` and
carries `expose_functions` in its place, and `frmtmb.sample` does not
re-export `expose_functions` (checked against its NAMESPACE), so with
only `frmtmb.sample` attached that name is not on the search path and
cannot be lost there. The review's 22 are exactly this lane's 23
without `loo_compare`. The count stood either way, and ACTIVE is 0 of
26 by both.

The last two rows were not obvious. Active-ness PROPAGATES through
`importIntoEnv()`, so the attached `package:frmtmb` binding is active
AND so is the binding in a package that re-exports the name from
frmtmb. That closes what round 1 filed as defect 5, which the swap
could not reach because both are copies taken at attach time.

**The three designs, run against the same unmodified test file.** The
swap build is still installed in the review's library, so this is one
comparison and not three:

    tests/testthat/test-generic-collision.R, 11 blocks, 43 assertions
      BASE   (rellib-r3)        PASS 24  FAIL 19
      SWAP   (genrev-lib)       PASS 36  FAIL  7
      ACTIVE (generics-lib)     PASS 43  FAIL  0

The swap build's seven failures are in four blocks and all four are
real: the stale namespace (`STALE TRUE`, `loo.brmsfit` unreachable);
the binding not being active in either the namespace or the attached
environment, with `unlockBinding` present in the sources; adopting a
non-generic where the active binding refuses to; and `hypothesis()`
carrying work in its generic, which the next section is about.

### A regression both designs had, which only `R CMD check` could see

**R2, found during this round and not by the review.** The first
`R CMD check` of the active-binding build reported
`Status: 1 ERROR`: `test-naming-collisions.R` lost 8 assertions. The
one-file-per-process suite had passed the same file.

`hypothesis()`'s GENERIC did work. It armed the reserved-name shadowing
note around `UseMethod()`. Once frmtmb's binding resolves to brms's
generic, which is a bare `UseMethod()`, that work is simply not done
and the note never fires. `R CMD check` runs the suite in ONE process,
where an earlier file had already loaded brms's namespace; one file per
process never has brms loaded when that file runs. The SWAP build had
the identical defect: it swapped the binding the same way.

```
== the regression R CMD check found, and the fix ==
dev/generics-notice.R: the shadowing note, with and without brms
--- BASE, nobrms
brms loaded            FALSE
hypothesis resolves to frmtmb
shadowing notes        1
--- BASE, brms
brms loaded            TRUE
hypothesis resolves to frmtmb
shadowing notes        1
--- SWAP, nobrms
brms loaded            FALSE
hypothesis resolves to frmtmb
shadowing notes        1
--- SWAP, brms
brms loaded            TRUE
hypothesis resolves to brms
shadowing notes        0
--- ACTIVE, nobrms
brms loaded            FALSE
hypothesis resolves to frmtmb
shadowing notes        1
--- ACTIVE, brms
brms loaded            TRUE
hypothesis resolves to brms
shadowing notes        1
```

The fix is structural. The arming moves into
`hypothesis.frmtmb_fit()` and `hypothesis.frmtmb_multiple()`, where it
runs whichever generic dispatched, and the save-and-restore in
`hyp_shadow_arm()` already makes nesting safe. And a new block in the
collision test asserts that EVERY generic in the ownership table is a
bare `UseMethod()` when no owner is loaded, so work cannot be put back
into a replaceable generic without a test failing. That block was seen
failing on the SWAP build, naming `hypothesis`, and passes on ACTIVE
(`dev/generics-bodies.R` gives 0 of 25).

This is the lane rule about one file per process read the other way.
The rule exists because a shared process HIDES leakage; here the shared
process was the only run that could EXPOSE a dependency on what another
file had loaded. Both kinds of run were needed, and the checklist had
only one of them as a gate for correctness.

### BLOCKER E. Round 1's `hypothesis()` fix removed the note from draws

**R3.** Round 1 moved the arming of the reserved-name shadowing note
out of core's `hypothesis()` generic and into core's two methods,
`hypothesis.frmtmb_fit()` and `hypothesis.frmtmb_multiple()`. That was
right for core. But there was a THIRD method that had relied on the
generic's arming: `frmtmb.sample`'s `hypothesis.frmtmb_draws()`, which
calls `hyp_env_vals()`, never armed the note, and could not, because
`hyp_shadow_arm()` was not on core's internal export list. So on draws
the note was gone in EVERY session, brms loaded or not, and no test in
either package covered it. Review measured it with
`dev/genrev-r2-drawsnote.R`: BASE 1 note on the fit and 1 on its draws,
the round-1 tree 1 and 0.

**The reason my guard missed it is precise, and it generalizes.** The
"no shared generic carries work in its own body" test stops work
coming BACK into a generic. It cannot see work that LEFT the generic
and did not reach every method. A guard that checks where work is must
also check where work went.

The fix, as review described it:

* `hyp_shadow_arm()` and `hyp_shadow_disarm()` join core's internal
  export list on `?frmtmb-sampling-api`, beside the `hyp_*` helpers
  `frmtmb.sample` already uses, with a paragraph saying every
  `hypothesis()` method must arm the note itself and why the generic
  cannot. `frmtmb.sample` uses `import(frmtmb)`, so it needs no
  NAMESPACE change of its own.
* `hypothesis.frmtmb_draws()` arms on entry and restores on exit, the
  same two lines as core's methods.
* Two tests in `frmtmb.sample`'s `test-draws-methods.R`, both seen
  failing first. One asserts the note appears exactly once on draws,
  with the same count on the fit as a control, so a construction that
  shadows nothing cannot pass it. It uses `fake_draws()`, because the
  note is emitted while the hypothesis is parsed against the fit and
  before any draw is read, so it needs no Stan build. The other is the
  WHERE-THE-WORK-WENT guard: every `hypothesis.*` method in both
  packages' namespaces must call the arm and disarm pair.

```
== BLOCKER E: the shadowing note on draws, frmtmb.sample ==
extensions/frmtmb.sample/tests/testthat/test-draws-methods.R, one process
dev/generics-runtests-sample.R; libraries named are put FIRST
build                                                      assertions / pass / fail
BASE core + BASE sample (rellib-r3), note test only        99 / 99 / 0
FIX core + UNFIXED draws method, before the fix            99 / 98 / 1
    failing: hypothesis() on draws gives the reserved-name note, once 
FIX core + FIX sample, note test only                      99 / 99 / 0
BASE core + BASE sample, with the guard                    101 / 100 / 1
    failing: every hypothesis() method arms the note itself 
FIX core + UNFIXED draws method, with the guard            101 / 99 / 2
    failing: hypothesis() on draws gives the reserved-name note, once 
    failing: every hypothesis() method arms the note itself 
FIX core + FIX sample, with the guard                      101 / 101 / 0

the guard names, on the construction that matters (FIX core, UNFIXED draws):
    `actual`:   "hypothesis.frmtmb_draws"

the review's own probe, dev/genrev-r2-drawsnote.R, FIX build, 1 chain x 300:
  brms-loaded FALSE
frmtmb_fit   note lines: 1
frmtmb_draws note lines: 1
  brms-loaded TRUE
frmtmb_fit   note lines: 1
frmtmb_draws note lines: 1
```

Read the rows in order. On BASE the note test passes, because the
generic armed; the guard fails on all three methods, because on BASE
none of them armed, which is correct for a guard about where the work
lives. With core fixed and the draws method not, both tests fail and
the guard names exactly `hypothesis.frmtmb_draws`. With both fixed,
both pass, and the review's own probe on the real sampler gives 1 note
on the fit and 1 on draws, with and without brms.

`hypothesis` was the only borrowed generic with work in its body at
the base commit, 1 of 28 (`dev/genrev-r2-basebodies.R`), so this one
method is the whole blast radius.

**What it costs.** One closure call per ACCESS instead of one per
session, so the cost model changed and was re-measured:

```
== active-binding cost, dev/generics-adoptcost.R ==
proc.time() tick here: 0.0100 s
names in the table 25, distinct owners 6

-- NO owner loaded: the fallback path, arms interleaved --
active binding, get('loo')                   1048576 reps     1.29 us
active binding, get('ngrps'), 2 owners        262144 reps     4.35 us
CONTROL get('frm'), ordinary binding         4194304 reps     0.64 us
CONTROL identity(1)                         16777216 reps     0.16 us
fixef(fit), a real accessor call              262144 reps    14.61 us
<environment: namespace:posterior>
<environment: namespace:loo>

-- posterior and loo loaded: the adopted path --
active binding, get('loo')                   1048576 reps     3.18 us
active binding, get('as_draws_df')           1048576 reps     3.28 us
CONTROL get('frm'), ordinary binding         4194304 reps     0.80 us
fixef(fit), a real accessor call               65536 reps    12.36 us

fallback path 1.29 us, adopted path 3.18 us
adopted path is 4.0x an ordinary binding lookup (0.80 us)
and 0.26 of one fixef() call (12.36 us), same round
```

3.18 us on the adopted path, 1.29 us on the fallback, against 0.80 us
for an ordinary binding lookup: 4.0x, and 0.26 of one `fixef()` call.
These names are reached once per user action, not in a loop. Two
instrument notes, both of which cost me a wrong number first:

* The arms are INTERLEAVED. Run one after another, the ratio to
  `fixef()` came out 0.73 in one run and 0.36 in the next on identical
  code, because `fixef()` drifted from 12 to 22 us between arms.
* The first spelling measured 8.89 us, and about 7 us of that was a
  `tryCatch()` on the hot path rather than the work:
  `dev/generics-adoptparts.R` times the pieces at `isNamespaceLoaded`
  0.20 us, `asNamespace` 0.92, `identical` 0.27, `getExportedValue`
  0.37, and the whole body inline 1.60 us. `asNamespace()` cannot fail
  behind an `isNamespaceLoaded()` that has just returned TRUE, so the
  handler came off. The two-owner path (`ngrps`, `refit`) was 49.13 us
  because it ran `sub()` over every entry of `search()`; as two vector
  operations it is 4.35 us.

## 2.5a and 2.5b: one mechanism, not two

The plan filed 2.5a and 2.5b as separate problems. The audit collapses
them, because **where a brms method is registered decides the fix** and
`dev/generics-audit2.R` measures that directly:

- brms IMPORTS 20 of the 27 shared names from their owners (11 from
  posterior, 3 from loo, 3 from nlme, 2 from rstantools, 1 from
  bayesplot), so `loo.brmsfit` lives in **loo's** table and
  `as_draws_df.brmsfit` in **posterior's**.
- brms DEFINES 7: `conditional_effects`, `expose_functions`,
  `hypothesis`, `posterior_summary`, `LOO`, `WAIC` and `ngrps`.

Either way the rule is the same: stop owning the name and dispatch
through the generic the owner defined. **The fallback the plan
reserved for the user, a distinct spelling in frmtmb, is not needed,
and no decision on it is being asked for.**

For nlme (Recommended, ships with R) and generics (base-package
dependencies only), `importFrom()` plus `export()`, which is what lme4
and brms do with nlme. For the five optional owners, the active
binding.

**What is NOT done, and why.** `s3_register()` is not vendored. The
plan asked for it; the measurement says it would be dead code. R has
supported `S3method(pkg::generic, class)` since 3.6.0, frmtmb declares
`R (>= 4.1)`, and frmtmb already uses that directive 38 times. rlang's
helper exists for packages that must register a class not known at
build time, or that support R < 3.6. Review confirmed this refusal.

## BLOCKER C, first half. The base damage is 27 of 27, not 25 of 27

**R2, and it is a correction to the INSTRUMENT.** Round 1 decided "did
dispatch work" by matching the string "no applicable method" in the
error. `posterior_summary` and `loo_compare` are the two names frmtmb
registers a `.default` for in its own table, so on the unfixed build a
`brmsfit` reaching frmtmb's generic does not error at all: it falls
silently into frmtmb's `.default`. The string test scored those two as
healthy. They are the WORST two of the 27, because a silent wrong
answer beats an error.

`dev/generics-check.R` no longer reads error messages. It reproduces
`UseMethod()`'s lookup, the `.__S3MethodsTable__.` of
`environment(generic)`, for the class and then for `default`, and
reports which package supplies the method that would run.

```
== Phase 2.5 generic-collision measurement ==
dev/generics-check.R + dev/generics-summary.R; ONE PROCESS PER ROW
seed 20260915, R 4.6.1, brms 2.23.0, lme4 2.0.6
DAMAGE is the loss of the CLASS method, found by reproducing
UseMethod()'s lookup in the method table of environment(generic),
not by reading an error message. The denominator is the generics
for which the control resolves a class method.

BASE = base commit, SWAP = round 1, FIX = the active binding

BASE  brms then frmtmb     class-method lost 27/27  getS3method lost 27/27
    SILENTLY answered by a .default instead:  loo_compare, posterior_summary 
    no method at all:   as_draws, as_draws_array, as_draws_df, as_draws_list, as_draws_matrix, as_draws_rvars, bayes_R2, conditional_effects, expose_functions, fixef, hypothesis, loo, LOO, nchains, ndraws, ngrps, niterations, nvariables, pp_check, prior_summary, ranef, VarCorr, variables, waic, WAIC 
SWAP  brms then frmtmb     class-method lost  0/27  getS3method lost  0/27
FIX   brms then frmtmb     class-method lost  0/27  getS3method lost  0/27
BASE  frmtmb then brms     class-method lost  0/27  getS3method lost  0/27
SWAP  frmtmb then brms     class-method lost  0/27  getS3method lost  0/27
FIX   frmtmb then brms     class-method lost  0/27  getS3method lost  0/27
BASE  brms loaded only     class-method lost 27/27  getS3method lost 27/27
    SILENTLY answered by a .default instead:  loo_compare, posterior_summary 
    no method at all:   as_draws, as_draws_array, as_draws_df, as_draws_list, as_draws_matrix, as_draws_rvars, bayes_R2, conditional_effects, expose_functions, fixef, hypothesis, loo, LOO, nchains, ndraws, ngrps, niterations, nvariables, pp_check, prior_summary, ranef, VarCorr, variables, waic, WAIC 
SWAP  brms loaded only     class-method lost  0/27  getS3method lost  0/27
FIX   brms loaded only     class-method lost  0/27  getS3method lost  0/27
BASE  lme4 then frmtmb     class-method lost  5/5  getS3method lost  5/5
    no method at all:   fixef, ngrps, ranef, refit, VarCorr 
SWAP  lme4 then frmtmb     class-method lost  0/5  getS3method lost  0/5
FIX   lme4 then frmtmb     class-method lost  0/5  getS3method lost  0/5

frmtmb's OWN methods on a frmtmb_fit: does the method that
would run belong to frmtmb? A foreign .default is a FAILURE.
  base-A   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  swap-A   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  fix-A    FITMETHOD_NOT_OURS 0 of 26 | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  base-C   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  swap-C   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  fix-C    FITMETHOD_NOT_OURS 0 of 26 | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  base-D   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  swap-D   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  fix-D    FITMETHOD_NOT_OURS 0 of 26 | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  base-L   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  swap-L   FITMETHOD_NOT_OURS 10 of 26 | FITMETHOD_LIST as_draws,as_draws_df,as_draws_array,as_draws_list,as_draws_matrix,as_draws_rvars,ndraws,nchains,niterations,nvariables | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5
  fix-L    FITMETHOD_NOT_OURS 0 of 26 | REFIT_OK TRUE | NGRPS 6 | POSTSUMMARY Estimate,Est.Error,Q2.5,Q97.5

the message a frmtmb_fit gets from two of the six that regressed:
  base-A   MSG as_draws: no applicable method for 'as_draws' applied to an object of clas
  base-A   MSG as_draws_df: no applicable method for 'as_draws_df' applied to an object of c
  swap-A   MSG as_draws: All list elements must be lists themselves.
  swap-A   MSG as_draws_df: All list elements must be lists themselves.
  fix-A    MSG as_draws: as_draws() needs posterior draws and a frmtmb_fit has none: frm(
  fix-A    MSG as_draws_df: as_draws_df() needs posterior draws and a frmtmb_fit has none: f

which namespace each generic resolves to:
  base-A   frmtmb=28
  fix-A    bayesplot=1 brms=7 generics=1 loo=3 nlme=3 posterior=11 rstantools=2
  fix-B    bayesplot=1 brms=7 generics=1 loo=3 nlme=3 posterior=11 rstantools=2
  fix-C    frmtmb=24 generics=1 nlme=3
  fix-D    bayesplot=1 brms=7 generics=1 loo=3 nlme=3 posterior=11 rstantools=2
  fix-L    frmtmb=23 lme4=2 nlme=3
```

## BLOCKER B. A regression my own metric hid

**R2.** On the six `as_draws*` generics a `frmtmb_fit` went from R's
own "no applicable method" to posterior's "All list elements must be
lists themselves", because a fit is a bare list, posterior's generic
has a default that walks into one, and frmtmb had no method for the
class at all. The round-1 counter reported `BROKEN_FRMTMBFIT` going
1 of 17 to 0 of 17, and **it went to zero BECAUSE of that**: it tested
for a string the new error does not contain.

Both halves are fixed, and the second matters more.

* The regression. `as_draws()`, `as_draws_array()`, `as_draws_df()`,
  `as_draws_list()`, `as_draws_matrix()`, `as_draws_rvars()`,
  `ndraws()`, `nchains()`, `niterations()` and `nvariables()` now have
  `frmtmb_fit` methods that refuse by name and point at
  `frmtmb.sample::frm_sample()`. That also discharges round 1's
  defect 3 rather than leaving it filed: `as_draws(fit)` gives the
  refusal `?as_draws` promises.
* The metric. The counter now asks whether the method that WOULD RUN
  belongs to frmtmb; a foreign `.default` is a failure. On that
  measure the base build is **10 of 26**, not 1 of 17, and the six
  `as_draws*` plus the four size generics are exactly the ten. The
  fixed build is 0 of 26. The block above carries both, with the
  before and after error text.

## The signature mismatches, and the fourth one

| name | base core | fix core | owner |
|---|---|---|---|
| `nvariables` | `(x)` | `(x, ...)` | posterior `(x, ...)` |
| `posterior_summary` | `(object, ...)` | `(x, ...)` | brms `(x, ...)` |
| `VarCorr` | `(x, ...)` | `(x, sigma, ...)` | nlme `(x, sigma, ...)` |
| `refit` **R2** | `(object, newresp, ...)` | `(object, ...)`, or `(object, newresp, ...)` with lme4 loaded | generics / lme4 |

**Only ONE caller breaks**, and round 1's wording implied three:
`posterior_summary(object = m)`. `nvariables` and `VarCorr` only GAIN
an argument. `refit` is a fourth changed signature whose formals now
depend on what else is loaded; no caller breaks on it, but a downstream
method author has to be told, and the NEWS bullet says so.

A source comment above the `nvariables` block asserted that posterior's
four size generics all take `x` alone. It is right about three.

### The sibling edit, which is not a widening of scope

`frmtmb.sample` imports these generics from core and registers
`frmtmb_draws` methods on them. An S3 method must carry every argument
of its generic, so aligning core's formals would have left
`frmtmb.sample` reporting a generic/method mismatch under `R CMD check`
that this lane created. Three one-line signature changes plus one
`@param`. Note the direction: `frmtmb.sample` was ALREADY registering
`posterior_summary.frmtmb_draws` on `brms::posterior_summary`, whose
first formal is `x`, while declaring `object`. Aligning core made an
existing mismatch visible.

## 2.5c: `ngrps`, decided by measurement

brms DEFINES `ngrps` and lme4 DEFINES `ngrps`; neither imports the
other, so they already collide with each other whether or not frmtmb
exists. Review confirmed by reading both tables: the two generics are
different closures and each one's method sits only in its own table, so
one `importFrom()` reaches exactly one of them.

```
== load cost, dev/generics-loadcost.R, R 4.6.1 ==
30 replicates, interleaved, one fresh process each
proc.time() tick here: min 0.0100 s, median 0.0200 s
Sys.time()  tick here: min 0.000001 s, median 0.000002 s
The first is the instrument that produced the withdrawn 0.040 s.

arm                     min(pt)  med(pt) |   min(hr)   med(hr)
frmtmb FIX (active)       0.170    0.180 |    0.1669    0.1783
frmtmb SWAP (round 1)     0.170    0.185 |    0.1700    0.1797
frmtmb BASE               0.140    0.160 |    0.1459    0.1549
BASE + nlme + gen         0.150    0.170 |    0.1627    0.1680
control Matrix FIX        0.360    0.405 |    0.3595    0.3983
control Matrix BASE       0.350    0.385 |    0.3650    0.3824
control nothing FIX       0.000    0.000 |    0.0002    0.0002
control nothing BASE      0.000    0.000 |    0.0002    0.0002
lme4 (2.5c price)         0.610    0.750 |    0.6199    0.7529
nlme (2.5a price)         0.030    0.050 |    0.0386    0.0422
generics (2.5a price)     0.000    0.020 |    0.0077    0.0086

ACTIVE - BASE on minima, hi-res:  +0.0210 s  p = 0.0000
SWAP   - BASE, same:              +0.0241 s  p = 0.0000
ACTIVE - SWAP, same:              -0.0031 s  p = 0.1251
CONTROL Matrix, same:             -0.0056 s  p = 0.5043
CONTROL nothing, same:            -0.0000 s  p = 0.0293
ACTIVE - BASE on proc.time():      +0.0300 s (3.0 ticks)

ATTRIBUTION, hi-res minima:
  BASE + nlme + generics, minus BASE: +0.0168 s
  ACTIVE minus that arm:              +0.0042 s
  SWAP   minus that arm:              +0.0073 s
20000 permutations, seed 20260915
```

- **Import lme4.** 22 recursive dependencies, including Rcpp,
  RcppEigen, minqa, nloptr, boot, Rdpack and rbibutils, against 6 for
  nlme and 5 for generics, all base. And it fixes only half the
  problem.
- **The active binding.** Zero dependencies. frmtmb takes whichever
  generic the user would have reached without it: the search path
  first, then `c("brms", "lme4")`.

**Decision: the active binding.** It is strictly better than the
import, not merely cheaper, because it repairs `ngrps(merMod)` for an
lme4 user and `ngrps(brmsfit)` for a brms user where one hard import
could only ever have repaired one. `ngrps` is off the damage list in
both `FIX lme4 then frmtmb` and `FIX brms then frmtmb`.

### `refit`, which the plan's answer got half right

The plan said to import `refit` from `generics`. That is done, and **it
does not fix the lme4 collision**: lme4 DEFINES its own `refit`, so
`refit.merMod` is in lme4's table and generics' generic cannot see it.
So `refit` carries both routes, the static import as the floor and
`lme4` in the run-time table, and the frmtmb method is registered on
both.

## BLOCKER C, second half. The load cost was 2.4x too large

**R2, withdrawn and replaced.** Round 1 published +0.040 s from
`system.time()[["elapsed"]]`, which is `proc.time()`. That clock
advances on this box in steps whose minimum is 0.010 s, so 0.220
against 0.180 was one or two ticks, and the round-1 CONTROL moved a
full tick, half the claimed effect. **This is the third timing claim to
evaporate on that instrument in this project, and the rule is already
in `dev/lane-rules.md`.** The figure is withdrawn.

The block above is the replacement: both clocks timed on every
replicate, both ticks printed, 30 interleaved replicates one fresh
process each, all three builds as arms of ONE run, two controls, and a
permutation test on the paired minima. It was run with nothing else on
the machine, after both `R CMD check` runs had finished.

What it says, and what it does not:

* **The two new hard Imports cost +0.0168 s.** That is the BASE arm
  with `loadNamespace("nlme")` and `loadNamespace("generics")` added,
  minus BASE. An earlier run of this script, on a busier machine, gave
  +0.0158 s.
* **R3, the disagreement is settled.** Review re-measured with 80
  interleaved rounds, arm order shuffled within each round, a
  high-resolution clock and paired bootstrap intervals
  (`dev/genrev-r2-loadcost.R`): IMPORTS minus BASE +0.0213 s
  [+0.0164, +0.0251], ACTIVE minus SWAP -0.0008 s [-0.0074, +0.0032],
  null control -0.0000 s. It withdrew its round-1 figure. That +0.0167 s
  was SWAP minus BASE on the minima of 30 runs, and it falls below the
  5th percentile of the same statistic resampled: **a low draw, not the
  cost of the whole change.** The finding on this page stands as it
  was: the two Imports are the cost, the two designs cannot be told
  apart at load, and the whole-change headline is a range.
* **The whole change costs +0.0210 s on the active-binding build**, and
  +0.0241 s on the swap build in the same run. The residual after the
  Imports is +0.0042 s for ACTIVE and +0.0073 s for SWAP.
* **The two designs cannot be told apart at load.** `ACTIVE - SWAP` is
  -0.0031 s with p = 0.1251. Across this session's three runs the
  headline `ACTIVE - BASE` read +0.0280, +0.0229 and +0.0210 s, a
  spread of 7 ms that is larger than any difference between designs,
  so the page reports the Imports attribution as the finding and the
  headline as a range.
* **One control is not clean and it is reported rather than hidden.**
  `CONTROL nothing` gives p = 0.0293 at a magnitude of 0.0000 s. The
  arm is two `Sys.time()` reads, 0.0002 s, so its minima differ by a
  microsecond or two, and a permutation test on the minima of values
  that close ranks them by noise rather than by effect. The magnitude is the
  honest reading there; the p-value on that arm is not.
* `proc.time()` measured a minimum tick of 0.0100 s and a median of
  0.0200 s in this run, so the same difference on the old instrument
  reads as +0.0300 s, three ticks. That is the instrument that produced
  the withdrawn +0.040 s.

## 2.5d: the scale contract

`?frmtmb-scales` (`R/scales.R`) states, per method, which of three
scales it reports: link, response, or unitless. Measured with
`dev/generics-scale.R` on a 400-row lognormal fit, seed 2026.

```
== lognormal fit, seed 2026, n = 400, dev/generics-scale.R ==
sigma(fit)                      0.4014523891
predict(fit)[1]  (default)      8.7277667035
predict(type='link')[1]         8.7277667035
predict(type='response')[1]     6689.8632629972
fitted(fit)[1]                  6689.8632629972
exp(mu)[1]        (median)      6171.9289610356
exp(mu+s^2/2)[1]  (mean)        6689.8632629972

relative residuals at full precision (identities, not fits):
predict() vs predict(type='link')          0.000e+00
predict(type='response') vs fitted()       0.000e+00
fitted() vs exp(mu + sigma^2/2)            0.000e+00
fitted()/exp(mu), the plan's 1.079         1.083918
exp(sigma^2/2)                             1.083918

residuals:
residuals() vs y - fitted()                0.000e+00
residuals('pearson') vs resp/sd            8.857e-17

coefficient tables, link scale or response scale:
fixef()$mu['(Intercept)']       8.0658105190
fixef()$sigma['(Intercept)']    -0.9126663352
exp(that)                       0.4014523891
sigma(fit)                      0.4014523891

summary() coefficient blocks:  mu, sigma 
  mu Estimate[1] = 8.0658105190 
  sigma Estimate[1] = -0.9126663352 

VarCorr()[[1]] (variance of the mu linear predictor)  0.1042556532 
sqrt of it                                            0.3228864402 

simulate()[1:3] (response scale)  15459.63 1478.36 1750.83 
observed y[1:3]                   7975.90 2475.47 4096.99 

== poisson fit ==
predict()[1] (link)             1.4581725525
fitted()[1]  (response)         4.2980977990
exp(link)[1]                    4.2980977990
sigma(f2) (no dispersion param) 1.0000000000
```

The plan's identity is reproduced and IS an identity: `fitted()` on a
lognormal is `exp(mu + sigma^2/2)` by the family's own `mean_fn`, so
the residual is 0.000e+00 at full precision and `identical()` is TRUE.
The plan's 1.079 is the same identity at its own seed's sigma of 0.386.

### BLOCKER D. Three holes, all now measured and on the page

**R2.**

```
== the rows the contract was missing, dev/generics-scale2.R ==
same 400-row lognormal fit, seed 2026

-- conditional_effects(), the row that matters most --
effect plotted                  x
estimate__ range                1274.05 to 9689.08
predict() on the same fit       8.7278
fitted() range                  698.20 to 11504.04
so conditional_effects() is RESPONSE, the opposite of predict()

-- the identity, and the two conditions it is stated without --
constant sigma: identical(fitted, exp(mu+s^2/2))  TRUE
distributional sigma: sigma(fit) is              NA
  identity with sigma():   NA, the formula fails
  identity with predict(dpar = 'sigma'): identical TRUE
truncated y | trunc(lb = 2000): max relative error of
  exp(mu + sigma^2/2) against fitted()            0.1150
  fitted()[1] 6596.12 against the naive formula 6570.00

-- the four methods stated nowhere --
posterior_summary(fit)          ERR: is.atomic(x) is not TRUE
posterior_summary(matrix)       Estimate,Est.Error,Q2.5,Q97.5
pp_check(fit)                   returned
bayes_R2(fit)                   refuses: bayes_R2() is computed per posterior dra
hypothesis() returns frmtmb_hypothesis/data.frame with columns hypothesis,estimate,se,lwr,upr,z,p
hypothesis(fit, 'x = 0') Estimate 0.391238  vs fixef 0.391238
```

* `conditional_effects()` was absent from the table and it is the one
  method whose scale is the OPPOSITE of `predict()`'s default:
  `estimate__` runs 1274.05 to 9689.08 where `predict()` on the same
  fit is 8.7278. It is now the row after `fitted()`.
* The identity was stated flatly and has two conditions. With a
  distributional sigma, `sigma()` returns `NA` with a warning and the
  formula returns `NA` with it; `predict(dpar = "sigma", type =
  "response")` reproduces `fitted()` exactly. Under truncation
  `fitted()` is the TRUNCATED mean and the formula answers a different
  question.

  **R3, and round 2's number was mislabeled.** Round 2 called 0.1150
  "the relative error at the worst row". It is `max |diff| /
  max(fitted)`; the worst row is 0.5163. And it is not a measurement at
  all but an IDENTITY: with `a = (log(lb) - mu) / sigma`, the per-row
  relative shortfall is `1 - pnorm(-a) / pnorm(sigma - a)`, because the
  truncated lognormal mean is `exp(mu + sigma^2 / 2) * pnorm(sigma - a)
  / pnorm(-a)`. It has no single value; it depends on where the bound
  sits in each row's distribution. The page now gives the formula and
  no number, and the test asserts the identity row by row to machine
  precision instead of a threshold on its size. That assertion is seen
  failing under mutants M3 and M6, and M1 fails the same block on its
  untruncated half. (Review's earlier 0.4548 was on its own
  construction and is superseded by this.)

```
== lognormal under trunc(lb), dev/generics-trunc.R, seed 2026 ==
per-row relative error of exp(mu + s^2/2) against fitted():
  (fitted - naive) / fitted = 1 - Phi(-a) / Phi(s - a),
  a = (log(lb) - mu) / s

    lb  rows  old stat worst row    median        a range   identity
   500   400    0.0098    0.1429    0.0000  -7.60..-0.57   2.08e-16
  2000   308    0.1150    0.5163    0.0516  -3.76.. 1.41   2.50e-16
  5000   104    0.3650    0.8430    0.3513  -1.77.. 4.02   2.22e-16

`old stat` is max|diff|/max(fitted), the number the page used to
call the worst row. `identity` is the largest per-row difference
between the measured error and the closed form.
```
* `hypothesis()`, `posterior_summary()`, `pp_check()` and `bayes_R2()`
  are now rows. `hypothesis()`'s `estimate` is `fixef()`'s number to
  the digit, so it is the link scale per dpar.
* **The page is findable.** `@keywords internal` is off `R/scales.R`
  and `frmtmb-scales` is in `_pkgdown.yml`, where the other nine
  comparable concept pages already were.

### Against brms, on the response scale

```
== row 1, the four numbers ==
frmtmb predict()   default (link)   8.727767
frmtmb fitted()    (response)       6689.863
brms   fitted()    Estimate         6756.484
brms   predict()   Estimate         6677.817
brms predict() / exp(frmtmb link)   1.081966
exp(sigma^2/2), sigma = 0.401452       1.083918
== agreement on the response scale, n = 400 rows ==
max |frmtmb/brms - 1|                      0.010627
median relative posterior SD of brms epred 0.089755
max |difference| / its own posterior SD    0.1130
median |difference| / its own posterior SD 0.0584
== the DEFAULT of predict() is a different scale ==
max |frmtmb predict() / brms predict() - 1| 0.999207
ratio brms predict() to frmtmb predict()   441.3
== sigma: link scale in summary(), response in sigma() ==
frmtmb summary() sigma Estimate  -0.912666 (log scale)
frmtmb sigma()                   0.401452
brms   sigma posterior mean      0.403424
```

`fitted()` agrees with brms's `fitted()` to a maximum relative
disagreement of 0.010627 over 400 rows, which is 0.1130 of brms's own
posterior standard deviation for that row at the worst row and 0.0584
at the median. Those four figures repeat to the last digit across runs,
because both sides are deterministic given the cached fit.

The `predict()` ratio is NOT one of those. brms's `predict()`
summarizes fresh predictive draws, so it moved from 443.1 to 441.3 at
the median row between two runs on the SAME fit. The page now says
"roughly 440, and 765 at row 1" rather than four significant figures on
a Monte Carlo quantity.

### The tests, and seeing them fail

`tests/testthat/test-scale-contract.R`: 12 blocks, 36 assertions, 1.7 s,
asserted THROUGH THE GENERICS. Every block is paired with the inverse
assertion the wrong scale would satisfy. Where a relation is an
IDENTITY the yardstick is machine epsilon in ulps; where it is a
MEASUREMENT it is a standard error the run reports. No absolute number
is a tolerance.

Two things the page should not claim, and now does not:

* The replaced yardstick. The mean and the median of a lognormal sit
  **0.9134** of one prediction standard error apart on this design, so
  the prediction's own error cannot separate them. Kept on the page: it
  is a result about the instrument.
* **R2, N7.** The replacement, `(ratio - 1) / se_ratio > 5`, reduces to
  about `1 / (2 * se(log sigma))`, which is about `sqrt(n/2)` and does
  not depend on sigma. It is a sample-size bound, not a measurement of
  the separation. It still discriminates, which is the job: under M1
  the ratio is exactly 1 and the statistic is 0. The test comment says
  all of this.

```
== item 2.5d: the scale test seen FAILING, per mutant ==
dev/generics-mutants.R and dev/generics-mutants2.R build each
mutant from this worktree, install it into its own library, and
dev/generics-runtests.R runs the UNMODIFIED test file against it.

mut  what it breaks                                      pass  fail
M1   lognormal fitted() is the MEDIAN exp(mu)              31     5
       caught by: fitted() on a lognormal is the MEAN, not the median
       caught by: the lognormal identity holds only for a constant sigma
       caught by: the lognormal identity fails under truncation
M2   predict() defaults to the RESPONSE scale              30     6
       caught by: predict() reports the LINK scale and nothing else
       caught by: a log-link count fit puts predict() and fitted() an exp apart
       caught by: conditional_effects() is on the RESPONSE scale
M3   sigma() returns the LINK scale, log(sigma)            31     5
       caught by: fitted() on a lognormal is the MEAN, not the median
       caught by: summary() and fixef() report every dpar on its own LINK
       caught by: the lognormal identity fails under truncation
M4   residuals() default is formed on the LINK scale       33     3
       caught by: residuals() default is on the RESPONSE scale
       caught by: residuals(type = 'pearson') is unitless
M5   simulate() draws on the LINK scale                    34     2
       caught by: simulate() draws on the RESPONSE scale
M6   fitted() is the median, predict(response) the mean    28     8
       caught by: predict(type = 'response') and fitted() are the same call
       caught by: fitted() on a lognormal is the MEAN, not the median
       caught by: residuals() default is on the RESPONSE scale
       caught by: the lognormal identity holds only for a constant sigma
       caught by: the lognormal identity fails under truncation
M7   VarCorr() is on the RESPONSE scale                    34     2
       caught by: the LINK scale is the log scale, proved by the twin fit

distinct test blocks caught by at least one mutant: 12
```

Seven mutant builds, each a one-line change to this worktree installed
into its own library, with the UNMODIFIED test file run against it. All
12 blocks are caught by at least one mutant.

## The tests for 2.5a to 2.5c, and seeing them fail

`tests/testthat/test-generic-collision.R` launches CHILD R PROCESSES,
one per load order, because the search path is built once per session.
**R2**: the probe inside it no longer reads error strings either, and
four blocks were added for the orders review found missing.

| build | blocks | assertions | pass | fail |
|---|---|---|---|---|
| BASE `rellib-r3` | 11 | 43 | 24 | 19 |
| SWAP `genrev-lib`, round 1 | 11 | 43 | 36 | 7 |
| ACTIVE, this worktree | 11 | 43 | 43 | 0 |

The five added blocks are: an owner unloaded and reloaded, the
construction that broke the swap; detach and reattach in both packages;
that the binding IS active in the namespace and in the attached
environment and that `unlockBinding` appears nowhere in the namespace;
that an owner exporting a NON-generic under one of these names does not
take the binding; and that no generic in the table carries work in its
body.

That last one had to be rebuilt. The first spelling shadowed
`posterior` with a package exporting only `ndraws`, and frmtmb then
failed to LOAD, because R resolves `S3method(posterior::as_draws, ...)`
against the shadow during `loadNamespace`. That failure is
**pre-existing and identical on the base build**
(`dev/generics-shadowload.R`), so it is filed below rather than fixed,
and the guard is built on `bayesplot` instead, which frmtmb registers
exactly one delayed method on. Measured there: the shadow loads, its
`pp_check` is not a generic, and `pp_check` still resolves to frmtmb.

The owner-absent case is built by construction: a directory named for
each of the six owners is written into a temporary library, put FIRST
on the child's library path with an `R/` file that calls `stop()`, and
the test asserts `loadNamespace()` on each one actually fails before
checking that frmtmb still loads and works.

## What was run

The session driving this lane ended once, mid-round, with an
`R CMD check` of core half-built. Before anything else the resumed run
confirmed every R file this lane touched still parses
(`dev/generics-parsecheck.R`: 128 files, the one failure being the
review's own `dev/genrev-scale2.R`, which this lane never edited), that
the lane's library holds no hollow package directory, and that the
installed core is newer than its newest source edit. The check was
then re-run from scratch, alone on the machine.

**R3, punch round 2.** Round 2 changed core in comments, in two names
added to the internal export list with their `@aliases` and a
paragraph on `?frmtmb-sampling-api`, in one Suggests floor and in
prose, and changed `frmtmb.sample`'s code in one method and its tests.
As instructed, `R CMD check` was re-run on `frmtmb.sample`, whose code
changed, on a quiet machine with no other R process and without
`--no-manual`: `Status: OK`, tests OK. It was NOT re-run on core; the
core block below is the round-1 resumed run. Core's NAMESPACE did
change, by two exports, so the checks that change could break were run
directly on the installed package instead: `tools::undoc()`,
`tools::codoc()`, `tools::checkDocFiles()` and `tools::checkS3methods()`
all report nothing, on core and on `frmtmb.sample`. Both new Rd
changes were verified by rendering (`dev/generics-rdcheck.R`, and
`Rd2txt` on `?frmtmb-sampling-api`).

The core test files round 2 can reach, one process each:

```
== suite: punch round 2, the core test files this round can reach, one process each ==
test files with a result   6 of 6
test blocks                74
assertions                 287
pass                       287
FAIL                       0
ERROR                      0
skip                       0
warning                    0

no file reported a failure or an error
```

The full core suite below ran against the build with the round-1
`hypothesis()` fix in; round 2 changed no core code path it exercises
beyond the files just listed. The `frmtmb.sample` suite below is round
2's.

```
== suite: frmtmb core, active-binding build, 136 files, one process per file, NOT_CRAN and FRMTMB_BRMS_FIT_TESTS true ==
test files with a result   136 of 136
test blocks                1462
assertions                 10207
pass                       10207
FAIL                       0
ERROR                      0
skip                       1
warning                    2

no file reported a failure or an error
```

```
== suite: frmtmb.sample, 17 files, one process per file, punch round 2 ==
test files with a result   17 of 17
test blocks                195
assertions                 1145
pass                       1145
FAIL                       0
ERROR                      0
skip                       3
warning                    0

no file reported a failure or an error
```

```
== R CMD check --as-cran, built WITH vignettes, no --no-manual ==
dev/generics-cran.sh, mirroring dev/release/run-check.ps1

frmtmb         Status: 1 NOTE
   HTML version of manual ... [19s] NOTE
     Skipping checking math rendering: package 'V8' unavailable

frmtmb.sample  Status: OK

```

The first core check of the active-binding build was
`Status: 1 ERROR, 1 NOTE`: the `hypothesis()` regression described
under BLOCKER A, found here and nowhere else. After the fix, the block
above.

## Defects found and NOT fixed

### 1. `frmtmb.sample` has the same defect for the generics it defines

Core's fix now reaches **all 26** names `frmtmb.sample` re-exports from
core, in every order including the one round 1 did not run. It does not
reach the 28 `frmtmb.sample` defines itself.

```
== does core's fix reach frmtmb.sample? dev/generics-propagate.R ==
S = library(brms) then library(frmtmb.sample)
T = the reverse; U = frmtmb.sample attached, brms merely LOADED
BASE = the base commit, SWAP = this lane's round 1, ACTIVE = now
lost = the CLASS method is unreachable (method-table lookup)
--- BASE mode S
re-exported from core      lost 26 of 26
  of which answered by a .default  2  loo_compare,posterior_summary
defined by frmtmb.sample   lost 28 of 28
as_draws binding active in package:frmtmb.sample  FALSE
as_draws resolves to       frmtmb
--- BASE mode T
re-exported from core      lost 0 of 26
  of which answered by a .default  0  
defined by frmtmb.sample   lost 0 of 28
as_draws binding active in package:frmtmb.sample  FALSE
as_draws resolves to       posterior
--- BASE mode U
re-exported from core      lost 26 of 26
  of which answered by a .default  2  loo_compare,posterior_summary
defined by frmtmb.sample   lost 28 of 28
as_draws binding active in package:frmtmb.sample  FALSE
as_draws resolves to       frmtmb
--- SWAP mode S
re-exported from core      lost 0 of 26
  of which answered by a .default  0  
defined by frmtmb.sample   lost 28 of 28
as_draws binding active in package:frmtmb.sample  FALSE
as_draws resolves to       posterior
--- SWAP mode T
re-exported from core      lost 0 of 26
  of which answered by a .default  0  
defined by frmtmb.sample   lost 0 of 28
as_draws binding active in package:frmtmb.sample  FALSE
as_draws resolves to       posterior
--- SWAP mode U
re-exported from core      lost 23 of 26
  of which answered by a .default  2  loo_compare,posterior_summary
defined by frmtmb.sample   lost 28 of 28
as_draws binding active in package:frmtmb.sample  FALSE
as_draws resolves to       frmtmb
--- ACTIVE mode S
re-exported from core      lost 0 of 26
  of which answered by a .default  0  
defined by frmtmb.sample   lost 28 of 28
as_draws binding active in package:frmtmb.sample  TRUE
as_draws resolves to       posterior
--- ACTIVE mode T
re-exported from core      lost 0 of 26
  of which answered by a .default  0  
defined by frmtmb.sample   lost 0 of 28
as_draws binding active in package:frmtmb.sample  TRUE
as_draws resolves to       posterior
--- ACTIVE mode U
re-exported from core      lost 0 of 26
  of which answered by a .default  0  
defined by frmtmb.sample   lost 28 of 28
as_draws binding active in package:frmtmb.sample  TRUE
as_draws resolves to       posterior
```

```
== frmtmb.sample, dev/generics-audit4.R ==
exported S3 generics                       54
  re-exported from core, adopted with core 26
  defined by frmtmb.sample itself          28
  of those, a name another package owns    28

the names frmtmb.sample still owns that someone else owns:
  as.mcmc                brms, coda
  bayes_factor           bridgesampling, brms
  bridge_sampler         bridgesampling, brms
  kfold                  brms, loo
  log_lik                brms, rstantools
  log_posterior          bayesplot, brms
  loo_moment_match       brms, loo, rstan
  loo_subsample          brms, loo
  mcmc_plot              brms
  neff_ratio             bayesplot, brms
  nsamples               brms, rstantools
  nuts_params            bayesplot, brms
  parnames               bbmle, brms
  post_prob              bridgesampling, brms
  posterior_epred        brms, rstantools
  posterior_interval     brms, rstantools
  posterior_linpred      brms, rstantools
  posterior_predict      brms, rstantools
  posterior_samples      brms, gratia
  pp_mixture             brms
  predictive_error       brms, rstantools
  predictive_interval    brms, rstantools
  psis                   brms, loo
  reloo                  brms
  restructure            brms
  rhat                   bayesplot, brms, posterior
  stancode               brms
  standata               brms

how those re-exports resolved in THIS process, where no optional owner was attached:
   frmtmb=23 nlme=3 
```

Mode U, `library(frmtmb.sample)` with brms merely LOADED, is the order
review added and the one a user hits the moment any other attached
package Imports brms. BASE loses 26 of 26 there, the round-1 SWAP loses
23 of 26 (the three it saves are the static nlme names), and the active
binding loses 0 of 26, because the re-exported binding is itself
active. Every one of the 28 names `frmtmb.sample` defines belongs to
rstantools, loo, bridgesampling, bayesplot, posterior, coda or brms.
The same two routes would work there. Out of this lane's four items.

### 2. 40 names frmtmb and brms both export are NOT generics

Found by accident and it cost a Stan compile.
`dev/generics-scale-brms.R` called `brm(family = lognormal())` with
both attached and brms refused the object, because `lognormal()`
resolved to frmtmb's family constructor:

    Error in `validate_family()`: Argument 'family' is invalid.

A generic can be shared; a family constructor cannot, because the two
return genuinely different objects. R warns at attach time and
`brms::lognormal()` is the workaround, but a ported script will hit it
and the error names neither package. **This is a decision for the user,
not a lane**: accept it, document the prefix, or rename frmtmb's
constructors.

```
== names frmtmb and brms both export, dev/generics-audit3.R ==
brms exports                       306
frmtmb exports                     187
both                               67
  of which S3 generics             27
  of which already handed to owner 3
  of which NOT generics            40

the non-generic collisions:
   acat, asym_laplace, bernoulli, Beta, beta_binomial 
   bf, categorical, cox, cratio, cumulative 
   custom_family, exgaussian, exponential, geometric, get_prior 
   hurdle_gamma, hurdle_lognormal, hurdle_poisson, lf, lognormal 
   mixture, multinomial, mvbf, negbinomial, nlf 
   prior, prior_, prior_string, set_prior, set_rescor 
   shifted_lognormal, skew_normal, sratio, student, von_mises 
   weibull, zero_inflated_beta, zero_inflated_binomial, zero_inflated_negbinomial, zero_inflated_poisson 
```

### 3. An incomplete `posterior` stops frmtmb loading

Found while building the non-generic guard, PRE-EXISTING, identical on
the base build (`dev/generics-shadowload.R`). A `posterior` on the
library path that does not export `as_draws` makes frmtmb fail to load:

    Error: package or namespace load failed for 'frmtmb' in
    get(genname, envir = envir): object 'as_draws' not found

That is R's own `registerS3methods()` resolving
`S3method(posterior::as_draws, ...)`, which the base commit also
carries. In practice an older or partial `posterior` earlier on the
library path breaks frmtmb entirely rather than degrading.

**R3, narrowed: "identical on BASE" is true of THAT construction and
not in general.** Review built a `posterior` that exports every generic
frmtmb registers on except `as_draws_rvars`, loaded before frmtmb. BASE
loads, SWAP loads, and this lane's build FAILS. Reproduced here:

```
== a posterior missing ONLY as_draws_rvars, loaded before frmtmb ==
dev/generics-shadowload.R; the shadow exports the other ten names
--- BASE
SHADOW_LOADS: TRUE 
FRMTMB: loaded 
--- SWAP
SHADOW_LOADS: TRUE 
FRMTMB: loaded 
--- ACTIVE
SHADOW_LOADS: TRUE 
FRMTMB: FAILED: package or namespace load failed for 'frmtmb' in get(genname, envir =  
--- CRAN posterior on this machine exports as_draws_rvars:
1.7.0 TRUE 
```

The cause is not the active binding. It is BLOCKER B's new
`S3method(posterior::as_draws_rvars, frmtmb_fit)` refusal, which is a
name no earlier NAMESPACE asked `posterior` for; `as_draws_list` is new
the same way. So this lane DID widen the exposure, by two names. The
realistic exposure is small: `posterior` 1.0.0, its first CRAN release,
already exports `as_draws_rvars` (the v1.0.0 tag's NAMESPACE, fetched
from the stan-dev repository), and the 1.7.0 installed here does. What
breaks is a partial or pre-CRAN install. `posterior (>= 1.0.0)` is now
declared in core's Suggests. That states the requirement and lets
`R CMD check` and installers see it; it does NOT stop a partial
`posterior` already on a library path from breaking `loadNamespace`,
which only removing the directives could do. The base-build failure
for a `posterior` lacking `as_draws` itself is unchanged and still not
this lane's.

### 4. `posterior_summary()` on a `frmtmb_fit` dies inside `var()`

`is.atomic(x) is not TRUE`, identically on both builds, because it
falls to `posterior_summary.default`, which calls `as.matrix()` on the
fit. It is the brms spelling a porter tries first. Filed rather than
fixed: unlike the `as_draws*` family, this class is currently ANSWERED
by a default rather than unmatched, so the fix is a decision about what
the right answer is, refuse or summarize the coefficient table, not a
missing method.

### 5. There is no new `R CMD check` NOTE

Round 1 introduced one, `unlockBinding`, and put a choice to the user
between keeping it and giving back a load order. **That choice was
false and is withdrawn**: the active binding removes the NOTE and
repairs more than the swap did. Confirmed on the final build:
`checking R code for possible problems ... OK`, and core's one NOTE is
the pre-existing V8 math-rendering one.

## Version

This needs a bump in core and a matching floor in `frmtmb.sample`, and
I am not choosing the numbers. `NEWS.md` carries the bullet, under
`# frmtmb (development version)`.

- Core gains two Imports (`generics`, `nlme`) and no longer DEFINES
  `fixef`, `ranef`, `VarCorr` or `refit`.
- One caller-visible break: `posterior_summary(object = ...)`.
  `nvariables` and `VarCorr` only gain arguments. `refit`'s formals now
  depend on whether lme4 is loaded.
- A downstream package that registers a method on any of these generics
  must match the OWNER's formals. `frmtmb.sample` is edited here for
  that reason, so its `frmtmb (>= 0.55.1)` floor has to move to the
  bumped core version.
- **R3.** The floor is now a hard requirement, not only a formals one:
  `frmtmb.sample`'s draws method calls `hyp_shadow_arm()` and
  `hyp_shadow_disarm()`, which no released core exports. Against the
  current release it would fail at the first `hypothesis()` on draws.
- **R3.** Core's Suggests now declares `posterior (>= 1.0.0)`; see
  defect 3. That is a dependency floor, not a version choice for either
  package.
- Nothing in the active binding needs a floor of its own.

## A mistake of mine that is worth the space

Two of the round-1 edit scripts wrote source files with `writeLines()`,
which opens a TEXT connection, and on Windows a text connection turns
every `\n` into `\r\n`. Eight LF source files came out CRLF in a tree
that is LF. Nobody noticed until an edit failed to match its own
pattern because of the stray `\r`.

**The first repair was much worse than the defect.** It walked every
tracked file and stripped any CR that was followed by an LF. That
condition is true of BINARY files too, and it rewrote 123 `.rds`,
`.png`, `.ttf` and `.woff2` files with those bytes removed, corrupting
all of them, while reporting "converted 424 tracked files back to LF"
as though that were a result. They were restored with
`git restore --worktree` and re-read to confirm, `readRDS()` on three
of them returning a list, and the narrowed version in
`dev/generics-eol.R` takes an explicit list of the lane's own files and
refuses anything containing a NUL byte or a high byte.

This is the project's own rule about guards, in a place nobody thought
to look: the condition has to be true ONLY of the thing being guarded,
and "construct the case where the guarded thing is absent" means
running it against a file it must not touch. `dev/generics-sub.R` now
holds the one edit helper, and it writes with `writeBin`.

## Scripts, all prefixed `generics-`

| script | what it produced |
|---|---|
| `dev/generics-repro.R` | seeded; the original 15-of-17 reproduction |
| `dev/generics-audit.R` | seeded; who exports each of the 30 names |
| `dev/generics-audit2.R` | who DEFINES each, where each brms method is registered, and every owner's formals |
| `dev/generics-audit3.R` | the non-generic collisions with brms |
| `dev/generics-audit4.R` | what is left in `frmtmb.sample` |
| `dev/generics-check.R` | the per-load-order probe, by method-table lookup |
| `dev/generics-summary.R` | the damage table, read at fixed widths |
| `dev/generics-propagate.R` | whether core's fix reaches `frmtmb.sample`, three orders |
| `dev/generics-whichbuild.R` | which mechanism a given library carries |
| `dev/generics-shadowload.R` | the incomplete-owner and non-generic-owner shadows |
| `dev/generics-loadcost.R` | `library()` cost on two clocks, interleaved, with a permutation test |
| `dev/generics-adoptcost.R`, `dev/generics-adoptparts.R` | the active binding's per-access cost, and where it goes |
| `dev/generics-scale.R`, `dev/generics-scale2.R` | the scale contract and the rows it was missing |
| `dev/generics-scale-brms.R` | the same against a brms fit, cached in `dev/stan-cache` |
| `dev/generics-mutants.R`, `dev/generics-mutants2.R` | the seven wrong-scale builds |
| `dev/generics-mutsummary.R` | the mutant matrix |
| `dev/generics-runtests.R`, `dev/generics-runtests-sample.R` | one test file per process, errors counted too |
| `dev/generics-suite.sh`, `dev/generics-suitesummary.R` | the suite, three at a time, and its counts |
| `dev/generics-cran.sh`, `dev/generics-checksummary.R` | `R CMD check --as-cran` and its status lines |
| `dev/generics-build.R` | roxygenise and install into the lane's library |
| `dev/generics-sub.R`, `dev/generics-edit*.R` | the one edit helper and every source edit |
| `dev/generics-eol.R` | the line-ending repair, and the story above |
| `dev/generics-assemble.R` | splices the generated blocks into this file |
| `dev/generics-bodies.R`, `dev/generics-notice.R` | the bare-generic guard, and the `hypothesis()` note on three builds |
| `dev/generics-parsecheck.R` | after the interrupted session: every file parses, no hollow package |
| `dev/generics-rdcheck.R` | the scale page verified by rendering it |
| `dev/generics-trunc.R` | **R3** the truncation identity, swept over three bounds |
| `dev/generics-edit39.py` to `dev/generics-edit43.py` | **R3** BLOCKER E's test, its fix, the guard, the nits, this page |
| `dev/generics-edit*.py` | the round-2 text edits made while an R CMD check held the machine |
