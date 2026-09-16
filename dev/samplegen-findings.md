# Lane `samplegen`: frmtmb.sample must not break the owners of its names

Worktree `frmtmb-wt-samplegen`, branch `wt-samplegen`, based on
`1eff19e`, the merged Phase 2.5 commit. Nothing is committed.

Private library `C:/Users/adf44/source/r/samplegen-lib`, holding this
worktree's core (DESCRIPTION still says 0.55.2) and frmtmb.sample
(0.4.1). The base, read-only, is `C:/Users/adf44/source/r/rellib-r3`:
core 0.56.0 and frmtmb.sample 0.4.2. R 4.6.1, brms 2.23.0, posterior
1.7.0, loo 2.10.1, rstantools 2.7.1, bayesplot 1.16.0, bridgesampling
1.2.1, coda 0.19.4.1, gratia 0.11.2.

Every count below is emitted by a script into a marked block and pasted
by `dev/samplegen-assemble.R`; none is typed. `R CMD check` was NOT run:
it was held for the release. See "What was run".

## The defect, restated from the measurement

This is defect 1 of `dev/generics-findings.md`. frmtmb.sample DEFINES
28 S3 generics, and every one of those names is owned by another
package. `UseMethod()` reads the method table of the namespace where
the generic it reached was defined, so frmtmb.sample's generic, reached
first, sends a `brmsfit` into frmtmb.sample's table. The
`S3method(owner::generic, frmtmb_draws)` directives the package already
carried for every one of the 28 put its methods in the owners' tables.
They did nothing about its rival generics, which is core's decisive
finding again.

The measurement found two failures the plan did not name, both silent:
with posterior loaded AFTER frmtmb.sample, `rhat()` on draws fell into
posterior's `rhat.default`, and with gratia loaded after it,
`posterior_samples()` fell into gratia's. Both were missing twin
directives, not rival generics, and both sit in the reverse orders that
the generics lane's probe, which looked only at `brmsfit`, reported as
0 of 28.

## The ownership audit

`dev/samplegen-audit.R`. Step one is `parseNamespaceFile()` over every
installed package, never grep, for the candidates that EXPORT a name.
Step two loads each candidate and reads the environment of its exported
function, because exporting a name is not defining it: brms exports all
28 names and defines 8. The owner is the package that DEFINES the
generic, because that is where brms's own methods are registered.

```
== ownership audit, dev/samplegen-audit.R ==
installed packages scanned        407
candidate exporters                bayesplot, bbmle, bridgesampling, brms, coda, gratia, loo, posterior, rstan, rstantools 

frmtmb.sample exports 54 S3 generics, defines 28 itself

### per name: definers (generic formals), importers, brmsfit method

as.mcmc   frmtmb.sample(x, ...)
  defines:   coda(x, ...) 
  imports:   brms<-coda 
  brmsfit method in table of: coda 

bayes_factor   frmtmb.sample(x, ...)
  defines:   bridgesampling(x1, x2, log, ...) 
  imports:   brms<-bridgesampling 
  brmsfit method in table of: bridgesampling 

bridge_sampler   frmtmb.sample(x, ...)
  defines:   bridgesampling(samples, ...) 
  imports:   brms<-bridgesampling 
  brmsfit method in table of: bridgesampling 

kfold   frmtmb.sample(x, ...)
  defines:   loo(x, ...) 
  imports:   brms<-loo 
  brmsfit method in table of: loo 

log_lik   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

log_posterior   frmtmb.sample(object, ...)
  defines:   bayesplot(object, ...) 
  imports:   brms<-bayesplot 
  brmsfit method in table of: bayesplot 

loo_moment_match   frmtmb.sample(x, ...)
  defines:   loo(x, ...) 
  imports:   brms<-loo  rstan<-loo 
  brmsfit method in table of: loo 

loo_subsample   frmtmb.sample(x, ...)
  defines:   loo(x, ...) 
  imports:   brms<-loo 
  brmsfit method in table of: loo 

mcmc_plot   frmtmb.sample(object, ...)
  defines:   brms(object, ...) 
  brmsfit method in table of: brms 

neff_ratio   frmtmb.sample(object, ...)
  defines:   bayesplot(object, ...) 
  imports:   brms<-bayesplot 
  brmsfit method in table of: bayesplot 

nsamples   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

nuts_params   frmtmb.sample(object, ...)
  defines:   bayesplot(object, ...) 
  imports:   brms<-bayesplot 
  brmsfit method in table of: bayesplot 

parnames   frmtmb.sample(x, ...)
  defines:   brms(x, ...) 
  other:     bbmle(not a generic) 
  brmsfit method in table of: brms 

post_prob   frmtmb.sample(x, ...)
  defines:   bridgesampling(x, ..., prior_prob, model_names) 
  imports:   brms<-bridgesampling 
  brmsfit method in table of: bridgesampling 

posterior_epred   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

posterior_interval   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

posterior_linpred   frmtmb.sample(object, ...)
  defines:   rstantools(object, transform, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

posterior_predict   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

posterior_samples   frmtmb.sample(x, ...)
  defines:   brms(x, pars, ...)  gratia(model, ...) 
  brmsfit method in table of: brms 

pp_mixture   frmtmb.sample(x, ...)
  defines:   brms(x, ...) 
  brmsfit method in table of: brms 

predictive_error   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

predictive_interval   frmtmb.sample(object, ...)
  defines:   rstantools(object, ...) 
  imports:   brms<-rstantools 
  brmsfit method in table of: rstantools 

psis   frmtmb.sample(log_ratios, ...)
  defines:   loo(log_ratios, ...) 
  imports:   brms<-loo 
  brmsfit method in table of: loo 

reloo   frmtmb.sample(x, ...)
  defines:   brms(x, ...) 
  brmsfit method in table of: brms 

restructure   frmtmb.sample(x, ...)
  defines:   brms(x, ...) 
  brmsfit method in table of: brms 

rhat   frmtmb.sample(object, ...)
  defines:   bayesplot(object, ...)  posterior(x, ...) 
  imports:   brms<-posterior 
  brmsfit method in table of: posterior 

stancode   frmtmb.sample(object, ...)
  defines:   brms(object, ...) 
  brmsfit method in table of: brms 

standata   frmtmb.sample(object, ...)
  defines:   brms(object, ...) 
  brmsfit method in table of: brms 

### frmtmb.sample's own S3method directives on these names
  as.mcmc              (own) coda
  bayes_factor         (own) bridgesampling
  bridge_sampler       (own) bridgesampling
  kfold                (own) loo
  log_lik              (own) rstantools
  log_posterior        bayesplot (own)
  loo_moment_match     loo (own)
  loo_subsample        loo (own)
  mcmc_plot            brms (own)
  neff_ratio           bayesplot (own)
  nsamples             (own) rstantools
  nuts_params          bayesplot (own)
  parnames             brms (own)
  post_prob            bridgesampling (own)
  posterior_epred      (own) rstantools
  posterior_interval   (own) rstantools
  posterior_linpred    (own) rstantools
  posterior_predict    (own) rstantools
  posterior_samples    brms (own)
  pp_mixture           brms (own)
  predictive_error     (own) rstantools
  predictive_interval  (own) rstantools
  psis                 loo (own)
  reloo                brms (own)
  restructure          brms (own)
  rhat                 bayesplot (own)
  stancode             brms (own)
  standata             brms (own)

DONE
```

What the audit decided:

* 26 names have exactly one defining owner. brms imports 20 of the 28,
  so for those brms's method is in that owner's table.
* `rhat` has two rival definers: posterior `(x, ...)` and bayesplot
  `(object, ...)`. brms imports posterior's, and `rhat.brmsfit` is in
  posterior's table. Both are owners, posterior first.
* `posterior_samples` has two rival definers: brms `(x, pars = NA, ...)`
  and gratia `(model, ...)`. gratia is an owner the plan did not name.
  Both are owners, brms first.
* bbmle exports a `parnames` that is not a generic. It is not an owner.
* `hypothesis` was core's one generic with work in its body. None of
  these 28 has any: all 28 are a bare `UseMethod()` on the base build
  (the owner-absent test below reads them there and passes), so no work
  had to move and there is no "where did the work go" guard to write.
  The guard that stops work coming BACK is in, and was seen failing on
  a mutant.
* No other extension in this monorepo defines an S3 generic. That is a
  source search for `UseMethod(` under `extensions/*/R`, used as a
  screen for definitions, not as an ownership audit.

## The design: core's mechanism, called with a second table

**Reuse worked, and core needed one change.** `frm_install_generics()`
was already parameterized by `pkgname` but read the hard-coded
`frm_generic_owners` table. It now takes `owners = frm_generic_owners`
as an argument, refuses a malformed table, and is exported on
`?frmtmb-sampling-api` with a new section, "The borrowed-generic seam",
that states the three rules an extension takes on with it.
`frm_bind_generic()`, `frm_adopt_target()` and `frm_is_generic()` are
unchanged and stay internal: the extension never calls them.

frmtmb.sample gains `R/generic-owners.R`, which holds only the table
`sample_generic_owners` and the reason for it, and its `.onLoad()` calls
`frm_install_generics(pkgname, sample_generic_owners)` as its first
line. The bindings land in frmtmb.sample's namespace, not core's:
`makeActiveBinding(gen, fun, ns)` receives frmtmb.sample's namespace,
which is unsealed during its own `.onLoad()`. The test "the bindings
are ACTIVE, in this namespace and not frmtmb's" asserts all 28 are
active in `asNamespace("frmtmb.sample")`, active in the attached
environment, and absent from core's namespace.

What reuse bought beyond not copying: the stale-owner repair, the
non-generic refusal and the multi-owner resolution came with it, and
each is asserted on frmtmb.sample by its own test. Active-ness also
propagates through `importIntoEnv()` into a third package that imports
frmtmb.sample, which is the importer row below.

The other changes, each one required by the mechanism:

* Two twin directives, `S3method(posterior::rhat, frmtmb_draws)` and
  `S3method(gratia::posterior_samples, frmtmb_draws)`. Written with
  `@rawNamespace`, because roxygen2 8.1.0 refuses two `@exportS3Method`
  tags in one block ("Block must contain only one @exportS3Method");
  the package already does this for `lme4::ngrps`.
* Formals, below.
* `gratia (>= 0.9.0)` in Suggests, the first gratia with
  `posterior_samples()` (gratia's NEWS, fetched from its repository),
  and `posterior (>= 1.0.0)`, whose NAMESPACE at the v1.0.0 tag
  exports `rhat` (fetched from stan-dev/posterior).

## The formals

The standard is the first owner in the table, which is brms's choice.
`dev/samplegen-formals.R` prints every owner's generic beside brms's
method (`dev/samplegen-out/formals.txt`); the test "each generic's
formals are its first owner's" compares deparsed formals and was seen
failing on base for exactly these six.

| name | base | now | owner |
|---|---|---|---|
| `bridge_sampler` | `(x, ...)` | `(samples, ...)` | bridgesampling |
| `bayes_factor` | `(x, ...)` | `(x1, x2, log = FALSE, ...)` | bridgesampling |
| `rhat` | `(object, ...)` | `(x, ...)` | posterior, not bayesplot |
| `post_prob` | `(x, ...)` | `(x, ..., prior_prob = NULL, model_names = NULL)` | bridgesampling |
| `posterior_samples` | `(x, ...)` | `(x, pars = NA, ...)` | brms, not gratia |
| `posterior_linpred` | `(object, ...)` | `(object, transform = FALSE, ...)` | rstantools |

Three break a caller that NAMED the first argument: `bridge_sampler(x
= )`, `bayes_factor(x = )` and `rhat(object = )`. The other three only
gain arguments; `posterior_linpred.frmtmb_draws()` already had
`transform`. Each method takes the same formals as its generic. NEWS
says all of this.

## The collision table, before and after, in every load order

`owners` is measured against a CONTROL computed in the same process:
the generic the user would reach with `package:frmtmb.sample` removed
from the search path. It counts only methods that control reaches,
because posterior's and bayesplot's rival `rhat` generics hide each
other's methods with no help from this package. It reads `0/0` where no
owner is attached, because nothing is then reachable by name from the
global environment in either arm; the brmsfit column carries those
rows.

```
== frmtmb.sample collision table, dev/samplegen-check.R ==
one process per row; lookup in the method table of environment(generic)
brmsfit: brms's class method unreachable (rows where brms is loaded)
draws:   the method that would run for a frmtmb_draws is not this
         package's own (a foreign .default counts as a failure)
owners:  owner class methods the CONTROL reaches and the test does not

build mode order                                     brmsfit   draws   owners
BASE  S    brms, then frmtmb.sample                  28/28     0/28    77/77
  resolves: frmtmb.sample=28 | active in namespace: FALSE 
BASE  T    frmtmb.sample, then brms                  0/28      1/28    0/77
  draws to a foreign .default: rhat<-posterior.default
  resolves: bayesplot=3 bridgesampling=3 brms=8 coda=1 loo=4 posterior=1 rstantools=8 | active in namespace: FALSE 
BASE  U    frmtmb.sample; brms only loaded           28/28     0/28    0/0
BASE  I    brms; an importer loaded, probed inside   28/28     0/28    n/a
  resolves: frmtmb.sample=28 | active in namespace: FALSE | in importer: FALSE
BASE  D    brms, sample; detach and reattach both    0/28      1/28    0/77
BASE  R    loo, sample; loo unloaded; brms loaded    28/28     0/28    0/0
  psis resolves to frmtmb.sample  stale: TRUE
  resolves: frmtmb.sample=28 | active in namespace: FALSE 
BASE  N    frmtmb.sample alone                       n/a       0/28    0/0
BASE  P    7 other owners, then frmtmb.sample        n/a       0/28    41/41
BASE  Q    frmtmb.sample, then 7 other owners        n/a       1/28    0/41
  draws to a foreign .default: posterior_samples<-gratia.default
  resolves: bayesplot=4 bridgesampling=3 coda=1 frmtmb.sample=7 gratia=1 loo=4 rstantools=8 | active in namespace: FALSE 
BASE  G    gratia, then frmtmb.sample                n/a       0/28    3/3
  posterior_samples.gam reachable: FALSE
  resolves: frmtmb.sample=28 | active in namespace: FALSE 
FIX   S    brms, then frmtmb.sample                  0/28      0/28    0/77
  resolves: bayesplot=3 bridgesampling=3 brms=8 coda=1 loo=4 posterior=1 rstantools=8 | active in namespace: TRUE 
FIX   T    frmtmb.sample, then brms                  0/28      0/28    0/77
  resolves: bayesplot=3 bridgesampling=3 brms=8 coda=1 loo=4 posterior=1 rstantools=8 | active in namespace: TRUE 
FIX   U    frmtmb.sample; brms only loaded           0/28      0/28    0/0
FIX   I    brms; an importer loaded, probed inside   0/28      0/28    n/a
  resolves: bayesplot=3 bridgesampling=3 brms=8 coda=1 loo=4 posterior=1 rstantools=8 | active in namespace: TRUE | in importer: TRUE
FIX   D    brms, sample; detach and reattach both    0/28      0/28    0/77
FIX   R    loo, sample; loo unloaded; brms loaded    0/28      0/28    0/0
  psis resolves to loo  stale: FALSE
  resolves: bayesplot=3 bridgesampling=3 brms=8 coda=1 loo=4 posterior=1 rstantools=8 | active in namespace: TRUE 
FIX   N    frmtmb.sample alone                       n/a       0/28    0/0
FIX   P    7 other owners, then frmtmb.sample        n/a       0/28    0/41
FIX   Q    frmtmb.sample, then 7 other owners        n/a       0/28    0/41
  resolves: bayesplot=4 bridgesampling=3 coda=1 frmtmb.sample=7 gratia=1 loo=4 rstantools=8 | active in namespace: TRUE 
FIX   G    gratia, then frmtmb.sample                n/a       0/28    0/3
  posterior_samples.gam reachable: TRUE
  resolves: coda=1 frmtmb.sample=26 gratia=1 | active in namespace: TRUE 
```

### What the `owners` denominators count, and `getS3method()`

**R1, from the review.** My `owners` column counts every entry in every
owner's method table, so it reads 77 in mode S, 41 in P and 3 in G. The
review counts (name, class) pairs deduplicated across owners with
`default` excluded, so the same cells read 39, 29 and 2
(`dev/reviews/20260915-samplegen.md`, section 1). Both are internally
consistent, both are 100 percent lost on BASE and 0 on FIX, and the
difference is in what is counted, not in what was found.

**`getS3method()` is not a column in my table and should be**, as it is
in core's. The review measured it over 14 load orders rather than my
10, one process per cell, with a control process per order:

```
== dev/sgrev-probe.R, the review's table, getS3method() column ==
BASE  S  brms then frmtmb.sample                     67/67 lost
BASE  D  detach and reattach both                    67/67
BASE  CA, CA2  core attached as well, either order   67/67
BASE  T, R, RN, RN2, BR                              0/67
BASE  P  7 other owners then frmtmb.sample           28/29
BASE  G  gratia then frmtmb.sample                   2/2
FIX   every one of the 14 orders                     0 lost
```

`RN`, `RN2` and `BR` are the three orders I did not run: frmtmb.sample
reached by `requireNamespace()` only, in both orders, and brms unloaded
and reattached under a live frmtmb.sample. `BR` is the one that could
have broken the memo, and it does not, because the memo compares
namespace environments by identity.

## brms held to output

Reachability says the right function is dispatched. This says the
result is the same: 22 calls on a cached brmsfit, value and warnings
hashed, in a brms-only control and in each order. The six calls that
would refit or recompile Stan (`kfold`, `reloo`, `loo_moment_match`,
`bridge_sampler`, `bayes_factor`, `post_prob`) are covered by the
method-table probe only.

```
== brms output with and without frmtmb.sample, dev/samplegen-brmsout.R ==
cached 400-row lognormal brmsfit, seed 20260915 before every call
control: library(brms) alone, fix-none.txt
calls in the control: 22 (1 of them brms's own error)

build  arm                                    calls differing from control
base   library(brms); library(frmtmb.sample)  22 of 22
base   library(frmtmb.sample); brms loaded    22 of 22
fix    library(brms); library(frmtmb.sample)  0 of 22
fix    library(frmtmb.sample); library(brms)  0 of 22
fix    library(frmtmb.sample); brms loaded    0 of 22

base-S, the first three differing calls:
  as.mcmc ERROR   no applicable method for 'as.mcmc' applied to an object of c
  log_lik ERROR   no applicable method for 'log_lik' applied to an object of c
  log_posterior ERROR   no applicable method for 'log_posterior' applied to an objec
```

`0 of 22` in S and T also shows the hashes are deterministic across
processes, seeded `posterior_predict()` included, so a zero here is a
measured zero.

Output brms adds when it is loaded: brms's own `parnames()` and
`posterior_samples()` GENERICS warn that the name is deprecated. On
draws that warning now precedes this package's refusal whenever brms is
loaded, because brms's generic is the one that runs. That is the owner's
behavior and NEWS says so.

### The other side of the same standard: output on a frmtmb_draws

**R1.** I measured only the brmsfit side. The fix changes which generic
runs for a `frmtmb_draws` too, and a generic can warn or coerce before
it dispatches, so the standard governs this side as well. The review
ran 46 calls on one cached draws object, 4 chains of 500, on four
channels (value, warnings, messages, printed form), its reference arm
run twice agreeing on 0 of 46:

```
== dev/sgrev-drawsout.R, 46 calls on a frmtmb_draws ==
BASE-S   0 of 46
BASE-T   3 of 46   rhat ERROR "unimplemented type 'list' in 'greater'"
                   parnames, posterior_samples: brms's deprecation warning
BASE-Q   2 of 46   rhat, the same error
                   posterior_samples: gratia's error, not this package's
                     refusal: "Don't know how to sample from the
                     posterior of <frmtmb_draws>FALSE"
FIX-N    0 of 46
FIX-S    2 of 46   the two brms deprecation warnings
FIX-T    2 of 46   the same two
FIX-P    0 of 46
FIX-Q    0 of 46
```

Two things this adds. The fix changes nothing else on the draws side:
every difference from BASE is one of the two deprecation warnings NEWS
declares. And my two "silent" BASE failures have a user-visible form,
which a method-table probe cannot show: with posterior in front
`rhat(ds)` errors "unimplemented type 'list' in 'greater'", and with
gratia in front `posterior_samples(ds)` reports gratia's "Don't know
how to sample from the posterior of <frmtmb_draws>FALSE" instead of
this package's refusal. Silent means the dispatch says nothing about
being wrong, not that the call is quiet.

## Guards, and seeing them fail

`extensions/frmtmb.sample/tests/testthat/test-generic-collision.R`, 13
blocks, a child R process per load order as in core. Every probe
compares the method reached with `identical()` to the method that
should run, never an error string, and a foreign `.default` counts as
lost.

```
== the collision tests, seen failing, dev/samplegen-runtests*.R ==
build                                        blocks assert  pass  fail error
-- extensions/frmtmb.sample/tests/testthat/test-generic-collision.R
BASE core + BASE sample (rellib-r3)              13     58    41    17     0
    failing: the package installs the owner table this file tests 
    failing: attaching frmtmb.sample after brms loses no brms method 
    failing: attaching frmtmb.sample before brms loses no method either way 
    failing: brms loaded but never attached loses no brms method 
    failing: a package importing frmtmb.sample does not lose brms methods 
    failing: the repair survives detach and reattach in both packages 
    failing: an owner unloaded and reloaded does not leave a stale generic 
    failing: every other owner's methods survive, in both orders 
    failing: the bindings are ACTIVE, in this namespace and not frmtmb's 
    failing: each generic's formals are its first owner's 
    failing: every method on an owned name is in each owner's table 
FIX, stancode() generic given work (mutant)      13     58    57     1     0
    failing: with no owner usable, it loads and its generics are bare 
FIX core + FIX sample                            13     58    58     0     0
-- tests/testthat/test-generic-collision.R (core)
BASE core (rellib-r3)                            13     47    47     0     1
    failing: frm_install_generics() takes a table and refuses a bad one 
FIX core                                         13     56    56     0     0
```

Reading the rows:

* On base, 11 of 13 blocks fail. The two that pass on base are guards
  with nothing to catch there. "with no owner usable, it loads and its
  generics are bare" is the body guard; the mutant row installs a copy
  of the fixed package whose `stancode()` generic calls `message()`
  before dispatching, and the guard names exactly `stancode`
  (`dev/samplegen-mutant.R`). "an owner exporting a NON-generic does
  not take the name" passes on base vacuously, because base has no
  binding to take; the refusal it relies on is core's
  `frm_is_generic()`, which core's own test covers.
* Every guard carries its positive condition as an assertion: the eight
  owner shadows must fail to load, the brms control must find 28
  methods, the owner probe must find more than 20, the twin check must
  read 28 rows and must name a twin it drops.
* The core test's base failure is the weak form, an unused argument,
  because base `frm_install_generics()` has no `owners`. The behavioral
  failures are on the frmtmb.sample side.
* **R1. The behavioral failure for the core block was constructed by
  the review** (`dev/sgrev-coremutant.R`): a core that TAKES `owners`
  and does not validate it, installed into its own library, run against
  the unmodified core test file. FIX core 13 blocks, 56 pass, 0 fail;
  the mutant 47 pass, 5 fail and 1 error, the five being the malformed
  tables the block expects to be refused and the error `invalid first
  argument` from `exists()` on the unvalidated table. So the block
  tests the refusal and not the signature, and the weak form was a page
  defect rather than a test defect.
* **R1.** The review also built a second body-guard mutant, on a
  different generic in a different file (`log_lik` in `R/loo.R`,
  `dev/sgrev-samplemutant.R`): 57 pass, 1 fail, and the guard names
  `log_lik`. A guard that only catches the generic its own author
  mutated is not a guard, so that run matters more than mine.

### A trap the mechanism sets for tests: mocking a bound generic

The first per-file run of `test-loo.R` on the fixed build ERRORED in
"a group-unit matrix says so, because loo() cannot". The block mocked
`log_lik` with `testthat::local_mocked_bindings()`. An assignment to an
active binding CALLS the binding's function with the value, and the
binding takes no argument, so the mock failed with "unused argument".

The block now mocks the METHOD, `log_lik.frmtmb_draws`, and passes a
draws-classed object. `UseMethod()` looks for a method from the calling
environment before the method table, and the caller is this namespace,
so the mock is what runs whichever generic dispatched. The edited file
passes 79 of 79 on base and on the fix; on base the mock is equally
the thing that answers, since the empty draws object would otherwise
error. Core's suite has no test that mocks one of its 25 bound names. A
future test in either package that mocks any of the 53 bound names will
hit this, and the comment at the mock says why.

**R1, the half my comment was missing.** `local_mocked_bindings()`
patches the namespace binding AND this package's own method table, but
NOT the owner's table (`dev/sgrev-mock2.R`). So with the owner loaded,
which is the normal session, the mock is visible only to a caller
INSIDE frmtmb.sample's namespace; a caller at the prompt reaches the
owner's generic, whose table still holds the real method, and the
review's probe gets "attempt to set an attribute on NULL" there. The
test is sound because `loo_matrix()` is such an inside caller, and the
comment in `test-loo.R` now says all of it, because the next author
will read that comment as a recipe.

This is a RULE and not a guard, and it is now in `dev/lane-rules.md`. A
guard would have to grep test sources for a mock of one of 53 names and
would fire on a correct test that mocks a same-named local, while the
real failure is loud, names "unused argument" and is one test run
away.

## What it costs

```
== load cost, dev/samplegen-loadcost.R ==
Sys.time() tick here: min 0.0000010 s, median 0.0000019 s

25 rounds, 100 fresh processes, arms shuffled per round
arm                        min    median
sample_BASE             0.1829    0.2024
sample_FIX              0.1841    0.2072
control_core_BASE       0.1544    0.1687
control_core_FIX        0.1544    0.1684

paired per-round difference, median [95% bootstrap interval]
frmtmb.sample FIX - BASE   +0.0000 s [-0.0063, +0.0059]
CONTROL frmtmb FIX - BASE  -0.0021 s [-0.0056, +0.0031]
seed 20260915, B = 20000
```

Load cost is inside the noise. The paired median difference for
frmtmb.sample is +0.0000 s with a 95% interval of about 6 ms either
side, the same width as the control, which loads frmtmb alone and
reads -0.0021 s. The minima differ by about a millisecond. The arm for
the fixed frmtmb.sample also loads the fixed core, so the figure is the
whole change. 25 rounds, 100 fresh processes, arms shuffled within each
round, `Sys.time()` with a measured tick of about 2 us, run with no
other R process on the machine.

**R1. The interval width is a property of machine state, and the
instrument carries no positive control.** A null without one says the
effect is smaller than the interval and nothing about what the
instrument could have seen. The review re-ran the design on a busier
machine, medians about 60 percent larger than mine, and got
`sample_FIX - sample_BASE` of -0.0132 s [-0.0594, +0.0042]: so my plus
or minus 6 ms describes a quiet machine rather than the change. Its
positive-control arm is itself uncalibrated, because `Sys.sleep(0.010)`
on Windows returns on the system timer quantum of about 15.6 ms; the
next run of this script should burn a measured 10 ms of CPU in a loop
instead. The finding, "inside the noise", is not falsified by either
run.

The per-access cost is core's measurement and was not repeated: 1.29 us
on the fallback path and 3.18 us on the adopted path against 0.80 us for
an ordinary binding, and 4.35 us for a two-owner name such as `rhat`
(`dev/generics-adoptcost.R`). These names are reached once per user
action.

## What I decided NOT to do

* **Copy the mechanism into frmtmb.sample.** One argument in core made
  it callable, and a copy would have to repeat every review round core's
  lane went through.
* **Export `frm_bind_generic()` or core's owner table.** The extension
  needs one entry point; exporting more would make internals API.
* **Rename `frm_install_generics()` on export** to the
  `frmtmb_register_*` spelling of the other load-time seams. The name is
  referenced from core's findings, scripts and review; a rename is the
  user's call.
* **Guard the two new directives at run time** to remove the exposure
  in defect 1. That is `s3_register()` again, which core's lane refused
  on measurement, and the exposure needs a partial or old owner loaded
  first.
* **Match `rhat()`'s and `neff_ratio()`'s OUTPUT to brms**, relabel
  what they return, or change the positional order of method
  arguments. See defects 2, 3 and 4, which are the user's calls.
* **Rebuild the pkgdown sites.** `docs/` still shows the old usage for
  `rhat`, `bridge_sampler` and `bayes_factor`; the release rebuilds it.
* **Choose version numbers.**

## What the review could not falsify

**R1**, from `dev/reviews/20260915-samplegen.md`, which derived its own
name list off the BASE install and reran nothing of mine. Kept here
because a confirmation that is only in a review file is one the next
reader will not find: the 0 of 28 claim in 14 load orders on four
channels with a control PROCESS per order; the brms output claim,
widened from 22 calls to 60 and from two channels to four, 0 of 60 in
four FIX arms with a two-process control at 0 of 60; the three rules
the new Rd section states, checked against the installed build (30
(name, owner) pairs, 0 missing twins, 0 formals disagreeing with the
first owner); brms AGREEING with the owner on all six changed formals,
so the tiebreaker question is moot there; all 28 bodies bare on BASE;
the partial-owner load failure in both halves and both orders on both
builds, with both upstream floors verified; the body guard against a
mutant on a different generic in a different file; the split-suite
provenance by timestamp; and that frmtmb.sample cannot load against
core 0.56.0.

The review also corrected my model of the instrument in my favor. Its
first probe modelled `UseMethod()`'s caller-environment leg as reaching
the search path and reported a loss I had not; measured, that leg walks
the enclosure chain to the global environment and stops, so an owner
that EXPORTS a `.default` cannot capture dispatch from a rival generic,
and the caller leg fired 0 times in its 56 cells. My table-only probe is
right for calls from the prompt.

Two nits on the newly public `frm_install_generics()` that I have NOT
acted on, since the review filed them as low risk and the function is
one line from covering both: called outside `.onLoad()` it reports R's
own "cannot remove bindings from a locked environment" rather than
saying when it may be called, and `frm_is_generic()` screens the
OWNER's function but not the caller's own fallback.

## Defects found and NOT fixed

### 1. The two twin directives widen the partial-owner exposure

`dev/samplegen-shadowload.R` builds a posterior that exports everything
the real one does except `rhat`, and a gratia without
`posterior_samples`, and loads each before and after frmtmb.sample.

```
== partial owners, dev/samplegen-shadowload.R ==
rellib-r3 is BASE, samplegen-lib is FIX
BUILD rellib-r3 shadow first 
  SHADOW posterior without rhat loaded 
  FRMTMB.SAMPLE loaded 
BUILD rellib-r3 sample first 
  FRMTMB.SAMPLE loaded 
  SHADOW posterior without rhat loaded 
BUILD rellib-r3 shadow first 
  SHADOW gratia without posterior_samples loaded 
  FRMTMB.SAMPLE loaded 
BUILD rellib-r3 sample first 
  FRMTMB.SAMPLE loaded 
  SHADOW gratia without posterior_samples loaded 
BUILD samplegen-lib shadow first 
  SHADOW posterior without rhat loaded 
  FRMTMB.SAMPLE FAILED: object 'rhat' not found 
BUILD samplegen-lib sample first 
  FRMTMB.SAMPLE loaded 
  SHADOW posterior without rhat loaded 
BUILD samplegen-lib shadow first 
  SHADOW gratia without posterior_samples loaded 
  FRMTMB.SAMPLE FAILED: object 'posterior_samples' not found 
BUILD samplegen-lib sample first 
  FRMTMB.SAMPLE loaded 
  SHADOW gratia without posterior_samples loaded 
```

The same class as core's defect 3, for two more names. Loaded AFTER
frmtmb.sample, nothing breaks; loaded BEFORE, frmtmb.sample does not
load, and base did. posterior's exposure is academic, since 1.0.0
exported `rhat`. gratia's is not: `posterior_samples()` first shipped in
gratia 0.9.0, so a session that loads an older gratia first cannot load
frmtmb.sample. The Suggests floor states the requirement and does not
enforce it. **This is a decision for the user**: accept it, or drop
gratia as an owner and give a gratia user back the base defect, where
`library(gratia); library(frmtmb.sample)` loses all 3 of gratia's
`posterior_samples` methods.

### 2. `rhat()` AND `neff_ratio()` answer a different question from brms

**R1. This replaces the version of this entry that named only `rhat()`
and left the size unmeasured. `neff_ratio()` has the same defect and it
is the worse of the two.** Both are measured on one cached draws
object, 4 chains of 500, seed 20260915
(`dev/sgrev-rhat.R`, `dev/sgrev-rhat2.R`, draws at
`dev/stan-cache/sgrev-draws.rds`).

brms's `rhat.brmsfit()` is `posterior::summarise_draws(draws, rhat =
posterior::rhat)`, the rank-normalized split-Rhat, the maximum of bulk
and tail. `rhat.frmtmb_draws()` calls `bayesplot::rhat()` on the
`stanfit`, which is `rstan::summary(object)$summary[, "Rhat"]`, the
classic split-Rhat.

```
== dev/sgrev-rhat.R, 4 chains x 500, seed 20260915 ==
max relative difference                       0.00563108
max |posterior - 1|, the whole signal         0.0100318
difference as a fraction of (rhat - 1)        1.15351
max |rhat(ds) - rstan's Rhat| over 11 vars    0   (exactly)
```

The disagreement is 1.15 times the ENTIRE excess over 1 that the
diagnostic reports, and `rhat(ds)` equals rstan's number to the bit,
which names the definition in use rather than leaving it to be guessed.

`neff_ratio.brmsfit()` is `min(ess_bulk, ess_tail) / ndraws` from
posterior; `neff_ratio.frmtmb_draws()` is bayesplot on the `stanfit`,
which is rstan's `n_eff` over the total.

```
== dev/sgrev-rhat2.R, same draws ==
max relative difference:  0.412677
per variable, b[2]:  0.41789470 here against 0.29581769 from posterior
```

0.41 against 0.0056 is 73 times larger, and an effective sample size
is what a user reads to decide whether to run the sampler longer.

**The blast radius is exactly these two, and the other three are clean
BY CONSTRUCTION.** Five draws methods reach bayesplot.
`nuts_params()` and `log_posterior()` are right, because brms's own
`brmsfit` methods also delegate to `bayesplot::nuts_params(object$fit)`
and `bayesplot::log_posterior(object$fit)`: matching bayesplot there IS
matching brms. `mcmc_plot()` builds from `as.array()`, which carries
this package's own names. So the gap is two names and no more.

None of this was caused by this lane, none of it is a dispatch defect,
and `test-draws-methods.R` asserts the current behavior on purpose.
**Whether to change the two methods is the user's call**, and it is a
real choice: following brms changes two published diagnostics, and
keeping rstan's keeps them consistent with `ds$stanfit`.

### 3. `rhat()` and `neff_ratio()` return raw Stan parameter names

**R1, and nobody had filed it.** Every other draws accessor reports
frmtmb names. These two report the Stan ones, because they read the
`stanfit` directly:

```
== dev/sgrev-rhat2.R ==
variables(ds)        :  Intercept x sigma_Intercept b[1..6] theta_1 lp__
names(rhat(ds))      :  beta[1] beta[2] betad     b[1..6] theta   lp__
in rhat() but not variables():  beta[1] beta[2] betad theta
in variables() but not rhat():  Intercept x sigma_Intercept theta_1
```

4 of 11 entries are not addressable by the names the rest of the
package uses, so `rhat(ds)["x"]` is `NA` and a user lining the two up
gets a silent mismatch. The package's DESCRIPTION promises draws "under
'frmtmb' parameter names". `?draws-diagnostics` does say these three
read the `stanfit` and therefore show Stan's own names, so it is
documented, not hidden. It is the same two functions as defect 2 and
the same decision: **the user's call**, together with whether they
should go through posterior on the relabeled array, which would settle
both items at once.

### 4. Positional arguments differ from brms on 10 of 28 methods

**R1. This replaces two examples with the review's full audit**
(`dev/sgrev-positional.R`, `dev/sgrev-api2.R`). Each row is the first
position at which a positional caller is answering a different
question, counting only positions before brms's `...`:

```
generic              pos  brms's method      frmtmb.sample's
as.mcmc              2    pars               combine_chains
log_lik              2    newdata            ndraws
mcmc_plot            2    pars               type
posterior_epred      3    re_formula         resp
posterior_interval   2    pars               prob
posterior_linpred    4    re_formula         resp
posterior_predict    3    re_formula         resp
pp_mixture           2    newdata            summary
predictive_error     2    newdata            resp
psis                 2    newdata            ndraws

disagree at some positional slot: 10   agree as far as both go: 18
```

Eight of the ten ERROR when exercised. **Two return a value with
nothing said**: `as.mcmc(ds, TRUE)` gives an `mcmc` of length 22000 and
`posterior_interval(ds, 0.9)` a matrix of length 8. brms's
`re_formula = NA` idiom passed positionally is an error rather than a
silent answer: `posterior_epred(ds, nd, NA)` says `Unknown response:
'NA'`.

The GENERICS match their owners, which is this lane's subject; these
are METHOD formals and were left alone.

### 5. bbmle's `parnames` is still masked

bbmle exports a plain function under that name. The binding correctly
refuses to adopt it, so `library(bbmle); library(frmtmb.sample)` still
masks bbmle's function with this package's generic. That is a
non-generic name collision, core's defect 2, not a dispatch defect.
**R1.** The review priced it: `bbmle::parnames(x)` works, and the
replacement form `parnames<-` is not masked at all, so what a bbmle
user loses is the bare getter at the prompt.

## Version

I am not choosing numbers.

* **Core needs a bump.** It gains one export, `frm_install_generics()`,
  and one argument on it. Nothing existing changes behavior.
* **frmtmb.sample's floor on frmtmb must move to that bumped core.**
  It is a hard requirement: `.onLoad()` calls `frm_install_generics()`
  with two arguments, which no released core exports, so against core
  0.56.0 frmtmb.sample fails to LOAD. The worktree DESCRIPTION still
  reads `frmtmb (>= 0.55.1)`, because this worktree's base predates the
  release's 0.56.0 bump; the consolidator sets it. **R1: left at that
  number deliberately, on instruction. It must not merge with it.**
* **frmtmb.sample needs a bump**, with breaking changes to three
  generic signatures, listed in its NEWS.
* Suggests floors `gratia (>= 0.9.0)` and `posterior (>= 1.0.0)` are
  dependency floors, not version choices.

## What was run

Both full suites on the final build, one process per file,
`NOT_CRAN` and `FRMTMB_BRMS_FIT_TESTS` true, with this worktree's
`dev/stan-cache` (copied from main's). frmtmb.sample ran two files at a
time. Core's first 28 files ran two at a time before the driving session
ended; the other 108 ran three at a time after it resumed, on the same
installed build, which `dev/samplegen-parsecheck.R` confirmed was newer
than every source edit and whose 28 changed R files all parse.
`test-loo.R` was rerun alone after its mock was changed.

```
== suite: frmtmb.sample, FIX build, 18 files, one process per file, NOT_CRAN and FRMTMB_BRMS_FIT_TESTS true ==
test files with a result   18 of 18
blocks                     208
assert                     1213
pass                       1213
fail                       0
error                      0
skip                       1
warn                       0

per file: blocks / assertions / fail / error / skip
  test-bracket-access                         1     1   0   0   0
  test-compat-preflight                      16    81   0   0   0
  test-conditional-effects-draws             12    57   0   0   0
  test-draws-methods                         19   101   0   0   0
  test-draws-spellings                        5    48   0   0   0
  test-evidence-ratio                        10    33   0   0   0
  test-generic-collision                     13    58   0   0   0
  test-loo                                   20    79   0   0   0
  test-message-uniqueness                     1     6   0   0   0
  test-parallel-chains                        2     4   0   0   0
  test-prior-route                            3     9   0   0   0
  test-reparam                               28   260   0   0   0
  test-sample-direct                         17   136   0   0   0
  test-sampling-ported                       34   212   0   0   0
  test-scale                                  1     1   0   0   1
  test-simulators                             9    50   0   0   0
  test-stan-control                          11    54   0   0   0
  test-tmbstan-build-guard                    6    23   0   0   0

no file reported a failure or an error 
```

```
== suite: frmtmb core, FIX build, one process per file, NOT_CRAN and FRMTMB_BRMS_FIT_TESTS true ==
test files with a result   136 of 136
blocks                     1464
assert                     10221
pass                       10221
fail                       0
error                      0
skip                       1
warn                       2

per file: blocks / assertions / fail / error / skip
  test-aliased-grouping                       3    33   0   0   0
  test-api-spellings                         11    36   0   0   0
  test-autocor                               30   191   0   0   0
  test-autoscale                              6    46   0   0   0
  test-backlog                               11    35   0   0   0
  test-bcm-bart                               3    12   0   0   0
  test-bcm-binomial                          14    37   0   0   0
  test-bcm-data-analysis                     15    46   0   0   0
  test-bcm-esp                                4    17   0   0   0
  test-bcm-gaussian                           7    23   0   0   0
  test-bcm-gcm                                4    17   0   0   0
  test-bcm-latent-mixtures                   13    52   0   0   0
  test-bcm-model-selection                    9    29   0   0   0
  test-bcm-mpt                                5    27   0   0   0
  test-bcm-psychophysics                      4    28   0   0   0
  test-bcm-retention                          6    32   0   0   0
  test-bcm-signal-detection                   5    30   0   0   0
  test-bcm-simple                             3    13   0   0   0
  test-boot                                   7    44   0   0   0
  test-bracket-access                         6    33   0   0   0
  test-brms-agreement                        19   187   0   0   0
  test-brms-likelihood                       36   405   0   0   0
  test-brms-methods                          47   950   0   0   0
  test-brms-port                              6    18   0   0   0
  test-brms-priors                           12   103   0   0   0
  test-car-spde                              23   135   0   0   0
  test-case-studies                           9    30   0   0   0
  test-ce-bands                              22   167   0   0   0
  test-ce-facets                              6    92   0   0   0
  test-cens-lccdf                             6    27   0   0   0
  test-cens-trunc                            20    69   0   0   0
  test-ci-siblings                            1    11   0   0   0
  test-compat-register                       27   103   0   0   0
  test-compat                                34   573   0   0   0
  test-confint-anova                          9    35   0   0   0
  test-core-boundary                          2     4   0   0   0
  test-covstruct                              4     7   0   0   0
  test-custom-family                         38   126   0   0   0
  test-data2                                  8    30   0   0   0
  test-dates                                  7    13   0   0   0
  test-diagnostics-ux                        27   127   0   0   0
  test-distributional                         6    28   0   0   0
  test-dpar-eta-api                          14   123   0   0   0
  test-dry-run                                2    23   0   0   0
  test-edgecases                             17    60   0   0   0
  test-effects                                9    54   0   0   0
  test-famgaps                               21    99   0   0   0
  test-families                              17   227   0   0   0
  test-fit-end-hook                           4    18   0   0   0
  test-frame                                  5    22   0   0   0
  test-functional                             1     3   0   0   0
  test-fuzz                                   1     1   0   0   1
  test-generic-collision                     13    56   0   0   0
  test-get-prior-route                        8    75   0   0   0
  test-glmm-binomial                          2    10   0   0   0
  test-glmm-gaussian                          4    17   0   0   0
  test-glmm-poisson                           4    12   0   0   0
  test-gp-multidim                            5    42   0   0   0
  test-habit-prep                             3    16   0   0   0
  test-id-kron                               11    58   0   0   0
  test-importance                            38   226   0   0   0
  test-influence-plot                         4    21   0   0   0
  test-input-validation                       7    43   0   0   0
  test-interop                                6    44   0   0   0
  test-lkj                                   13   212   0   0   0
  test-map                                    2     7   0   0   0
  test-mean-fn                                1     4   0   0   0
  test-message-uniqueness                     1     6   0   0   0
  test-method-residue                        12    64   0   0   0
  test-methods-audit                         10    57   0   0   0
  test-methods                               11    59   0   0   0
  test-mixture-start                          2     5   0   0   0
  test-mm                                    27   116   0   0   0
  test-mo-terms                               7    74   0   0   0
  test-multinomial                            4    15   0   0   0
  test-multiple-pooling                       3    33   0   0   0
  test-multivariate                           6    23   0   0   0
  test-mvn-mixture                            9   277   0   0   0
  test-naming-collisions                      6    31   0   0   0
  test-nl-body-vars                           6    23   0   0   0
  test-nl-lexical                             7    31   0   0   0
  test-nl-rtmb-scope                         15   103   0   0   0
  test-nl                                    11    44   0   0   0
  test-nlf                                   11    74   0   0   0
  test-numerical-robustness                  37   679   0   0   0
  test-open-issues                           13    43   0   0   0
  test-ordinal-fitted                        13    87   0   0   0
  test-ordinal                               11   109   0   0   0
  test-osa-inference                          6    34   0   0   0
  test-par-template                          11    51   0   0   0
  test-parse                                 14    43   0   0   0
  test-perf                                   2     3   0   0   0
  test-pooled-anova                           6    64   0   0   0
  test-portability                           13    91   0   0   0
  test-predict-lp-basis                       5    30   0   0   0
  test-predict-newdata                        6    12   0   0   0
  test-prior-compat                          28   196   0   0   0
  test-priors-autocor-classes                10    63   0   0   0
  test-priors-bounds-grcov                   10    49   0   0   0
  test-ps                                    12    71   0   0   0
  test-quadrature-defects                     7    59   0   0   0
  test-reml-nongaussian                       4     9   0   0   0
  test-review-fixes                           6    17   0   0   0
  test-review-v25                            13    60   0   0   0
  test-review-v28                            13    73   0   0   0
  test-review-v29                            18   136   0   0   0
  test-rl-example                            14   102   0   0   0
  test-sandwich                               9    57   0   0   0
  test-scale-contract                        12    37   0   0   0
  test-setprior                               5    27   0   0   0
  test-simulate-density                      50   435   0   0   0
  test-simulate-ergonomics                    9    50   0   0   0
  test-smooth-population                      7    33   0   0   0
  test-smooths                                8    27   0   0   0
  test-sparsear1                              8    24   0   0   0
  test-sparsex                                5    64   0   0   0
  test-spectral                              21   113   0   0   0
  test-structure                             28   135   0   0   0
  test-sugar                                  8    38   0   0   0
  test-tabular-inputs                         4     9   0   0   0
  test-tmb-examples                          15    28   0   0   0
  test-tre                                   18   101   0   0   0
  test-trunc-postfit                          7    35   0   0   0
  test-unpinned-seams                         4    28   0   0   0
  test-v07                                    6    25   0   0   0
  test-v11                                    7    22   0   0   0
  test-v12                                    5    26   0   0   0
  test-v14                                   14    81   0   0   0
  test-v15                                    9    41   0   0   0
  test-v16                                    3    24   0   0   0
  test-v17                                    8    35   0   0   0
  test-v18                                    4    27   0   0   0
  test-v19                                    4    17   0   0   0
  test-verbose                                9    40   0   0   0
  test-vignette-wiener                        5    21   0   0   0
  test-zi                                     5    11   0   0   0

no file reported a failure or an error 
```

The same frmtmb.sample suite in ONE process, the way `R CMD check`
runs it, which is the run that caught core's `hypothesis()` defect:

```
== dev/samplegen-oneproc.R ==
FILES   18 
BLOCKS  208 
ASSERT  1213 
FAIL    0 
ERROR   0 
SKIP    1 
brms loaded at the end: TRUE 
```

`R CMD check --as-cran`, once each, built with vignettes and without
`--no-manual`, on a quiet machine (`dev/samplegen-cran.sh`):

```
== dev/samplegen-cran.sh ==
-- cran-core.txt
* checking HTML version of manual ... [19s] NOTE
Skipping checking math rendering: package 'V8' unavailable
Status: 1 NOTE
See
-- cran-sample.txt
Status: OK
```

Before the check was allowed, the sections a NAMESPACE, signature or
Rd change can break were run on the installed packages directly:

```
==== frmtmb from C:/Users/adf44/source/r/samplegen-lib/frmtmb 
-- undoc            clean
-- codoc            clean
-- checkDocFiles    clean
-- checkS3methods   clean
-- checkReplaceFuns clean
==== frmtmb.sample from C:/Users/adf44/source/r/samplegen-lib/frmtmb.sample 
-- undoc            clean
-- codoc            clean
-- checkDocFiles    clean
-- checkS3methods   clean
-- checkReplaceFuns clean
```

Every changed Rd page, rendered with `Rd2txt` and grepped:

```
ok   frmtmb-sampling-api.Rd                        The borrowed-generic seam
ok   frmtmb-sampling-api.Rd                        frm_install_generics(pkgname, owners)
ok   frmtmb-sampling-api.Rd                        S3method(owner::generic, class)
ok   frmtmb-sampling-api.Rd                        bare UseMethod()
ok   frmtmb-loo-refusals.Rd                        bridge_sampler(samples, ...)
ok   frmtmb-loo-refusals.Rd                        bayes_factor(x1, x2, log = FALSE, ...)
ok   frmtmb-loo-refusals.Rd                        post_prob(x, ..., prior_prob = NULL, model_names = NULL)
ok   frmtmb-loo-refusals.Rd                        x, x1, x2, samples, log, prior_prob, model_names, ...
ok   draws-diagnostics.Rd                          rhat(x, ...)
ok   frmtmb-draws-refusals.Rd                      posterior_samples(x, pars = NA, ...)
ok   frmtmb-draws-refusals.Rd                      object, x, pars, ...
ok   posterior_epred.Rd                            posterior_linpred(object, transform = FALSE, ...)
missing: 0 
```

## Scripts, all prefixed `samplegen-`

| script | what it produced |
|---|---|
| `dev/samplegen-audit.R` | the ownership audit |
| `dev/samplegen-formals.R` | every owner's formals beside brms's method |
| `dev/samplegen-check.R`, `dev/samplegen-checkall.sh` | the per-order probe, one process per mode |
| `dev/samplegen-summary.R` | the collision table |
| `dev/samplegen-brmsout.R`, `dev/samplegen-brmsout-summary.sh` | brms output hashes, per arm |
| `dev/samplegen-runtests.R`, `dev/samplegen-runtests-sample.R` | one test file per process, errors counted |
| `dev/samplegen-suite.sh`, `dev/samplegen-suitesummary.R` | the per-file suite, two at a time, and its counts |
| `dev/samplegen-oneproc.R` | the frmtmb.sample suite in ONE process |
| `dev/samplegen-testsummary.sh` | the seen-failing table |
| `dev/samplegen-mutant.R` | the work-in-the-generic mutant |
| `dev/samplegen-shadowload.R` | the partial-owner exposure |
| `dev/samplegen-loadcost.R` | load cost, interleaved, hi-res clock |
| `dev/samplegen-cran.sh` | `R CMD check --as-cran`, mirroring `dev/generics-cran.sh` |
| `dev/samplegen-parsecheck.R` | after the interrupted session: every changed R file parses, installs not stale |
| `dev/samplegen-toolscheck.R` | `undoc`, `codoc`, `checkDocFiles`, `checkS3methods` |
| `dev/samplegen-rdcheck.R` | the changed Rd pages, rendered |
| `dev/samplegen-build.R` | roxygenise and install into the private library |
| `dev/samplegen-whoscoda.R` | which owner frmtmb.sample's own load pulls in (coda) |
| `dev/samplegen-assemble.R` | splices the blocks into this page |

Outputs are under `dev/samplegen-out/`. The brmsfit is
`dev/stan-cache/samplegen-brmsfit.rds` (gitignored), a copy of the
generics lane's `generics-scale-brmsfit.rds`.
