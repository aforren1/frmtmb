# Review: lane `samplegen`, frmtmb.sample and the borrowed generics

Reviewer pass on `wt-samplegen`, 2026-09-15. Adversarial pass, because
this lane changes CORE, exports a new function, ships three breaking
signatures and moves a hard dependency floor. Nothing in the worktree
was edited by this review except this file and my own `sgrev-`
scripts. No git operation was run.

Reviewer library `C:/Users/adf44/source/r/sgrev-lib`, built by me from
this worktree (core and frmtmb.sample). BASE read-only from
`C:/Users/adf44/source/r/rellib-r3` (core 0.56.0, frmtmb.sample 0.4.2).
Two mutant libraries of my own, `sgrev-mut-novalidate-lib` (core) and
`sgrev-mut-loglik-lib` (frmtmb.sample), and a shadow library
`sgrev-shadowlib`. Scripts are in `dev/`, prefixed `sgrev-`, with
output under `dev/sgrev-out/`. Every number below came from a script
of mine. The lane's scripts were read, never rerun.

**Verdict: mergeable with named fixes.** Two BLOCKERs, both on the
record rather than on the mechanism, and nine nits. The dispatch
repair is correct and survives every order I could build, including
four the lane did not run. What does not survive is the size of
defect 2: `rhat()` is neither the only diagnostic that answers a
different question from brms, nor the worst one.

## Summary of the attack

| claim | outcome |
|---|---|
| 28 generics defined, all owned elsewhere | reproduced, derived independently |
| all 28 bodies are a bare `UseMethod()` on BASE | reproduced, 0 of 28 carry work |
| 0 of 28 lost in every load order | reproduced, and in 4 more orders |
| `getS3method()` checked separately | done here; 67 of 67 lost on BASE-S, 0 on FIX |
| the collision probe's lookup model | CORRECTED: the caller-env leg does not reach the search path |
| brms output, 22 calls, 0 of 22 on FIX | reproduced and widened to 60 calls, with messages and printed form |
| the draws side | NOT MEASURED BY THE LANE; measured here, 0 of 46 beyond the two warnings NEWS declares |
| brms dispatches `rhat` to posterior's | CONFIRMED by reading brms's own binding and table |
| `rhat()` on draws does not match brms | CONFIRMED and SIZED; and `neff_ratio()` is worse |
| the three breaking signatures follow the owner | confirmed, and brms agrees with the owner in all six |
| "some method positional arguments differ from brms" | enumerated in full: 10 of 28, 2 of them silent |
| the active-binding mocking trap | mechanism confirmed; the fix is right, its comment is incomplete |
| a posterior without `rhat` or an old gratia stops the load | reproduced, both halves, both orders, both builds |
| `gratia (>= 0.9.0)` and `posterior (>= 1.0.0)` | both floors confirmed from upstream |
| the body guard's mutant is real | confirmed, and it catches a mutant of mine on a different generic |
| the core block fails on BASE only in the weak form | confirmed; the behavioural failure is constructed here |
| the core suite ran against one unchanged build | confirmed by install and source timestamps |
| load cost inside the noise | not falsified; the instrument carries no positive control |
| bbmle's `parnames` is still masked, `::` suffices | confirmed, including the replacement form |
| frmtmb.sample cannot load against core 0.56.0 | confirmed by construction |

## 0. The instrument, before any of its readings

The lane's probe reads only the `.__S3MethodsTable__.` of
`environment(generic)`. R's documented dispatch also looks in the
environment the generic was CALLED from, and the two legs are
interleaved per class, not run one after the other. My first probe
modelled the caller leg as the search path and promptly reported a
loss the lane had not: `loo_moment_match` on draws falling into loo's
EXPORTED `loo_moment_match.default`, in modes P and Q, on BASE **and**
on FIX. It was my bug, and the measurement that settles it is worth
keeping because it justifies the lane's simpler probe
(`dev/sgrev-leg1.R`, `dev/sgrev-leg1b.R`):

```
loo exports loo_moment_match.default:  TRUE
visible from globalenv:  TRUE
rival generic, class zzz, caller = globalenv:   "no applicable method"
rival generic, caller = a function frame:       "no applicable method"
with a .default in the caller's own frame:      "LOCAL DEFAULT RAN"

1 top-level method, top-level call:              TOPLEVEL
2 top-level method, called inside a function:    TOPLEVEL
3 method on the search path, top-level call:     not found
3 method on the search path, called in a function: not found
```

So the caller leg walks the enclosure chain as far as the global
environment and stops there. An owner that EXPORTS a `.default`
cannot capture dispatch away from a rival generic. **The lane's
table-only probe is therefore right for calls from the prompt**, and
my probe now runs both legs interleaved and reports that the caller
leg fired 0 times in 56 cells.

The probe itself was then checked against real dispatch rather than
against my reading of R (`dev/sgrev-dispatchcheck.R`): the generic is
CALLED on an object of the class and the function that actually got
the frame is read out of the call stack. Two spellings of that script
were wrong before it worked, and both are instructive. A `tryCatch()`
established inside a `withCallingHandlers()` unwinds first, so the
calling handler never runs and the script reported 56 of 56
undetermined. And `max()` over the frames where the generic appears
reads a NESTED call, because `rhat.frmtmb_draws()` and the three
bayesplot delegations call the owner's generic again from inside.
With both fixed, prediction and observation agree on every cell except
four, and those four are an artefact of my "no applicable method"
shortcut firing on the nested call, which `dev/sgrev-stackdump.R`
shows frame by frame.

## 1. The collision table. Reproduced, and extended

`dev/sgrev-probe.R`, one fresh process per cell, 28 test cells and 28
control cells. The control is a SEPARATE process with the same load
order and frmtmb.sample never loaded, not an in-process reconstruction.
The 28 names are derived in their own process off the BASE install
(`dev/sgrev-names.R`: exported by frmtmb.sample, defined in its
namespace, body a dispatch) and not read from the lane's table. They
come to exactly 28 and match the lane's list name for name.

```
build mode order                                      brmsfit   draws   owners    gS3
BASE  S    brms, then frmtmb.sample                   28/28     0/28    39/39     67/67
BASE  T    frmtmb.sample, then brms                   0/28      1/28    0/39      0/67
      draws lost:  rhat<-posterior.default
BASE  U    frmtmb.sample; brms only loaded            0/0       0/28    0/0       0/0
BASE  N    frmtmb.sample alone                        0/0       0/28    0/0       0/0
BASE  P    7 other owners, then frmtmb.sample         0/0       0/28    29/29     28/29
BASE  Q    frmtmb.sample, then 7 other owners         0/0       2/28    0/29      0/29
      draws lost:  posterior_samples<-gratia.default  rhat<-posterior.default
BASE  G    gratia, then frmtmb.sample                 0/0       0/28    2/2       2/2
BASE  D    brms, sample; detach and reattach both     28/28     0/28    39/39     67/67
BASE  R    loo, sample; loo unloaded; brms loaded     0/28      1/28    0/39      0/67
BASE  RN   brms attached; sample requireNamespace ONLY 0/28      1/28    0/39      0/67
BASE  RN2  sample requireNamespace ONLY; then brms    0/28      1/28    0/39      0/67
BASE  BR   both attached; brms unloaded and RELOADED  0/28      1/28    0/39      0/67
BASE  CA   brms, frmtmb, then frmtmb.sample           28/28     0/28    39/39     67/67
BASE  CA2  brms, frmtmb.sample, then frmtmb           28/28     0/28    39/39     67/67
FIX   S    brms, then frmtmb.sample                   0/28      0/28    0/39      0/67
FIX   T    frmtmb.sample, then brms                   0/28      0/28    0/39      0/67
FIX   U    frmtmb.sample; brms only loaded            0/0       0/28    0/0       0/0
FIX   N    frmtmb.sample alone                        0/0       0/28    0/0       0/0
FIX   P    7 other owners, then frmtmb.sample         0/0       0/28    0/29      0/29
FIX   Q    frmtmb.sample, then 7 other owners         0/0       0/28    0/29      0/29
FIX   G    gratia, then frmtmb.sample                 0/0       0/28    0/2       0/2
FIX   D    brms, sample; detach and reattach both     0/28      0/28    0/39      0/67
FIX   R    loo, sample; loo unloaded; brms loaded     0/28      0/28    0/39      0/67
FIX   RN   brms attached; sample requireNamespace ONLY 0/28      0/28    0/39      0/67
FIX   RN2  sample requireNamespace ONLY; then brms    0/28      0/28    0/39      0/67
FIX   BR   both attached; brms unloaded and RELOADED  0/28      0/28    0/39      0/67
FIX   CA   brms, frmtmb, then frmtmb.sample           0/28      0/28    0/39      0/67
FIX   CA2  brms, frmtmb.sample, then frmtmb           0/28      0/28    0/39      0/67
active in namespace: 0/28 on every BASE row, 28/28 on every FIX row
resolutions from the caller-env leg, over all 56 cells: 0
```

`RN`, `RN2` and `BR` are the three orders the task named as untried,
and `CA`/`CA2` are two more of mine, with CORE attached as well so
that the 26 re-exported names have two attached copies. `BR` is the
one that could have broken the memo: brms is detached, its namespace
unloaded, and brms attached again under a live frmtmb.sample. The memo
holds a namespace ENVIRONMENT and compares by identity, and a reloaded
namespace is a new environment, so it re-resolves. 0 of 28.

Three notes on the numbers rather than on the finding.

* **The denominators are not the lane's and should not be read as
  disagreement.** My `owners` column counts (name, class) pairs
  deduplicated across owners with `default` excluded; the lane's
  counts every entry in every owner's table. 39 against 77 in mode S,
  29 against 41 in P, 2 against 3 in G. Both are 100 percent lost on
  BASE and 0 on FIX.
* **`getS3method()` is a column here, which the lane's frmtmb.sample
  table does not have** although core's did. 67 of 67 lost in mode S
  on BASE, 0 of 67 on FIX. Nit: add it, so the two packages' tables
  read alike.
* BASE mode P is the one row where dispatch and `getS3method()`
  disagree, 29 of 29 against 28 of 29. It is a BASE-only artefact and
  has no bearing on the fix.

The two silent BASE failures the lane reported as the measurement's
own finding both reproduce, and I can give them their user-visible
form, which the lane's table cannot (section 2).

## 2. brms held to output. Reproduced, widened, and turned around

The lane hashes value and warnings on 22 calls and **muffles messages
entirely**, so a difference in messages was invisible to it. It also
reduces `mcmc_plot()` to `nrow(...$data)` and `restructure()` to a
class and a dimension.

`dev/sgrev-brmsout.R` runs 60 calls on the same cached brmsfit and
hashes four things: value, warnings, MESSAGES and the PRINTED form.
The extra 38 are the forms a brms user actually types: positional
second and third arguments, `newdata`, `prob`, `version = FALSE`,
`ggplot_build()` on the two plotting methods rather than a row count,
`print()` and `summary()` on the fit, and the core names
frmtmb.sample re-exports, which the lane's list omitted although they
sit on the search path in every one of its arms.

```
CONTROL brms-only vs brms-only, two processes: 0 of 60 differ
deterministic calls: 60 of 60

arm                                   value     warnings  messages  printed
BASE-S                                12/60     3/60      1/60      12/60
BASE-T                                0/60      0/60      0/60      0/60
FIX-S                                 0/60      0/60      0/60      0/60
FIX-T                                 0/60      0/60      0/60      0/60
FIX-U                                 0/60      0/60      0/60      0/60
FIX-SC                                0/60      0/60      0/60      0/60
```

`FIX-SC` attaches core as well. The zero is a measured zero: the
control is two independent processes and they agree on all 60,
seeded `posterior_predict()` included. BASE-S differs on 12 rather
than the lane's 22 of 22 because my arms attach the other seven owners
AFTER frmtmb.sample, so only the seven brms-owned names are in front;
the calls that differ are exactly those seven names' calls.

**A message the lane's harness could not have seen.** On BASE-S,
`mcmc_plot(bfit, type = "hist")` loses ggplot2's `stat_bin()` message
along with its value. That is one message in 60, it is on BASE, and it
disappears on the fix; the point is that the instrument now covers the
channel.

**The side the lane did not measure.** The fix changes which generic
runs for a `frmtmb_draws` as well, and a generic can warn or coerce
before it dispatches. `dev/sgrev-drawsout.R` runs 46 calls on one
cached `frmtmb_draws` (4 chains, 500 iterations), BASE against FIX, in
five load orders, on the same four channels. The reference arm run
twice differs on 0 of 46.

```
BASE-S   differs on  0 of 46 calls
BASE-T   differs on  3 of 46 calls
    rhat            ERROR "unimplemented type 'list' in 'greater'"
    parnames            + brms's deprecation warning
    posterior_samples   + brms's deprecation warning
BASE-Q   differs on  2 of 46 calls
    rhat            ERROR "unimplemented type 'list' in 'greater'"
    posterior_samples   gratia's error, not this package's refusal:
      "Don't know how to sample from the posterior of <frmtmb_draws>FALSE"
FIX-N    differs on  0 of 46 calls
FIX-S    differs on  2 of 46 calls  (the two deprecation warnings)
FIX-T    differs on  2 of 46 calls  (the two deprecation warnings)
FIX-P    differs on  0 of 46 calls
FIX-Q    differs on  0 of 46 calls
```

Two things fall out. The lane's two "silent" BASE failures are not
silent at the prompt, they are a cryptic error and a foreign package's
error, and those strings belong on the page because they are what a
user would have reported. And the fix changes NOTHING else on the
draws side: the only difference from BASE is the two brms deprecation
warnings, which NEWS already declares in as many words.

Nit: put the draws-side run on the findings page. Under a standard
that holds frmtmb.sample to brms "down to output", the output on ITS
OWN class is the half the standard governs, and the page currently
carries only the brmsfit half.

## 3. BLOCKER 1. Defect 2 is understated, and `rhat` is not the worst

Two halves, both measured on the same cached draws
(`dev/sgrev-rhat.R`, `dev/sgrev-rhat2.R`).

**The tiebreaker question is settled in the lane's favour.** brms's own
`rhat()` is posterior's generic, read out of brms's namespace, and
`rhat.brmsfit` is in posterior's table and not in bayesplot's. Taking
posterior's `(x, ...)` is right.

**The number is not.** brms's `rhat.brmsfit()` computes
`posterior::summarise_draws(draws, rhat = posterior::rhat)`, the
rank-normalized split-Rhat, maximum of bulk and tail.
`rhat.frmtmb_draws()` calls `bayesplot::rhat()` on the `stanfit`,
which is `rstan::summary(object)$summary[, "Rhat"]`, the classic
split-Rhat.

```
variable                   rhat(ds)    posterior  rel.diff
b[1]                     1.00772300   1.00530366  2.407e-03
b[2]                     1.00333508   1.00533295  1.987e-03
b[3]                     1.00544315   1.00588723  4.415e-04
b[4]                     1.00696662   1.00457556  2.380e-03
b[5]                     1.00052987   1.00619584  5.631e-03
b[6]                     0.99941111   1.00383618  4.408e-03
lp__                     1.01059816   1.01003176  5.608e-04

max relative difference                       0.00563108
max |posterior - 1|, the whole signal         0.0100318
difference as a fraction of (rhat - 1)        1.15351
max |rhat(ds) - rstan's Rhat| over 11 vars    0   (exactly)
```

The disagreement is **1.15 times the entire excess over 1** that the
diagnostic is reporting, and `rhat(ds)` equals rstan's number to the
bit, which names the definition being used.

**And `neff_ratio()` is the same defect, unfiled and four times
larger.** brms's `neff_ratio.brmsfit()` is
`min(ess_bulk, ess_tail) / ndraws` from posterior;
`neff_ratio.frmtmb_draws()` is bayesplot on the stanfit, which is
rstan's `n_eff` over the total.

```
variable       neff_ratio(ds) posterior bulk   rel.diff
b[1]               0.28588095     0.24788810  1.533e-01
b[2]               0.41789470     0.29581769  4.127e-01
b[3]               0.41717875     0.29958122  3.925e-01
b[4]               0.33033134     0.27561938  1.985e-01
b[5]               0.38534959     0.31367754  2.285e-01
b[6]               0.33663946     0.26107820  2.894e-01
lp__               0.22158448     0.22394504  1.054e-02
max relative difference:  0.412677
```

**The blast radius is exactly two names, and I checked the other
three.** Five draws methods reach bayesplot. `nuts_params()` and
`log_posterior()` are FINE, because brms's own methods also delegate
to `bayesplot::nuts_params(object$fit)` and
`bayesplot::log_posterior(object$fit)`; matching bayesplot there IS
matching brms. `mcmc_plot()` builds from `as.array()`. So the gap is
`rhat` and `neff_ratio`, and no more.

**A third thing nobody has filed.** Those two return RAW STAN
parameter names, where every other draws accessor returns frmtmb
names:

```
variables(ds)        :  Intercept x sigma_Intercept b[1..6] theta_1 lp__
names(rhat(ds))      :  beta[1] beta[2] betad     b[1..6] theta   lp__
names(neff_ratio(ds)):  beta[1] beta[2] betad     b[1..6] theta   lp__
in rhat() but not variables():  beta[1] beta[2] betad theta
in variables() but not rhat():  Intercept x sigma_Intercept theta_1
```

4 of 11 entries are not addressable by the names the rest of the
package uses, so `rhat(ds)["x"]` is `NA` and a user lining the two up
gets a silent mismatch. The package's own DESCRIPTION promises draws
"under 'frmtmb' parameter names".

**BLOCKER, and it is a record fix, not a code fix.** None of this was
caused by this lane and none of it is a dispatch defect. But defect 2
currently says "the size of the difference was not measured here" and
names only `rhat`. Leaving it at that means the next reader will not
reopen it, and this is the silent-wrong-answer class the project ranks
highest. Defect 2 must name `neff_ratio` as well, carry these numbers
and their construction, state that `nuts_params` and `log_posterior`
are clean and why, and file the raw-name mismatch as its own item.

## 4. The three breaking signatures, and the positional audit in full

The tiebreaker is brms, not the owner, so I compared against brms
directly: the value of each name INSIDE brms's namespace, and brms's
own `brmsfit` method wherever it is registered
(`dev/sgrev-owners.R`, `dev/sgrev-out/owners.txt`).

| name | brms's generic comes from | brms's formals | the lane's |
|---|---|---|---|
| `bridge_sampler` | bridgesampling | `(samples, ...)` | `(samples, ...)` |
| `bayes_factor` | bridgesampling | `(x1, x2, log = FALSE, ...)` | same |
| `rhat` | posterior | `(x, ...)` | `(x, ...)` |
| `post_prob` | bridgesampling | `(x, ..., prior_prob, model_names)` | same |
| `posterior_samples` | brms | `(x, pars = NA, ...)` | same |
| `posterior_linpred` | rstantools | `(object, transform = FALSE, ...)` | same |

brms and the owner AGREE in all six, so following the owner and
following brms are the same act here and the lane picked a defensible
authority either way. The three rules the new Rd section states are
also true of the build (`dev/sgrev-twins.R`): 30 (name, owner) pairs
checked, 0 missing twins, 28 of 28 in this package's own table, 0
generics whose formals differ from their first owner's, 0 methods
missing a formal their generic has.

**The lane's own flag, enumerated.** `dev/sgrev-positional.R` lines
every `frmtmb_draws` method up against brms's `brmsfit` method and
reports the first position at which a positional caller is answering a
different question, counting only positions before brms's `...`:

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

Eight of the ten ERROR when exercised; two return a value with nothing
said (`dev/sgrev-api2.R`):

```
  as.mcmc(ds, TRUE)                        RETURNED mcmc of length 22000
  posterior_interval(ds, 0.9)              RETURNED matrix of length 8
  log_lik(ds, nd)                          ERROR
  mcmc_plot(ds, 'b_x')                     ERROR
  posterior_epred(ds, nd, NA)              ERROR: Unknown response: 'NA'
  posterior_linpred(ds, TRUE, nd, NA)      ERROR
  posterior_predict(ds, nd, NA)            ERROR
  pp_mixture(ds, nd)                       ERROR
  predictive_error(ds, nd)                 ERROR
  psis(ds, nd)                             ERROR
```

I went looking for the worst case, brms's `re_formula = NA` idiom
passed positionally, and it is an error and not a silent answer:
`posterior_epred(ds, nd, NA)` says `Unknown response: 'NA'`, while
`re_formula = NA` named does the right thing and differs from the
default. Nit: replace defect 3's two examples with this table. Ten of
twenty-eight with the position named is a finding; "several, not
audited in full" is not.

## 5. The mocking trap. Mechanism confirmed, comment incomplete

`dev/sgrev-mock.R`, `dev/sgrev-mock2.R`.

The mechanism is exactly as described. An assignment to an active
binding calls the binding function with the value, and
`frm_bind_generic()` installs a function of no arguments:

```
assign to a NO-ARG active binding: unused argument (base::quote(1))
local_mocked_bindings(log_lik = ) on the FIX build:
    unused argument (base::quote(function (x, ...) 1))
```

The surface is 25 active bindings in frmtmb and 28 in frmtmb.sample,
53 as the lane says, and I listed both sets.

The fix is correct, and correct for the reason the comment gives, but
the comment is not the whole mechanism and the missing half is the one
that could mislead the next author. `local_mocked_bindings()` patches
the namespace binding AND this package's own method table, but NOT the
owner's:

```
                                   rstantools NOT loaded   loaded
namespace binding changed               TRUE               TRUE
frmtmb.sample's method table changed    TRUE               TRUE
rstantools' table entry is the mock      -                 FALSE
called from globalenv               MOCK RAN    "attempt to set an attribute on NULL"
called from inside the namespace    MOCK RAN    MOCK RAN
```

So with the owner loaded, which is the normal session, the mock is
visible ONLY to a caller inside frmtmb.sample's namespace. That is the
case the test is in, so the test is sound. Nit: the comment should say
that the caller-environment leg is what carries it when the owner is
loaded, and that a mock of this shape does not cover a caller outside
the namespace, because the next author will read this comment as a
recipe.

**Guard or rule?** A rule. A guard would have to grep test sources for
`local_mocked_bindings(<one of 53 names> =` and would fire on a
correctly written test that mocks a same-named LOCAL. The failure mode
is loud, immediate and has a distinctive message ("unused argument"
naming a quoted function), so the cost of discovering it is one test
run. `dev/lane-rules.md` is the right home, in one line beside the
existing active-binding entries.

## 6. The new load failure, and the floors

Reproduced independently (`dev/sgrev-shadow-build.R`,
`dev/sgrev-shadow-load.R`). Building and loading must be in separate
processes: my first spelling read the real package's exports in the
loading process, so the real owner was already loaded and the guard
"the shadow really lacks the name" read FALSE. The guard caught it,
which is the rule about constructing the absent case working as
intended.

```
LIB rellib-r3  owner posterior (no rhat)  shadow-first   sample: loaded
LIB rellib-r3  owner posterior (no rhat)  sample-first   sample: loaded
LIB rellib-r3  owner gratia (no posterior_samples) shadow-first  loaded
LIB rellib-r3  owner gratia (no posterior_samples) sample-first  loaded
LIB sgrev-lib  owner posterior (no rhat)  shadow-first
    frmtmb.sample: FAILED: object 'rhat' not found
LIB sgrev-lib  owner posterior (no rhat)  sample-first   sample: loaded
LIB sgrev-lib  owner gratia (no posterior_samples) shadow-first
    frmtmb.sample: FAILED: object 'posterior_samples' not found
LIB sgrev-lib  owner gratia (no posterior_samples) sample-first  loaded
```

Every cell asserts that the SHADOW is what got loaded before reading
the result. Both halves of the lane's defect 1 hold exactly.

Both floors check out upstream. gratia's NEWS puts `posterior_samples()`
in **0.9.0** as a new function; posterior's NAMESPACE at the **v1.0.0**
tag carries `export(rhat)`. So `gratia (>= 0.9.0)` and
`posterior (>= 1.0.0)` are the right numbers, and the lane's reading
that posterior's exposure is academic and gratia's is not is correct.

**Pricing the alternative.** Dropping gratia as an owner costs a
gratia user, in the order `library(gratia); library(frmtmb.sample)`,
the methods gratia has for this name: `posterior_samples.gam`,
`posterior_samples.scam` and `posterior_samples.default`. The lane's
"3 of 3" counts the default; my table says 2 of 2 because it excludes
defaults by construction. Either way it is every method gratia has. In
the other direction it also costs the twin, so `posterior_samples(ds)`
falls into gratia's default and reports
"Don't know how to sample from the posterior of <frmtmb_draws>FALSE",
which I measured in section 2. Keeping gratia is the better trade and
the lane is right to put the choice to the user rather than to make it.

## 7. The two weak-form claims

**The core block's behavioural failure, constructed.**
`dev/sgrev-coremutant.R` builds a core that TAKES `owners` and does
not validate it, one substitution, asserted to match exactly once,
installed into its own library. `dev/sgrev-runtest.R` runs the
unmodified core test file against it:

```
FIX core         BLOCKS 13  PASS 56  FAIL 0  ERROR 0
MUTANT core      BLOCKS 13  PASS 47  FAIL 5  ERROR 1
  failing: frm_install_generics() takes a table and refuses a bad one
     Expected frm_install_generics("frmtmb", owners = b) to throw a error.  (x5)
     Error in exists(gen, envir = ns, inherits = FALSE): invalid first argument
```

So the block does test the refusal and not the signature. Nit: record
this mutant on the page; the weak form is a page defect, not a test
defect, and the fix is two sentences and a script.

**The split suite ran against one build.** Timestamps, and no install
falls inside the window:

```
installed core in samplegen-lib   R/frmtmb.rdb      16:53:16
installed frmtmb.sample           R/frmtmb.sample.rdb 16:54:03
newest core source edit                             16:52:53.33
newest frmtmb.sample source edit                    16:54:01.34
core suite outputs                        17:14:51 .. 18:33:52
    one gap of 4032 s ending 18:27:26  (the session restart)
frmtmb.sample suite outputs               16:58:45 .. 17:07:22
test-loo.R edited                                   17:05:30.80
test-loo.txt written                                17:07:22.57
```

Both installs precede every suite output, every source edit precedes
its install, and the one file edited mid-suite has a result written
after its edit. The claim stands as recorded.

## 8. The body guard's mutant, and one of mine

The lane's mutant is real: `dev/samplegen-mutant.R` copies the
package, asserts the target string occurs exactly once, rewrites
`stancode()`'s generic to `message()` before dispatching, and installs
it. A guard that only catches the generic its author mutated is not a
guard, so I built a second one on a DIFFERENT generic in a DIFFERENT
source file (`dev/sgrev-samplemutant.R`, `log_lik` in `R/loo.R`):

```
FIX sample                BLOCKS 13  PASS 58  FAIL 0  ERROR 0
MUTANT log_lik            BLOCKS 13  PASS 57  FAIL 1  ERROR 0
  failing: with no owner usable, it loads and its generics are bare
     Expected parse_field(out, "WITHWORK") to equal character().
     actual: "log_lik"
```

The guard names it. The block also carries its own positive
conditions, `OWNERSLOADED` empty and `NOTOURS` empty, so it cannot
pass by reading nothing.

Incidentally, an early misrouted run of mine put the shared user
library in front and so ran the sample collision file against the
RELEASED frmtmb.sample: 13 blocks, 58 assertions, 41 pass, 17 fail.
That is the lane's BASE row to the assertion, reproduced by accident
and from a different library.

Two affected test files rerun on my own build, one process each:
`test-loo.R` 20 blocks, 79 pass, 0 fail, 0 error; `test-draws-methods.R`
19 blocks, 101 pass, 0 fail, 0 error. Both match the lane's per-file
counts exactly.

## 9. The new export, the Rd page, and bbmle

`dev/sgrev-api.R`, `dev/sgrev-api2.R`.

```
  frm_install_generics   exported: TRUE   in namespace: TRUE
  frm_bind_generic       exported: FALSE  in namespace: TRUE
  frm_adopt_target       exported: FALSE  in namespace: TRUE
  frm_is_generic         exported: FALSE  in namespace: TRUE
  frm_generic_owners     exported: FALSE  in namespace: TRUE
  frmtmb.sample re-exports it: FALSE
  formals: function (pkgname = "frmtmb", owners = frm_generic_owners)
```

Only the one entry point is public, the helpers stayed internal, and
the extension reaches it through `import(frmtmb)` rather than by
re-exporting it. The Rd page RENDERS, the "borrowed-generic seam"
section is in the rendered text in full, it carries no unescaped
percent sign, it names `frm_install_generics(pkgname, owners)`, and it
does not leak the internal helpers. I read the rendered section against
the code and it is accurate on all three rules and on why the call has
to be first in `.onLoad()`.

Two nits on the newly public function.

* Misuse off `.onLoad()` gives a raw R message with no direction:
  `frm_install_generics("stats", list(sd = "stats"))` says
  "cannot remove bindings from a locked environment". The page says
  when to call it; the function does not.
* The `frm_is_generic()` refusal is applied to the OWNER's function
  and not to the caller's own fallback, so an extension that puts a
  non-generic of its own in the table gets a binding installed around
  it with nothing said. Low risk, since the caller writes the table,
  but the validation block is one line from covering it.

**bbmle.** Confirmed, in the order defect 4 names:

```
search head: .GlobalEnv package:frmtmb.sample package:bbmle
parnames resolves to: frmtmb.sample   formals (x, ...)
plain parnames(x):        ERROR "no applicable method"
bbmle::parnames(x):       works
parnames(x) <- c(...):    works (the replacement name is not masked)
bbmle::parnames(x):       [1] "a" "b" "c"
```

`::` suffices, as it did for core's non-generic collisions, and the
replacement form `parnames<-` is not masked at all, so the only thing
a bbmle user loses is the bare getter at the prompt. The binding
correctly refuses to adopt a non-generic, so this is unchanged from
BASE rather than caused by the lane. No action.

## 10. Load cost. Not falsified, and the instrument has one gap

The design is right: interleaved arms, order shuffled within a round,
one fresh process per measurement, a high-resolution clock with its
tick reported, paired bootstrap intervals, and a control. What it does
not carry is a POSITIVE control, an arm with a known non-zero effect
that must come back at its known size. A null with no positive control
says the effect is smaller than the interval, and nothing about
whether the instrument could have seen it.

I added one (`dev/sgrev-loadcost.R`, 30 rounds, 150 fresh processes,
seed 20260915):

```
Sys.time() tick: min 9.54e-07 s  median 2.15e-06 s
arm                           min   median
sample_BASE                0.2231   0.3480
sample_FIX                 0.2161   0.3268
core_BASE                  0.1896   0.2849
core_FIX                   0.1860   0.2613
POSCTRL_BASE_plus_10ms     0.2440   0.3454

paired per-round difference, median [95% bootstrap interval]
sample_FIX - sample_BASE              -0.0132 s [-0.0594, +0.0042]
core_FIX - core_BASE                  -0.0076 s [-0.0192, +0.0039]
POSCTRL_BASE_plus_10ms - sample_BASE  +0.0260 s [-0.0166, +0.0438]

on the MINIMA over 30 rounds
sample_FIX - sample_BASE              -0.0069 s
core_FIX - core_BASE                  -0.0036 s
POSCTRL_BASE_plus_10ms - sample_BASE  +0.0209 s
```

**My positive control is itself not calibrated, and I say so rather
than quoting it.** `Sys.sleep(0.010)` on Windows returns on the system
timer quantum, about 15.6 ms, so the arm's true effect is not 10 ms
and +0.0209 on the minima is consistent with the quantum plus noise.
What the run does establish is that an effect of roughly 20 ms is
visible in the minima on a machine this noisy, while
`sample_FIX - sample_BASE` comes back NEGATIVE at -0.0069 s. My
medians are 60 percent larger than the lane's on the same arms, so my
machine was busier; the lane's interval of about 6 ms either side is a
claim about a quiet machine and I cannot reproduce that width.

So: the finding "inside the noise" is not falsified, and the residual
after core's own +0.0210 s (the generics lane's number) is not
separable from run-to-run spread. Nit: say on the page that the
interval width is a property of machine state, and add a positive
control arm the next time the script runs. A `for` loop burning a
measured 10 ms of CPU is better than `Sys.sleep()` on this box.

## 11. The version claim, measured

```
released core version: 0.56.0
released core exports frm_install_generics: FALSE
lane core exports it:                       TRUE
lane frmtmb.sample against released core:
  FAILED: .onLoad failed in loadNamespace() for 'frmtmb.sample':
  call: frm_install_generics(pkgname, sample_generic_owners)
  error: could not find function "frm_install_generics"
```

Confirmed by construction. The requirement is hard, not advisory.

## Fixes

### BLOCKER 1. Widen defect 2 and put the numbers on it

`dev/samplegen-findings.md`, defect 2. It currently names `rhat` only
and says the size was not measured. Replace it with: `rhat` and
`neff_ratio` both delegate to bayesplot on the `stanfit` where brms
uses posterior on the draws; the measured relative difference is up to
0.0056 for `rhat`, which is 1.15 times the whole excess over 1, and up
to 0.41 for `neff_ratio`; `rhat(ds)` equals rstan's classic split-Rhat
exactly, which names the definition; `nuts_params()` and
`log_posterior()` are clean because brms's own methods delegate to
bayesplot too; and both functions return raw Stan parameter names, so
4 of 11 entries do not match `variables(ds)`. Construction:
`dev/sgrev-rhat.R`, `dev/sgrev-rhat2.R`, draws cached at
`dev/stan-cache/sgrev-draws.rds`, seed 20260915.

This is a record fix. Whether to change the two methods is the user's
call and belongs in the same entry, not this lane's work.

### BLOCKER 2. The floor in DESCRIPTION

`extensions/frmtmb.sample/DESCRIPTION` still reads
`frmtmb (>= 0.55.1)` while `.onLoad()` calls an export no released
core has. The lane says the consolidator sets it, and that is fine,
but it must not merge with the current number: the failure is a hard
load failure and nothing in the tree states the requirement.

### Nits

1. Add `getS3method()` as a column to the frmtmb.sample collision
   table, as core's has. Measured: 67 of 67 lost in BASE-S, 0 on FIX.
2. Put the draws-side output run on the page. 46 calls, five orders,
   BASE against FIX, 0 of 46 differ beyond the two brms deprecation
   warnings NEWS declares.
3. Give the two BASE silent failures their user-visible strings:
   `rhat(ds)` errors with "unimplemented type 'list' in 'greater'" and
   `posterior_samples(ds)` reports gratia's
   "Don't know how to sample from the posterior of <frmtmb_draws>FALSE".
4. Replace defect 3's two examples with the full ten-row table, with
   the disagreeing position named, and say that eight error and two
   return quietly.
5. Say that the `owners` denominators (77, 41, 3) count per-owner
   table entries, so a reviewer counting (name, class) pairs gets 39,
   29 and 2 and has not found a disagreement.
6. Finish the mock comment in `test-loo.R`: `local_mocked_bindings()`
   also patches this package's own method table but not the owner's,
   so with the owner loaded only a caller inside this namespace sees
   the mock.
7. Add the mocking trap to `dev/lane-rules.md` as a rule, not a guard.
   A guard over 53 names would fire on a correct test that mocks a
   same-named local, and the real failure is loud and one run away.
8. Record the core mutant that produces the behavioural failure for
   the new core block, so "the weak form" stops being the last word.
9. Load cost: say the interval width is a property of machine state,
   and carry a positive control next time. Do not use `Sys.sleep()`
   for it on Windows; the timer quantum is about 15.6 ms.

## What I could not falsify

* The 0 of 28 claim, in fourteen load orders including four the lane
  never ran, on four channels, with an independently derived name list
  and a per-order control process.
* The output claim, widened from 22 calls to 60 and from two channels
  to four, in four FIX arms.
* The three rules the Rd section states, checked against the installed
  build rather than against the lane's test.
* The formals table. brms agrees with the owner on all six.
* The partial-owner load failure, both halves, both orders, both
  builds, with the shadow's absence asserted in every cell.
* The body guard, against a mutant of my own on a different generic in
  a different file.
* The split-suite provenance, by timestamp.
* The `frm_install_generics()` export surface and the Rd rendering.
* The bbmle finding, including that `::` and `parnames<-` both work.
* The load-cost null, subject to nit 9.

## Scripts

| script | what it produced |
|---|---|
| `dev/sgrev-names.R` | the 28 names and their bodies, off the BASE install |
| `dev/sgrev-probe.R`, `dev/sgrev-runprobe.sh`, `dev/sgrev-summary.R` | the collision table, 56 processes |
| `dev/sgrev-dispatchcheck.R`, `dev/sgrev-stackdump.R` | the probe checked against real dispatch |
| `dev/sgrev-leg1.R`, `dev/sgrev-leg1b.R` | what the caller-environment leg reaches |
| `dev/sgrev-brmsout.R`, `dev/sgrev-brmsout-summary.R` | 60 calls on a brmsfit, four channels |
| `dev/sgrev-drawsout.R`, `dev/sgrev-drawsout-summary.R` | 46 calls on draws, BASE against FIX |
| `dev/sgrev-rhat.R`, `dev/sgrev-rhat2.R` | the R-hat and ESS definitions and numbers |
| `dev/sgrev-owners.R` | brms's generic and method beside every definer's |
| `dev/sgrev-positional.R` | the positional audit, all 28 |
| `dev/sgrev-twins.R` | the three rules against the installed build |
| `dev/sgrev-mock.R`, `dev/sgrev-mock2.R` | the active-binding trap and the mock's real mechanism |
| `dev/sgrev-shadow-build.R`, `dev/sgrev-shadow-load.R` | the partial-owner load failure |
| `dev/sgrev-coremutant.R` | core that takes `owners` and does not validate it |
| `dev/sgrev-samplemutant.R` | work put into `log_lik()`'s generic |
| `dev/sgrev-runtest.R` | one test file per process, errors counted |
| `dev/sgrev-api.R`, `dev/sgrev-api2.R` | exports, the rendered Rd, bbmle, positional outcomes |
| `dev/sgrev-loadcost.R`, `dev/sgrev-loadcost2.R` | load cost with a positive control |
| `dev/sgrev-floor.R` | the load failure against core 0.56.0 |

Outputs are under `dev/sgrev-out/`. The draws object is
`dev/stan-cache/sgrev-draws.rds` (gitignored), 4 chains of 500 on a
120-row Gaussian mixed model, seed 20260915.
