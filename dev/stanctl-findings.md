# stanctl lane: items 1.1 and 1.2 of Phase 1

Date: 2026-09-08. Worktree `frmtmb-wt-stanctl`, branch `wt-stanctl`,
based on 780dec1. Nothing committed.

Versions under test: frmtmb 0.55.0, frmtmb.sample 0.3.2,
frmtmb.ode 0.1.1, RTMB 1.9, RTMBode 1.0 (commit 5242257), tmbstan 1.2.0,
rstan 2.32.7, deSolve 1.42, R 4.6.1 on Windows 11.

---

## 1. Item 1.1: the sampler control list

### 1.1 What the plan asked for, and why it was not built that way

The plan's row asked for a second argument, `stan_control`, so that
`control` could stay `frmtmb_control()` and a brms user's
`control = list(adapt_delta = 0.99)` would keep being refused by name.
That was overruled during the lane and the row in
`dev/extension-gaps-plan.md` is edited to say so. The reasons, each
checked rather than taken on trust:

1. **brms compatibility is a locked goal.** SPEC.md section 5: "brms
   code ports mechanically, priors included", and `frm_sample()`
   already spells `prior` as brms does for exactly this reason. In
   brms, `control` IS the Stan control list.
2. **The package already documented the collision as a limitation.**
   `vignettes/brms-posterior.Rmd` line 54 carried the row
   "`control = list(adapt_delta = )` | not through `frm_sample()`,
   whose `control` is `frmtmb_control()`; see below", and
   `vignettes/posterior-diagnostics.Rmd` sent a user with divergences
   to `as_tmbstan()`. A `stan_control` argument would have left both
   rows standing with a footnote. Taking brms's spelling removes them.
3. **`as_tmbstan()` already means `control` the brms way.** It has no
   `control` formal at all; `...` goes straight to tmbstan. So the
   rename makes the two entry points of this package agree, rather
   than introducing a third convention.

The other option the coordinator raised, putting `adapt_delta` and
`max_treedepth` into `frmtmb_control()`, was rejected for reasons that
also check out: SPEC.md section 3 fixes core's Imports as "No Stan, no
compiled code of our own, no brms", `frmtmb_control()` validates every
field eagerly in core (`check_flag`, `check_count`, `check_positive`),
so core would be validating options for a sampler it does not know
about and tracking rstan's option set; and `frm()` would then accept
`adapt_delta` and silently do nothing with it. The two also apply at
different times: `frmtmb_control()` is fit-time, the sampler's control
is sample-time.

Recorded here so nobody re-proposes `stan_control`.

### 1.2 What was built

`frm_sample()` now takes:

- `control = NULL`, the SAMPLER's control list, passed to
  `tmbstan::tmbstan()` and so to `rstan::sampling()` unchanged.
- `fit_control = frmtmb_control()`, what `frm()` calls `control`. Used
  only on the formula route, as before.

`as_tmbstan()` is unchanged: it never had the collision.

### 1.3 What it used to do, measured

`frm_sample(fit, control = list(adapt_delta = 0.99))` on the FIT route
was **silently ignored**. The fit route never reads `control` at all
(it is forwarded only by `sample_assemble()`, which the fit route does
not call), so the call sampled with rstan's defaults and said nothing.
Measured on the unfixed build: `ds$stanfit@stan_args[[1]]$control` was
`NULL` where 0.97 had been asked for.

On the FORMULA route the same call did error, but with `frm()`'s
message, "`control` must come from frmtmb_control(); this list is
missing optimizer, optCtrl, restarts, grad_tol", which names neither
`adapt_delta` nor any route to it.

### 1.4 How the two lists are told apart

`frmtmb_control()` returns a PLAIN unclassed list, so there is no class
to dispatch on. The refusal keys on field names instead, and that is
decidable because the two vocabularies are disjoint:

| list | names |
|---|---|
| `frmtmb_control()` (13) | optimizer, optCtrl, restarts, grad_tol, profile, sparse_x, autoscale, check_nlev_1, check_olre, importance_seed, importance_rounds, importance_ess, verbose |
| rstan `control` (15) | adapt_engaged, adapt_gamma, adapt_delta, adapt_kappa, adapt_t0, adapt_init_buffer, adapt_term_buffer, adapt_window, stepsize, stepsize_jitter, metric, int_time, max_treedepth, epsilon, error |

Intersection: empty. `test-stan-control.R` asserts that it stays empty,
since both vocabularies belong to other packages and can grow.

The rstan list is copied from `rstan:::config_argss()`, found by
searching the rstan namespace for the function carrying the string
"unknown members".

### 1.5 Why the unknown-option check is ours and not rstan's

rstan does validate `control` names, but it does not raise. It prints
`'control' list contains unknown members of names: nonsense_option`,
then `error in specifying arguments; sampling not done`, and RETURNS AN
EMPTY stanfit. Measured directly through `tmbstan::tmbstan()`. An empty
stanfit lands on `frm_sample()`'s existing "the sampler returned no
draws" refusal, whose text is about a tape calling an external solver,
so a typo in `control` would have been reported to the user as a solver
failure. Hence `check_stan_control()` refuses an unknown name here.

### 1.6 Can `control` arrive twice, through `...` and by name?

Not twice, but the guard was bypassable and the first version of this
section drew the wrong conclusion from a true premise. It is true that
`control` follows `...` in the formals, so R requires exact name
matching there and `control =` always binds to the formal, never to
`list(...)`. The reasoning stopped one step early: an ABBREVIATION
lands in `...`, and `...` is handed to `tmbstan::tmbstan()` and then to
`rstan::sampling()`, whose `control` argument comes BEFORE that
method's own `...`, so R partial-matches it there.

Measured on the build before this fix:

| call | what reached the sampler |
|---|---|
| `frm_sample(fit, control = list(adapt_delta = 0.97))` | adapt_delta 0.97, after `check_stan_control()` |
| `frm_sample(fit, contro = list(adapt_delta = 0.97))` | adapt_delta 0.97, checked by nothing |
| `frm_sample(fit, contro = list(nonsense = 1))` | rstan declined, and the call reported "the sampler returned no draws ... a tape that calls an external solver" |
| `frm_sample(fit, contro = frmtmb_control())` | the same solver mis-report |

The last two are the exact mis-report `check_stan_control()` exists to
prevent, reached by mistyping the argument it guards.
`refuse_partial_control()` now refuses any dots name that `pmatch()`es
`"control"`, with `duplicates.ok = TRUE` so a second abbreviation in
the same call is not silently dropped by pmatch's own de-duplication.
Verified after the fix: `contro`, `cont` and `c` are each refused by
name, and `control` still reaches rstan.

The hole was only ever the partial-match set of rstan's own formals.
rstan already rejects a dots name that is not a prefix of one of them
("passing unknown arguments: fit_contro"), and `control` is the one
member of that set this package validates itself.

### 1.8 The arguments that were accepted and discarded

`fit_control`, `start`, `data2`, `na.action` and `REML` are formula
route arguments: they assemble a model, and a fitted object is already
assembled. On the fit route they were read by nothing, and nothing was
validated either, so complete nonsense was accepted in silence. `data`
and `family` were already refused there by name.

Four of the five predate this lane. `fit_control` does not, and that
is why all five are now closed rather than reported: the rename fixed
"a sampler list bound to the fit-time argument is discarded" and, in
the same edit, would have created "a sampler list bound to the NEW
fit-time argument is discarded", under the very name the `control`
refusal points a confused user at. `frm_sample(fit, fit_control =
list(adapt_delta = 0.99))` sampled four chains with rstan's defaults
and said nothing.

All five now refuse by name on the fit route, using `missing()` so that
passing the default explicitly is still passing it, and the message
names `control` for the caller who meant the sampler. The formula
route is untouched.

### 1.7 The check the plan asked for: adapt_delta on a centered funnel

Construction: `y ~ 1 + (1 | g)`, gaussian, 6 groups of 3 observations,
true group sd 0.05 against residual sd 1, data seed 2026;
`reparameterize = FALSE` so the posterior stays centered (the default
non-centering removes the funnel and with it the thing being measured);
brms default priors; 1 chain, 1000 iterations, chain seeds 11 to 18.

adapt_delta 0.80 against 0.99, per chain seed:

| seed | div @0.80 | div @0.99 | mean stepsize @0.80 | @0.99 |
|---|---|---|---|---|
| 11 | 0 | 0 | | |
| 12 | 0 | 0 | | |
| 13 | 1 | 0 | | |
| 14 | 4 | 0 | | |
| 15 | 0 | 0 | | |
| 16 | 42 | 0 | | |
| 17 | 13 | 0 | | |
| 18 | 15 | 0 | | |
| **total** | **75** | **0** | mean 0.2517 | 0.0809 |

The step size falls for every one of the eight seeds individually; the
divergence count is heavy-tailed across seeds (42 of the 75 come from
one), so it is compared as a sum. Both comparisons are between the two
arms of the SAME runs, so no threshold from another machine enters.
`nuts_params()` is what reads both quantities back, which is the
plan's second check.

The construction was chosen by measurement, not by eye. Four
alternatives were run over the same seeds; every one moved the counts
in the same direction, and this one has the widest margin:

| construction | div @0.80 | div @0.99 | stepsize | per-seed monotone |
|---|---|---|---|---|
| 6 groups x 3, tau .05, sigma 1, data seed 2026 | 75 | 0 | 0.252 to 0.081 | yes |
| 6 groups x 3, tau .05, sigma 1, data seed 7 | 5 | 1 | 0.296 to 0.068 | yes |
| 5 groups x 3, tau .05, sigma 1 | 3 | 0 | 0.124 to 0.0015 | yes |
| 6 groups x 2, tau .05, sigma 1 | 4 | 0 | 0.101 to 0.0009 | yes |

Two constructions with a FLAT prior were tried first and rejected: they
produce more divergences (24 and 23 over four seeds) but the adaptation
becomes erratic and the step size is no longer monotone per seed, which
is the mechanism the test is really about. A flat prior leaves the
variance parameter's tail open at zero, which is what makes both the
funnel worse and the adaptation less reliable.

The first construction tried, the one the test carried before this
search, gave 2 divergences against 1 over four seeds. That is a real
direction but not an assertion worth shipping, and it was replaced.

The assertion is gated behind `sampler_gates_on()`, the switch this
package already uses for chain-agreement gates, because a seeded Stan
chain is not platform-deterministic.

---

## 2. Item 1.2: the pre-flight refusal

### 2.1 The ODE breakage still reproduces on 0.55.0

Confirmed before building anything on top of it, with the worktree's
0.55.0 / 0.3.2 / 0.1.1 installed into the lane's private library.

Model: one-compartment oral PK, `frm_ode()` in a nonlinear body, 6
subjects x 8 timepoints, **no random effects at all** (the simplest ODE
tape there is; `dev/feature-gaps.md` recorded that "whether single-node
ODE fits also fail is untested", and they do).

| call | result |
|---|---|
| `frm()` | converges, logLik -47.1102 |
| `obj$fn(mode)` | 47.1102 |
| `obj$gr(mode)` | 6.9e-06, -3.5e-05, -6.5e-05, -4.0e-06 |
| `obj$fn(mode + 1)` | 156.3939 |
| `obj$fn(mode + 10)` | 503.1102 |
| `obj$fn(mode + 100)` | 4823.11 |
| `frm_sample(fit, chains = 1, iter = 30)` | **fails in 2.96 s** |

The failure output is the recorded signature verbatim: DLSODA warns
`Internal T (=R1) and H (=R2) are such that in the machine, T + H = T`,
then `DINTDY- T (=R1) illegal` at 0.25 and 0.5, then
`Error in lsoda(...): illegal input detected before taking any
integration steps`, then `error occurred during calling the sampler;
sampling not done`, and `frm_sample()` reports its "returned no draws"
refusal.

So the item does not change shape. Repro script:
`scratchpad/stanctl-ode-repro.R`.

### 2.2 The upstream patches are still current

`RTMBode` installed here is version 1.0 with
`RemoteSha: 52422571529ba0a02341b6fbbbf0f9b01813df3e`, which is the
exact commit `dev/upstream/rtmbode-issues.md` names (`5242257`), built
by r-universe on 2026-08-31, two days before the patch work. It is the
only build available. Checked function by function against the
installed namespace:

| patch | what it adds | present in the installed build? |
|---|---|---|
| 01 solver-failure-nan | `tryCatch` around the deSolve call in `ODEadjoint()`, NaN fill, row mapping | no: no `tryCatch` anywhere in `ODEadjoint`, no NaN fill |
| 02 events-ad-path | `names(y)` carried through `addInfo()`/`augment()`, `replace`/`multiply` refused | no: no state names carried, no method refusal |
| 03 augmented-workspace-guard | `checkWorkspace()` | no: no such object, no `lrw` mentioned anywhere in the namespace |

The namespace's object list is unchanged from the write-up:
`addInfo`, `augment`, `desolve_derivs`, `func2tape`, `ode`,
`ODEadjoint`, `set_parms`, `set_pointers`, `setTape`. All three patches
apply and none is superseded. Probe: `scratchpad/stanctl-rtmbode-check.R`.

**Nothing was filed and nobody was contacted.** Section 6 below is the
filing-ready text, for the user to send or not.

### 2.3 What was built

`frmtmb.ode/R/zzz.R` registers, from `.onLoad()`:

- the feature `c("frm_ode()" = "special")`;
- one rule, `frm_ode() x frm_sample` at status `refused`, with a note
  that says what the failure is and what to do instead;
- `expects = "frm_sample"`, because `frm_sample` is frmtmb.sample's
  feature and frmtmb.ode neither depends on nor suggests it.

`frmtmb.sample/R/sample.R` gains `sample_preflight()`, called from
`frm_sample()` immediately after the route is determined and BEFORE
`sample_assemble()` on the formula route or any retape on the fit
route. It reads `frm_compat("frm_sample", status = "refused")` and, for
each row, decides whether the model in front of it uses the other side.

### 2.4 What the pre-flight can decide, and what it will not guess at

This is the part worth arguing with, so it is stated plainly. It has
now been wrong twice and corrected twice, both times by construction
rather than by argument, and both corrections are kept below because
the shape of the error is the useful part.

A registry row is a PAIR. Reading the refused rows is easy. Deciding
whether the model uses the other side is the hard half, because **what
identifies a feature in a model is WHERE its call sits in the formula.**
`se` inside the bar on the response's left is the `se()` addition term;
`se` in a nonlinear body is whatever function the user has of that name.

**The first error: matching the registry KEY.** Three keys are shared
by two display names each, deliberately (`R/compat.R`: "The special
vocabulary collides with the covariance vocabulary: gp, cs and mi name
both a bar-term structure and a predictor special"). A refused
`mi_pred()` row fired on `bf(y | mi() ~ x)`, which writes the ATERM.
Constructed and reproduced.

**The second error: claiming the repair was SOUND.** The repair added
two tests, name uniqueness and a base R test, and I wrote that display
name plus those two gives the largest set decidable soundly without a
core seam. That is false, and the counterexample is cheap: `se()` has a
unique name and is not base R, so it took the name route; with the
ordinary helper `se <- function(v) v / 2`, the model
`bf(y ~ a * se(x), a ~ 1, nl = TRUE)` FITS, at logLik -4.438, and a
refused `se()` row refused it. The two tests remove two ENUMERABLE
classes of collision, a name two features share and a name base R also
uses. They do not touch the class the thesis itself names. `trunc()`
was caught only because it happens to live in base.

**What it does now: three routes, and only two of them are sound.**

1. **family**, matched by name against the response families. Exact: a
   family is named, never called. It does not decompose `mixture()`,
   whose family name is the whole string
   `"mixture(gaussian, gaussian)"`, so a refused `gaussian` row misses
   a mixture of gaussians. A miss, not a false alarm.
2. **aterm**, matched by POSITION. An addition term is written inside
   the bar on the left-hand side of the response formula and nowhere
   else, and the `bf()` object carries that, so this needs nothing from
   core and is sound. Twelve lines. It separates `bf(y | se(v) ~ x)`
   from `bf(y ~ a * se(x), ...)`, which is exactly the counterexample
   above, and it covers 9 of the 20 call-route rows: `trials`, `cens`,
   `se`, `vint`, `vreal`, `dec`, `reward`, `payoff`, `stage2`.
3. **everything else the vocabulary writes as a call**, matched by
   NAME, on the display name, over calls OUTSIDE aterm position, and
   only when the name is not a base R function.

Route 3 is the largest set justifiable FROM NAMES ALONE, and that is
the honest claim; it is not sound. `mo()`, `ma()`, `cosy()` and `mm()`
pass every test it applies and would each fire on a user's function of
that name in a nonlinear body, and a qualified `somepkg::s(x)` matches
`s()` because the walker drops the namespace. Route 2 removed four of
the eight names review found this way; these four are what is left.
Exposure today is zero, because a refusal has to be registered on one
of them first and the only refused row is `frm_ode()`.

**Why the display name and not the key, and why that is what makes the
guard un-revocable.** Keying on the parser name forced a downgrade to a
warning whenever a key was shared, and any package could manufacture
that: registering a bare `frm_ode` display name, which core accepts,
downgraded frmtmb.ode's refusal to a warning and the model then
sampled, which is the thing that leaves rstan unusable for the rest of
the session. A guard a third party can switch off is not a guard.

Matching the display name closes it, and the closure is a PROOF, which
is better than the measurement I would otherwise have shipped:

- core refuses a second registration of a display name, verified:
  `frmtmb_register_compat(features = c("frm_ode()" = "grammar"))` gives
  "the registry already holds that name under the kind 'special'. One
  display name carries one KIND". The kind is as load-bearing as the
  name, because the route is chosen from both, and core's refusal
  covers both. Review checked this against `R/compat.R` rather than one
  message: `compat_new_features()` sits on both public entry points and
  there is a read-time net besides, since `compat_features_build()`
  ends in `out[!duplicated(out$name), ]`, so display names stay unique
  even when the registration guard is bypassed;
- dropping a fixed `()` suffix is injective over the names that carry
  it;
- a BARE name cannot reach the call route. **This is a property of the
  CODE, not of the vocabulary**, and the first version of this document
  had it the wrong way round. A formula DOES write bare vocabulary rows
  as calls: `ar1` and `gr_cov` are bare, and `formula_calls()` on
  `y ~ x + (1 | ar1(t, g))` returns `ar1`. What stops them is that
  `preflight_route()` returns the without-parentheses sentence before a
  bare name can reach the call route at all. As a claim about the
  vocabulary it is falsifiable in one line; as a claim about the code
  it holds.

Therefore two call-shaped rows cannot want the same call name, the
uniqueness test is unreachable, and it was DELETED rather than shipped
as dead code. Measured alongside the proof: 0 duplicate display names
over 111 vocabulary rows, 0 duplicate stripped names over the 23
call-shaped ones, and the bare-`frm_ode` attack, performed for real
against the live registry rather than a hand-built table, now leaves
the route at `call` and the model still refused end to end.

What premise 1 costs is the other direction, recorded rather than
fixed: because the registry cannot be silently taken over, it can be
loudly DEADLOCKED. A package that claims `frm_ode()` under a different
kind BEFORE frmtmb.ode loads makes frmtmb.ode fail to load at all
(`.onLoad failed in loadNamespace()`). That is pre-existing core
behavior and it is the loud direction, which is the right one.

What matching the display name costs is `mi_pred()`, `gp_pred()` and
`cs_pred()`, whose display names no formula writes. They do NOT miss in
silence, which was the shape of the first version of this fix and is
the failing-open shape the project keeps closing: the registry itself
says which rows are affected, since its key is the parser name, so a
call-shaped row whose display name and key differ is routed to the
warning. The three are exactly those rows, asserted as a set.

**The seam, now scoped to three kinds not four.** `frm_features_used(spec)`
is filed in the plan's Core seams table. Route 2 closes `aterm` without
it, so the seam is needed for `special`, `autocor` and `grammar`. Core
owns the parser and can say exactly which vocabulary rows a model uses;
with it, every consulting caller could refuse exactly and warn never.

### 2.5 The checks the plan asked for

**"An ODE fit refuses in under a second with no Stan call."** Met, and
met unqualified once the pre-flight was moved.

Timed with `Sys.time()`, whose resolution here is sub-millisecond (an
empty closure times at 7 to 16 microseconds, the control), not
`proc.time()`, whose 10 ms tick would have been too coarse.

The first version ran the pre-flight AFTER `requireNamespace("tmbstan")`
and `requireNamespace("rstan")`, so a refused model paid the namespace
loads before being refused, and the first call of a session cost 2.48 s
here and 3.05 s on the reviewer's machine, which is SLOWER than the
2.88 to 2.96 s failure it replaces. The sentence "so a refused model
costs one registry read" was therefore not earned. The pre-flight needs
neither namespace, so it now runs before both, and a model that cannot
be sampled at all is a more useful thing to be told than which package
is missing.

Measured after the move, one fresh process per row, three processes for
the first-call figure:

| measurement | seconds |
|---|---|
| first `frm_sample()` of a fresh session, refused ODE fit, two extensions loaded | 0.355, 0.405, 0.420 |
| the same, seven extensions loaded (review's library) | 0.580, 0.582, 0.582 |
| the same, warm, 10 repeats per process | min 0.081, median 0.093 to 0.173 |
| `sample_preflight()` alone, 30 repeats | min 0.080, median 0.091 to 0.127 |
| `base_r_function()`, one call, 200 repeats | 0.0001 |
| CONTROL, empty closure, 200 repeats | 0.000007 to 0.000009 |
| the unguarded failure it replaces, fresh session | 2.96, plus a session that must be restarted |

The two first-call rows differ by which extensions the session had
loaded, not by machine noise: the pre-flight's cost is the registry
read and the registry grows with the ecosystem. Section 2.6 has that
dependence measured. Both are under a second, so the plan's check holds
either way, and the refusal is 5x to 8x faster than the failure it
replaces rather than 5 percent slower.

`rstan` is verifiably not in `loadedNamespaces()` after the refusal,
which is the check that the move actually took. The refusal is now
about 7x faster than the failure it replaces rather than 5 percent
slower, and "under a second" holds on the first call and not only warm.

"No Stan call" is asserted STRUCTURALLY rather than by a clock. The
end-to-end test calls `frm_sample(form)` on the formula route with
`data` NOT supplied. The formula route refuses a missing `data` the
moment it tries to assemble, so receiving the ODE refusal instead
proves the pre-flight ran before assembly, and therefore before any
tape was built. That test failed on the unfixed build with exactly the
"needs data" message, which is recorded below. The reviewer's caveat is
right and worth repeating: that proves ordering relative to assembly,
not that nothing before it calls Stan; the latter follows from reading
the statements in between, and from `rstan` not being loaded at all.

**"The eam and latent rows stay `conditional`."** After the change the
only refused row is `frm_ode()`, and `wiener`, `hmm` and `lca` are
unchanged; `test-compat-preflight.R` asserts both.

The before-and-after counts depend on which extensions the session has
loaded, which the first version of this document failed to say. With
frmtmb, frmtmb.sample, frmtmb.ode, frmtmb.latent and frmtmb.eam loaded,
`frm_compat("frm_sample")` had 112 conditional, 4 works, 1 untested and
ZERO refused. With all seven extensions loaded the same query gives 131
conditional, 4 works, 1 untested and, after this change, 1 refused.
Either way the point stands: this is the registry's first refusal of
any kind for `frm_sample`.

### 2.6 The cost of the pre-flight, and why the round-1 judgment was
measured at the wrong end

Every `frm_sample()` and every `as_tmbstan()` pays the registry read.
Round 1 attributed it and called 0.12 s immaterial. The attribution was
right and the judgment was measured at the small end of something that
GROWS WITH THE ECOSYSTEM, which is the part that matters.

Attribution, warm, median of 30 repeats, in a session with frmtmb.sample
and frmtmb.ode loaded and nothing else, empty-closure control at 8
microseconds:

| | seconds |
|---|---|
| `frm_compat("frm_sample", status = "refused")` | 0.0897 |
| `frm_compat_features()` (cached) | 0.00001 |
| the formula walk on a bform | 0.00007 |
| `base_r_function()`, one name, memoized after the first | 0.0001 |
| `sample_preflight()` whole | 0.117 |

So the cost is the registry read and nothing else. But the registry
grows, and the pair table it resolves grows faster than the vocabulary.
Interleaved, five rounds, blocks of 20 calls, with a control arm built
from the same code that should report 1.0:

| session | vocabulary | rules | per call | control |
|---|---|---|---|---|
| frmtmb.sample + frmtmb.ode | 111 | 447 | 0.117 s | 0.84 |
| + latent, eam | 119 | 614 | 0.177 s | 1.05 |
| all seven extensions | 138 | 991 | 0.324 s | 1.04 |
| seven + 4 synthetic packages | 158 | 1091 | 0.508 s | 0.80 |
| seven + 8 synthetic packages | 178 | 1191 | 0.581 s | 1.00 |

The control arm sits at 1.00 to 1.05 in three rows and at 0.80 to 0.84
in two, so the instrument carries about 20 percent of noise on this
machine; the effect is 5x across the table and is not in doubt. The
synthetic packages register five features and twenty-five rules each,
which is well under the ~5 features and ~109 rules the seven real ones
average, so the bottom two rows UNDERSTATE what a fifteen-extension
ecosystem would cost. Even understated, the answer at fifteen
extensions is **at least 0.58 s per call**, five times the figure the
round-1 judgment was made on.

Restated: 0.12 s is what a TWO-extension session pays. At seven it is
0.32 s, which is 2.7 times the figure the round-1 judgment rested on,
and at a fifteen-extension scale it is 0.58 s or more, paid by every
call whether or not anything is refused.

I first wrote that the short-circuit should be taken before the eighth
extension lands. Review disagreed and is right on both counts. The
trigger is already met, since the seventh extension has landed. And a
package count is the wrong variable: the cost tracks the RULE count and
the pair table it resolves, so one extension registering 200 rules
would trip the real threshold without tripping a package-count one,
which is an adjective standing in for the quantity that matters.

The short-circuit is still declined, for the reason round 1 gave: a
second code path that can disagree with the first is a bad trade for a
refusal. The right fix is to make the ONE path cheap, and it is core's
to make. `frm_compat()` does not cache its resolved answer, so five
rounds of twenty identical calls each cost the full 0.32 s. Core
already has the discipline this needs and uses it for the vocabulary:
`compat_cache_store()` validates against `frmtmb_compat_contrib$features`
rather than a dirty flag, precisely so a caller that restores the list
cannot leave a stale flag behind. Caching the resolved per-feature
answer under that same key keeps one code path, keeps the call-time
semantics that make a package's refusal appear the moment it loads, and
takes the cost to near zero after the first call in a session.

Filed in the plan's Core seams table with these measurements and that
pointer. Not built here: it is a change under core's `R/`, which this
lane does not make.

### 2.7 The warn branch's false-alarm behavior

The branch is reachable through the real registry, not only through the
injected `rows =` the tests use: registering any `refused` row the
pre-flight cannot decide makes `frm_sample()` warn.

The property worth recording is the cost of choosing loud over silent.
The warning fires on EVERY `frm_sample()` and `as_tmbstan()` call for
the rest of the session, including on models that could not possibly
use the feature, because not being able to tell is exactly the
condition that triggers it. False-alarm rate on the registry as it
stands: zero, because the only refused row is one the pre-flight can
decide. False-alarm rate the day someone registers a refusal on a
covariance structure, a fitting mode, a base R name or one of the three
`_pred` specials: one warning per call, forever, on correct models.

That is a real cost of "loud rather than silent" and it belongs next to
the choice. It is also the strongest argument for the
`frm_features_used(spec)` seam: with it, none of these rows would need
a warning at all.

### 2.8 The core contract question

`?frmtmb_register_compat`'s status vocabulary said only that a
`refused` pair "is rejected with a message". It did not say who
rejects it, and until now nothing in either package read the registry
at run time: this pre-flight is the first consumer.

The first version of the added section opened "the pair raises an
error, and a row that says `refused` while the pair runs is a defect of
the same size as a missing guard", which is stronger than what the code
does: the warn branch is a constructed counterexample, a `refused` row
whose pair runs, by design, in the very caller the section holds up as
the worked example. The promise now reads: **a refused row is enforced
wherever the enforcing code can see the feature, and that code says
what it can see.** Hand-written guards see everything, so for them it
is unqualified; a consulting caller documents what it cannot decide and
warns there. The registry's behavior is unchanged.

### 2.9 `as_tmbstan()` no longer returns an empty stanfit

Recorded first as a defect left unfixed, then fixed on instruction. It
is a silent wrong answer, which this project ranks above a refusal and
above a missing feature.

rstan answers some failures by declining to sample, printing its reason
to the console, and returning the SHELL of a `stanfit`: an object of
the right class with `length(sf@sim) == 0`. `frm_sample()` has always
refused that. `as_tmbstan()` returned it, so a failed run looked like a
successful one until an accessor came back empty, and the caller had to
notice for themselves.

Both now call one helper, `check_stan_draws(sf, what)`, so there is one
spelling of the check and one message template, which is also what
`test-message-uniqueness.R` wants. The message names the caller and
says which of the two failures it is:

- pre-flight: `as_tmbstan() is refused on a model using frm_ode(): ...`
- no draws: `as_tmbstan(): the sampler returned no draws ... This is a
  failure of THIS RUN, not a model that cannot be sampled at all:
  those are refused before any chain starts, and
  frm_compat("frm_sample") lists them.`

They want different fixes, which is why the sentences differ rather
than sharing a generic one.

**The reachability question, answered.** `frm_sample()`'s check had
neither the ordering nor the partial-matching hazard that `control`
had. Measured rather than asserted
(`scratchpad/stanctl-reach-probe.R`):

| question | answer |
|---|---|
| statements between the tmbstan call and the check, in `frm_sample()` | 0 |
| any `return()`, `stop()` or `if` between them | no |
| the same in `as_tmbstan()` | `sf <- tmbstan::tmbstan(...)`, `check_stan_draws(...)`, `sf` |
| can `@` partial-match a slot name | no: `o@si` raises `no slot of name "si"` |
| names in rstan's `sim` | 12 |
| any of them that `pmatch("samples", .)` hits besides `samples` | none |

So the check is on every path that reaches a `stanfit`, and nothing was
reachable through the old `$samples`. The read was still converted to
`[[ ]]`, by the house rule rather than as a fix, and the reason it was
never a hazard is written next to it so the next reader does not have
to re-derive it.

**What depends on the old contract.** Nothing. Grepped the package's
`R/`, its tests and its four vignettes: the two `as_tmbstan()` test
blocks both assert on a real fit, and `vignette("sampling")` mentions
the function only for the tmbstan-build check. The change is announced
in `?as_tmbstan`'s `@return` and in NEWS.

**How it is tested.** With `control = list(nonsense_option = 1)`, which
makes rstan decline and return the shell deterministically, in about a
second, with no solver anywhere near it. That routes through
`as_tmbstan()` rather than `frm_sample()` on purpose: `frm_sample()`
screens an unknown option out before rstan sees it, and the escape
hatch does not, which is exactly the asymmetry that makes this the
right place to test the guard. Seen failing first on the pre-change
build: `Expected suppressWarnings(...) to throw a error.`

---

## 3. Tests, and seeing them fail first

Every new assertion was run against the UNFIXED build first, in two
rounds: the original three files against the pre-change 0.55.0 / 0.3.2
/ 0.1.1, and the blocks added for the review's findings against the
build that stood before those fixes.

**A capped tally is not a measurement.** The first version of this
document reported before-counts of 1/9/3, 1/4/5 and 1/8. Those do not
reproduce, and the reason is the instrument: testthat's summary
reporter stops after ten failures and prints "Maximum number of 10
failures reached, some test results may be missing", so the tally under
it is a floor, not a count. The reviewer got 2/18/3, 1/4/6 and 2/9/0
from the same files. The counts are dropped here. What reproduces, and
what the lane rules actually ask for, is the FAILURES, quoted verbatim
below; the after-counts come from a silent reporter that does not cap
and are reported per block.

The failures that matter, quoted from runs against the unfixed builds:

- `Expected rec$adapt_delta to equal 0.97. actual is NULL`. The sampler
  control was dropped on the floor on the fit route.
- `Expected frm_sample(fit, control = frmtmb_control()) to throw a
  error.` The old shape was accepted.
- `Expected conditionMessage(err) to match "frm_ode\(\)". Actual
  text: frm_sample() from a formula needs data =: there is no fitted
  model to take the design from`. The pre-flight did not exist, so
  assembly was reached first.
- `Expected "frm_ode()" %in% rows$feature_b to be TRUE. actual FALSE`.
  No refused row existed.
- `Expected frm_sample(fit, contro = list(adapt_delta = 0.97)) to throw
  a error.` The abbreviation reached rstan unchecked (F1).
- `Error: frm_sample(): the sampler returned no draws ... A known case:
  a tape that calls an external solver can fail inside tmbstan`, raised
  by `frm_sample(fit, contro = frmtmb_control())`. This is the exact
  mis-report the `control` guard exists to prevent, reached by
  mistyping the argument it guards (F1).
- `Expected frm_sample(fit, fit_control = list(adapt_delta = 0.99)) to
  throw a error.`, and the same for `start`, `data2`, `na.action` and
  `REML` (F2).
- `Error: frm_sample() is refused on a model using mi_pred(): No.`,
  raised by `bf(y | mi() ~ x)`, which writes the ATERM. The false
  refusal the key collision produced (F3).
- `Expected as_tmbstan(fit, chains = 1, iter = 20, refresh = 0) to
  throw a error.` The second door was unguarded (F5).
- `Expected suppressWarnings(...) to throw a error.`, from
  `as_tmbstan(fit, control = list(nonsense_option = 1))`. The escape
  hatch returned the empty stanfit instead of refusing (section 2.9).

After the fixes, per block, from a reporter that does not cap:

| file | blocks | pass | fail | error | skip |
|---|---|---|---|---|---|
| `frmtmb.sample/tests/testthat/test-stan-control.R` | 11 | 54 | 0 | 0 | 0 |
| `frmtmb.sample/tests/testthat/test-compat-preflight.R` | 16 | 81 | 0 | 0 | 0 |
| `frmtmb.ode/tests/testthat/test-compat.R` | 3 | 11 | 0 | 0 | 0 |

### Full suites, one file per R process

`NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`, `FRMTMB_STAN_CACHE` left
at its default so it reads `dev/stan-cache` in this worktree, and
`R_MAKEVARS_USER` pointed at that directory's `makevars-cxx17.mk`.

frmtmb.ode, all six files:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-bracket-access.R | 1 | 0 | 0 | 0 |
| test-compat.R | 11 | 0 | 0 | 0 |
| test-ode.R | 72 | 0 | 0 | 0 |
| test-ode-events.R | 108 | 0 | 0 | 0 |
| test-ode-nlf.R | 4 | 0 | 0 | 0 |
| test-ode-tv.R | 27 | 0 | 0 | 0 |

frmtmb.sample, all fifteen files. Every file was rerun after the
review fixes, because the `fit_control`/`start`/`data2`/`na.action`/
`REML` refusal is a behavior change that any existing call could have
tripped; none did:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-bracket-access.R | 1 | 0 | 0 | 0 |
| test-compat-preflight.R | 81 | 0 | 0 | 0 |
| test-conditional-effects-draws.R | 57 | 0 | 0 | 0 |
| test-draws-methods.R | 97 | 0 | 0 | 0 |
| test-draws-spellings.R | 48 | 0 | 0 | 0 |
| test-evidence-ratio.R | 33 | 0 | 0 | 0 |
| test-loo.R | 79 | 0 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| test-parallel-chains.R | 4 | 0 | 0 | 0 |
| test-prior-route.R | 9 | 0 | 0 | 0 |
| test-reparam.R | 260 | 0 | 0 | 0 |
| test-sample-direct.R | 136 | 0 | 0 | 0 |
| test-sampling-ported.R | 215 | 0 | 0 | 0 |
| test-simulators.R | 50 | 0 | 0 | 0 |
| test-stan-control.R | 54 | 0 | 0 | 0 |

frmtmb, the three files that read the registry (the only core file
changed is `R/compat.R`, and only its roxygen):

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-compat.R | 573 | 0 | 0 | 0 |
| test-compat-register.R | 103 | 0 | 0 | 0 |
| test-prior-compat.R | 195 | 0 | 0 | 0 |

One trap worth recording for the next lane: `dev/stan-cache/` is
gitignored, so a fresh worktree does not have it, and
`R_MAKEVARS_USER` then points at a file that does not exist. On the
first pass `test-loo.R` and `test-sampling-ported.R` each reported one
ERROR, both inside `brms::brm()` compiling a Stan program from scratch
and failing in the C++ compile. Copying the 112 cached programs plus
`makevars-cxx17.mk` from the main checkout into
`<worktree>/dev/stan-cache/` fixed both with no code change. Neither
error was related to this lane's work.

### R CMD check --as-cran

With `RSTUDIO_PANDOC` and the quarto and TinyTeX tool directories on
PATH, as the lane rules describe.

| package | status |
|---|---|
| frmtmb.ode 0.1.1 | 1 WARNING |
| frmtmb.sample 0.3.2 | 1 WARNING |

The WARNING is `checking CRAN incoming feasibility` in both cases, and
it is the pre-existing "frmtmb is not on CRAN" one. Verified by
checking the SAME two packages from the main checkout: identical text,
differing only in that this worktree's frmtmb.sample lists
`frmtmb.latent, frmtmb.ode` under "Suggests or Enhances not in
mainstream repositories" where main lists `frmtmb.latent` alone, which
is this lane's added Suggests. Everything else is OK, tests included
(frmtmb.ode 16 s, frmtmb.sample 18 s), and vignettes rebuild.

The V8 NOTE the lane rules expect did not appear for either extension;
it appears on core, whose HTML manual renders math. Core was checked
with tests, examples and vignettes off, since the change there is
documentation: every Rd check is OK and the notes are the V8 one and
CRAN incoming feasibility.

---

## 4. Defects found and NOT fixed

**Re-registering a package's compat rules duplicates them, and an
EDITED rule then makes core's own validator report a defect that is not
in the sources.** `frmtmb_register_compat(features =)` is a no-op for a
name it already carries, but `rules =` is appended unconditionally.
Running `frmtmb.sample:::.onLoad()` a second time in one session took
the rule table from 446 rows to 452 and left two identical
`rescor x frm_sample` rules.

For IDENTICAL rules that is harmless: they tie on specificity, agree on
the status, and `frm_compat()`'s answer does not move. The review
sharpened this, and the sharper version is the one worth acting on.
Re-register the SAME pair with a DIFFERENT status, which is what
`pkgload::load_all()` produces after you edit a rule, and the stale and
fresh rules tie on precedence and DISAGREE. `frm_compat()` still gives
the right answer, because the documented tie-break is later-wins. But
that disagreement is exactly the registry defect
`frmtmb_compat_validate()` exists to report, and it reports it: one
row. Core's `test-compat.R` asserts that table is empty. So a developer
who edits a rule's status and reloads gets a core test failure
describing a defect that is not in their sources.

The wrongness is a false positive rather than a missed refusal, which
is the right direction, but it is a real cost of a dev workflow the
project uses. A registry that deduplicated identical rules, or keyed a
contribution by registrant so a re-registration replaced rather than
appended, would close it. Probes:
`scratchpad/stanctl-reload-probe.R` and the review's
`rvsc-p5-reload2.R`. Not fixed here: changing the registry's behavior
is out of this lane's scope, and it is now the strongest item in the
core-side queue behind `frm_features_used()`.

**`?frm_ode`'s example gradient.** Not touched, but
`dev/upstream/rtmbode-issues.md` section 4 records that the `?ode`
example's gradient disagrees with central differences by 10 percent at
the starting value on both the stock and the patched build. Still open
upstream.

---

## 5. What was NOT done

- **Nothing was filed upstream and nobody was contacted.** Section 6 is
  the text, prepared for the user to send.
- **The pre-flight was not made to decide a covariance structure, a
  fitting mode, a shared parser name or a base R name.** It cannot be
  done from a display name, because the discriminator is formula
  position; those cases warn and sample. The seam that would close it,
  `frm_features_used(spec)`, is filed in the plan's Core seams table
  and argued in section 2.4.
- **The registry read was not short-circuited.** Measured at 0.09 s to
  0.13 s per call and judged immaterial; see 2.6.
- **`frm_ode()` was given no blanket `untested` rule.** Core's default
  for an undeclared pair is already `untested` with the note "No rule
  covers this pair", so a blanket would have added 117 rows and said
  the same thing.
- **The warn branch was not made once-per-session.** It fires per call,
  which is per model fit and not per iteration, and a once-per-session
  warning can be missed by whoever reads only the tail of a log. The
  cost is recorded in section 2.7 instead.

---

---

## 6. Filing-ready text for the RTMBode patches

Three issues against `kaskr/RTMB`, subdirectory `RTMBode/`, and one
against `deSolve`. The full analysis, the thirteen reproduction
scripts and the three-commit patch series are in
`extensions/frmtmb.ode/dev/upstream/`. The series applies to commit
`5242257`, which is what r-universe is currently building, so a filing
today needs no rebase.

Suggested order, because issue 2's second defect is only reachable once
issue 1 is fixed: file 1, then 2, then 3, then the deSolve one.

Each issue body below is the text from
`dev/upstream/rtmbode-issues.md` section 5, which was written for this
purpose and is still accurate. Attach the corresponding patch file.

### 6.1 Issue 1 (attach `patches/01-solver-failure-nan.patch`)

Title: **`RTMBode`: a failed `deSolve` call escapes the adjoint node as
an R error**

Body: section 5.1 of `dev/upstream/rtmbode-issues.md`, unchanged. Its
claims re-verified today: `ODEadjoint()` in RTMBode 1.0 at commit
5242257 contains no `tryCatch`, and the frmtmb-side reproducer in
`scratchpad/stanctl-ode-repro.R` still produces the
`illegal input detected before taking any integration steps` signature
from `frm_sample()` on a fixed-effect-only ODE model.

One addition worth making to the body when filing, because it
strengthens the report and was measured in this lane rather than the
last: the abort reproduces on an ODE model with **no random effects
and one solve group per subject**, not only on the two-node
Lotka-Volterra of the original write-up. The earlier note that
"whether single-node ODE fits also fail is untested" is now answered:
they do.

### 6.2 Issue 2 (attach `patches/02-events-ad-path.patch`)

Title: **`RTMBode`: `events` are unusable on the AD path, and two of
the three methods would be silently wrong**

Body: section 5.2, unchanged. Re-verified: `addInfo()`, `setTape()` and
`augment()` in the installed build carry no state names and refuse no
event method.

This is the highest-severity of the three, because `replace` and
`multiply` return correct SOLUTION values and wrong GRADIENTS (42 and
59 percent off central differences), `replace` is deSolve's fallback
when no method is given, and nothing warns. An optimum found with
those gradients is not the maximum likelihood estimate.

### 6.3 Issue 3 (attach `patches/03-augmented-workspace-guard.patch`)

Title: **`RTMBode`: silent `NaN` gradients above ~9 states under
Laplace; `lsoda`'s work array overflows R's integer range**

Body: section 5.3, unchanged. Re-verified: no `checkWorkspace()` in the
installed namespace, and no function in it mentions `lrw`.

### 6.4 Issue 4, against deSolve rather than RTMBode

Title: **`lsoda`: validate `lrw` before it overflows integer range**

Body: section 5.4, unchanged. This one is a message request, not a bug
report: the limit is structural (`lrw` is a Fortran `INTEGER`), and the
ask is that `lrw = 22 + neq * max(16, neq + 9)` be computed as a double
and refused with a message naming the limit and the banded and sparse
routes, instead of overflowing silently at `neq = 46337` and printing
`cannot allocate memory block of size 134217728 Tb`.

### 6.5 What to caveat in the covering note

`dev/upstream/rtmbode-issues.md` section 4 lists FIVE things the
write-up does not explain, not four as the first version of this
document said. Three belong in the covering note:

1. **`checkWorkspace()`'s warning frequency was observed, not proved.**
   This is the one that matters most and the one the first version
   missed. `ODEadjoint()` is reconstructed on every gradient, and the
   claim that `Df` is cached per node, so the warning fires once per
   node per order rather than once per evaluation, was read off the
   observed behavior and not off the code. That caveat is about patch
   03, the patch being offered, so it belongs in the filing more than
   the other two. The write-up itself suggests `warning(call. = FALSE)`
   plus a `once` guard if a caller builds many nodes; offer that in the
   issue.
2. **The `ode45`, `rk4` and `euler` crashes at 32 states** reported
   earlier were NOT reproduced and may be a separate allocation
   failure in the fixed-step methods. Say so rather than implying the
   ceiling story covers them.
3. **`forcings` is unusable through RTMBode** for a reason nothing here
   patches: `desolve_derivs` has no forcing hook, and `approxfun()` in
   the dynamics is silently frozen at t = 0 because `func2tape()` tapes
   the derivative function once at an all-zero point. This is a fourth
   report with no patch. File it separately or hold it back; do not
   fold it into the three that have patches.

The other two items in that section (the deSolve console diagnostics
that patch 01 does not suppress, and the `?ode` example's 10 percent
gradient disagreement, which predates every patch) are fine as written.

---

---

## 7. Files touched

Changed:

- `R/compat.R` (roxygen only: the "What a refused row promises"
  section) and the regenerated `man/frmtmb_register_compat.Rd`
- `NEWS.md`, under `# frmtmb (development version)`
- `extensions/frmtmb.sample/R/sample.R` and the regenerated
  `man/frm_sample.Rd` and `man/as_tmbstan.Rd`
- `extensions/frmtmb.sample/DESCRIPTION` (Suggests: frmtmb.ode)
- `extensions/frmtmb.sample/NEWS.md`, under
  `# frmtmb.sample (development version)`
- `extensions/frmtmb.sample/vignettes/brms-posterior.Rmd`
- `extensions/frmtmb.sample/vignettes/posterior-diagnostics.Rmd`
- `extensions/frmtmb.ode/R/zzz.R`
- `extensions/frmtmb.ode/R/ode.R` (the "Sampling an ODE fit" section)
  and the regenerated `man/frm_ode.Rd`
- `extensions/frmtmb.ode/NEWS.md`, under
  `# frmtmb.ode (development version)`
- `.github/workflows/check-frmtmb-sample.yaml`: the sibling-install
  step covers frmtmb.ode and frmtmb.eam as well as frmtmb.latent, and
  the `paths:` filter triggers on all three
- `dev/extension-gaps-plan.md`, rows 1.1 and 1.2 edited in place, and
  two rows added to Core seams (`frm_features_used(spec)` and a
  resolved-answer cache for `frm_compat()`)

Added:

- `extensions/frmtmb.sample/tests/testthat/test-stan-control.R`
- `extensions/frmtmb.sample/tests/testthat/test-compat-preflight.R`
- `extensions/frmtmb.ode/tests/testthat/test-compat.R`
- this file

Not committed. `dev/stan-cache/` was populated from the main checkout
and is gitignored.

## 8. What the review changed, and what it did not

The review (`dev/reviews/2026-09-08-stanctl.md`) reproduced every
headline number exactly, including the divergence counts seed for seed,
and named eleven findings. All eleven are addressed above:

| finding | what it was | where |
|---|---|---|
| F1 | `contro =` bypassed `check_stan_control()` through the dots and rstan's own partial matching | closed, section 1.6 |
| F2 | the rename created a new silent-ignore under `fit_control` | closed, all five arguments, section 1.8 |
| F3 | matching on `key` refused the wrong feature where two share a parser name | closed, section 2.4 |
| F4 | "parentheses keep base R out" was false; `trunc()` is the counterexample | claim corrected and the code now checks, section 2.4 |
| F5 | `as_tmbstan()` was unguarded and still bricked sessions | closed, the pre-flight is on both doors, and section 2.9 closes the empty-stanfit silence behind it |
| F6 | the pre-flight ran after the namespace loads, so a first-call refusal was slower than the failure it replaced | closed by moving it, section 2.5 |
| F7 | the warn branch's false-alarm behavior was not recorded | recorded, section 2.7 |
| F8 | the core promise was stronger than the code | reworded, section 2.8 |
| F9 | no NEWS bullet for the silent discard | the reviewer added it; kept |
| F10 | four record slips, including capped before-tallies | all four corrected, sections 1.3, 2.5, 3 and 6.5 |
| F11 | `control = list()` was refused with the wrong message | closed, section 1.6's guard rewritten |

Two of the review's own conclusions are adopted rather than argued
with: that the `control`/`fit_control` redirection was right, on
`docs/SPEC.md:456` and `:149` rather than on the grounds first given;
and that the re-registration defect's real cost is a false positive in
`frmtmb_compat_validate()` after an EDITED rule is reloaded, which is
sharper than the version first recorded here and is now in section 4.

Nothing in the review was disputed.

### Round 2

The re-check reproduced every round-1 closure and named four more.

| finding | what it was | where |
|---|---|---|
| R1 | "display name plus uniqueness plus base R is the largest SOUND set" is false; `se()` in a nonlinear body is the counterexample, on a model that fits | claim corrected in the code comment, in `?frm_sample` and in section 2.4; route 2 removes four of the eight names it names |
| R2 | an aterm-position test is available without the core seam | BUILT, not just recorded: 9 of the 20 call-route rows move from name-matching to position-matching, and the filed seam is rescoped to three kinds |
| R3 | the uniqueness test made a shipped refusal revocable by a third party's registration | closed by matching the display name, which is un-revocable by proof rather than by measurement; the now-unreachable uniqueness branch was deleted |
| R4 | the `wiener`/`hmm`/`lca` block ran zero assertions | closed: it loads the owning packages, skips with a reason, and counts what it checked, so a future no-op fails loudly. 10 assertions where there was 1 |
| R5 | the no-draws guard had no test independent of rstan's behavior | added, on a hand-made `new("stanfit", sim = list())`; the vacuous end-of-string assertion dropped |
| R6 | `base_r_function()` loaded four namespaces and listed 13 of 14 base packages | memoized, and `tcltk` added to make it 14 |
| R7 | the family branch does not decompose `mixture()` | said so, in the code comment and in `?frm_sample` |
| numbers | 1085 should be 1081; the first-call figure was 0.580 s not 0.355 to 0.420 | corrected, and the DIFFERENCE is now the finding: see section 2.6 |

Two of the review's judgments are adopted rather than argued with: that
the no-draws test should stay routed through `as_tmbstan()` because the
asymmetry is the point, and that its round-1 suggestion of matching the
display name was right after all. My round-1 reason for rejecting it,
that `mi_pred` is a name no formula writes, was correct about the miss
and wrong about which failure to prefer: a miss is safe and the key
route was both a false-alarm source and revocable.

### Round 3

| finding | what it was | where |
|---|---|---|
| F-D | `frmtmb.eam` and `frmtmb.ode` in Suggests break the frmtmb.sample CI job: `R CMD check` gives `checking package dependencies ... ERROR` for a Suggests it cannot resolve, before any test runs, and neither exists in any repository | fixed in `.github/workflows/check-frmtmb-sample.yaml`; see below |
| F-A | the deletion left the three `_pred` rows routing to `call` on a name no formula writes, so they could never match and said nothing, while every other undecidable case warned | closed: a call-shaped row whose display name differs from the registry key routes to the warning, and the three are asserted as a set |
| F-B | `preflight_route()`'s `ft` was a dead parameter, so the revocation test's key assertion passed because the argument was ignored, and the block never performed the registration it described | `ft` is now live (F-A uses it), and the block registers the bare name for real against the live registry |
| premise 3 | false as a claim about the vocabulary (`ar1` and `gr_cov` are bare rows a formula writes as calls) and true as a claim about the code | reworded in section 2.4 and in the code comment |
| the kind clause | the route is chosen from the row's KIND as well as its name, and premise 1 covers both | said, in section 2.4 and in the code |
| F-C | premise 1's other half: the registry cannot be silently taken over, so it can be loudly deadlocked instead | recorded in section 2.4 |
| the cost trigger | an eighth-extension trigger is the wrong variable and the trigger is already met | replaced by a core filing, section 2.6 |
| the total | 1117, not 1081 | corrected below |

**F-D, and what can and cannot be verified from here.** The R-side
semantics I verified; the runner I could not.

Verified locally: a copy of frmtmb.sample carrying one unresolvable
Suggests gives `checking package dependencies ... ERROR`, `Package
suggested but not available`, `Status: 1 ERROR`, and the check stops
before a single test. So the break is real and it is fatal rather than
degrading. Also verified: core, then frmtmb.latent, frmtmb.ode and
frmtmb.eam each install cleanly from the checkout into an empty
library, in that order, which is exactly what the new workflow step
does; none of the three needs anything but core, RTMB and base
packages to install.

NOT verifiable from here: whether the GitHub runner behaves as the
local check does. `_R_CHECK_FORCE_SUGGESTS_` appears nowhere under
`.github/`, so it is at its default of true, and if
`check-r-package@v2` set it internally the ERROR would degrade to a
NOTE. The project's own precedent says it does not: the
`frmtmb.latent` install step exists for exactly this reason and its
comment says so. The fix follows `check-frmtmb-learn.yaml`, which the
coordinator identifies as the pattern, and it also widens the `paths:`
filter so a change to either sibling retriggers this job. **This one
needs a push to confirm.**

F-D also decides R4. `skip_if_not_installed("frmtmb.eam")` means the
eam-and-latent block currently SKIPS on CI, so the plan's second check
for item 1.2 was not running there. The same workflow step fixes both,
and that is why the block's skip is loud rather than silent.
