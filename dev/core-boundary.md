# What goes in core, and what goes in an extension

Status: adopted 2026-10-23. In 0.74.0 (2026-10-23), all five proposed
moves were made: curve inference into core; `cox()`, `royston_parmar()`
and `rp_floored()` into the new `frmtmb.survival`; `mixture_mvn()` into
`frmtmb.latent`; `car()` and `spde()` into the new `frmtmb.spatial`.
`frmtmb.spline` is retired. The section "Moves made in 0.74.0" has the
numbers. 0.75.0 (2026-10-24) moved no feature. It added the seams
`frmtmb_register_aterm(option =)` and `frmtmb_register_autocor()`,
through which `frmtmb.survival` supplies brms's `bhaz()` and
`frmtmb.spatial` brms's `sar()`; `fcor()` is in core.

This document states the rule that decides where a feature lives, the
guardrails that a move must satisfy, and the moves that the rule
proposes. It applies to core `frmtmb` and to the `extensions/frmtmb.*`
packages.

## The rule

Core holds what most users need and what every extension builds on.
Extensions hold model classes and workflows for one domain.

brms parity is measured on the whole ecosystem: core plus every
extension. A brms feature does not have to be in core. brms is a
general-purpose package that puts almost everything in one place, and
frmtmb does not have to copy that layout to match brms's behavior.

### What core holds

- The formula grammar: `bf()`, `mvbf()`, `nlf()`, addition terms, and
  the parser.
- The general families: the GLM families, the distributional families,
  ordinal families, mixtures, zero-inflated and hurdle families.
- General structure: smooths, random effects, structured covariance,
  autocorrelation, missing data, and measurement error.
- Inference: maximum likelihood, REML, the Laplace approximation,
  standard errors, profiles, bootstrap, and the SE checks.
- Post-fit methods: prediction, simulation, residuals, diagnostics,
  model comparison, and `conditional_effects()`.
- The seams that extensions use (see "The seams").

### What an extension holds

- A model class for one domain, for example decision models
  (`frmtmb.eam`), ODE models (`frmtmb.ode`), learning models
  (`frmtmb.learn`), or latent-state models (`frmtmb.latent`).
- A second inference engine with its own dependencies
  (`frmtmb.sample`).
- Workflows and post-fit tools that only that domain uses.

### How to decide

Ask these questions in order. The first "yes" decides.

1. Do extensions build on it? Then it is core.
2. Would most users of frmtmb meet it, whatever their field? Then it
   is core.
3. Is it a model class or a workflow for one domain? Then it is an
   extension.

"brms has it in core" is not a reason to keep a feature in core. It is
a reason to make sure that a brms user finds it (see "No silent gaps").

## Guardrails

### 1. No silent gaps

Every brms spelling must resolve to one of these:

- a core feature;
- a feature of a named extension;
- a refusal that says why frmtmb does not support it.

When the extension is not loaded, the refusal names the package and
tells the user how to load it. The Wiener family already works this
way. brms has `wiener` in core; frmtmb has it in `frmtmb.eam`. Without
that package, `frm(..., family = "wiener")` and `brmsfamily("wiener")`
stop with "The wiener family is in the extension package frmtmb.eam,
which is not loaded." The mechanism is `extension_family_homes` in
`R/families.R`.

Since 0.74.0 the same mechanism covers formula specials and addition
terms: `extension_special_homes` and `extension_aterm_homes` in
`R/brms-vocabulary.R`. Without `frmtmb.spatial`, `car()` and `spde()`
stop with "The car() special is in the extension package
frmtmb.spatial, which is not loaded". The brms names that frmtmb does
not support at all are in three tables of refusals in the same file
(`custom` since 0.75.0), each with its reason. `bhaz()` is in
`extension_aterm_homes`, and frmtmb.survival refuses `bhaz(gr = )` by
name. `fcor()` is in core, and `sar()` is in `frmtmb.spatial`
(`extension_special_homes`).

A test enforces this guardrail: `tests/testthat/test-guardrail1.R`.
It probes every family, addition term and formula special of brms
2.23.0 in a child process that holds core alone. Each name must
resolve in core, be refused with the name of the extension that holds
it, or be refused by name with a reason. Then, for each extension
that a refusal names and that is installed, a child with that
extension loaded must resolve the same calls. The compat registry lets
an extension contribute its features and rules
(`frmtmb_register_compat()`), and `test-compat-extensions.R` checks the
contributed rows with all eight extensions loaded.

A saved fit is a second way to meet a moved feature. Since 0.74.0 each
fit records the package of its family and of each covariance
structure, and since 0.75.0 of each residual correlation term. Read
where such a package is not installed, the fit stops with the
package's name (class `frmtmb_extension_missing`) at the first method
that needs it (`R/extension-home.R`). Since 0.75.0 that one refusal
names every package the fit lacks.

### 2. Interface before move

No feature leaves core until the public seam that it needs exists, and
until an extension uses that seam with no `frmtmb:::` calls. All eight
extensions meet this today: no code in their `R/` reads a core
internal. Some of their tests do.

If a move needs a seam that does not exist, the seam is the first piece
of work. The seam must be useful on its own, so that a third party can
use it without the move.

### 3. Measure before move

A move must have a reason that is measured, not only a cleaner layout.
For each candidate, record:

- its lines of code in `R/`;
- its test time, and that time as a share of the core suite;
- how often it changed in the last five rounds;
- the dependencies that it alone needs;
- the ported brms rows that depend on it.

Record the numbers in the round's findings. A feature that is small,
stable, and needs no extra dependency is a weak candidate, whatever its
domain.

### 4. Parity is measured on the ecosystem

The ported brms suite (`dev/brmsport-*`) and the vignette remeasure run
against core plus every extension. Each row that needs an extension is
labeled with that extension, so that a regression in an extension shows
in the ledger. The families count already works this way: at 0.73.0,
39 of 40 brms families were in core, and 40 of 40 with the extensions.
Since 0.74.0 `cox` is in an extension too.

## The seams

The public interfaces that an extension can use today:

| Seam | What it registers |
|---|---|
| `frmtmb_family()`, `custom_family()`, `as_frmtmb_family()` | a family |
| `frmtmb_register_aterm()` | an addition term, such as `dec()` |
| `frmtmb_register_aterm(option =)` (0.75.0) | an addition term that sets an argument of the family, such as `bhaz()` of `cox()` |
| `frmtmb_register_compat()` | compat features and rules |
| `frmtmb_register_frame_check()` | a check on the model frame |
| `frmtmb-sampling-api` exports | what `frmtmb.sample` reads from a fit |
| `frm_hazard_reads()` | a test that a package reads no hazard container |
| `frmtmb_register_covstruct()` (0.74.0) | a covariance structure: a bar structure `name(terms \| g)`, or a field written as a formula special, such as `car()` |
| `frmtmb_family(simplex_pars =)` (0.74.0) | which extra parameters of a family are simplexes, such as the baseline of `cox()` |
| `frmtmb_family(aterm_data =)` (0.74.0) | data that a family builds once from the observed response |
| `post$curve_check` (0.74.0) | a family hook that lets `frm_curve()` and its companions refuse a fit |
| `us_sigma()`, `us_theta()` (0.74.0) | the map between an unstructured covariance and its parameter vector, on the tape |
| `frmtmb_local_registries()` (0.74.0) | a test helper that saves every registry and restores it |
| `frmtmb_register_autocor()` (0.75.0) | a residual correlation term: one affine map from the response to independent innovations, such as `sar()`; its `bounds` slot gives a bounded parameter brms's class prior and frmtmb.sample's flat default |
| `frmtmb_autocor_pars()` (0.75.0) | lists a fit's registered residual-term parameters with their intervals, for a package that would otherwise keep its own list of names (frmtmb.sample's `lagsar`, `errorsar`) |
| `frmtmb_covstruct_info()` (0.75.0) | reads a structure's properties (field, coupled levels, full rank) from the registry, for a package that would otherwise keep its own list of names |
| `frmtmb_register_covstruct(full_rank =)` (0.75.0) | whether a row's own effect can reach any value with a finite density, which the weights check reads; a structure that states nothing leaves the check open |

Until 0.74.0 a structured covariance or a precision matrix that
depends on parameters had no seam: the parser built every covariance
structure in core. `frmtmb_register_covstruct()` is that seam. It is
how `frmtmb.spatial` registers `car()` and `spde()`.

## Proposed moves

These were the proposals under the rule, each to be measured under
guardrail 3 before it was planned. All five were made in 0.74.0
("Moves made in 0.74.0").

| Feature | From | To | Prerequisite |
|---|---|---|---|
| curve inference: `frm_curve()`, `frm_curve_deriv()`, `frm_curve_feature()` | `frmtmb.spline` | core | none; it uses only public interfaces |
| `cox()` | core | `frmtmb.survival` | an `extension_family_homes` entry |
| `royston_parmar()`, `rp_floored()` | `frmtmb.spline` | `frmtmb.survival` | none |
| `mixture_mvn()` | core | `frmtmb.latent` | none; no string spelling reaches it, and brms has no such family |
| `car()`, `spde()` | core | `frmtmb.spatial` | the covariance-structure seam, and a refusal by name for a formula special |

Curve inference moves the other way: it applies to any fitted curve,
not only to splines.

After these moves `frmtmb.spline` held nothing, and it was retired in
0.74.0. Its NEWS is an article of `frmtmb.survival`, "History of
frmtmb.spline", and the pages of its site redirect to their new
homes.

These stay in core: `gp()` and HSGP terms, functional regression, the
general families, `mixture()`, smooths, missing data, and measurement
error. Each is general, and extensions build on several of them.

### A suggested order

1. Curve inference into core. It needs no new seam.
2. `frmtmb.survival`: `cox()` and Royston-Parmar. It needs only the
   family pointer. Its later contents are listed in "Not decided".
3. `mixture_mvn()` into `frmtmb.latent`.
4. The covariance-structure seam, then `frmtmb.spatial`.

## Moves made in 0.74.0

The user kept all five moves at the consolidation of 2026-10-23,
on the numbers below. Each row is from the lane record named. The
range for commits is `4f5ea39f..19b162b9` (0.68.1 to 0.73.1, the 0.69
to 0.73 rounds, 37 commits). A test time is the sum of per-test times
from one run with eight processes on a shared machine, so it is an
order of magnitude only: the same 262 core files summed to 2045.6 s,
1678.9 s and 1471.4 s in three runs (`dev/mvnlatent-findings.md`).

| Move | Lines in `R/` (all / code) | Test time, share of the source suite | Commits that changed it | Dependencies only it needs | Ported brms rows | Record |
|---|---|---|---|---|---|---|
| curve inference, `frmtmb.spline` to core | 1846 / 805, and 43 / 40 of compat rules | 62.7 s, 72.9% of `frmtmb.spline`'s 86.0 s | 1 (f5770ca0) | gratia (tests) | 0 | `dev/splitsurv-findings.md` |
| `royston_parmar()`, `rp_floored()`, `frmtmb.spline` to `frmtmb.survival` | 1139 / 411, and 37 / 36 of compat rules | 22.1 s, 25.7% of `frmtmb.spline`'s suite | 0 | flexsurv, rstpm2 (tests, vignette), survival (tests) | 0 | `dev/splitsurv-findings.md` |
| `cox()`, `cox_baseline()`, core to `frmtmb.survival` | 324 / 174 in `R/families.R`, and 12 of compat rules | 1.2 s, 0.09% of core's 1428.7 s | 1 (58b97a2c) | splines2 (tests) | 0 | `dev/splitsurv-findings.md` |
| `mixture_mvn()`, core to `frmtmb.latent` | 355 / 235 in `R/families.R` (0.41% of `R/`), and 13 of compat rules | 6.3 s, 0.31% of core's 2045.6 s | 0 | mclust (tests) | 0 | `dev/mvnlatent-findings.md` |
| `car()`, core to `frmtmb.spatial` | 446 of its own, and 66 that name it elsewhere | 7.33 s, 0.43% of core's 1711.23 s | 1 (ed824ffb, shared with `spde()`) | none | 0 | `dev/spatial-findings.md` |
| `spde()`, core to `frmtmb.spatial` | 157 of its own, and 49 that name it elsewhere | 2.19 s, 0.13% | 1 (ed824ffb) | fmesher (a test, a vignette) | 0 | `dev/spatial-findings.md` |

`car()` and `spde()` also share 2.64 s (0.15%) of tests. For
`mixture_mvn()`, core lost 235 lines net, because the seam
`us_sigma()` and its help came back (`git diff --numstat`). By the
rule's words in guardrail 3, every move but curve inference is a weak
candidate: small, stable, and with no dependency of its own outside
the tests. The lanes and their reviews said so. The user kept them
for the rule's question 3 (a model class for one domain) and for the
seams that the moves made.

What each move added to core:

- curve inference: three exports, two print methods, the family hook
  `post$curve_check`, and the article "How to report a fitted curve".
- `cox()`: the seams `frmtmb_family(simplex_pars =)` and
  `frmtmb_family(aterm_data =)`, which replace the branches that named
  `cox()` in four files, an `extension_family_homes` entry, and
  `R/extension-home.R`.
- `mixture_mvn()`: the exports `us_sigma()` and `us_theta()`. No string
  spelling reached the family, so it has no `extension_family_homes`
  entry.
- `car()`, `spde()`: `frmtmb_register_covstruct()`, the covariance
  structure registry as an environment, `R/brms-vocabulary.R` with its
  tables of homes and refusals, `test-guardrail1.R`, and the export
  `frmtmb_local_registries()`.

Costs measured at the consolidation (`dev/round-20261023.md`): the
release is in lockstep, because `frmtmb.latent` 0.6.2 does not load
with frmtmb 0.74.0, and the three packages that moved code need frmtmb
0.74.0. The compat table with every extension loaded keeps every
status; its notes differ from 0.73.0's on the pairs the lanes
accepted.

## Costs that a move adds

1. **Install friction.** A brms user whose code uses `car()` or `cox()`
   must install one more package. The refusal of guardrail 1 tells the
   user which one, but it does not remove the step. Since 0.74.0 the
   README and `vignette("brms-migration")` document one way to install
   every extension; there is no meta-package ("Not decided").
2. **Version floors.** Each move adds a floor from an extension on core.
   Floors slip easily: in the 0.72.0 and 0.73.0 rounds,
   `frmtmb.sample`'s floor had to rise in each. Each release must
   measure every floor, as `dev/rel073-floor.sh` does.
3. **Cross-package defects.** An interaction between two features is
   harder to see when the features are in different packages. In the
   0.73.0 round the weights check and the new `cens()` families
   interacted. Both were in core, and only the release review found
   the defect. Each release must test the merged build of all packages,
   in one session, as the release gate already does
   (`dev/release/run-sessions.sh`, three orders).
4. **CI cost.** Each extension adds a job to the coverage matrix and a
   check workflow.

## Not decided

- What else goes in `frmtmb.survival`. The candidates:
  - more parametric families: Gompertz, log-logistic, generalized
    gamma, generalized F, piecewise exponential;
  - a `Surv()` front end, delayed entry, and counting-process data;
  - cure models, competing risks, and multistate models;
  - relative survival;
  - joint longitudinal and survival models;
  - survival and hazard curves with bands, RMST, marginal survival,
    and prediction metrics.

  The reference packages for these are survival, flexsurv, rstpm2,
  mstate, and cmprsk, not brms.
- Decided on 2026-10-23: no meta-package installs core and every
  extension. Users rarely fit models from two domains in one analysis.
  The guardrail-1 refusal already names the package to install. A
  meta-package adds release, floor and check cost. The README and
  `vignette("brms-migration")` document one way to install every
  extension, which covers the rest.
- Decided in 0.74.0: `spde()` moved with `car()`, though the
  measurements show that it is small and stable. fmesher went with it,
  from core's Suggests to `frmtmb.spatial`'s.
