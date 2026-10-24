## The 0.74.0 round, 2026-10-23

`dev/round-20261023.md` is this round's record. Each lane has
`dev/<lane>-findings.md` and `dev/reviews/2026-10-23-<lane>.md`. It is
the first round under the core boundary rule (`dev/core-boundary.md`,
adopted the same day): every lane moved a feature between packages,
measured its guardrail-3 numbers first, and the user decided at the
consolidation to keep all five moves. Three lanes off `19b162b9`
(0.73.1), merged in this order:

- `splitsurv`: curve inference (`frm_curve()`, `frm_curve_deriv()`,
  `frm_curve_feature()`) from frmtmb.spline into core, with the family
  hook `post$curve_check`; the new extension frmtmb.survival with
  `cox()`, `cox_baseline()`, `royston_parmar()` and `rp_floored()`; the
  seams `frmtmb_family(simplex_pars =, aterm_data =)`;
  `R/extension-home.R`, so a saved fit names the extension its family
  needs; frmtmb.spline retired.
- `mvnlatent`: `mixture_mvn()` from core into frmtmb.latent; the
  exports `us_sigma()` and `us_theta()`; the mixture component
  interface documented in `?frmtmb-extension-api`.
- `spatial`: `frmtmb_register_covstruct()` and the covariance structure
  registry; `R/brms-vocabulary.R` and `test-guardrail1.R` (every brms
  name resolves in core, in a named extension, or is refused by name);
  the new extension frmtmb.spatial with `car()` and `spde()`;
  `frmtmb_local_registries()` exported.

Versions: frmtmb **0.74.0** (minor, breaking); frmtmb.survival
**0.1.0** and frmtmb.spatial **0.1.0** (new); frmtmb.latent **0.7.0**;
frmtmb.sample **0.21.1**; frmtmb.learn **0.7.3**; frmtmb.spline
retired; eam 0.13.0, coupling 0.6.1 and ode 0.7.2 unchanged. Every
floor was measured (`dev/rel074-floor.sh`): survival, spatial, latent
and sample need frmtmb 0.74.0; learn keeps 0.68.0. The release is in
lockstep: frmtmb.latent 0.6.2 does not load with frmtmb 0.74.0. The
release library is `C:/Users/adf44/source/r/rellib-r12`; `rellib-r11`
keeps 0.73.0.

**A release fix.** frmtmb.survival's first R CMD check --as-cran was
`Status: 1 ERROR, 1 NOTE`: `test-saved-fits.R` could not build its
"not installed" library stack under R CMD check. The coordinator
approved a release fix (`dev/rel074-release-fix.R`): one helper both
packages copy (`helper-absent-package.R`, frmtmb.spatial's two stacks
and probe), and an "installed but unloadable" arm for frmtmb.survival
that runs under R CMD check. The check is now `Status: 1 NOTE` (V8);
`dev/round-20261023.md`, "R CMD check --as-cran".

**The release review** (`dev/reviews/2026-10-23-release.md`, FIX, no
blocker) is applied as release fixes (`dev/rel074-release-review-fix.R`;
the record's "After the release review"): the README names
frmtmb.spatial and says `car()` and `spde()` are there; the README and
`vignette("brms-migration")` give one way to install every extension
(the user decided on 2026-10-23: no meta-package); the README's first
example runs; a `cox()` fit saved by frmtmb 0.73.0 is refused by name
with "Refit the model with library(frmtmb.survival) loaded"
(`old_core_family_home()`); NEWS tells users to remove frmtmb.spline;
the one-session runner watches the covariance-structure registry; a CI
job keeps the two helper copies identical; and the nits. After them:
441 of 441 files, 29752 assertions, 0 failures; frmtmb `Status: 1
NOTE`, frmtmb.sample `Status: OK`.

### Decisions made at the consolidation, accepted by the user

- **frmtmb.sample's floor rises to frmtmb 0.74.0** although its code
  loads on 0.73.0: four of its test files fail there (they read
  `frmtmb:::cs_entry()`, `frmtmb_family(simplex_pars =)` and
  frmtmb.latent 0.7.0's `mixture_mvn()`), and the rule is that a file
  that fails only on the old core names the floor.
- **The ported brms suite ran with every extension loaded**
  (`FRMTMB_PORT_EXTENSIONS=all`, a new switch of `dev/brmsport-run.R`
  and `dev/brmsport-bin2-run.R`), as guardrail 4 asks.
- **`dev/ciharden-scan.sh` has a plain tier** (`CIHARDEN_UNGATED=true`)
  and runs the nine suites.
- **The NEWS sentences over 25 words are split** (19), as at 0.73.0,
  with four edits of fact (`dev/rel074-news-split.R`).
- **`mixture_mvn()`'s row of `dev/core-boundary.md`** keeps the user's
  wording ("none; no string spelling reaches it").

### What the 0.74.0 round is evidence for

**A lane measures what it can load.** splitsurv compared the compat
table with core and the one extension it changed; the merged table with
every extension loaded also differs on 7 notes of frmtmb.learn, the
note splitsurv reworded there (`dev/rel074-compat.R`). A move's compat
check loads every extension, as `test-compat-extensions.R` does.

**Two absent-package checks compose, but one at a time.** A fit with a
family from frmtmb.survival and a structure from frmtmb.spatial is
refused by a package's name in every method that needs one, and never
fails on an internal (`dev/rel074-saved-mixed.R`); the message names
one missing package per call.

**R CMD check is a different library stack.** A test that builds its
own `.libPaths()` for a child process passed one file per process and
in one session, and failed only under R CMD check --as-cran, where the
tests see one temporary library holding the package and every
dependency. A lane that writes such a test runs R CMD check after it,
and probes the stack before it asserts, as frmtmb.spatial's does.

**A gate can hide a test from the skip count.** A file that defines its
gated tests only when the gate is set reports no skip in the plain
tier, so "rerun the files with a skip" missed 634 assertions; the
release reran every file with every gate on.

**A new package's first roxygen pass links to another package.** In a
fresh `man/`, `[cox()]` resolved to `brms::cox`, because the package's
own topic did not exist yet; the second pass links locally. Run
roxygen twice on a new package and compare.
