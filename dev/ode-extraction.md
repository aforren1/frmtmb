# Extracting the ODE seam: frmtmb.ode

Status: EXECUTED 2026-09-03 on branch wt-odeext against v0.46.0.
extensions/frmtmb.ode exists and holds frm_ode(); core keeps neither
R/ode.R nor any RTMBode dependency. Written as design against v0.45.0,
to execute after protocol steps 6 through 9 landed (the ratchet and
frame.R would otherwise have been contested by two lanes).

Corrections found in execution, kept here because the note is the
record: frame.R reached ode.R ONCE, not twice. The only core reference
to any R/ode.R symbol was check_ode_constancy() at frame.R:2221. The
supposed second hook, the nl-body auto-wrap riding ad_overload_fn, is
not ODE-specific at all: drop_nl_lexical_datavars() accepts any
function, list, environment or language object, and frm_ode appeared in
it only as prose. Core also carried NO ODE compat rows, so the
compat-row migration this note anticipated had nothing to move.

## Layout

The repo root STAYS the core package. Extensions live under
`extensions/<pkg>/`, excluded from the core build by one
`.Rbuildignore` line (`^extensions$`). Rationale:

- Moving core into a subdirectory would churn every CI path, the
  pkgdown config, and the ergonomics of four days of git history, for
  no functional gain.
- `R CMD build` of the root is unaffected once the directory is
  build-ignored; CRAN never sees the extensions.
- r-universe discovers packages by DESCRIPTION anywhere in the repo,
  so `extensions/frmtmb.ode` builds and serves without configuration.
- Development installs are
  `remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.ode")`.
- Precedent: kaskr/adcomp is a non-package root with the TMB package
  in a subdirectory; this is the same shape with the dominant package
  kept at the root.

## What moves into extensions/frmtmb.ode

- `R/ode.R` wholesale (frm_ode, the deSolve bridge, the failure log,
  `frm_ode_failures()`).
- `tests/testthat/test-ode.R` and every RTMBode install guard.
- `vignettes/ode.Rmd`.
- DESCRIPTION lines: `RTMBode` leaves core Suggests;
  `Additional_repositories: https://kaskr.r-universe.dev` leaves core
  entirely and lands in frmtmb.ode's DESCRIPTION (where it is an
  Imports-adjacent Suggests, same guard pattern).
- `dev/upstream/` RTMBode patches, if they still exist when this
  executes.

## The two core hooks

frame.R reaches ode.R twice (the nl-body auto-wrap for frm_ode
dynamics rides `ad_overload_fn`). The extraction turns that into a
registration: core exports a small "nl-body transformer" seam (the
same shape as `frmtmb_structure()`: an extension registers a
predicate plus a body rewriter at load time, core consults the
registry where it consults ode.R today). The registry is the export;
frm_ode is its first customer and the worked example. If on
inspection the two hooks are trivial enough that a plain generic
suffices, prefer the generic; the registry earns its keep only if the
hooks carry state.

## The tmbstan refusal string

interop.R's sampler-failure message names frm_ode and RTMBode (the
known ODE-under-tmbstan failure case). After extraction core must not
name the extension. Options, in order of preference: (a) the message
keeps a generic "a tape that calls external solvers" sentence and the
ode package DOCUMENTS the specific failure in its own vignette;
(b) a registered "known sampler failure causes" list extensions
append to, consulted when the refusal fires. Take (a) unless the lane
finds the specific naming is load-bearing for users; (b) is
machinery for one string.

## Ratchet endgame

After extraction, `interop.R|frm_ode` and `interop.R|RTMBode` leave
the pinned inventory. Combined with steps 6 and 7 and the compat-row
decision, the inventory should then be empty or compat-only, and the
boundary test's comment flips from ratchet language to zero-tolerance
language.

## CI

Core workflows unchanged (they run at the root). One added workflow
for the extension with a path filter on `extensions/frmtmb.ode/**`,
installing core from the checkout first (`R CMD INSTALL .`), then
checking the extension against it. The pkgcheck docker gate stays
core-only.

## The draws package (sequencing revised 2026-09-03)

STATUS 2026-09-03: executed. `extensions/frmtmb.sample` exists and
carries the whole sampling surface; the inventory that preceded the
move is `dev/draws-extraction.md`. The tmbstan-refusal-string question
this note left to the draws lane is resolved there (option (a): a
generic "external solver" sentence, with no package named), so core
names neither `frm_ode` nor RTMBode any more. One consequence for the
ratchet is recorded at the end of that document.

Maintainer decision: frmtmb.sample separates SOONER than the
original "at CRAN time" plan, in the same round as the ODE
extraction, because the sampling tests are the heaviest part of the
suite and the split buys iteration speed immediately; development
continues from the monorepo top level either way. Coordination rule
for the round: the draws lane owns R/interop.R outright (its content
mostly leaves), so this note's tmbstan-refusal-string question
resolves in the draws package, not in core; the ODE lane does not
touch interop.R. The draws extraction starts with a call-graph
inventory of every core internal the sampling surface reaches (the
export list is much larger than the structured-family API:
build_objective on the frame, the covstruct factor accessors for
non-centering, the log_lik row-density machinery, sdreport paths for
check_laplace), and that inventory is its first deliverable. The
compat rows for sampling features move onto the registration seam
steps 6-9 built, exactly as hmm's and lca's did. hmm and lca move
out under the structured-family protocol, not this note, and may
share one extension package with a name to be chosen.
