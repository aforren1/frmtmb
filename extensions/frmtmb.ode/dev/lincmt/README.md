# The frm_lincmt() measurement scripts

Every number in `dev/lincmt-findings.md` (repository root) was produced
by a script in this directory. They are here rather than in a session
scratchpad because a scratchpad has been cleaned mid-round before and
took a process document with it.

They are RUN scripts, not tests: they print tables rather than assert,
and several take minutes. The assertions that guard the same behavior
in the shipped package are in `tests/testthat/test-lincmt.R`.

## What each one does

| script | what it measures |
|---|---|
| `lincmt-src.R` | sourced by the others: sets `.libPaths()` and reads `R/ode.R` and `R/lincmt.R` directly, for iteration without an install |
| `lincmt-proto.R` | the primitives on their own, against `Matrix::expm` |
| `lincmt-degen.R` | the disposition coefficients as rate constants converge, against `Matrix::expm` |
| `lincmt-ident.R` | 118 schedules against `frm_ode()` at `atol = rtol = 1e-12` |
| `lincmt-edge.R` | nine schedule-semantics corners, and how far `frm_ode()`'s `n_ss` run-in is from the exact steady state |
| `lincmt-mpfr.R` | against a 240-bit `Rmpfr` reference, through every coalescence down to exact equality, plus a random parameter box |
| `lincmt-ad.R` | value, gradient and Hessian on the tape, against `frm_ode()`'s adjoint and against central differences |
| `lincmt-seenfail.R` | the gradient test run against the UNFIXED primitive, to record the failure it pins |
| `lincmt-coverage.R` | which NONMEM-shaped records `frm_lincmt()` accepts and which it refuses |
| `lincmt-time-small.R` | interleaved gradient timing with a control arm, at a subject count the solver can afford |
| `lincmt-phase0.R` | the Phase 0 design, before and after, in one process |
| `lincmt-postfit.R` | `fitted()`, `predict()`, `simulate()`, `diagnose()` through a `frm_lincmt()` predictor |
| `lincmt-sample.R` | `frmtmb.sample::frm_sample()` on a `frm_lincmt()` model |

## Running them

They hold absolute library paths from the lane that wrote them:

    .libPaths(c("<a private library holding frmtmb.ode>",
                "<a library holding frmtmb, RTMB, RTMBode, deSolve>",
                "<the user library>"))

Edit that line, or set `R_LIBS`. `lincmt-mpfr.R` needs **Rmpfr**;
`lincmt-sample.R` needs **frmtmb.sample** and its Stan toolchain; the
rest need only what the package already suggests.
