# RTMB pitfalls met while building frmtmb

Each item below cost real time once. Read this before you write or
review objective code, a covariance block, a custom family, or a
hand-rolled reference in a test. Items 1 to 13 and 16 to 18 are about
the tape. Items 14, 15, 19 and 20 are about the toolchain around it.

This file was moved into the repository on 2026-09-28 from a
machine-local memory so that it survives a change of machine.

## Class stripping

1. `base::matrix()` strips the advector class. In a script,
   `library(RTMB)` masks it, so code works there and breaks inside a
   package namespace. Use `RTMB::matrix()` and `RTMB::diag()` in AD
   code. `t()`, `rep(x, k)`, `as.vector()` and `%*%` are generics with
   registered methods and dispatch correctly.
2. `base::c()` on mixed numeric and advector input strips the class the
   same way. Any closure that concatenates advectors, including a test
   reference, needs a local `"c" <- RTMB::ADoverload("c")`.
3. `rep(adv, each = n)` and `rowSums(adv_matrix)` strip the class. Item
   1 holds only for plain `rep(x, k)`. Broadcast with
   `RTMB::matrix(1, n, 1) %*% RTMB::matrix(v, 1, K)`. Row-sum with
   `as.vector(M %*% rep(1, K))`. `cumsum(adv)` works.
4. `ADoverload("[<-")` is lexically scoped. Every helper that
   sub-assigns advectors needs its own local
   `"[<-" <- RTMB::ADoverload("[<-")`.
5. `drop()` does not collapse a `dgeMatrix` in numeric mode.
   `as.vector()` collapses both Matrix and advector results, which keeps
   the objective evaluable in plain-numeric mode. That matters: a
   numeric `nll(template)` against the taped `obj$fn` is how a tape bug
   is localized.

## R semantics that bite on the tape

6. `$` partial matching: `pars$b` silently matches `pars$beta` when `b`
   is absent, which is the GLM case. Use `[[ ]]` for every
   parameter-list access.
7. A `for` loop over an empty sequence sets its variable to `NULL`,
   which clobbers a same-named function argument. This bit `spec`.
8. `RTMB::OBS(x)` registers observations under the DEPARSED ARGUMENT
   EXPRESSION. `OBS(y[[r]])` in a loop over responses collides on the
   identical expression and silently swaps data between responses. The
   tape value is wrong and nothing errors. frmtmb applies `OBS()` only
   to univariate, non-matrix responses.
9. `reformulas::splitForm` silently DROPS any fixed term whose
   expression tree mentions a special name, such as `exp(x)` or even
   `I(exp(x))`, not only `s()` and `t2()`. frmtmb protects against this
   with `.frm_<name>` alias substitution in `parse_linpred`
   (`sub_specials`), with aliases bound in per-linpred and
   response-shared environments.
10. `RTMB::qgamma` takes a RATE positionally where a reader expects a
    scale. A "scale" coefficient in an nl body fits as 1/scale. The
    likelihood is identical, so only a cross-check against an external
    parameterization has to invert it.

## Numerics

11. Tape-safe floors. The maximum form `0.5 * (x + lo + abs(x - lo))`
    with a tiny `lo` (1e-300) is annihilated by rounding whenever
    |x| >> lo. For x = -1e-17 it returns exactly 0, not `lo`, so the
    floor fails in the one regime it exists for, and `log(0)` then
    poisons a mixture's gradient. Write floors as the positive part
    plus the floor, `0.5 * (x + abs(x)) + lo`, or as
    `lo + 0.5 * ((x - lo) + abs(x - lo))`. Both keep `lo` outside the
    cancelling sum. `0.5 * (e + abs(e))` for max(e, 0) and
    `0.5 * (x + cap - abs(cap - x))` for min(x, cap) are exact. A
    difference of two `smin()` evaluations of the same cap rounds
    NEGATIVE in about 9 percent of draws, so floors must tolerate tiny
    negatives, not only exact zero. Found independently by two lanes on
    2026-09-04.
12. Matern `besselK`: AD works in both arguments, but raw `exp()`
    transforms for range and shape let the optimizer reach `0 * Inf`
    and a NaN gradient, which kills `nlminb`. Bound the transforms with
    a logistic shape cap and a range floor taken from the distance
    scale. glmmTMB's `mat` dies the same way.
13. `TapeConfig(comparison = "tape")` is the documented RTMB mechanism
    for indicator functions on the tape, as in
    `f(x) * (x < 0) + g(x) * (x >= 0)`. The default `"forbid"` errors
    on AD comparisons. frmtmb has never adopted it and writes everything
    branch-free with `abs()` and `tanh()`. A lane that says "no
    tape-safe conditional exists" describes the convention, not RTMB.
    It is a global setting, so adopting it is an ecosystem policy
    decision, set and restored around tape builds. It would unlock
    regime switches such as the LBA pnorm-tail fix.

## Performance

14. Elementwise `a[i] <- x` in a taping loop is slow at TAPE BUILD, not
    at run time, and the penalty grows with n because each
    sub-assignment copies the vector: about 1.3x at n = 1000 and about
    25x at n = 100000 on the perf canary's Poisson GLMM. Gradient cost
    is the same either way. Observation-length loops are forbidden in
    `objective.R` and `covstruct.R`, enforced by
    `tests/testthat/test-perf.R` as a wall-clock bound. Vectorize across
    units and loop over the short dimension, as the
    reinforcement-learning vignette does.
15. `mkReTrms` `Zt` row ordering is level-major, so `matrix(b, nrow = d)`
    has one level per column. Pass `reorder.terms = FALSE` to keep the
    user's term order, as glmmTMB does.

## Simulation and sampling

16. `obj$simulate()` replays the objective with simref objects. The
    `dim<-` reshape is simref-safe and covariance blocks use it, but a
    sparse `Z %*% simref` product is not supported, so `obj$simulate` is
    dead for this architecture. The numeric `simulate()` is final.
17. `tmbstan(cores > 1)` on Windows works with frmtmb objectives.
    tmbstan writes a temporary `R_PROFILE` for the PSOCK workers and
    retapes from the serialized closure, and frmtmb's generated closures
    are self-contained, so every route works, seed-identical to
    sequential, with about 9 s of fixed startup. The one real failure: a
    `devtools::load_all()` namespace cannot be rebuilt on a worker, and
    chains die with "Stan model 'tmb_generic' does not contain samples".
    Parallel anything needs the INSTALLED package.
    `test-parallel-chains.R` skips under pkgload for that reason.

## Toolchain

18. PowerShell `-e` with `$` inside double quotes mangles R code even
    with backslash escapes. Debug scripts go in files, always.
19. `knitr:::tangle_block()` extracts a `read_chunk()` call with the
    regex `read_chunk\(([^)]+)\)`, so a nested call such as
    `knitr::read_chunk(bcm_file("marginal.R"))` is cut at the first `)`
    and fails to parse. pkgdown 2.2 tangles every article, so the whole
    article build dies with "unexpected end of input" while
    `rmarkdown::render()` and `R CMD check` pass. Fix: `purl = FALSE` on
    the chunk that calls `read_chunk()`.
20. The session scratchpad is shared by every lane and reviewer, and
    their R processes all carry the same session GUID in their paths. A
    reviewer that stopped "its" jobs by matching the GUID killed another
    lane's as-cran run and truncated a shared suite log twice. Kill
    only processes whose command line carries your own lane prefix.
    Never truncate or delete a file you did not create. Before you
    declare an agent dead, list processes by every prefix it may have
    used and by its worktree path.
