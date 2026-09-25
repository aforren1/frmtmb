# brms parity round, 2026-09-25 (frmtmb 0.64.0)

This round closes the gaps that the model menu listed after 0.63.0. It
is a record of what was done, how it was verified, and what is left.
Each lane's `dev/<lane>-findings.md` has its own validation, numbers
and scripts.

## Lanes

| lane | gap | findings |
|---|---|---|
| me | `me()` noise-free terms, `set_mecor()` | `dev/me-findings.md` |
| thres | `thres()` for the ordinal families | `dev/thres-findings.md` |
| grby | `gr(g, by = )`, `mm(g1, g2, by = )` | `dev/grby-findings.md` |
| arcov | `ar()`, `ma()`, `arma()` without `cov = TRUE` | `dev/arcov-findings.md` |
| mv | three or more responses with `+`; `student` with `rescor`; ordinal families in a multivariate model | `dev/mv-findings.md` |
| icpt0 | `0 + Intercept`, `bf(center = FALSE)` | `dev/icpt0-findings.md` |
| fams | `hurdle_negbinomial()`, `zero_one_inflated_beta()` | `dev/fams-findings.md` |
| emm | emmeans on nonlinear and multivariate fits | `dev/emm-findings.md` |
| sratio | `sratio()` thresholds unordered (found by lane thres) | `dev/sratio-findings.md` |

Each lane worked in its own worktree on the branch `lane/<name>`. The
consolidating session reviewed each diff, reran the lane's test files
on the lane build, did one independent check per lane, and merged the
lanes with `--no-ff`.

Independent checks by the consolidating session:

- **emm:** a multivariate fit without `resp` matches the univariate
  fit to 6.8e-6 relative.
- **thres:** a grouped fit with per-group slopes equals `MASS::polr()`
  per group to 2.7e-12. The groups have 3 and 2 thresholds.
- **fams:** a hurdle negative binomial with `(1 | g)` equals glmmTMB's
  `truncated_nbinom2` hurdle to 1.3e-11.
- **grby:** `(1 + x | gr(g, by = f))` equals lme4 fitted with one
  indicator term per by-level, to 11 digits.

## Defects the merge created, fixed on the branch

A lane can only test its own feature. These defects appeared where two
lanes met, or in test files a lane did not run.

| defect | cause | pinned by |
|---|---|---|
| `ordinal_ncat()` read the first response's family | thres and mv, a clean textual merge | `test-parity-integration.R` |
| `cs()` in a multivariate ordinal response had zero coefficients | thres changed the count to `extras[["tau_raw"]]`, and mv renamed that per response | `test-parity-integration.R` |
| `lp$center` read with `$` | icpt0 | `test-bracket-access.R` |
| two refusals shared one message | arcov | `test-message-uniqueness.R` |
| a nonlinear body with `m[, 1]` stopped with 'argument "ei" is missing' | me first called the older `calls_function()` on a body | `test-nl-body-vars.R`, `test-nl-lexical.R` |
| nine stale manual port verdicts, and three reasons that regeneration reverted | emm and icpt0 edited only the generated verdicts file | `dev/brmsport-ledger.R` |

Both `test-parity-integration.R` pins were seen failing on a build
without the fixes.

## Verification on the release head

The machine is a Linux container. R 4.5.3 comes from conda-forge,
because the proxy blocks CRAN. RTMB 2.0, RTMBdist and tmbstan come from
their GitHub sources. Stan compiles against RcppParallel 5.1.9 in a
separate library, because the conda build of RcppParallel lacks
`tbb/tbb_stddef.h`. The reference log-likelihoods of
`dev/machine-library.md` reproduce exactly: -295.602189818 and
-332.137876329.

- **Suite**, core and all seven extensions, one file per process: 305
  files, 17,815 assertions. Two files fail, and both fail identically
  on the unmodified base build here.
  - `test-pp-check-types.R` (4): conda's bayesplot 1.15.0 has 42 ppc
    types and no `ppc_dots()`.
  - frmtmb.sample `test-reparam.R` (1): a sampler comparison on the
    correlated model.
- **Gated tier**: 54 files, 3,767 assertions, 0 fail, 0 error. This
  includes `test-brms-likelihood.R` (428) and `test-brms-methods.R`
  (979) with Stan compiling. Run separately with their own gates:
  `test-drmtmb-agreement.R` (131) and `test-fuzz.R`, 0 skip. The
  remaining skips are the scale tier, which the release harness does
  not run, and reference packages that are not installed here:
  brokenstick, fmesher, rtdists, hmmTMB, depmixS4, gratia and rxode2.
  None of those reach code changed in this round.
- **Ported brms suite**: recorded on the release build, then
  regenerated with `dev/brmsport-ledger.R` and the generator. The
  regenerated verdicts file and every generated test file are
  byte-identical to the committed ones. Bin 1 passes 275 of 494, up
  from 252.
- **R CMD build**: OK, including the vignettes.
- **R CMD check --as-cran** runs with `--no-manual`, because there is
  no LaTeX here. The in-check tests fail only on
  `test-pp-check-types.R`, as above: 9,098 pass, 4 fail. The warnings
  and notes come from the container: the locale, a missing `qpdf`, no
  access to the CRAN index, and the clock check. One warning was real,
  an undeclared `ordinal::` in `test-thres.R`, and `ordinal` is now
  suggested.

## Left open

- **Adding a family to a multivariate formula.** In
  `bf(o ~ x) + cumulative() + bf(y ~ x) + gaussian()`, frmtmb fills
  only the responses that have no family, so `o` stays ordinal. brms
  gives the last family to every response. frmtmb already behaved this
  way before this round, and the behavior is now documented. Decide
  whether it should follow brms.
- **Priors on `me()` hyperparameters.** `meanme`, `sdme` and `corme`
  priors are refused, not implemented (lane me).
- **REML with `me()` in `mu`** is approximate, and is registered as
  conditional. The same argument applies to the existing `mi()`
  predictor, which is registered as "works" (lane me).
- **The convergence check fires on correct fits.** "Large maximum
  absolute gradient" appears on ordinal and multivariate fits whose
  likelihood identities hold to 1e-12, and on fits with an active
  bound. The 1e-3 threshold is absolute (lanes thres, mv, arcov).
- **Refits can lose a threshold.** `frm_bootstrap()` and other refits
  recount the thresholds on each simulated data set (lane thres).
- **`y ~ x + cs(x)` is not identified and is not refused**
  (lane sratio).
- **`draw_prior_entry()` on an unordered threshold vector** may copy
  one draw into every threshold. It is not known to be reachable
  (lane sratio).
- **Not yet implemented:** frmtmb.sample `log_lik()` and `loo()` for
  `cov = FALSE` ARMA, and brms's latent-residual AR for non-gaussian
  families (lane arcov). Also `gr(g, by = f, cov = A)` (lane grby).
- **Older record inconsistencies:** frmtmb.eam's tests call
  `frmtmb.sample::` without declaring it. `codemeta.json` still says
  0.50.0.
- **Upstream brms defect:** `posterior_predict_hurdle_negbinomial()`
  does not draw the zero-truncated negative binomial (lane fams).
