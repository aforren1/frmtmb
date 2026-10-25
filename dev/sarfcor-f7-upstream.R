# Fix round 1, F7: brms-36 to brms-38 and two notes in
# dev/upstream-bugs.md. Rscript <this>
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/upstream-bugs.md"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
entries <- "### brms-36. `sar()` passes complex eigenvalues of a row-standardized matrix to Stan

- Seen: 2.23.0, 2026-10-24 (`dev/sarfcor-upstream-eigen.R`, log
  `dev/sarfcor-out/upstream-eigen.txt`).
- Repro: an 8 x 8 rook lattice, `W <- A / rowSums(A)`,
  `brm(y ~ x + sar(W), data2 = list(W = W))`. `eigen(W)` returns its
  values as complex, max |Im| 7.97e-17, and `standata()` keeps
  `eigenMsar` complex.
- Expected: the real parts. `W` is similar to the symmetric
  `D^-1/2 A D^-1/2`, so its eigenvalues are real. Actual: Stan stops in
  data initialization, \"variable does not exist ... variable
  name=eigenMsar; base type=double\", and `brm()` returns a fit with no
  draws. A row-standardized matrix is spdep's default and the usual
  `sar()` input.
- Workaround: frmtmb.spatial takes the real part when the imaginary
  part is below 1e-8 of the largest modulus, and refuses a matrix whose
  eigenvalues are complex beyond that. The brms rows of frmtmb.spatial's
  tests use the symmetric normalization.
- Record: `dev/sarfcor-findings.md`; `dev/reviews/2026-10-24-sarfcor.md`,
  finding 7.
- Status: not reported.

### brms-37. `sar()` with `weights()`, `se()`, `cens()` or `trunc()` writes Stan code that does not parse

- Seen: 2.23.0, 2026-10-24 (`dev/sarfcor-upstream-stanc.R`, log
  `dev/sarfcor-out/upstream-stanc.txt`).
- Repro: `stancode(y | weights(w) ~ x + sar(W), ...)` and the same with
  `se(se, sigma = TRUE)`, `cens(cen)` or `trunc(lb = lb)`, on a
  row-standardized 4 x 4 lattice, then
  `rstan::stanc(model_code = code)`.
- Expected: a refusal from brms, as for the other residual correlation
  terms, or a program that compiles. Actual: stanc gives \"Semantic
  error\" at the `normal_lagsar_lpdf` call for each of the four; the
  plain `y ~ x + sar(W)` type-checks. The user sees a Stan compiler
  error, not a message about the model.
- Workaround: frmtmb refuses each of the four with the term and the
  reason (the density is joint over the rows).
- Record: `dev/sarfcor-findings.md`; `dev/reviews/2026-10-24-sarfcor.md`,
  finding 7.
- Status: not reported.

### brms-38. `posterior_epred()` and `fitted()` fail on every `sar(type = \"error\")` fit

- Seen: 2.23.0, 2026-10-24 (`dev/sarfcor-upstream-epred.R`, log
  `dev/sarfcor-out/upstream-epred.txt`).
- Repro: a 5 x 5 lattice with the symmetric normalization,
  `fit <- brm(y ~ x + sar(Ws, type = \"error\"), ...)`, then
  `posterior_epred(fit)` or `fitted(fit)`.
- Expected: `mu`, since the error form does not move the mean. Actual:
  \"'a' is 0-diml\". `prepare_predictions_ac()` reads `lagsar` with
  `prepare_draws()`, which gives a matrix with no columns (200 x 0)
  rather than `NULL` when the fit has no `lagsar`, so
  `posterior_epred_gaussian()` takes the lag branch and `solve()` stops.
  `posterior_predict()` and `log_lik()` work.
- Workaround: none needed; frmtmb's `fitted()` gives `mu` for the error
  form. frmtmb.spatial's `posterior_epred()` comparison runs on the lag
  form only.
- Record: `dev/sarfcor-findings.md`, fix round 1 (F3).
- Status: not reported.

## RTMB
"
rep1("\n## RTMB\n", paste0("\n", entries))
notes <- "- brms `fcor()` checks that `M` is positive definite with
  `min(eigen(M)$values <= 0)`, which passes a matrix with a negative
  eigenvalue (`dev/review-sarfcor-brms-nocompile.R`). Not a wrong
  model: the Stan program takes `cholesky_decompose(Mfcor)` in
  transformed data, and the fit fails there. brms fails late, where
  frmtmb refuses at the start.
- brms `sar()` bounds `rho` by the eigenvalues of `M` themselves. For a
  matrix that is not row-standardized (a 0/1 adjacency), that interval
  reaches past `1 / max eigenvalue`, where `log1m(rho * eigenMsar)` is
  not a number and Stan rejects the proposal, so brms's support is
  frmtmb's interval, at the cost of rejected proposals and initial
  values. For a row-standardized `M` the interval is narrower than the
  valid one: a deliberate divergence, recorded in frmtmb.spatial's
  NEWS and `?sar`.
- System GSL for rtdists on CI: a system requirement."
rep1("- System GSL for rtdists on CI: a system requirement.", notes)
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
