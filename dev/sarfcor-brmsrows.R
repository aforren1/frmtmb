# The ported brms-suite rows that use sar() or fcor(): bin 2,
# dev/brmsport-bin2-ledger.tsv rows log_lik:150, 152, 157, 159, 174,
# 176, 230 and posterior_predict:111, 113, 118, 120, 133, 135. All 13
# were "cannot transfer / absent / gap" at 0.74.0 because frm() refused
# both terms. This replays what each row asserts through frmtmb's public
# API on a fit of the same shape (10 observations; 3 draws), and prints
# the outcome each row would now get. It edits no ledger.
# Run: Rscript dev/sarfcor-brmsrows.R > dev/sarfcor-out/brmsrows.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.spatial"),
    "|", find.package("frmtmb.sample"), "\n")
set.seed(1)
n <- 10
d <- data.frame(y = rnorm(n), x = rnorm(n))
# brms's rows use Msar = diag(10), whose eigenvalues are all 1
I10 <- diag(n)
# a path graph, a weights matrix with a zero diagonal
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
Wp <- P / rowSums(P)
row <- function(id, outcome, note) {
  cat(sprintf("%-22s %-10s %s\n", id, outcome, note))
}
try_fit <- function(f, fam, d2) {
  tryCatch(frm(f, family = fam, data = d, data2 = d2),
           error = function(e) e)
}
# posterior_predict rows: dim(pred) == c(ndraws, nobs) with 3 draws
pp_dims <- function(fit) {
  s <- suppressMessages(suppressWarnings(
    frm_sample(fit, chains = 1, iter = 200, warmup = 100, seed = 1,
               refresh = 0)))
  dim(posterior_predict(s, ndraws = 3))
}
for (cs in list(
  list(id = "posterior_predict:111", f = y ~ x + sar(M),
       fam = gaussian()),
  list(id = "posterior_predict:113", f = y ~ x + sar(M),
       fam = student()),
  list(id = "posterior_predict:118", f = y ~ x + sar(M, type = "error"),
       fam = gaussian()),
  list(id = "posterior_predict:120", f = y ~ x + sar(M, type = "error"),
       fam = student()))) {
  lit <- try_fit(cs$f, cs$fam, list(M = I10))
  if (inherits(lit, "error")) {
    row(cs$id, "literal", paste("M = diag(10) refused:",
                                substr(conditionMessage(lit), 1, 90)))
  }
  fit <- try_fit(cs$f, cs$fam, list(M = Wp))
  dm <- tryCatch(pp_dims(fit), error = function(e) conditionMessage(e))
  row(cs$id, if (is.numeric(dm) && all(dm == c(3, n))) "pass" else "fail",
      paste("path-graph W: dim(posterior_predict) =",
            paste(dm, collapse = " x ")))
}
for (cs in list(
  list(id = "posterior_predict:133", fam = gaussian()),
  list(id = "posterior_predict:135", fam = student()))) {
  fit <- try_fit(y ~ x + fcor(M), cs$fam, list(M = I10))
  dm <- tryCatch(pp_dims(fit), error = function(e) conditionMessage(e))
  row(cs$id, if (is.numeric(dm) && all(dm == c(3, n))) "pass" else "fail",
      paste("Mfcor = diag(10): dim(posterior_predict) =",
            paste(dm, collapse = " x ")))
}
# log_lik rows: brms's pointwise log_lik is a leave-one-out conditional
# density; frmtmb gives none for a residual covariance term
ll_try <- function(fit) {
  s <- suppressMessages(suppressWarnings(
    frm_sample(fit, chains = 1, iter = 200, warmup = 100, seed = 1,
               refresh = 0)))
  tryCatch({
    dim(log_lik(s))
  }, error = function(e) paste("refused:", substr(conditionMessage(e), 1, 110)))
}
for (cs in list(
  list(id = "log_lik:150", f = y ~ x + sar(M), fam = gaussian()),
  list(id = "log_lik:152", f = y ~ x + sar(M), fam = student()),
  list(id = "log_lik:157", f = y ~ x + sar(M, type = "error"),
       fam = gaussian()),
  list(id = "log_lik:159", f = y ~ x + sar(M, type = "error"),
       fam = student()),
  list(id = "log_lik:174", f = y ~ x + fcor(M), fam = gaussian()),
  list(id = "log_lik:176", f = y ~ x + fcor(M), fam = student()),
  list(id = "log_lik:230", f = y ~ x + fcor(M), fam = student()))) {
  m <- if (grepl("sar", deparse(cs$f))) Wp else I10
  fit <- try_fit(cs$f, cs$fam, list(M = m))
  r <- ll_try(fit)
  row(cs$id, if (is.character(r)) "refused" else "pass",
      paste(r, collapse = " x "))
  # the ML fit's own log_lik
  r2 <- tryCatch(dim(log_lik(fit)), error = function(e) {
    paste("refused:", substr(conditionMessage(e), 1, 110))
  })
  row("", "core", paste(r2, collapse = " x "))
}
