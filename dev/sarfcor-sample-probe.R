# frmtmb.sample on fcor() and sar() fits: sampling, draws names, the
# pointwise log_lik refusal, posterior_predict and posterior_epred.
# Seed 3. Run: Rscript dev/sarfcor-sample-probe.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.sample"),
    "|", find.package("frmtmb.spatial"), "\n")
set.seed(3)
g <- expand.grid(r = 1:6, c = 1:6)
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
V <- 0.5^abs(outer(1:n, 1:n, "-"))
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n)))
for (f in list(y ~ x + sar(W), y ~ x + fcor(V))) {
  fit <- frm(f, data = d, data2 = list(W = W, V = V))
  cat("====", deparse(f), "\n")
  s <- tryCatch(frm_sample(fit, chains = 2, iter = 600, warmup = 300,
                           seed = 3),
                error = function(e) {
                  cat("frm_sample ERROR:", conditionMessage(e), "\n")
                  NULL
                })
  if (is.null(s)) next
  print(utils::head(colnames(as.matrix(s)), 12))
  print(utils::head(posterior::summarise_draws(posterior::as_draws_df(as.matrix(s))), 6))
  r <- tryCatch(dim(log_lik(s)), error = function(e) {
    paste("log_lik refused:", substr(conditionMessage(e), 1, 160))
  })
  print(r)
  r <- tryCatch(dim(posterior_predict(s)), error = function(e) {
    paste("posterior_predict ERROR:", conditionMessage(e))
  })
  print(r)
  r <- tryCatch({
    ep <- posterior_epred(s)
    c(dim(ep), max(abs(colMeans(ep) - fitted(fit)[, 1])))
  }, error = function(e) paste("posterior_epred ERROR:", conditionMessage(e)))
  print(r)
}
