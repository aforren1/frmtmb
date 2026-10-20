# Probe 4: the exact blockers of the G3 and G4 rows left "not ported
# (cost)": is each model one frmtmb fits, and does log_lik() of draws
# answer for it? Log: dev/bin2port-out/probe4.txt. Seed 4.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.sample",
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(4)
n <- 10
dd <- data.frame(y = rnorm(n), x = rnorm(n), t = 1:n, g = 1)
W <- matrix(0, n, n)
W[cbind(1:(n - 1), 2:n)] <- 1
W <- W + t(W)
show <- function(label, expr) {
  r <- tryCatch({
    force(expr)
    "fits"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(label, ":", r, "\n")
}
show("sar(W) lag", frm(bf(y ~ x + sar(W, type = "lag")), data = dd,
                        data2 = list(W = W)))
show("sar(W) error", frm(bf(y ~ x + sar(W, type = "error")), data = dd,
                          data2 = list(W = W)))
show("fcor(M)", frm(bf(y ~ x + fcor(M)), data = dd,
                     data2 = list(M = diag(n))))
# ARMA with cov = TRUE fits; log_lik() of draws takes its refusal from
# draws_loglik_factors(), which reads only the fit, before any draw
fit <- frm(bf(y ~ x + arma(time = t, gr = g, p = 1, q = 1, cov = TRUE)),
           data = dd)
cat("arma(cov = TRUE) fit converged:", fit$opt$convergence == 0, "\n")
r <- tryCatch(frmtmb.sample:::draws_loglik_factors(fit, "log_lik()"),
              error = conditionMessage)
cat("log_lik() of draws on arma(cov = TRUE):", r, "\n")
fit0 <- frm(bf(y ~ x + arma(time = t, gr = g, p = 1, q = 1)), data = dd)
r <- tryCatch({
  frmtmb.sample:::draws_loglik_factors(fit0, "log_lik()")
  "no refusal"
}, error = conditionMessage)
cat("log_lik() of draws on arma(cov = FALSE):", r, "\n")
# a multivariate gaussian with residual correlation fits
dm <- data.frame(y1 = rnorm(30), y2 = rnorm(30), y3 = rnorm(30))
fm <- frm(bf(mvbind(y1, y2, y3) ~ 1) + set_rescor(TRUE), data = dm)
cat("mvbind rescor fit converged:", fm$opt$convergence == 0, "\n")
r <- tryCatch({
  frmtmb.sample:::draws_loglik_factors(fm, "log_lik()")
  "no refusal"
}, error = conditionMessage)
cat("log_lik() of draws on a rescor model:", r, "\n")
# what frmtmb.eam's wiener simulator returns through sim_draw()
suppressMessages(library(frmtmb.eam))
w <- brmsfamily("wiener")
s <- sim_draw(list(family = w, dpars = list(mu = rnorm(5), bs = rep(1.5, 5),
                                            ndt = rep(0.3, 5),
                                            bias = rep(0.5, 5)),
                   aterms = list(), n = 5, autocor = NULL, extra = NULL,
                   max_iter = NULL))
cat("wiener sim_draw class:", class(s), " dim:", dim(s), "\n")
print(s)
