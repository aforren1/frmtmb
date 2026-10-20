# Probe 3: can the thin adapter reach frmtmb.eam's wiener family the way
# brms's log_lik_wiener(i, prep) reads its prep (tests.log_lik.R:555)?
# Log: dev/bin2port-out/probe3.txt. Seed 555.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.eam))
cat("frmtmb.eam", format(packageVersion("frmtmb.eam")), "\n")
f <- tryCatch(brmsfamily("wiener"), error = conditionMessage)
if (is.character(f)) cat("brmsfamily('wiener'):", f, "\n") else {
  cat("brmsfamily('wiener') dpars:", f$dpars, " accepts:", f$accepts_aterms,
      " finalize:",
      !is.null(f$family_finalize), "\n")
}
w <- wiener(link_ndt = "log")
cat("wiener(link_ndt = 'log') dpars:", w$dpars, " accepts:",
    w$accepts_aterms, "\n")
print(names(formals(w$lpdf)))
set.seed(555)
ns <- 5
nobs <- 3
prep <- list(ndraws = ns, nobs = nobs)
prep$dpars <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                   bs = rchisq(ns, 3), ndt = rep(0.5, ns),
                   bias = rbeta(ns, 1, 1))
prep$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
i <- 2
b <- brms:::log_lik_wiener(i, structure(prep, class = "brmsprep"))
dp <- list(mu = prep$dpars$mu[, i], bs = prep$dpars$bs,
           ndt = prep$dpars$ndt, bias = prep$dpars$bias)
y <- rep(prep$data$Y[i], ns)
for (fam in list(w, if (!is.character(f)) f)) {
  if (is.null(fam)) next
  a <- tryCatch(as.numeric(row_lpdf(fam, y, y, dp,
                                    list(dec = rep(prep$data$dec[i], ns)),
                                    NULL)),
                error = conditionMessage)
  print(a)
}
print(b)
