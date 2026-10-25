# Reviewer: the help page example of frmtmb_register_autocor(), verbatim,
# then fitted.
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb)})
frmtmb_local_registries(environment())
myfix_build <- function(spec, rows, n_all, lookup) {
  V <- lookup(spec$V, "myfix()")[rows, rows, drop = FALSE]
  L <- t(chol(V))
  list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
       logdet = sum(log(diag(L))))
}
myfix_whiten <- function(y, mu, sigma, theta, aux) {
  u <- as.vector(aux$Linv %*% ((y - mu) / sigma))
  list(u = u, logjac = -aux$logdet - length(y) * log(sigma))
}
myfix_unwhiten <- function(u, mu, sigma, theta, aux) {
  mu + sigma * as.vector(aux$L %*% u)
}
frmtmb_register_autocor("myfix", args = "V", build = myfix_build,
                        whiten = myfix_whiten,
                        unwhiten = myfix_unwhiten, sigma_varies = FALSE)
set.seed(5); n <- 24
V <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
d <- data.frame(x = rnorm(n)); d$y <- 1 + d$x + as.vector(t(chol(V)) %*% rnorm(n))
r <- tryCatch({f <- frm(y ~ x + myfix(V), data = d, data2 = list(V = V))
  f0 <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  sprintf("FIT logLik %.8f, fcor %.8f", logLik(f), logLik(f0))},
  error = function(e) paste("ERROR:", conditionMessage(e)))
cat("help-page example, sigma_varies = FALSE:", r, "\n")
