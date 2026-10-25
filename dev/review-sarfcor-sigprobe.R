.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb)})
frmtmb_local_registries(environment())
V <- 0.5^abs(outer(1:20, 1:20, "-"))
b <- function(spec, rows, n_all, lookup) {
  L <- t(chol(V)); list(L = L, Linv = solve(L), logdet = sum(log(diag(L))))
}
wh <- function(y, mu, sigma, theta, aux) {
  message("whiten: length(sigma) = ", length(sigma), " length(mu) = ", length(mu))
  u <- as.vector(aux$Linv %*% ((y - mu) / sigma))
  list(u = u, logjac = -aux$logdet - length(y) * log(sigma[1]))
}
uw <- function(u, mu, sigma, theta, aux) {
  message("unwhiten: length(sigma) = ", length(sigma))
  mu + sigma * as.vector(aux$L %*% u)
}
frmtmb_register_autocor("pf", args = "M", build = b, whiten = wh,
                        unwhiten = uw, sigma_varies = FALSE)
set.seed(1); d <- data.frame(x = rnorm(20)); d$y <- d$x + rnorm(20)
f <- frm(y ~ x + pf(1), data = d)
invisible(fitted(f))
