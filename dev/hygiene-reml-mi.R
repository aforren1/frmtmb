# Item 1222: does REML run with mi(), with and without a predictor mi(x)
# reading the latent values, and does the latent-value case behave as an
# approximation? Seed 7, n = 120, 12 groups, 20 missing x.
#   Rscript dev/hygiene-reml-mi.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(7)
n <- 120L
d <- data.frame(g = factor(rep(1:12, each = 10)), z = rnorm(n))
d$x <- 0.5 + 0.8 * d$z + rnorm(n, 0, 0.6)
d$y <- 1 + 0.7 * d$x + rnorm(12, 0, 0.5)[d$g] + rnorm(n, 0, 0.8)
d$x[sample(n, 20)] <- NA
fm <- bf(y ~ mi(x) + (1 | g)) + bf(x | mi() ~ z) + set_rescor(FALSE)
r <- frm(fm, data = d, REML = TRUE)
cat("REML with mi(x) reading 20 latent values: code",
    r$opt$convergence, "logLik", format(as.numeric(logLik(r)), digits = 10),
    "\n")
# the same with the missing rows' x filled in: no latent value, REML is
# the exact restricted likelihood of a linear gaussian model
d2 <- d
d2$x[is.na(d2$x)] <- 0.5 + 0.8 * d2$z[is.na(d2$x)]
r2 <- frm(fm, data = d2, REML = TRUE)
cat("REML with mi(x), nothing missing: code", r2$opt$convergence, "\n")
# a missing response no density reads leaves the objective
d3 <- d2
d3$y[1:5] <- NA
r3 <- frm(bf(y | mi() ~ x + (1 | g)), data = d3, REML = TRUE)
r3c <- frm(bf(y ~ x + (1 | g)), data = d3[-(1:5), ], REML = TRUE)
cat("REML y | mi() with no reader minus complete case logLik:",
    format(as.numeric(logLik(r3)) - as.numeric(logLik(r3c)), digits = 3),
    "\n")
print(frm_compat("REML", "mi()")$status)
print(frm_compat("REML", "mi_pred()")$status)
# Laplace is exact only when the integrand is gaussian in what it
# integrates, that is when the inner Hessian does not move with the
# random vector. Compare it at the inner mode and at a moved point.
hess_move <- function(fit) {
  env <- fit$obj$env
  p <- env$last.par.best
  # spHess() writes into one sparse matrix in place, so copy it out
  h0 <- as.matrix(env$spHess(p, random = TRUE))
  q <- p
  q[env$random] <- q[env$random] + 0.3
  h1 <- as.matrix(env$spHess(q, random = TRUE))
  max(abs(h1 - h0)) / max(abs(h0))
}
cat("inner Hessian relative move, 20 latent x:", format(hess_move(r), digits = 3),
    "\n")
cat("inner Hessian relative move, nothing missing:",
    format(hess_move(r2), digits = 3), "\n")
for (fit in list(r, r2)) {
  env <- fit$obj$env
  cat("random:", paste(names(table(names(env$last.par.best)[env$random])),
                       table(names(env$last.par.best)[env$random]),
                       collapse = " "), "\n")
  cat("outer:", paste(names(fit$opt$par), collapse = " "), "\n")
}
# the instrument's positive control: a poisson GLMM is not gaussian in b
dp <- data.frame(g = factor(rep(1:12, each = 10)), x = rnorm(n))
dp$k <- rpois(n, exp(0.3 + 0.4 * dp$x + rnorm(12, 0, 0.5)[dp$g]))
rp <- frm(bf(k ~ x + (1 | g)), family = poisson(), data = dp)
cat("control, poisson GLMM:", format(hess_move(rp), digits = 3), "\n")
