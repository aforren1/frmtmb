# Upstream repro (dev/upstream-bugs.md): brms 2.23.0 passes eigen()'s
# values of a row-standardized sar() matrix to Stan as they come. Such a
# matrix is similar to a symmetric one, so its eigenvalues are real, but
# eigen() of the unsymmetric matrix returns them as complex with a
# rounding imaginary part, and the fit then fails. Brms only, no frmtmb.
# Rscript dev/sarfcor-upstream-eigen.R  (compiles one Stan program)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
k <- 8
g <- expand.grid(r = 1:k, c = 1:k)
n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
ev <- eigen(W)$values
cat(sprintf("eigen(W): complex %s, max |Im| %.3g\n", is.complex(ev),
            max(abs(Im(ev)))))
set.seed(5)
d <- data.frame(y = rnorm(n), x = rnorm(n))
s <- standata(y ~ x + sar(W), data = d, data2 = list(W = W))
cat("standata eigenMsar class:", class(s$eigenMsar), "\n")
# brm() prints Stan's error and returns a fit with no draws
fit <- brm(y ~ x + sar(W), data = d, data2 = list(W = W), chains = 1,
           iter = 200, refresh = 0, seed = 1)
cat("brm(): chains with draws in the fit:", length(fit$fit@sim$samples),
    "\n")
