# brms 2.23.0 at a small Dirichlet concentration: its rdirichlet(), and
# posterior_predict_dirichlet_multinomial() and
# posterior_predict_dirichlet() on a prep with a small phi.
#   Rscript dev/nb2dm-brms-rdirichlet.R > dev/nb2dm-out/brms-rdirichlet.txt 2>&1
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
p <- exp(c(0, 0.3, -0.2)) / sum(exp(c(0, 0.3, -0.2)))
n <- 20000
for (phi in c(0.01, 1e-3)) {
  set.seed(5)
  x <- brms::rdirichlet(n, alpha = phi * p)
  cat(sprintf("rdirichlet phi %g: NaN rows %d of %d\n", phi,
              sum(!is.finite(rowSums(x))), n))
}
ns <- 2000
prep <- structure(list(ndraws = ns, nobs = 1), class = "brmsprep")
prep$dpars <- list(mu2 = rep(0.3, ns), mu3 = rep(-0.2, ns),
                   phi = rep(1e-3, ns))
prep$data <- list(Y = matrix(c(4L, 4L, 4L), 1), trials = 12L, ncat = 3L)
prep$refcat <- 1
prep$family <- dirichlet_multinomial()
set.seed(6)
r <- tryCatch(dim(brms:::posterior_predict_dirichlet_multinomial(1, prep)),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("posterior_predict_dirichlet_multinomial at phi 1e-3:", r, "\n")
prep$family <- dirichlet()
prep$data$Y <- matrix(1 / 3, 1, 3)
set.seed(6)
r <- tryCatch({
  x <- brms:::posterior_predict_dirichlet(1, prep)
  sum(!is.finite(rowSums(x)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("posterior_predict_dirichlet at phi 1e-3: NaN rows", r, "of", ns, "\n")
