# Prior classes of an exact gp(): names against brms's get_prior(); and what each
# package does with class = "zgp". Also the MAP's Jacobian against the closed form:
# frm() at normal(0, 1) on sdgp must maximize the marginal plus log prior plus
# log sd (brms's Stan chart). Seed 20261102.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261102); n <- 60; x <- round(runif(n, 0, 6), 3); d <- data.frame(x = x, y = sin(2 * x) + 0.1 * rnorm(n))
bp <- brms::get_prior(brms::bf(y ~ gp(x)), data = d)
fp <- frmtmb::default_prior(frmtmb::bf(y ~ gp(x)), data = d, route = "sample")
print(as.data.frame(bp)[, c("prior", "class", "coef", "group")])
print(as.data.frame(fp)[, c("prior", "class", "coef", "group")])
r <- tryCatch(brms::make_stancode(brms::bf(y ~ gp(x)), data = d, prior = brms::set_prior("normal(0, 2)", class = "zgp")), error = function(e) conditionMessage(e))
cat("brms class zgp:", substr(r, 1, 200), "\n")
r <- tryCatch({frm(bf(y ~ gp(x)), data = d, prior = set_prior("normal(0, 2)", class = "zgp")); "accepted"}, error = function(e) conditionMessage(e))
cat("frm class zgp:", substr(r, 1, 200), "\n")
fm <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, prior = set_prior("normal(0, 1)", class = "sdgp")))
th <- fm$estimates$theta
nll <- function(p) {
  K <- exp(2 * p[2]) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + 1e-12) + exp(2 * p[4])
  R <- chol(K); r <- backsolve(R, d$y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
}
post <- function(p) nll(p) - dnorm(exp(p[2]), 0, 1, log = TRUE) - log(2) - p[2]
o <- nlminb(c(fm$estimates$beta[[1]], th, log(sigma(fm)[1])), post, control = list(rel.tol = 1e-12, iter.max = 500, eval.max = 800))
cat(sprintf("MAP frm sdgp %.6f rho %.6f sigma %.6f | closed form with the Jacobian sdgp %.6f rho %.6f sigma %.6f\n",
            exp(th[1]), exp(th[2]), sigma(fm)[1], exp(o$par[2]), exp(o$par[3]), exp(o$par[4])))
