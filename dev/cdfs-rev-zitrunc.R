# Reviewer: zero_inflated_poisson with trunc(ub = 30) stops at a NaN
# gradient. Where does the NaN come from?
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
ns <- asNamespace("frmtmb")
set.seed(5)
n <- 60
x <- rnorm(n)
d <- data.frame(x = x, yp = rpois(n, exp(1 + 0.3 * x)))
cat("max y", max(d$yp), " zeros", sum(d$yp == 0), "\n")
fit <- frm(yp | trunc(ub = 30) ~ x, family = zero_inflated_poisson(),
           data = d, dry_run = "objective")
obj <- if (!is.null(fit$obj)) fit$obj else fit
p0 <- obj$par
cat("start:", p0, "\n")
cat("fn:", obj$fn(p0), "\n")
cat("gr:", obj$gr(p0), "\n")
f <- ns$as_frmtmb_family(zero_inflated_poisson())
for (mu in c(0.5, 2.7, 10)) for (zi in c(1e-4, 0.2, 0.5)) {
  cat(sprintf("mu %g zi %g: lcdf(30) %.17g lccdf(30) %.6g lcdf(0) %.6g\n", mu, zi,
              f$lcdf(30, list(mu = mu, zi = zi), list()),
              f$lccdf(30, list(mu = mu, zi = zi), list()),
              f$lcdf(0, list(mu = mu, zi = zi), list())))
}
cat("== the tape: d/d(log mu, logit zi) of log(lcdf(30)) at mu = 2.7, zi = 0.2\n")
g <- function(p) {
  mu <- exp(p[1])
  zi <- 1 / (1 + exp(-p[2]))
  log(f$lcdf(30, list(mu = mu, zi = zi), list()))
}
tp <- RTMB::MakeTape(g, c(log(2.7), qlogis(0.2)))
print(tp$jacobian(c(log(2.7), qlogis(0.2))))
cat("== poisson_cdf_log(30, mu) on the tape\n")
h <- function(p) ns$poisson_cdf_log(30, exp(p[1]))$lF
th <- RTMB::MakeTape(h, log(2.7))
print(th$jacobian(log(2.7)))
cat("numeric:", (ns$poisson_cdf_log(30, 2.7 * exp(1e-6))$lF -
                 ns$poisson_cdf_log(30, 2.7 * exp(-1e-6))$lF) / 2e-6, "\n")
cat("== same with ub = 12 and ub = 60\n")
for (ub in c(12, 20, 25, 60)) {
  r <- tryCatch({
    ft <- frm(yp | trunc(ub = ub) ~ x, family = zero_inflated_poisson(),
              data = transform(d, ubv = ub))
    sprintf("fits logLik %.6f", as.numeric(logLik(ft)))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 60)))
  cat("ub", ub, ":", r, "\n")
  h2 <- function(p) ns$poisson_cdf_log(ub, exp(p[1]))$lF
  cat("   d lF / d log mu at mu = 2.7:",
      RTMB::MakeTape(h2, log(2.7))$jacobian(log(2.7)), "\n")
}
cat("== poisson() itself with trunc(ub = 30)\n")
r <- tryCatch({
  ft <- frm(yp | trunc(ub = 30) ~ x, family = poisson(), data = d)
  sprintf("fits logLik %.6f", as.numeric(logLik(ft)))
}, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 60)))
cat(r, "\n")
cat("== hurdle_poisson, negbinomial, zinb with trunc(ub = 30), y >= 0\n")
d$yp1 <- pmax(d$yp, 1)
for (fm in list(negbinomial(), zero_inflated_negbinomial(), hurdle_poisson(),
                hurdle_negbinomial(), geometric(), brmsfamily("com_poisson"))) {
  r <- tryCatch({
    ft <- suppressWarnings(frm(yp | trunc(ub = 30) ~ x, family = fm, data = d))
    sprintf("fits logLik %.6f", as.numeric(logLik(ft)))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 60)))
  cat(sprintf("%-28s %s\n", ns$as_frmtmb_family(fm)$family, r))
}
