# Reviewer: zero_inflated_poisson at the optimizer's bad point, with
# and without trunc(ub = 30), on the tape; and the slot pieces on the
# tape at eta_zi = 1802.
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
bad <- c(-6.231083, 6.055766, 1802.885)
o1 <- frm(yp ~ x, family = zero_inflated_poisson(), data = d,
          dry_run = "objective")$obj
o2 <- frm(yp | trunc(ub = 30) ~ x, family = zero_inflated_poisson(), data = d,
          dry_run = "objective")$obj
for (lz in c(5, 20, 36, 37, 40, 100, 745, 1802.885)) {
  p <- bad
  p[3] <- lz
  cat(sprintf("eta_zi %8g: no trunc fn %-14s trunc fn %-14s | trunc gr %s\n", lz,
              format(o1$fn(p)), format(o2$fn(p)),
              paste(format(o2$gr(p), digits = 4), collapse = " ")))
}
cat("== the cdf pair on the tape, one row y = 3, mu = 2, eta_zi = t\n")
f <- ns$as_frmtmb_family(zero_inflated_poisson())
lk <- f$links$zi
tapefun <- function(p) {
  zi <- lk$linkinv(p[1])
  attr(zi, "eta") <- p[1]
  dp <- list(mu = 2 + 0 * p[1], zi = zi)
  c(log(f$lcdf(30, dp, list())), f$lccdf(30, dp, list()),
    f$lpdf(3, dp, list()))
}
print(names(f$links))
r <- tryCatch({
  tp <- RTMB::MakeTape(tapefun, 1)
  for (t in c(5, 37, 40, 745, 1802)) {
    cat("eta", t, ": value", format(tp(t)), " jac", format(tp$jacobian(t)), "\n")
  }
}, error = function(e) cat("tape error:", conditionMessage(e), "\n"))
cat("== dpar_log_complement source\n")
print(ns$dpar_log_complement)
