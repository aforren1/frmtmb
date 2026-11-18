# Reviewer of lane drawdir, check 1: what "beyond the parameter's own
# slice" measures. A 2-D objective, a quadratic with correlation rho
# plus a cubic in parameter a ALONE (r(u) = s * max(u_a, 0)^3 / 6), so
# every excess is one parameter's own skew and none is joint. Under
# that objective the normal law's marginal of a is off by exactly
# r(k sd_a) at k marginal sds (the others integrate out of the
# quadratic), so eA is the marginal's error and eB = r(k csd_a) the
# same function nearer in; draw_check_dirs() reads eA - max(eB, 0).
# Usage: Rscript dev/reviews/drawdir-review-skew.R
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), " tol", ns$draw_dir_tol, "\n")
for (rho in c(0, 0.6, 0.8, 0.9)) {
  for (s in c(2, 5, 8, 12)) {
    V <- matrix(c(1, rho, rho, 1), 2)
    Vi <- solve(V)
    fn <- function(p) sum(p * (Vi %*% p)) / 2 + s * max(p[1], 0)^3 / 6
    dw <- list(L = chol(V), law = NULL, v0 = c(0, 0), held = integer(0),
               off = 0L, labels = c("a", "b"))
    object <- list(obj = list(fn = fn, env = new.env()),
                   opt = list(par = c(0, 0)))
    dd <- ns$draw_check_dirs(object, dw)
    csd <- sqrt(1 - rho^2)
    eB <- s * (2 * csd)^3 / 6
    cat(sprintf(paste0("rho %.1f sratio %.2f s %.1f | slice excess eB %6.3f",
                       " (side law steep at 1 csd: %s) | stat %7.3f (%s)",
                       " fire %s | eB*(sratio^3-1) %7.3f\n"),
                rho, 1 / csd, s, eB, s * csd^3 / 6 >= 1.5, dd$stat, dd$par,
                dd$stat > ns$draw_dir_tol, eB * (1 / csd^3 - 1)))
  }
}
