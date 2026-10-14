# Lane optstart, item R1: step rules for the finite-difference outer
# Hessian (fd_hessian()), on brms_missings' fit_imp2 and fit_imp3, with
# 0.70.0's units (8 8 8 8 128 128 1 1, which the lane's mi() fix
# removes) and with the fit's own. Per rule: the eigenvalue range and
# the largest relative error of the standard errors against a
# Richardson Jacobian of the AD gradient (numDeriv).
# Rules: fixed (a step of 1e-3 in par / u, 0.70.0's), and curv(tau)
# (that pass's diagonal sets each raw step h_i = sqrt(2 tau / H_ii), the
# step along which the objective rises by tau; then a second pass), and
# the build's own fd_hessian() in par / u.
#   Rscript dev/optstart-r1-steps.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
source("dev/brms-port/shim.R")
data("nhanes", package = "mice")
set.seed(port_seed("brms_missings.10.1"))
nhanes$se <- rexp(nrow(nhanes), 2)
forms <- list(
  fit_imp2 = bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
    set_rescor(FALSE),
  fit_imp3 = bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi(se) ~ age) +
    set_rescor(FALSE))
cdiff <- function(p, gr, h) {
  n <- length(p)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    a <- p
    b <- p
    a[i] <- a[i] + h[i]
    b[i] <- b[i] - h[i]
    H[, i] <- (as.numeric(gr(a)) - as.numeric(gr(b))) / (2 * h[i])
  }
  (H + t(H)) / 2
}
for (nm in names(forms)) {
  f <- suppressWarnings(suppressMessages(frm(forms[[nm]], data = nhanes)))
  p <- f$opt$par
  gr <- function(x) f$obj$gr(x)
  Hr <- numDeriv::jacobian(function(x) as.numeric(gr(x)), p)
  ref <- sqrt(diag(solve((Hr + t(Hr)) / 2)))
  units <- list(r8 = c(8, 8, 8, 8, 128, 128, 1, 1),
                lane = f$par_units %||% rep(1, length(p)),
                none = rep(1, length(p)))
  cat("\n", nm, "lane units:", units$lane, "\n")
  for (un in names(units)) {
    u <- units[[un]]
    H1 <- cdiff(p, gr, 1e-3 * u)
    rules <- list(fixed = H1)
    fh <- frmtmb:::fd_hessian(p / u, function(q) as.numeric(gr(q * u)) * u)
    rules[["fd_hessian"]] <- fh$H / outer(u, u)
    for (tau in c(1e-2, 1e-3, 1e-4)) {
      d <- diag(H1)
      h <- ifelse(is.finite(d) & d > 0, sqrt(2 * tau / d), 1e-3 * u)
      rules[[sprintf("curv(%g)", tau)]] <- cdiff(p, gr, h)
    }
    for (rn in names(rules)) {
      H <- rules[[rn]]
      ev <- range(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
      se <- tryCatch(sqrt(diag(solve(H))), error = function(e) NA,
                     warning = function(w) NA)
      cat(sprintf("  units %-4s %-10s eigen %10.4g to %10.4g | max SE rel err %s\n",
                  un, rn, ev[1], ev[2],
                  if (all(is.finite(se))) signif(max(abs(se / ref - 1)), 3)
                  else "NaN"))
    }
  }
}
