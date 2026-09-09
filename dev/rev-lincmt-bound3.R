# A hill climb on the scale-relative error, started from the worst
# point the random search over seven decades found, to answer the one
# question the item's 1e-8 bar turns on: can the error relative to the
# trajectory's own maximum be pushed past 1e-8?
#
# Seed 8675. Script path: dev/rev-lincmt-bound3.R.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

lags <- c(0.02, 0.2, 1, 5, 20, 60, 180)
errS <- function(nc, p) {
  pl <- c(p[c("ke", "k12", "k21", "k13", "k31", "ka")], list(V = 1))
  a <- tryCatch(frm_lincmt(parms = pl, times = lags, ncmt = nc,
                           depot = TRUE, init = list(depot = 1),
                           output = "central"),
                error = function(e) NULL)
  if (is.null(a)) return(-1)
  b <- vapply(lags, function(u) as.numeric(bolus_ss(nc, TRUE, p, u)),
              0)
  pk <- max(abs(b))
  if (!(pk > 0) || anyNA(a)) return(-1)
  max(abs(a - b)) / pk
}

x0 <- list(ke = 0.0005058, k12 = 43.43, k21 = 0.4567,
           k31 = 2.904e-05, k13 = 0.0008687, ka = 0.01594)
x0 <- x0[c("ke", "k12", "k21", "k13", "k31", "ka")]
set.seed(8675)
cur <- x0
best <- errS(3L, cur)
cat("start:", format(best, digits = 4), "\n")
for (step in c(1, 0.5, 0.25, 0.1, 0.05)) {
  for (it in 1:220) {
    cand <- lapply(cur, function(v) v * exp(rnorm(1, 0, step)))
    names(cand) <- names(cur)
    e <- errS(3L, cand)
    if (e > best) { best <- e; cur <- cand }
  }
  cat(sprintf("after step %.2f: worst error / peak %10.3e\n", step,
              best))
}
cat("\nthe worst point found:\n")
for (nm in names(cur)) cat(sprintf("   %-4s %.6g\n", nm, cur[[nm]]))
pl <- c(cur, list(V = 1))
a <- frm_lincmt(parms = pl, times = lags, ncmt = 3, depot = TRUE,
                init = list(depot = 1), output = "central")
b <- vapply(lags, function(u) as.numeric(bolus_ss(3L, TRUE, cur, u)),
            0)
cat(sprintf("\n%8s %16s %16s %12s %12s\n", "lag", "frm_lincmt",
            "300-bit", "pointwise", "/ peak"))
for (i in seq_along(lags))
  cat(sprintf("%8g %16.8e %16.8e %12.3e %12.3e\n", lags[[i]], a[[i]],
              b[[i]], abs(a[[i]] - b[[i]]) / abs(b[[i]]),
              abs(a[[i]] - b[[i]]) / max(abs(b))))
d <- frmtmb.ode:::lincmt_disp(3L, cur$ke, cur$k12, cur$k21, cur$k13,
                              cur$k31)
cat("\neigenvalues", paste(format(vapply(d$lam, as.numeric, 0),
                                  digits = 6), collapse = "  "), "\n")
cat("coefficients", paste(format(vapply(d$coef, as.numeric, 0),
                                 digits = 6), collapse = "  "), "\n")
cat("\nverdict:", if (best < 1e-8)
  "the item's 1e-8 bar was NOT crossed" else
  "the item's 1e-8 bar WAS crossed", "\n")
