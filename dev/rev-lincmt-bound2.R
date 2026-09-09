# Push the scale-relative error until it breaks or until it is bounded.
#
# dev/rev-lincmt-bound.R found 3.18e-11 of the trajectory's own scale
# for three compartments on a six-decade box, which is five decades
# above the "eps of its own scale" the findings claim. This searches
# for the worst it can be, and reports the disposition coefficients
# beside it, because a large coefficient of either sign is the only way
# a positive combination can cancel.
#
# Seed 5150. Script path: dev/rev-lincmt-bound2.R.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

lags <- c(0.02, 0.2, 1, 5, 20, 60, 180)
coefmax <- function(nc, p) {
  d <- frmtmb.ode:::lincmt_disp(
    nc, p$ke, if (nc >= 2) p$k12 else 0, if (nc >= 2) p$k21 else 0,
    if (nc == 3) p$k13 else 0, if (nc == 3) p$k31 else 0)
  max(abs(vapply(d$coef, as.numeric, 0)))
}
score <- function(nc, p) {
  pl <- c(p[c("ke", "k12", "k21", "k13", "k31",
              "ka")[c(TRUE, nc >= 2, nc >= 2, nc == 3, nc == 3,
                      TRUE)]], list(V = 1))
  a <- frm_lincmt(parms = pl, times = lags, ncmt = nc, depot = TRUE,
                  init = list(depot = 1), output = "central")
  b <- vapply(lags, function(u) as.numeric(bolus_ss(nc, TRUE, p, u)),
              0)
  pk <- max(abs(b))
  if (!(pk > 0)) return(NULL)
  list(errS = max(abs(a - b)) / pk,
       pw = max(abs(a - b) / pmax(abs(b), 1e-300)), cf = coefmax(nc, p))
}

cat("\n=== A. random search, rates over seven decades ===\n")
set.seed(5150)
for (nc in 2:3) {
  best <- NULL; bs <- 0
  for (i in 1:90) {
    p <- as.list(exp(runif(6, log(1e-5), log(50))))
    names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
    z <- score(nc, p)
    if (!is.null(z) && z$errS > bs) { bs <- z$errS; best <- c(p, z) }
  }
  cat(sprintf("%d cmt: worst error / peak %10.3e   pointwise %10.3e",
              nc, bs, best$pw),
      sprintf("  max|coef| %8.3g\n", best$cf))
  cat("   at ", paste(sprintf("%s=%.4g", names(best)[1:6],
                              unlist(best[1:6])), collapse = " "),
      "\n")
}

cat("\n=== B. the near-coalescent direction, driven hard ===\n")
cat("k21 and k31 driven together with k12, k13 spread far apart:",
    "\nthis is where a three-compartment coefficient can grow.\n\n")
cat(sprintf("%10s %12s %12s %12s %12s\n", "k31-k21", "max|coef|",
            "err/peak", "pointwise", "verdict"))
for (d in c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9, 1e-12, 0)) {
  p <- list(ke = 3, k12 = 8, k21 = 0.05, k13 = 0.004,
            k31 = 0.05 * (1 + d), ka = 1.3)
  z <- score(3L, p)
  cat(sprintf("%10.1e %12.4g %12.3e %12.3e %12s\n", d * 0.05, z$cf,
              z$errS, z$pw, if (z$errS < 1e-8) "under 1e-8" else
                "OVER 1e-8"))
}

cat("\n=== C. the two-compartment decoupling, driven hard ===\n")
cat(sprintf("%10s %12s %12s %12s %12s\n", "k12 k21/g^2", "c_small",
            "err/peak", "pointwise", "verdict"))
for (r in 10^-(1:10)) {
  p <- list(ke = 30, k12 = r, k21 = r, ka = 0.3)
  g <- (p$ke + p$k12 - p$k21) / 2
  dd <- sqrt(g * g + p$k12 * p$k21)
  z <- score(2L, p)
  cat(sprintf("%10.1e %12.3e %12.3e %12.3e %12s\n",
              p$k12 * p$k21 / (g * g), (dd - g) / (2 * dd), z$errS,
              z$pw, if (z$errS < 1e-8) "under 1e-8" else "OVER 1e-8"))
}
