# rel084-kt copy of dev/reviews/squant-short-tape.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Short check of lane squant fix round 1, item 1: the taped derivatives of
# quantile_elf()'s cdf_log in mu and sigma against the analytic -f/F, f/S
# (and the sigma forms), with log F and log S from numerical integration
# of an independent density, over t in [-2000, 2000]; continuity across
# t = +-690 and t = 0 in steps of 1e-6.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-tape.R > dev/reviews/squant-short-tape.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "| elf_pos0 in namespace:",
    exists("elf_pos0", ns), "| BLAS", extSoftVersion()[["BLAS"]], "\n")
src <- parse("R/families-quantile.R", keep.source = FALSE)
env <- new.env()
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                   as.character(e[[2]]) == "elf_cdf_log") eval(e, env)
cat("installed elf_cdf_log identical to worktree source:",
    identical(deparse(body(env$elf_cdf_log)), deparse(body(ns$elf_cdf_log))),
    "\n")

# independent log density (y units)
lpdf_ref <- function(y, mu, sigma, p, lam) {
  e <- (y - mu) / sigma
  z <- e / lam
  sp <- ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  (1 - p) * e - lam * sp - log(lam * sigma) - lbeta(lam * (1 - p), lam * p)
}
# log of the integral of exp(lpdf(q + s v) - lpdf(q)) over v in [0, Inf)
# (s = -1 integrates to the left), in geometric pieces
tail_int <- function(q, s, p, lam) {
  l0 <- lpdf_ref(q, 0, 1, p, lam)
  br <- c(0, lam * 2^seq(-12, 22, by = 1))
  tot <- 0
  for (k in seq_len(length(br) - 1L)) {
    v <- stats::integrate(function(v) exp(lpdf_ref(q + s * v, 0, 1, p, lam) - l0),
                          br[k], br[k + 1L], rel.tol = 1e-13,
                          subdivisions = 500L)$value
    tot <- tot + v
    if (v < 1e-17 * tot) break
  }
  l0 + log(tot)
}
log1mexp <- function(l) ifelse(l > -log(2), log(-expm1(l)), log1p(-exp(l)))
ref_pair <- function(q, p, lam) {
  mode <- lam * stats::qlogis(1 - p)
  if (q <= mode) {
    lF <- tail_int(q, -1, p, lam)
    lS <- log1mexp(lF)
  } else {
    lS <- tail_int(q, 1, p, lam)
    lF <- log1mexp(lS)
  }
  c(lF = lF, lS = lS)
}
# q with q / lam == t exactly where one exists near lam * t
q_of_t <- function(t, lam) {
  q <- lam * t
  if (q / lam == t) return(q)
  for (k in 1:8) for (s in c(-1, 1)) {
    qq <- q + s * k * .Machine$double.eps * abs(q)
    if (qq / lam == t) return(qq)
  }
  q
}

tgrid <- sort(unique(c(
  seq(-2000, 2000, by = 50), -1000, -800, -700, -500, -300, -200, -100,
  -40, -10, 10, 40, 100, 200, 300, 500, 700, 800, 1000,
  -690 + c(-10, -1, -0.1, -1e-3, -1e-6, 0, 1e-6, 1e-3, 0.1, 1, 10),
  690 + c(-10, -1, -0.1, -1e-3, -1e-6, 0, 1e-6, 1e-3, 0.1, 1, 10),
  c(-1, -0.1, -1e-3, -1e-6, 0, 1e-6, 1e-3, 0.1, 1))))
cat("t grid:", length(tgrid), "points\n\n")

cat("== A. tape d/dmu and d/dsigma at mu 0, sigma 1 against analytic, ",
    "reference log F, log S by integration ==\n", sep = "")
cat("value err: max |lX - ref| on the small side (log scale = relative in",
    "prob), max |lX / ref - 1| on the large side\n")
cat("deriv err: relative where |analytic| > 1e-12 x slope scale; else",
    "absolute / scale (the 'tiny' column)\n")
worst <- list()
for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) for (lam in c(0.05, 0.25, 1, 4)) {
  fam <- quantile_elf(smoothing = lam)
  q <- vapply(tgrid, q_of_t, 0, lam = lam)
  t <- q / lam
  ref <- t(vapply(q, ref_pair, c(lF = 0, lS = 0), p = p, lam = lam))
  lf <- lpdf_ref(q, 0, 1, p, lam)
  dp <- function(m) list(mu = m[1] + 0 * q, sigma = m[2], quantile = p)
  v <- fam$cdf_log(q, dp(c(0, 1)), list())
  scale <- max(p, 1 - p)
  for (side in c("lF", "lS")) {
    sgn <- if (side == "lF") -1 else 1
    tp <- RTMB::MakeTape(function(m) fam$cdf_log(q, dp(m), list())[[side]],
                         c(0, 1))
    J <- tp$jacobian(c(0, 1))
    an_mu <- sgn * exp(lf - ref[, side])
    an_sg <- sgn * q * exp(lf - ref[, side])
    small <- ref[, side] < log(0.5)
    verr <- ifelse(small, abs(v[[side]] - ref[, side]),
                   abs(v[[side]] / ref[, side] - 1))
    verr[ref[, side] == 0 & v[[side]] == 0] <- 0
    meas_mu <- abs(an_mu) > 1e-12 * scale
    meas_sg <- abs(an_sg) > 1e-12 * scale * pmax(1, abs(q))
    re_mu <- abs(J[, 1] / an_mu - 1)
    re_sg <- abs(J[, 2] / an_sg - 1)
    ab_mu <- abs(J[, 1] - an_mu) / scale
    ab_sg <- abs(J[, 2] - an_sg) / (scale * pmax(1, abs(q)))
    nonfin <- sum(!is.finite(J)) + sum(!is.finite(v[[side]]))
    wmu <- if (any(meas_mu)) which(meas_mu)[which.max(re_mu[meas_mu])] else NA
    wsg <- if (any(meas_sg)) which(meas_sg)[which.max(re_sg[meas_sg])] else NA
    cat(sprintf(paste0("p %.2f lam %.2f %s: value err %.2g (t %g) | mu: ",
                       "rel %.2g (t %g, n %d) tiny %.2g | sigma: rel %.2g ",
                       "(t %g, n %d) tiny %.2g | n rel > 1e-8: %d | ",
                       "nonfinite %d\n"),
                p, lam, side, max(verr), t[which.max(verr)],
                max(re_mu[meas_mu]), t[wmu], sum(meas_mu),
                max(c(0, ab_mu[!meas_mu])),
                max(re_sg[meas_sg]), t[wsg], sum(meas_sg),
                max(c(0, ab_sg[!meas_sg])),
                sum(re_mu[meas_mu] > 1e-8) + sum(re_sg[meas_sg] > 1e-8),
                nonfin))
    worst[[length(worst) + 1L]] <- c(max(verr), max(re_mu[meas_mu]),
                                     max(re_sg[meas_sg]))
  }
}
w <- do.call(rbind, worst)
cat(sprintf("SUMMARY A: worst value err %.3g, worst rel d/dmu %.3g, worst rel d/dsigma %.3g over %d (p, lam, side)\n\n",
            max(w[, 1]), max(w[, 2]), max(w[, 3]), nrow(w)))

cat("== B. continuity across t0 = -690, 0, 690 in steps of 1e-6 ==\n")
cat("per side: dt = tape d/dt (= -lam d/dmu); jump = max |dt_{k+1} - dt_k|",
    "/ max|dt| over the 11 points; smooth = same away from t0 (k != 0);",
    "vslope = range of (l_{k+1} - l_k) / (h dt_k) (1 = consistent)\n")
jumps <- numeric()
for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) for (lam in c(0.05, 0.25, 1, 4)) {
  fam <- quantile_elf(smoothing = lam)
  for (t0 in c(-690, 0, 690)) {
    tt <- t0 + (-5:5) * 1e-6
    q <- vapply(tt, q_of_t, 0, lam = lam)
    tt <- q / lam
    dp <- function(m) list(mu = m[1] + 0 * q, sigma = m[2], quantile = p)
    v <- fam$cdf_log(q, dp(c(0, 1)), list())
    out <- character()
    for (side in c("lF", "lS")) {
      tp <- RTMB::MakeTape(function(m) fam$cdf_log(q, dp(m), list())[[side]],
                           c(0, 1))
      J <- tp$jacobian(c(0, 1))
      dt <- -lam * J[, 1]
      d1 <- abs(diff(dt)) / max(abs(dt))
      cross <- c(5, 6)
      jmp <- max(d1[cross])
      sm <- max(d1[-cross])
      h <- diff(tt)
      vs <- diff(v[[side]]) / (h * dt[-11])
      vs <- vs[is.finite(vs) & abs(dt[-11]) > 1e-200]
      # d sigma: compare against q * d mu identity (dF/dsigma = q dF/dmu
      # at mu 0, sigma 1)
      sid <- max(abs(J[, 2] - q * J[, 1])) /
        max(abs(J[, 1]) * pmax(1, abs(q)), 1e-300)
      jumps <- c(jumps, jmp)
      out <- c(out, sprintf(
        "%s dt[t0] %.6g jump %.2g smooth %.2g vslope %s sigma-id %.2g",
        side, dt[6], jmp, sm,
        if (length(vs)) sprintf("%.6f..%.6f", min(vs), max(vs)) else "n/a",
        sid))
    }
    cat(sprintf("p %.2f lam %.2f t0 %5g: %s\n", p, lam, t0,
                paste(out, collapse = " | ")))
  }
}
cat(sprintf("SUMMARY B: max jump at t0 relative to max|dt| %.3g\n", max(jumps)))

cat("\n== C. single-point tapes at exactly t = -690, 0, 690 (scalar q) ==\n")
for (p in c(0.01, 0.5, 0.99)) for (lam in c(0.05, 0.25, 4)) {
  fam <- quantile_elf(smoothing = lam)
  for (t0 in c(-690, 0, 690)) {
    q <- q_of_t(t0, lam)
    ref <- ref_pair(q, p, lam)
    lf <- lpdf_ref(q, 0, 1, p, lam)
    o <- character()
    for (side in c("lF", "lS")) {
      tp <- RTMB::MakeTape(function(m) {
        fam$cdf_log(q, list(mu = m[1], sigma = m[2], quantile = p),
                    list())[[side]]
      }, c(0, 1))
      J <- tp$jacobian(c(0, 1))
      an <- (if (side == "lF") -1 else 1) * exp(lf - ref[[side]])
      o <- c(o, sprintf("%s mu %.10g (an %.10g) sigma %.4g", side, J[1], an,
                        J[2]))
    }
    cat(sprintf("p %.2f lam %.2f t %g (q/lam==t %s): %s\n", p, lam, t0,
                q / lam == t0, paste(o, collapse = " | ")))
  }
}
cat("DONE\n")
