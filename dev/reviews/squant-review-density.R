# Review of lane squant, items 1 and 2: the ELF density against qgam's
# own elf() family, normalization, the asym_laplace limit, the mean and
# variance, the log CDF pair against numerical integration of the
# density, the quantile function, draws, and the tape's derivatives at
# the switch points (|t| = 690 and the cap).
#
# Run (worktree root):
#   Rscript dev/reviews/squant-review-density.R > dev/reviews/squant-review-density.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(qgam)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "qgam", format(packageVersion("qgam")), "\n")
ns <- asNamespace("frmtmb")
dp <- function(mu, sigma, p) list(mu = mu, sigma = sigma, quantile = p)
lpdf <- function(lam, y, mu, sigma, p) {
  quantile_elf(smoothing = lam)$lpdf(y, dp(mu, sigma, p), list())
}
lcdf <- function(lam, y, mu, sigma, p) {
  quantile_elf(smoothing = lam)$cdf_log(y, dp(mu, sigma, p), list())
}

## 1a. against qgam's elf(): -aic/2 per row, co = lambda sigma ----------
cat("\n== 1a. log density against qgam::elf()$aic, per row\n")
worst <- 0
for (tau in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
  for (co in c(0.01, 0.3, 2)) {
    for (th in c(-1.5, 0, 1)) {
      fq <- qgam::elf(theta = th, link = "identity", qu = tau, co = co)
      y <- c(-30, -3, -0.2, 0, 0.1, 0.7, 4, 25)
      mu <- 0.15
      ll_q <- vapply(y, function(yy) {
        -fq$aic(yy, mu, th, 1, NULL) / 2
      }, 0)
      sig <- exp(th)
      ll_f <- lpdf(co / sig, y, mu, sig, tau)
      rel <- max(abs(ll_q - ll_f) / pmax(1, abs(ll_q)))
      worst <- max(worst, rel)
      if (rel > 1e-10) {
        cat(sprintf("  tau %.2f co %.2f theta %.1f: max rel diff %.3g\n",
                    tau, co, th, rel))
      }
    }
  }
}
cat(sprintf("1a. worst relative difference over 45 settings x 8 y: %.3g\n",
            worst))

## 1b. normalization over a grid of tau, lambda, sigma -------------------
cat("\n== 1b. integral of the density\n")
int_dens <- function(lam, sigma, p, mu = 0.3) {
  s <- lam * sigma
  # split at mu and at several kernel widths and sigma scales either side
  br <- sort(unique(mu + c(-1, 1) %o% c(0, s, 5 * s, 30 * s, sigma,
                                         10 * sigma / min(p, 1 - p))))
  f <- function(z) exp(lpdf(lam, z, mu, sigma, p))
  pieces <- c(stats::integrate(f, -Inf, br[1], rel.tol = 1e-11,
                               subdivisions = 2000)$value,
              vapply(seq_len(length(br) - 1), function(i) {
                stats::integrate(f, br[i], br[i + 1], rel.tol = 1e-11,
                                 subdivisions = 2000)$value
              }, 0),
              stats::integrate(f, br[length(br)], Inf, rel.tol = 1e-11,
                               subdivisions = 2000)$value)
  sum(pieces)
}
nb <- 0
wb <- 0
for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
  for (lam in c(0.005, 0.05, 0.25, 1, 5, 20)) {
    for (sigma in c(0.05, 1, 20)) {
      tot <- int_dens(lam, sigma, p)
      nb <- nb + 1
      wb <- max(wb, abs(tot - 1))
      if (abs(tot - 1) > 1e-7) {
        cat(sprintf("  p %.2f lam %.3f sigma %.2f: integral %.10f\n",
                    p, lam, sigma, tot))
      }
    }
  }
}
cat(sprintf("1b. %d settings: max |integral - 1| = %.3g\n", nb, wb))

## 1c. the limit to asym_laplace() --------------------------------------
cat("\n== 1c. limit as smoothing -> 0, against asym_laplace()$lpdf\n")
for (p in c(0.01, 0.3, 0.99)) {
  y <- c(-5, -0.4, 0.09, 0.11, 0.9, 6)
  ald <- asym_laplace()$lpdf(y, dp(0.1, 0.8, p), list())
  gaps <- vapply(10^-(1:5), function(lam) {
    max(abs(lpdf(lam, y, 0.1, 0.8, p) - ald))
  }, 0)
  cat(sprintf("  p %.2f: max |ELF - ALD| at lambda 1e-1..1e-5: %s\n", p,
              paste(sprintf("%.3g", gaps), collapse = " ")))
}

## 1d. mean and variance -------------------------------------------------
cat("\n== 1d. mean and variance against numerical integration\n")
fam_post <- function(lam, mu, sigma, p) {
  f <- quantile_elf(smoothing = lam)
  d <- dp(mu, sigma, p)
  c(f$post$mean_fn(d, list()), f$post$var_fn(d, list()))
}
wm <- 0
wv <- 0
for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
  for (lam in c(0.02, 0.25, 1, 5)) {
    mu <- 0.3
    sigma <- 0.7
    s <- lam * sigma
    f <- function(z, k) {
      v <- (z - mu)^k * exp(lpdf(lam, z, mu, sigma, p))
      ifelse(is.finite(v), v, 0)
    }
    br <- mu + c(-1, 1) %o% c(0, s, 30 * s, sigma, 50 * sigma / min(p, 1 - p))
    br <- sort(unique(as.numeric(br)))
    m1 <- 0
    m2 <- 0
    for (k in 1:2) {
      v <- stats::integrate(f, -Inf, br[1], k = k, rel.tol = 1e-12,
                            subdivisions = 2000)$value +
        sum(vapply(seq_len(length(br) - 1), function(i) {
          stats::integrate(f, br[i], br[i + 1], k = k, rel.tol = 1e-12,
                           subdivisions = 2000)$value
        }, 0)) +
        stats::integrate(f, br[length(br)], Inf, k = k, rel.tol = 1e-12,
                         subdivisions = 2000)$value
      if (k == 1) m1 <- v else m2 <- v
    }
    num <- c(mu + m1, m2 - m1^2)
    got <- fam_post(lam, mu, sigma, p)
    rm <- abs(got[1] - num[1]) / sqrt(num[2])
    rv <- abs(got[2] / num[2] - 1)
    wm <- max(wm, rm)
    wv <- max(wv, rv)
    if (rm > 1e-6 || rv > 1e-6) {
      cat(sprintf("  p %.2f lam %.2f: mean %.8g vs %.8g, var %.8g vs %.8g\n",
                  p, lam, got[1], num[1], got[2], num[2]))
    }
  }
}
cat(sprintf(paste0("1d. 20 settings: max |mean - integral| / sd %.3g, max",
                   " |var / integral - 1| %.3g\n"), wm, wv))

## 2a. log CDF pair against integration of the density -------------------
# log F(q) = lpdf(q) + log int_0^Inf exp(lpdf(q - s) - lpdf(q)) ds, which
# stays finite in either tail; the same for S with q + s
cat("\n== 2a. log CDF pair against integration of the density\n")
log_tail <- function(lam, q, mu, sigma, p, lower) {
  l0 <- lpdf(lam, q, mu, sigma, p)
  sg <- if (lower) -1 else 1
  f <- function(s) exp(lpdf(lam, q + sg * s, mu, sigma, p) - l0)
  sc <- lam * sigma
  d0 <- sg * (mu - q)   # distance to mu in the direction of integration
  br <- c(0, sc, 10 * sc, 100 * sc)
  if (d0 > 0) br <- c(br, d0 + c(-30, -3, 0, 3, 30) * sc)
  br <- sort(unique(pmax(br, 0)))
  v <- sum(vapply(seq_len(length(br) - 1), function(i) {
    stats::integrate(f, br[i], br[i + 1], rel.tol = 1e-12,
                     subdivisions = 2000)$value
  }, 0)) + stats::integrate(f, br[length(br)], Inf, rel.tol = 1e-12,
                            subdivisions = 2000)$value
  l0 + log(v)
}
wc <- 0
cnt <- 0
for (p in c(0.01, 0.5, 0.99)) {
  for (lam in c(0.001, 0.02, 0.25, 5)) {
    mu <- 0.3
    sigma <- 0.7
    for (t in c(-5000, -800, -690, -689.9, -50, -2, -0.1, 0, 0.1, 2, 50,
                689.9, 690, 800, 5000)) {
      q <- mu + t * lam * sigma
      got <- lcdf(lam, q, mu, sigma, p)
      if (t <= 0) {
        lf <- log_tail(lam, q, mu, sigma, p, TRUE)
        e1 <- abs(got$lF - lf) / max(1, abs(lf))
        e2 <- abs(exp(got$lS) - (-expm1(lf))) / max(1e-300, -expm1(lf))
      } else {
        ls <- log_tail(lam, q, mu, sigma, p, FALSE)
        e1 <- abs(got$lS - ls) / max(1, abs(ls))
        e2 <- abs(exp(got$lF) - (-expm1(ls))) / max(1e-300, -expm1(ls))
      }
      cnt <- cnt + 1
      e <- max(e1, e2)
      wc <- max(wc, e)
      if (!is.finite(e) || e > 1e-8) {
        cat(sprintf("  p %.2f lam %.3f t %g: lF %.12g lS %.12g err %.3g %.3g\n",
                    p, lam, t, got$lF, got$lS, e1, e2))
      }
    }
  }
}
cat(sprintf("2a. %d points: max relative error of the small tail %.3g\n",
            cnt, wc))

## 2b. the quantile function at extreme tau and tiny smoothing ------------
cat("\n== 2b. elf_quantile(): lF(q(pr)) against log(pr)\n")
wq <- 0
nq <- 0
for (p in c(0.01, 0.5, 0.99)) {
  for (lam in c(0.001, 0.02, 0.25, 5)) {
    pr <- c(1e-300, 1e-100, 1e-12, 1e-3, 0.5, 1 - 1e-3, 1 - 1e-12)
    ws <- character()
    q <- withCallingHandlers(ns$elf_quantile(pr, 0.3, 0.7, p, lam),
                             warning = function(w) {
                               ws <<- c(ws, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             })
    cp <- lcdf(lam, q, 0.3, 0.7, p)
    # compare the smaller tail on the log scale
    err <- ifelse(pr <= 0.5, abs(cp$lF - log(pr)) / pmax(1, abs(log(pr))),
                  abs(cp$lS - log1p(-pr)) / pmax(1, abs(log1p(-pr))))
    nq <- nq + length(pr)
    wq <- max(wq, err, na.rm = TRUE)
    if (any(!is.finite(q)) || any(err > 1e-8) || length(ws)) {
      cat(sprintf("  p %.2f lam %.3f: q %s | rel err %s | warnings %d %s\n",
                  p, lam, paste(sprintf("%.6g", q), collapse = " "),
                  paste(sprintf("%.2g", err), collapse = " "), length(ws),
                  substr(paste(unique(ws), collapse = " | "), 1, 120)))
    }
  }
}
cat(sprintf("2b. %d points: max relative error of the log tail %.3g\n", nq,
            wq))

## 2c. draws: KS against the family's CDF at extreme settings -------------
cat("\n== 2c. draws (sim) against the CDF, KS on 1e5 draws\n")
set.seed(20261117)
for (p in c(0.01, 0.5, 0.99)) {
  for (lam in c(0.001, 0.02, 0.25, 5)) {
    f <- quantile_elf(smoothing = lam)
    x <- f$sim(dp(0.3, 0.7, p), list(), 1e5)
    cdf <- function(z) exp(lcdf(lam, z, 0.3, 0.7, p)$lF)
    ks <- suppressWarnings(stats::ks.test(x, cdf))
    mv <- fam_post(lam, 0.3, 0.7, p)
    cat(sprintf(paste0("  p %.2f lam %.3f: KS p-value %.3f | mean z %.2f",
                       " | finite %d of 1e5\n"), p, lam, ks$p.value,
                (mean(x) - mv[1]) / sqrt(mv[2] / 1e5), sum(is.finite(x))))
  }
}

## 2d. the tape at the switch points -------------------------------------
cat("\n== 2d. tape gradient in mu, sigma, quantile at the switch points\n")
tape_grad <- function(lam, q, mu, sigma, p, which) {
  fam <- quantile_elf(smoothing = lam)
  g <- function(th) {
    r <- fam$cdf_log(q, dp(th[1], th[2], th[3]), list())
    r[[which]]
  }
  tp <- RTMB::MakeTape(g, c(mu, sigma, p))
  list(val = tp(c(mu, sigma, p)), grad = as.numeric(tp$jacobian(c(mu, sigma,
                                                                  p))))
}
fd_grad <- function(lam, q, mu, sigma, p, which) {
  th <- c(mu, sigma, p)
  vapply(1:3, function(k) {
    h <- 1e-6 * max(abs(th[k]), 1e-3)
    up <- th
    dn <- th
    up[k] <- up[k] + h
    dn[k] <- dn[k] - h
    (lcdf(lam, q, up[1], up[2], up[3])[[which]] -
       lcdf(lam, q, dn[1], dn[2], dn[3])[[which]]) / (2 * h)
  }, 0)
}
wg <- 0
for (p in c(0.01, 0.5, 0.99)) {
  for (lam in c(0.02, 0.25, 5)) {
    mu <- 0.3
    sigma <- 0.7
    for (t in c(-700, -690, -689, -0.5, 0, 0.5, 689, 690, 700)) {
      q <- mu + t * lam * sigma
      for (which in c("lF", "lS")) {
        tg <- tape_grad(lam, q, mu, sigma, p, which)
        fg <- fd_grad(lam, q, mu, sigma, p, which)
        dv <- lcdf(lam, q, mu, sigma, p)[[which]]
        # the analytic mu-derivative: -f/F for lF, f/S for lS
        an <- exp(lpdf(lam, q, mu, sigma, p) - dv) * (if (which == "lF") -1
                                                      else 1)
        e <- max(abs(tg$grad - fg) / pmax(1e-8, abs(fg)))
        ea <- abs(tg$grad[1] - an) / max(1e-300, abs(an))
        if (!is.finite(tg$val) || abs(tg$val - dv) > 1e-12 * max(1, abs(dv)) ||
              !all(is.finite(tg$grad)) || e > 1e-4 || ea > 1e-6) {
          cat(sprintf(paste0("  p %.2f lam %.2f t %g %s: tape %.10g double",
                             " %.10g | grad %s | fd %s | analytic d/dmu",
                             " %.6g\n"), p, lam, t, which, tg$val, dv,
                      paste(sprintf("%.6g", tg$grad), collapse = " "),
                      paste(sprintf("%.6g", fg), collapse = " "), an))
        }
        wg <- max(wg, ea)
      }
    }
  }
}
cat(sprintf("2d. max relative error of the tape's d/dmu against -f/F, f/S: %.3g\n",
            wg))
cat("DONE\n")
