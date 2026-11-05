# The exact gp() gap on brms's own nugget, in extended precision.
#
#   Rscript dev/gpgap-mpfr.R "<spec key>"       (gaussian or lognormal specs)
#
# frmtmb's estimate is fixed, in brms's order of positions. With 256 bits
# (Rmpfr) the factors of frmtmb's and brms's kernels are read exactly, so
# double precision cannot blur these:
#   A1  ll(L_brms z) - ll(f), z = L_frm^-1 f: the joint gap in the chart
#       brms samples, against the same reading through Stan in double;
#   C   log N(f | K_brms) - log N(f | K_frm), the same gap at the same
#       FIELD, against R's and Stan's double readings of it (B).
# A spec's prior is dropped: the prior of sdgp is not what these read.
args <- commandArgs(trailingOnly = TRUE)
key <- if (length(args)) args[1] else "lognormal x gp"
source("dev/gpgap-lib.R")
library(Rmpfr)
E$key <- key
local(envir = E, {
  spec <- eq_specs[[key]]()
  spec$prior <- NULL
  s <- gpgap_setup(spec)
  prec <- 256L
  n <- nrow(s$kb)
  xb <- mpfr(as.numeric(s$xb), prec)
  sd2 <- mpfr(s$sdgp, prec)^2
  rho2 <- mpfr(s$rho, prec)^2
  # the kernel's columns, computed once from double inputs
  ker <- lapply(seq_len(n), function(j) sd2 * exp(-(xb - xb[j])^2 / (2 * rho2)))
  mchol <- function(nug) {
    L <- vector("list", n)
    for (j in seq_len(n)) {
      idx <- j:n
      v <- ker[[j]][idx]
      v[1] <- v[1] + nug
      for (k in seq_len(j - 1L)) {
        v <- v - L[[k]][j] * L[[k]][idx]
      }
      lc <- mpfrArray(0, prec, dim = n)
      lc <- mpfr(rep(0, n), prec)
      d <- sqrt(v[1])
      lc[idx] <- v / d
      L[[j]] <- lc
    }
    L
  }
  msolve <- function(L, f) {
    z <- mpfr(rep(0, n), prec)
    for (i in seq_len(n)) {
      s0 <- mpfr(f[i], prec)
      if (i > 1L) {
        ri <- new("mpfr", unlist(lapply(seq_len(i - 1L), function(k) {
          L[[k]][i] * z[k]
        })))
        s0 <- s0 - sum(ri)
      }
      z[i] <- s0 / L[[i]][i]
    }
    z
  }
  mmult <- function(L, z) {
    out <- mpfr(rep(0, n), prec)
    for (k in seq_len(n)) out <- out + L[[k]] * z[k]
    out
  }
  logdet_half <- function(L) {
    sum(new("mpfr", unlist(lapply(seq_len(n), function(k) log(L[[k]][k])))))
  }
  t0 <- proc.time()[[3]]
  Lf <- mchol(mpfr(1e-6, prec) * sd2)
  Lb <- mchol(mpfr(1e-12, prec))
  cat(sprintf("factors in %.0f s\n", proc.time()[[3]] - t0))

  fm <- mpfr(s$fb, prec)
  zf <- msolve(Lf, fm)
  zb <- msolve(Lb, fm)
  cat(sprintf("solves done at %.0f s\n", proc.time()[[3]] - t0))
  fam <- spec$family %||% "gaussian()"
  y <- spec$data[[if (grepl("y1", spec$bf)) "y1" else "y"]]
  ib <- s$sdat[[sub("^zgp_", "Jgp_", s$zn)]]
  icpt <- s$pars[[grep("Intercept$", names(s$pars), value = TRUE)[1]]]
  sg <- s$pars[[grep("^sigma", names(s$pars), value = TRUE)[1]]]
  ll <- function(f) {
    mu <- icpt + f[ib]
    if (grepl("lognormal", fam)) sum(stats::dlnorm(y, mu, sg, log = TRUE)) else
      sum(stats::dnorm(y, mu, sg, log = TRUE))
  }

  # A1, exactly and through Stan, at the z the instrument hands Stan
  z1 <- as.numeric(solve(t(chol(s$kb)), s$fb))
  fprime <- as.numeric(mmult(Lb, mpfr(z1, prec)))
  A1x <- ll(fprime) - ll(s$fb)
  sf <- gpgap_sf(s$code, s$sdat)
  p1 <- s$pars
  p1[[s$zn]] <- array(z1, n)
  lz <- function(z) sum(stats::dnorm(z, log = TRUE))
  A1s <- gpgap_lp(sf, p1) - lz(z1) - ll(s$fb)
  cat(sprintf(paste0("GPMPFR %s\n",
                     "  A1 exact %.8f, through Stan in double %.8f, ",
                     "difference %.3e\n",
                     "  max |L_brms z - f| %.3e (sd of the field's residual ",
                     "scale: sigma %.4g)\n"),
              key, A1x, A1s, A1x - A1s, max(abs(fprime - s$fb)), sg))

  # C and B at the same field
  Cx <- as.numeric(-0.5 * sum(zb * zb) - logdet_half(Lb) +
                     0.5 * sum(zf * zf) + logdet_half(Lf))
  kb2 <- s$kb
  diag(kb2) <- diag(kb2) - 1e-6 * s$sdgp^2 + 1e-12
  Rb <- t(chol(kb2))
  Rf <- t(chol(s$kb))
  zbd <- as.numeric(solve(Rb, s$fb))
  zfd <- as.numeric(solve(Rf, s$fb))
  Cd <- -0.5 * sum(zbd^2) - sum(log(diag(Rb))) + 0.5 * sum(zfd^2) +
    sum(log(diag(Rf)))
  p2 <- s$pars
  p2[[s$zn]] <- array(zbd, n)
  # Stan at z_b in double: the likelihood at Stan's rebuilt field
  llb <- gpgap_lp(sf, p2) - lz(zbd)
  zbm <- as.numeric(zb)
  p3 <- s$pars
  p3[[s$zn]] <- array(zbm, n)
  llb3 <- gpgap_lp(sf, p3) - lz(zbm)
  cat(sprintf(paste0("  C exact %.6f, in R's double %.6f\n",
                     "  sum z_b^2 exact %.6f, double %.6f; sum z_f^2 exact ",
                     "%.6f, double %.6f\n",
                     "  ll at Stan's L z: z_b from R double %.6f; z_b exact ",
                     "%.6f; ll(f) %.6f\n"),
              Cx, Cd, as.numeric(sum(zb * zb)), sum(zbd^2),
              as.numeric(sum(zf * zf)), sum(zfd^2), llb, llb3, ll(s$fb)))
  mins <- vapply(seq_len(n), function(k) as.numeric(Lb[[k]][k]), 0)
  cat(sprintf("  smallest pivot of L_brms %.3e, of L_frm %.3e\n", min(mins),
              min(vapply(seq_len(n), function(k) as.numeric(Lf[[k]][k]), 0))))
  # the double factor against the exact one, rows of the small pivots
  cat(sprintf("  R double factor of K_brms against the exact one: max %.3e\n",
              max(abs(Rb - vapply(seq_len(n), function(k) as.numeric(Lb[[k]]),
                                  numeric(n))))))
  cat(sprintf("total %.0f s\n", proc.time()[[3]] - t0))
})
