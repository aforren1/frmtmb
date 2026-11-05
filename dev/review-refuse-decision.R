# Item 3: the open decision, with replicates. For each family and missing
# share, R replicates of a design I chose (three shapes per family where a
# shape matters), joint Laplace against the complete-case fit (the exact ML of
# an unread response). Largest and mean |shift| in complete-case SE per
# parameter class, the logLik gap, and the error of each route against the
# truth (rmse of the dpar on its internal scale).
# Usage: Rscript dev/review-refuse-decision.R <family>. Seed 9400 + family.
fam <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/review-refuse-lib.R")
fams <- c("student", "exgaussian", "skew_normal", "gen_extreme_value", "beta",
          "Gamma", "inverse.gaussian", "von_mises")
fi <- match(fam, fams)
ff <- switch(fam,
  beta = frmtmb::Beta(),
  Gamma = stats::Gamma(link = "log"),
  inverse.gaussian = stats::inverse.gaussian(link = "log"),
  get(fam, asNamespace("frmtmb"))())
shapes <- switch(fam,
  student = c(3, 6, 15),
  exgaussian = c(0.5, 1.5, 4),   # beta / sigma = shape, sigma 0.5
  skew_normal = c(1, 3, 8),      # alpha
  gen_extreme_value = c(-0.2, 0.1, 0.3),  # xi
  beta = c(1, 4, 20),            # phi
  Gamma = c(0.7, 2, 8),          # shape
  inverse.gaussian = c(0.5, 2, 10),  # shape
  von_mises = c(1, 3, 10))       # kappa
gen <- function(n, s) {
  x <- rnorm(n)
  eta <- 0.3 + 0.5 * x
  y <- switch(fam,
    student = eta + rt(n, s),
    exgaussian = eta + rnorm(n, 0, 0.5) + rexp(n, 1 / (0.5 * s)),
    skew_normal = {
      del <- s / sqrt(1 + s^2)
      eta + del * abs(rnorm(n)) + sqrt(1 - del^2) * rnorm(n)
    },
    gen_extreme_value = eta + ((rexp(n))^(-s) - 1) / s,
    beta = rbeta(n, plogis(eta) * s, (1 - plogis(eta)) * s),
    Gamma = rgamma(n, s, s / exp(eta)),
    inverse.gaussian = {
      m <- exp(eta); v <- rnorm(n)^2; lam <- s
      yy <- m + m^2 * v / (2 * lam) -
        m / (2 * lam) * sqrt(4 * m * lam * v + m^2 * v^2)
      ifelse(runif(n) < m / (m + yy), yy, m^2 / yy)
    },
    von_mises = {
      m <- 2 * atan(eta)
      # by rejection, wrapped into (-pi, pi)
      out <- numeric(n)
      for (i in seq_len(n)) {
        repeat {
          v <- runif(1, -pi, pi)
          if (runif(1) < exp(s * (cos(v) - 1))) break
        }
        out[i] <- ((v + m[i] + pi) %% (2 * pi)) - pi
      }
      out
    })
  data.frame(x = x, y = y)
}
R <- 20L
set.seed(9400 + fi)
for (s in shapes) for (cfg in list(c(200, 15), c(200, 60))) {
  n <- cfg[1]; m <- cfg[2]
  gaps <- numeric(R); loc <- numeric(R); dp <- numeric(R)
  dpj <- numeric(R); dpc <- numeric(R); fails <- 0L; nbad <- 0L
  for (r in seq_len(R)) {
    d <- gen(n, s)
    d$y[sample(n, m)] <- NA
    a <- tr(frm(bf(y | mi() ~ x), family = ff, data = d))
    b <- tr(frm(bf(y ~ x), family = ff, data = d[!is.na(d$y), ]))
    if (is_err(a) || is_err(b)) { fails <- fails + 1L; next }
    if (a$opt$convergence != 0L) nbad <- nbad + 1L
    gaps[r] <- (as.numeric(logLik(a)) - as.numeric(logLik(b))) / m
    ea <- a$estimates; eb <- b$estimates
    H <- tr(stats::optimHess(b$opt$par, b$obj$fn, b$obj$gr))
    sb <- tr(sqrt(diag(solve(H))))
    if (is_err(sb)) { fails <- fails + 1L; next }
    nb <- length(eb$beta)
    loc[r] <- max(abs((ea$beta - eb$beta) / sb[seq_len(nb)]))
    if (length(eb$betad)) {
      dp[r] <- max(abs((ea$betad - eb$betad) / sb[nb + seq_along(eb$betad)]))
      dpj[r] <- ea$betad[[length(ea$betad)]]
      dpc[r] <- eb$betad[[length(eb$betad)]]
    }
  }
  pr("%-18s shape=%-5g n=%d m=%2d R=%d fail=%d nonconv=%d gap/miss mean %+.4f [%+.4f,%+.4f]; location |shift| SE mean %.3f max %.3f; dpar |shift| SE mean %.3f max %.3f; last dpar med joint %.3f cc %.3f",
     fam, s, n, m, R, fails, nbad, mean(gaps), min(gaps), max(gaps), mean(loc),
     max(loc), mean(dp), max(dp), median(dpj), median(dpc))
}
