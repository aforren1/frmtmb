# Fix round 2, item 4: smallpois's 26 and t_outlier's 1 are called non-collapses
# because Laplace reads high. Against my exact marginal at each fit's own estimates,
# data seeds 1:40.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
gh <- gh_rule(81L)
for (nm in c("smallpois", "t_outlier")) {
  n <- 0L; flagged <- 0L; lap_high <- 0L
  zi <- numeric(0); gl <- numeric(0); gi <- numeric(0); dth <- numeric(0)
  for (s in 1:40) {
    de <- design_list[[nm]](s)
    dl <- split(de$data, de$data$g)
    ex <- function(par) {
      sd <- de$sdof(par)
      tot <- 0
      for (d in dl) {
        fn <- function(u) de$ll(par, d, u) + dnorm(u, 0, sd, log = TRUE)
        o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
        m <- o$minimum
        cur <- -(fn(m + 1e-3) - 2 * fn(m) + fn(m - 1e-3)) / 1e-6
        w <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
        nd <- m + sqrt(2) * w * gh$x
        a <- fn(nd) + gh$x^2 + log(gh$w)
        tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * w))
      }
      tot
    }
    l <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
    i <- suppressWarnings(frm(de$formula, family = de$family, data = de$data, importance = de$draws))
    pl <- unname(l$opt$par); pi <- unname(i$opt$par)
    if (pl[length(pl)] < -4) next
    n <- n + 1L
    gain <- as.numeric(logLik(l)) - as.numeric(logLik(i))
    e_l <- -ex(pl); e_i <- -ex(pi)
    fl <- gain > 0.05 + 3 * i$importance$mcse
    if (fl) flagged <- flagged + 1L
    # Laplace's value against exact at its own estimates; importance's against exact at its own
    gl <- c(gl, as.numeric(logLik(l)) - e_l)
    gi <- c(gi, (as.numeric(logLik(i)) - e_i) / i$importance$mcse)
    dth <- c(dth, e_i - e_l)
    if (fl) {
      lap_high <- lap_high + as.integer(as.numeric(logLik(l)) - e_l > 0)
      cat(sprintf("  %s s %d gain %.3f | Laplace ll - exact at its estimate %+.3f | imp err %+.3f (%.2f MCSE) | exact at imp - exact at Laplace %+.4f | log sd %.3f %.3f\n",
                  nm, s, gain, as.numeric(logLik(l)) - e_l, as.numeric(logLik(i)) - e_i, (as.numeric(logLik(i)) - e_i) / i$importance$mcse, e_i - e_l, pl[length(pl)], pi[length(pi)]))
    }
  }
  cat(sprintf("%s: %d fits off the boundary, %d flagged, of which Laplace reads above exact on %d. All fits: Laplace - exact median %+.3f range %+.3f to %+.3f; imp err in MCSE range %+.2f to %+.2f; exact(imp) - exact(Laplace) min %+.5f\n",
              nm, n, flagged, lap_high, median(gl), min(gl), max(gl), min(gi), max(gi), min(dth)))
}
