# Review of lane wcheck: right censoring with a bounded intercept prior
# that keeps the mean below the censoring bound, so no constant mean puts
# the right-censored rows at probability 1. The lane's rule still takes
# them as free (it reads only that the design spans the constant).
#   Rscript dev/reviews/wcheck-review-cens2.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n)
q8 <- quantile(d$y, 0.8)
d$cc <- ifelse(d$y > q8, "right", "none")
d$yc <- pmin(d$y, q8)
nu <- sum(d$cc == "none")
fo <- bf(yc | weights(w) + cens(cc) ~ 1 + (1 | obs))
pr <- set_prior("uniform(-5, 0.5)", class = "Intercept", lb = -5, ub = 0.5)
for (E in c(2, 6)) {
  d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
  r <- tryCatch({f <- q(frm(fo, family = gaussian(), data = d, prior = pr))
    sprintf("fit logLik %.6g beta %s betad %s conv %d", as.numeric(logLik(f)),
            format(f$estimates$beta, digits = 4),
            format(f$estimates$betad, digits = 4), f$opt$convergence)},
    error = function(e) substr(gsub("\\s+", " ", conditionMessage(e)), 1, 120))
  cat(sprintf("E %g, bound %.3f, 16 right censored: %s\n", E, q8, r))
  ob <- q(frm(fo, family = gaussian(), data = d, prior = pr,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  nm <- names(p)
  cat("  par:", nm, " start:", format(p, digits = 4), "\n")
  ib <- which(nm == "beta"); isd <- which(nm == "betad")
  free <- setdiff(seq_along(p), isd)
  p[ib] <- 0
  cat("  fn at the start with beta 0:", ob$fn(p), "\n")
  gs <- c(-2, -4, -6, -8, -10, -12, -14)
  v <- vapply(gs, function(t) {
    p[isd] <- t
    fr <- function(r) {
      q2 <- p; q2[free] <- r
      vv <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(vv)) 1e10 else vv
    }
    op <- nlminb(p[free], fr)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
  cat("  log sigma", gs, "\n  -fn", format(v, digits = 7), "\n  slope per -log sigma",
      format(-diff(v) / diff(gs), digits = 3), "\n  beta at the end",
      format(p[ib], digits = 4), "\n")
}
cat("DONE\n")
