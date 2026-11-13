# Lane wcheck: is the Laplace objective of a censored model with an
# observation-level effect unbounded when the uncensored excess is below
# the 0.81.0 allowance of 1/2 per censored row? The construction: hold
# the mean's slope at 0 and its intercept M a distance above every
# right-censoring bound, so that each right-censored row's probability
# is 1 at its own effect's prior mean, and send log sigma to -Inf with
# the effects' log sd re-optimized at each step (a lower bound on the
# profile). dry_run = "objective" skips the weights check.
#   Rscript dev/wcheck-cens-path.R <lib|base>
# Data: dev/rel073-wcens-probe.R (seed 20261021, n = 80).
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
q <- quantile(d$y, 0.8)
d$cc <- ifelse(d$y > q, "right", "none")
d$yc <- pmin(d$y, q)
nu <- sum(d$cc == "none")
ql <- quantile(d$y, 0.1)
d$c2 <- ifelse(d$y > q, "right", ifelse(d$y < ql, "left", "none"))
d$y2c <- pmin(pmax(d$y, ql), q)
nu2 <- sum(d$c2 == "none")
d$ya <- 1 + 0.5 * d$x + rnorm(n)
o <- order(d$ya, decreasing = TRUE)
d$ca <- rep("none", n); d$ca[o[1:72]] <- "right"
qa <- d$ya[o[72]]
d$yac <- pmin(d$ya, qa)

lsig <- c(-2, -4, -6, -8, -10, -14, -18, -22, -26, -30)
path <- function(label, f, family, M, fix_quantile = FALSE) {
  ob <- suppressWarnings(suppressMessages(
    frm(f, family = family, data = d, dry_run = "objective")))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  nm <- names(p)
  ib <- which(nm == "beta")
  isd <- which(nm == "betad")
  # the sigma coefficient comes first (asym_laplace: sigma, quantile)
  isg <- isd[1L]
  # asym_laplace's quantile held at logit 0 (the median) when asked
  iq <- if (fix_quantile) isd[2L] else integer(0)
  ith <- which(nm == "theta")
  free <- setdiff(seq_along(p), c(ib, isg, iq))
  p[ib] <- c(M, 0)
  if (length(iq)) p[iq] <- 0
  out <- vapply(lsig, function(ls) {
    p[isg] <- ls
    fr <- function(r) {
      q2 <- p; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- nlminb(p[free], fr)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
  cat(sprintf("%s\n  par: %s\n  -fn at log sigma %s:\n  %s\n  slope per unit of -log sigma: %s\n  log sd of the effects at the last step: %s\n",
              label, paste(nm, collapse = " "), paste(lsig, collapse = ","),
              paste(format(out, digits = 7), collapse = " "),
              paste(format(-diff(out) / diff(lsig), digits = 3),
                    collapse = " "),
              paste(format(p[ith], digits = 5), collapse = " ")))
}
G <- gaussian()
for (E in c(0.5, 2, 4, 7.5)) {
  d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
  path(sprintf("gaussian, 16 right censored, E_unc %g, M = q + 2", E),
       bf(yc | weights(w) + cens(cc) ~ x + (1 | obs)), G, q + 2)
}
# the control: E_unc 0 has no growth along the same path
d$w <- 1
path("gaussian, 16 right censored, E_unc 0 (control), M = q + 2",
     bf(yc | weights(w) + cens(cc) ~ x + (1 | obs)), G, q + 2)
for (E in c(3.5, 4.5, 8)) {
  d$w2 <- ifelse(d$c2 == "none", 1 + E / nu2, 1)
  path(sprintf("gaussian, 16 right + 8 left, E_unc %g, M = q + 2", E),
       bf(y2c | weights(w2) + cens(c2) ~ x + (1 | obs)), G, q + 2)
}
fa <- suppressWarnings(frm(yac ~ x, family = asym_laplace(), data = d))
cat("asym_laplace betad names:", names(fa$estimates$betad), "\n")
path("asym_laplace, 72 of 80 right censored, M = qa + 2",
     bf(yac | cens(ca) ~ x + (1 | obs)), asym_laplace(), qa + 2)
path("asym_laplace, 72 of 80 right censored, M = qa + 2, quantile 0.5",
     bf(yac | cens(ca) ~ x + (1 | obs)), asym_laplace(), qa + 2,
     fix_quantile = TRUE)
# the control at the 0.81.0 allowance's own construction: the fitted
# mean, where the censored rows' probability is not 1
path("asym_laplace, 72 of 80 right censored, M = qa - 2, quantile 0.5",
     bf(yac | cens(ca) ~ x + (1 | obs)), asym_laplace(), qa - 2,
     fix_quantile = TRUE)
cat("DONE\n")
