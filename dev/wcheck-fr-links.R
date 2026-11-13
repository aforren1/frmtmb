# Lane wcheck, fix round 2: the Jacobian power of each link, confirmed by
# a profile on both sides of each boundary. The objective's power at the
# limit is excess - order - j toward 0 (per unit of log(1 / sigma)) and
# rate * excess + j - p toward infinity (per unit of log shape), with j
# from natural_jacobian_power(): log 1/1, identity 0/0, softplus 1/0,
# squareplus 2/0. Each arm: the check's verdict, then the slope of the
# profile (the rest re-optimized) along a grid even in log x, mapped to
# the link scale, against that prediction.
#   Rscript dev/wcheck-fr-links.R <lib|base>
# Data: seed 20261111, n = 80; gaussian y and Gamma yg with (1 | obs).
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
d$y <- 0.5 + 0.4 * d$x + rnorm(n)
mu <- exp(0.5 + 0.3 * d$x + rnorm(n, 0, 0.3))
d$yg <- rgamma(n, shape = 4, rate = 4 / mu)
# the link scale of x, for a grid even in log x
to_eta <- list(log = log, identity = identity,
               softplus = function(x) x + log(-expm1(-x)),
               squareplus = function(x) (x^2 - 1) / x)
verdict <- function(fo, fam, pr, data) {
  tryCatch({
    f <- q(frm(fo, family = fam, data = data, prior = pr))
    sprintf("fit logLik %.6g conv %d", as.numeric(logLik(f)),
            f$opt$convergence)
  }, error = function(e) {
    m <- gsub("\\s+", " ", conditionMessage(e))
    if (grepl("has no maximum|is improper", m)) "REFUSED" else
      paste("ERROR", substr(m, 1, 80))
  })
}
prof <- function(fo, fam, pr, data, eta) {
  ob <- q(frm(fo, family = fam, data = data, prior = pr,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  k <- which(names(p) == "betad")[1L]
  free <- setdiff(seq_along(p), k)
  vapply(eta, function(t) {
    p[k] <- t
    fr <- function(r) {
      q2 <- p; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- tryCatch(nlminb(p[free], fr), error = function(e) NULL)
    if (is.null(op)) return(NA_real_)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
}
arm <- function(label, fo, fam, pr, data, lx, ln, pred) {
  eta <- to_eta[[ln]](exp(lx))
  v <- prof(fo, fam, pr, data, eta)
  s <- diff(v) / abs(diff(lx))
  cat(sprintf("%s | check: %s\n  slope %s | predicted %g\n", label,
              verdict(fo, fam, pr, data),
              paste(format(s, digits = 3), collapse = " "), pred))
}
# toward 0: gaussian sigma, normal(0, 1) (order 0) or gamma(2, 1)
# (order 1), excess E; slope per unit of log(1 / sigma) = E - order - j
fo <- bf(y | weights(w) ~ x + (1 | obs))
lx0 <- -c(2, 4, 6, 8, 10)
j0 <- c(log = 1, identity = 0, softplus = 1, squareplus = 2)
for (ln in names(j0)) {
  fam <- brmsfamily("gaussian", link_sigma = ln)
  # order 0 unless the link needs order 1 for a bounded side (identity)
  ord <- if (j0[[ln]] == 0) 1 else 0
  pr <- set_prior(if (ord == 0) "normal(0, 1)" else "gamma(2, 1)",
                  class = "sigma")
  for (E in ord + j0[[ln]] + c(-0.5, 0.5)) {
    dd <- d
    dd$w <- 1 + E / n
    arm(sprintf("gaussian link_sigma %s, order %g, E %g", ln, ord, E),
        fo, fam, pr, dd, lx0, ln, E - ord - j0[[ln]])
  }
}
# toward infinity: Gamma shape, rate 1/2, E 12 (growth 6); student_t(nu)
# has p = nu + 1; slope per unit of log shape = 6 + j - p
fg <- bf(yg | weights(w) ~ x + (1 | obs))
lxi <- c(4, 6, 8, 10, 12)
ji <- c(log = 1, identity = 0, softplus = 0, squareplus = 0)
dg <- d
dg$w <- 1 + 12 / n
for (ln in names(ji)) {
  fam <- brmsfamily("Gamma", link = "log", link_shape = ln)
  for (p in 6 + ji[[ln]] + c(-0.5, 0.5)) {
    pr <- set_prior(sprintf("student_t(%g, 0, 2.5)", p - 1), class = "shape")
    arm(sprintf("Gamma link_shape %s, student_t(%g), growth 6", ln, p - 1),
        fg, fam, pr, dg, lxi, ln, 6 + ji[[ln]] - p)
  }
}
cat("DONE\n")
