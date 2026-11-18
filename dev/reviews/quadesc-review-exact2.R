# Review of lane quadesc: SYMSMALL and WEAK at roundings where a restart of
# the rule's escape raised (quadesc-review-gainhole-lib-ref-*.txt), scored
# by the review's own exact marginal log-likelihood (own density, piecewise
# stats::integrate() split at the cliffs): the lane's fit, 0.83.0's fit
# (from the base library in a second process), the Laplace fit.
#   Rscript dev/reviews/quadesc-review-exact2.R <lib|base> <case> <k-expr>
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-quadesc-lib"
.libPaths(c(if (a[1] == "lib") lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
case <- a[2]
ks <- eval(parse(text = a[3]))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("dev/quadesc-probe.R")
i0 <- grep("^sn <- function", src)
i1 <- grep("^fit1 <- function", src) - 1L
eval(parse(text = src[i0:i1]))
ldsn <- function(y, m, s, al) {
  dl <- al / sqrt(1 + al^2)
  om <- s / sqrt(1 - 2 * dl^2 / pi)
  xi <- m - om * dl * sqrt(2 / pi)
  z <- (y - xi) / om
  log(2) - log(om) + stats::dnorm(z, log = TRUE) +
    stats::pnorm(al * z, log.p = TRUE)
}
exact <- function(par, d) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  sdg <- exp(par[names(par) == "theta"][1])
  mu <- b[1] + b[2] * d$x
  sig <- exp(bd[1])
  al <- if (length(bd) == 3L) bd[2] + bd[3] * d$x else rep(bd[2], nrow(d))
  L <- 12 * sdg
  tot <- 0
  for (g in levels(d$g)) {
    i <- which(d$g == g)
    lf <- function(u) {
      out <- stats::dnorm(u, 0, sdg, log = TRUE)
      for (j in i) out <- out + ldsn(d$y[j], mu[j] + u, sig, al[j])
      out
    }
    dl <- al[i] / sqrt(1 + al[i]^2)
    om <- sig / sqrt(1 - 2 * dl^2 / pi)
    edge <- d$y[i] - mu[i] + om * dl * sqrt(2 / pi)
    edge <- edge[edge > -L & edge < L]
    br <- sort(unique(c(-L, edge, L)))
    lv <- lf(sort(c(seq(-L, L, length.out = 200001), edge)))
    mx <- max(lv[is.finite(lv)])
    s <- 0
    for (j in seq_len(length(br) - 1L)) {
      s <- s + stats::integrate(function(u) exp(lf(u) - mx), br[j],
                                br[j + 1], rel.tol = 1e-11, abs.tol = 0,
                                subdivisions = 2000L)$value
    }
    tot <- tot + mx + log(s)
  }
  tot
}
m <- make(case)
q <- function(expr) suppressMessages(suppressWarnings(expr))
for (k in ks) {
  d <- m$d
  d$y <- d$y * (1 + 2e-15 * k)
  f <- q(frm(m$fo, family = skew_normal(), data = d, quadrature = TRUE,
             start = list(betad = m$start)))
  l <- q(frm(m$fo, family = skew_normal(), data = d,
             start = list(betad = m$start)))
  al <- f$opt$par[names(f$opt$par) == "betad"][-1]
  cat(sprintf(paste0("EX %s k %d quad ll %.6f alpha %s code %d exact %.6f",
                     " | laplace fit ll %.6f exact %.6f | quad-laplace",
                     " exact %.4f\n"),
              case, k, as.numeric(logLik(f)),
              paste(format(al, digits = 4), collapse = ","),
              f$opt$convergence, ef <- exact(f$opt$par, d),
              as.numeric(logLik(l)), el <- exact(l$opt$par, d), ef - el))
}
cat("END\n")
