# Review of lane wcheck: (1) the width rates of unbounded_spike_rate,
# measured as the slope of the Laplace objective (no prior) along the log
# dispersion with every other parameter re-optimized; (2) each prior rule
# on both sides of its boundary: the check's verdict against the slope of
# the objective with the prior.
#   Rscript dev/reviews/wcheck-review-rates.R <lib|base> <part: rate|prior>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
part <- a[2]
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "| part", part, "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
mu <- exp(0.5 + 0.3 * d$x + rnorm(n, 0, 0.3))
d$yg <- rgamma(n, shape = 4, rate = 4 / mu)
d$yw <- rweibull(n, shape = 3, scale = mu)
d$yn <- log(mu) + rnorm(n, 0, 0.5)
d$yb <- rbeta(n, plogis(log(mu) - 0.5) * 8, (1 - plogis(log(mu) - 0.5)) * 8)
d$yi <- 1 / rgamma(n, 3, 3 * mu)   # positive, roughly inverse.gaussian-like
d$yv <- atan2(sin(0.3 * d$x + rnorm(n, 0, 0.4)), cos(0.3 * d$x +
                                                    rnorm(n, 0, 0.4)))
d$yf <- mu * (-log(runif(n)))^(-1 / 4)  # frechet, nu = 4
d$yl <- exp(d$yn)

verdict <- function(fo, family, prior = NULL, data = d) {
  tryCatch({q(frm(fo, family = family, data = data, prior = prior)); "fit"},
           error = function(e) {
             m <- gsub("\\s+", " ", conditionMessage(e))
             if (grepl("has no maximum|is improper", m)) "REFUSED" else
               paste("ERROR", substr(m, 1, 100))
           })
}

# the objective along the dispersion coefficient j over grid g, the rest
# re-optimized; returns -fn at each point
profile_disp <- function(fo, family, prior = NULL, g, j = 1L, data = d) {
  ob <- q(frm(fo, family = family, data = data, prior = prior,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  isd <- which(names(p) == "betad")[j]
  free <- setdiff(seq_along(p), isd)
  vapply(g, function(t) {
    p[isd] <- t
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

show <- function(label, g, v, expect) {
  s <- diff(v) / diff(g)
  cat(sprintf("%s\n  grid %s\n  -fn  %s\n  slope %s\n  predicted %s\n",
              label, paste(g, collapse = " "),
              paste(format(v, digits = 7), collapse = " "),
              paste(format(s, digits = 3), collapse = " "),
              format(expect, digits = 4)))
}

if (part == "rate") {
  E <- 8
  d$w <- 1 + E / n
  cases <- list(
    list("gaussian sigma", bf(yn | weights(w) ~ x + (1 | obs)), gaussian(),
         -c(2, 4, 6, 8, 10, 12), -1, 1),
    list("lognormal sigma", bf(yl | weights(w) ~ x + (1 | obs)), lognormal(),
         -c(2, 4, 6, 8, 10, 12), -1, 1),
    list("student sigma", bf(yn | weights(w) ~ x + (1 | obs)), student(),
         -c(2, 4, 6, 8, 10, 12), -1, 1),
    list("Gamma shape", bf(yg | weights(w) ~ x + (1 | obs)), Gamma("log"),
         c(4, 6, 8, 10, 12, 14), 1, 0.5),
    list("weibull shape", bf(yw | weights(w) ~ x + (1 | obs)), weibull(),
         c(3, 5, 7, 9, 11, 13), 1, 1),
    list("Beta phi", bf(yb | weights(w) ~ x + (1 | obs)), Beta(),
         c(4, 6, 8, 10, 12, 14), 1, 0.5),
    list("inverse.gaussian shape", bf(yi | weights(w) ~ x + (1 | obs)),
         inverse.gaussian("log"), c(4, 6, 8, 10, 12, 14), 1, 0.5),
    list("von_mises kappa", bf(yv | weights(w) ~ x + (1 | obs)),
         von_mises(), c(4, 6, 8, 10, 12, 14), 1, 0.5),
    list("frechet nu", bf(yf | weights(w) ~ x + (1 | obs)), frechet(),
         c(3, 5, 7, 9, 11, 13), 1, 1))
  for (cs in cases) {
    v <- tryCatch(profile_disp(cs[[2]], cs[[3]], g = cs[[4]]),
                  error = function(e) {
                    cat("  ERROR", conditionMessage(e), "\n"); NULL })
    if (is.null(v)) next
    # slope per unit of log(1 / width): sign by the direction of the limit
    show(sprintf("%s (E = %g, rate in the table %g)", cs[[1]], E, cs[[6]]),
         cs[[4]] * cs[[5]], v, cs[[6]] * E)
  }
}

if (part == "prior") {
  # Gamma: growth E/2 per unit log shape. E = 12 -> 6.
  E <- 12
  d$w <- 1 + E / n
  fg <- bf(yg | weights(w) ~ x + (1 | obs))
  GA <- Gamma("log")
  sh <- function(p) set_prior(p, class = "shape")
  gg <- c(4, 6, 8, 10, 12, 14, 16)
  for (nu in c(4.5, 5.5)) {
    pr <- sh(sprintf("student_t(%g, 0, 2.5)", nu))
    show(sprintf("Gamma, E %g, student_t(%g) on shape: p = %g vs E/2 = %g; check: %s",
                 E, nu, nu + 1, E / 2, verdict(fg, GA, pr)),
         gg, profile_disp(fg, GA, pr, gg), E / 2 - (nu + 1))
  }
  for (al in c(4.5, 5.5)) {
    pr <- sh(sprintf("inv_gamma(%g, 1)", al))
    show(sprintf("Gamma, E %g, inv_gamma(%g, 1) on shape: p = %g; check: %s",
                 E, al, al + 1, verdict(fg, GA, pr)),
         gg, profile_disp(fg, GA, pr, gg), E / 2 - (al + 1))
  }
  # weibull: growth E per unit log shape
  E2 <- 6
  d$w2 <- 1 + E2 / n
  fw <- bf(yw | weights(w2) ~ x + (1 | obs))
  gw <- c(3, 5, 7, 9, 11, 13)
  for (nu in c(4.5, 5.5)) {
    pr <- sh(sprintf("student_t(%g, 0, 2.5)", nu))
    show(sprintf("weibull, E %g, student_t(%g) on shape: p = %g vs E = %g; check: %s",
                 E2, nu, nu + 1, E2, verdict(fw, weibull(), pr)),
         gw, profile_disp(fw, weibull(), pr, gw), E2 - (nu + 1))
  }
  # gaussian sigma toward 0: gamma(a, 1) is order a - 1; E = 2
  E3 <- 2
  d$w3 <- 1 + E3 / n
  fn <- bf(yn | weights(w3) ~ x + (1 | obs))
  gs <- -c(2, 4, 6, 8, 10, 12)
  for (aa in c(2.5, 3.5)) {
    pr <- set_prior(sprintf("gamma(%g, 1)", aa), class = "sigma")
    show(sprintf("gaussian, E %g, gamma(%g, 1) on sigma: order %g; check: %s",
                 E3, aa, aa - 1, verdict(fn, gaussian(), pr)),
         -gs, profile_disp(fn, gaussian(), pr, gs), E3 - (aa - 1))
  }
  # log link, Gamma shape ~ 1: exponential(r) on the Intercept; E/2 = 6
  fg1 <- bf(yg | weights(w) ~ x + (1 | obs), shape ~ 1)
  pi <- function(p) set_prior(p, class = "Intercept", dpar = "shape")
  for (r in c(5.5, 6.5)) {
    pr <- pi(sprintf("exponential(%g)", r))
    show(sprintf("Gamma shape ~ 1, E %g, exponential(%g) on Intercept; check: %s",
                 E, r, verdict(fg1, GA, pr)),
         gg, profile_disp(fg1, GA, pr, gg), E / 2 - r)
  }
  for (r in c(5.5, 6.5)) {
    pr <- pi(sprintf("gamma(2, %g)", r))
    show(sprintf("Gamma shape ~ 1, E %g, gamma(2, %g) on Intercept; check: %s",
                 E, r, verdict(fg1, GA, pr)),
         gg, profile_disp(fg1, GA, pr, gg), E / 2 - r)
  }
  # lognormal on a log-scale coefficient toward +Inf: tail exp(-(log x)^2)
  pr <- pi("lognormal(0, 1)")
  show(sprintf("Gamma shape ~ 1, E %g, lognormal(0, 1) on Intercept; check: %s",
               E, verdict(fg1, GA, pr)),
       gg, profile_disp(fg1, GA, pr, gg), E / 2)
  # logistic on a log-scale sigma intercept: rate 1 / s vs E = 2
  fn1 <- bf(yn | weights(w3) ~ x + (1 | obs), sigma ~ 1)
  for (s in c(0.6, 0.4)) {
    pr <- set_prior(sprintf("logistic(0, %g)", s), class = "Intercept",
                    dpar = "sigma")
    show(sprintf("gaussian sigma ~ 1, E %g, logistic(0, %g) on Intercept: 1/s = %.3g; check: %s",
                 E3, s, 1 / s, verdict(fn1, gaussian(), pr)),
         -gs, profile_disp(fn1, gaussian(), pr, gs), E3 - 1 / s)
  }
  # exponential on a log-scale sigma intercept: the check says "bounds"
  # because the support is (0, Inf); is the support enforced?
  pr <- set_prior("exponential(1)", class = "Intercept", dpar = "sigma")
  show(sprintf("gaussian sigma ~ 1, E %g, exponential(1) on Intercept; check: %s",
               E3, verdict(fn1, gaussian(), pr)),
       -gs, profile_disp(fn1, gaussian(), pr, gs), NA)
}
cat("DONE\n")
