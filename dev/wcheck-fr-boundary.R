# Lane wcheck, fix round: both sides of every prior-tail boundary of the
# corrected rule (MAP objective with the log-Jacobian). For each arm, the
# check verdict and the profile of the objective along the dispersion,
# the rest re-optimized. Preamble from dev/reviews/wcheck-review-rates.R
# (seed 20261111, n = 80).
#   Rscript dev/wcheck-fr-boundary.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
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


# predicted: the objective's slope at the limit, per unit of log k (or
# of -log sigma): growth + 1 (the Jacobian) - tail power
GA <- Gamma("log")
sh <- function(p) set_prior(p, class = "shape")
E <- 12
d$w <- 1 + E / n
fg <- bf(yg | weights(w) ~ x + (1 | obs))
gg <- c(4, 6, 8, 10, 12, 14, 16)
for (pr in c("student_t(5.5, 0, 2.5)", "student_t(6.5, 0, 2.5)",
             "inv_gamma(5.5, 1)", "inv_gamma(6.5, 1)")) {
  # the tail power: nu + 1 or shape + 1, the first number of the string
  p <- as.numeric(regmatches(pr, regexpr("[0-9.]+", pr))) + 1
  show(sprintf("Gamma, E %g (growth 6), %s on shape; check: %s", E, pr,
               verdict(fg, GA, sh(pr))),
       gg, profile_disp(fg, GA, sh(pr), gg), E / 2 + 1 - p)
}
E2 <- 6
d$w2 <- 1 + E2 / n
fw <- bf(yw | weights(w2) ~ x + (1 | obs))
gw <- c(3, 5, 7, 9, 11, 13)
for (nu in c(5.5, 6.5)) {
  pr <- sh(sprintf("student_t(%g, 0, 2.5)", nu))
  show(sprintf("weibull, E %g (growth 6), student_t(%g) on shape; check: %s",
               E2, nu, verdict(fw, weibull(), pr)),
       gw, profile_disp(fw, weibull(), pr, gw), E2 + 1 - (nu + 1))
}
gs <- -c(2, 4, 6, 8, 10, 12)
for (E3 in c(0.5, 1.5)) {
  d$w3 <- 1 + E3 / n
  fn <- bf(yn | weights(w3) ~ x + (1 | obs))
  pr <- set_prior("normal(0, 1)", class = "sigma")
  show(sprintf("gaussian, E %g, normal(0, 1) on sigma (order 0); check: %s",
               E3, verdict(fn, gaussian(), pr)),
       -gs, profile_disp(fn, gaussian(), pr, gs), E3 - 1)
}
E3 <- 2
d$w3 <- 1 + E3 / n
fn <- bf(yn | weights(w3) ~ x + (1 | obs))
for (aa in c(1.5, 2.5)) {
  pr <- set_prior(sprintf("gamma(%g, 1)", aa), class = "sigma")
  show(sprintf("gaussian, E %g, gamma(%g, 1) on sigma (order %g); check: %s",
               E3, aa, aa - 1, verdict(fn, gaussian(), pr)),
       -gs, profile_disp(fn, gaussian(), pr, gs), E3 - aa)
}
fg1 <- bf(yg | weights(w) ~ x + (1 | obs), shape ~ 1)
pr <- set_prior("lognormal(0, 1)", class = "Intercept", dpar = "shape")
show(sprintf("Gamma shape ~ 1, E %g, lognormal(0, 1) on Intercept; check: %s",
             E, verdict(fg1, GA, pr)),
     gg, profile_disp(fg1, GA, pr, gg), E / 2)
cat("DONE\n")
