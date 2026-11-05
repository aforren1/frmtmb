# B1's neighbors: what else moves a predictor's scale. Each row is an affine
# map of something the model does not care about, against the same fit on the
# standard form: SE ratio of the slope and of theta, log-likelihood difference
# in MCSE, rounds, warnings. Data seeds 1:3. Log dev/impsamp-out/scales.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ses <- function(f) sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed))
fit1 <- function(form, fam, d, N = 300L) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(suppressMessages(frm(form, family = fam, data = d, importance = N)),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(fit = f, warn = w)
}
row <- function(lab, a, b, se_mult = 1, jb = 2L, jt = NULL) {
  if (inherits(b$fit, "error")) {
    cat(sprintf("%-34s STOP: %s\n", lab,
                substr(conditionMessage(b$fit), 1, 120)))
    return(invisible())
  }
  sa <- ses(a$fit)
  sb <- ses(b$fit)
  th <- which(names(a$fit$opt$par) == "theta")[[1L]]
  cat(sprintf("%-34s slope SE %.4f theta SE %.4f dll/mcse %+.3f rounds %d capped %s warns %d %s\n",
              lab, sb[[jb]] * se_mult / sa[[jb]], sb[[th]] / sa[[th]],
              (as.numeric(logLik(b$fit)) - as.numeric(logLik(a$fit))) /
                a$fit$importance$mcse, b$fit$importance$rounds,
              b$fit$importance$capped, length(b$warn),
              if (length(b$warn)) substr(b$warn[[1L]], 1, 40) else ""))
}
for (s in 1:3) {
  cat("--- data seed", s, "\n")
  set.seed(s)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10L))
  z <- rnorm(150L)
  u <- rnorm(G, 0, 0.8)
  # poisson with an offset: the predictor's level moves, the columns do not
  y <- rpois(150L, exp(0.5 + 0.3 * z + u[g]))
  d <- data.frame(y = y, z = z, g = g, o = 0)
  a <- fit1(bf(y ~ z + (1 | g)), poisson(), d)
  d2 <- d
  d2$y <- y
  d2$o <- log(1000)
  a2 <- fit1(bf(y ~ z + offset(o) + (1 | g)), poisson(), d2)
  row("poisson offset log(1000)", a, a2)
  # a column with a huge shift and a small sd in a factor-coded neighbour
  d3 <- d
  d3$w <- z + 1e5
  row("poisson z + 1e5", a, fit1(bf(y ~ w + (1 | g)), poisson(), d3))
  d3$w <- z * 1e5
  row("poisson 1e5 z", a, fit1(bf(y ~ w + (1 | g)), poisson(), d3), 1e5)
  d3$w <- z * 1e-4
  row("poisson 1e-4 z", a, fit1(bf(y ~ w + (1 | g)), poisson(), d3), 1e-4)
  d3$w <- 100 * z + 1e4
  row("poisson 100 z + 1e4", a, fit1(bf(y ~ w + (1 | g)), poisson(), d3), 100)
  # a response's units, on an identity link
  yt <- 0.3 + 0.5 * z + u[g] + 0.5 * rt(150L, 5L)
  dt <- data.frame(y = yt, z = z, g = g)
  a4 <- fit1(bf(y ~ z + (1 | g), nu = 5), student(), dt)
  for (k in c(1e6, 1e-6)) {
    dk <- dt
    dk$y <- dt$y * k
    row(sprintf("student y x %g", k), a4,
        fit1(bf(y ~ z + (1 | g), nu = 5), student(), dk), 1 / k)
  }
  # a dispersion's column, shifted and scaled
  v <- rnorm(150L)
  yn <- rnbinom(150L, mu = exp(1 + 0.3 * z + u[g]), size = exp(0.7 + 0.3 * v))
  dn <- data.frame(y = yn, z = z, v = v, g = g)
  a5 <- fit1(bf(y ~ z + (1 | g), shape ~ v), negbinomial(), dn)
  jd <- which(names(a5$fit$opt$par) == "betad")[[2L]]
  for (m in list(c(1, 2000), c(1000, 0), c(100, 1e4))) {
    dk <- dn
    dk$v <- m[[1L]] * v + m[[2L]]
    row(sprintf("shape ~ %g v + %g", m[[1L]], m[[2L]]), a5,
        fit1(bf(y ~ z + (1 | g), shape ~ v), negbinomial(), dk), m[[1L]], jd)
  }
  # sigma's regression on a gaussian
  yg <- 0.3 + 0.5 * z + u[g] + rnorm(150L, 0, exp(-0.3 + 0.3 * v))
  dg <- data.frame(y = yg, z = z, v = v, g = g)
  a6 <- fit1(bf(y ~ z + (1 | g), sigma ~ v), gaussian(), dg)
  jd <- which(names(a6$fit$opt$par) == "betad")[[2L]]
  dk <- dg
  dk$v <- v + 2000
  row("gaussian sigma ~ v + 2000", a6,
      fit1(bf(y ~ z + (1 | g), sigma ~ v), gaussian(), dk), 1, jd)
}
