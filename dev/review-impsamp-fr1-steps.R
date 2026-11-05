# Fix round 1, item 1: try to break imp_outer_steps(). Reference: my own exact
# marginal, optimized over coefficients of the column divided by its largest
# absolute value, and mapped back. Student, gaussian and dpar cases use the same
# route with base R densities. Usage: REV_ARM=after|before Rscript ... <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
seed <- as.integer(commandArgs(TRUE)[[1]])
gh <- gh_rule(61L)
ng_f <- function(f, p, h = 1e-4) {
  vapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, h)
    (f(p + e) - f(p - e)) / (2 * h)
  }, 0)
}
# X: fixed design of mu; rows: per-group list; dens(eta, extra, i): row log density
ref_fit <- function(X, g, dens, nextra, start, off = 0) {
  sc <- apply(abs(X), 2, max)
  sc[sc == 0] <- 1
  Xs <- sweep(X, 2, sc, "/")
  idx <- split(seq_len(nrow(X)), g)
  k <- ncol(X)
  nll <- function(p) {
    b <- p[seq_len(k)]
    ex <- if (nextra) p[k + seq_len(nextra)] else NULL
    sd <- exp(p[k + nextra + 1L])
    eta0 <- as.vector(Xs %*% b) + off
    tot <- 0
    for (ii in idx) {
      fn <- function(u) {
        vapply(u, function(ui) sum(dens(eta0[ii] + ui, ex, ii)), 0) +
          dnorm(u, 0, sd, log = TRUE)
      }
      o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
      m <- o$minimum
      cur <- -(fn(m + 1e-3) - 2 * fn(m) + fn(m - 1e-3)) / 1e-6
      s <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
      nodes <- m + sqrt(2) * s * gh$x
      a <- fn(nodes) + gh$x^2 + log(gh$w)
      tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * s))
    }
    tot
  }
  st <- start
  st[seq_len(k)] <- st[seq_len(k)] * sc
  o <- nlminb(st, nll, gradient = function(p) ng_f(nll, p),
              control = list(rel.tol = 1e-10, iter.max = 200))
  H <- numDeriv::jacobian(function(p) ng_f(nll, p), o$par, method = "simple")
  se <- sqrt(abs(diag(solve(0.5 * (H + t(H))))))
  list(par = o$par, se = se, nll = o$objective, sc = c(sc, rep(1, nextra + 1L)), fn = nll)
}
run <- function(lab, f, family, d, X, g, dens, nextra = 0L, draws = 300L, off = 0) {
  w <- character(0)
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    tryCatch(frm(f, family = family, data = d, importance = draws), error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 70))
      invokeRestart("muffleWarning")
    })
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  lap <- suppressWarnings(tryCatch(frm(f, family = family, data = d), error = function(e) NULL))
  if (inherits(fit, "error")) {
    cat(sprintf("%-20s STOP %.0fs: %s\n", lab, secs, substr(conditionMessage(fit), 1, 170)))
    return(invisible(NULL))
  }
  r <- ref_fit(X, g, dens, nextra, unname(if (is.null(lap)) fit$opt$par else lap$opt$par), off)
  p <- unname(fit$opt$par) * r$sc
  z <- (p - r$par) / r$se
  se <- tryCatch(sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed)) * r$sc / r$se,
                 error = function(e) rep(NA, length(p)))
  zl <- if (is.null(lap)) NA else max(abs((unname(lap$opt$par) * r$sc - r$par) / r$se))
  cat(sprintf(
    "%-20s rounds %d capped %s conv %d grad %.3g maxz %.3f (lap %.3f) SE %.3f..%.3f loss %.5f dll %+.4f mcse %.4f ess %.3f warns %d %s\n",
    lab, fit$importance$rounds, fit$importance$capped, fit$opt$convergence,
    fit$importance$grad, max(abs(z)), zl, min(se), max(se), r$fn(p) - r$nll,
    r$nll - fit$opt$objective, fit$importance$mcse, fit$importance$ess_min,
    length(w), paste(w, collapse = "; ")))
}
pdens <- function(y) function(eta, ex, ii) dpois(y[ii], exp(eta), log = TRUE)
bdens <- function(y) function(eta, ex, ii) dbinom(y[ii], 1, plogis(eta), log = TRUE)
G <- 15L
n <- 150L
mk <- function() {
  set.seed(seed)
  g <- factor(rep(seq_len(G), each = 10))
  list(g = g, z = rnorm(n), u = rnorm(G, 0, 0.7)[g])
}
one <- function(lab, x, beta_eff, family = "pois") {
  b <- mk()
  eta <- 0.3 + beta_eff * b$z + b$u
  # the covariate named x is what reaches the design; z drives the data
  y <- if (family == "pois") rpois(n, exp(eta)) else rbinom(n, 1, plogis(eta))
  d <- data.frame(y, x, g = b$g)
  run(paste0(lab, "-", family), bf(y ~ x + (1 | g)),
      if (family == "pois") poisson() else bernoulli(), d, cbind(1, x), b$g,
      if (family == "pois") pdens(y) else bdens(y))
}
b <- mk()
# a column of zeros with one large value
x0 <- numeric(n); x0[1] <- 1e4
one("spike_one_1e4", x0, 0)
x0 <- numeric(n); x0[1:3] <- c(1e6, 1, -1)
one("spike_mixed_1e6", x0, 0)
set.seed(seed + 100)
one("cauchy", rcauchy(n), 0.1)
one("lognorm_sd3", exp(3 * b$z), 0.2)
one("lognorm_sd3", exp(3 * b$z), 0.2, "bern")
xb <- b$z
xb[which.max(abs(xb))] <- 1e8
one("outlier_1e8", xb, 0.3)
xb[which.max(abs(xb))] <- 1e12
one("outlier_1e12", xb, 0.3)
set.seed(seed + 200)
one("dummy_rare", rbinom(n, 1, 0.03), 0.2)
one("dummy_half", rbinom(n, 1, 0.5), 0.2)
# an interaction of two scaled columns, with main effects
set.seed(seed + 300)
a1 <- 1000 * rnorm(n) + 5000
a2 <- 1000 * rnorm(n)
b2 <- mk()
y <- rpois(n, exp(0.5 + 2e-4 * a1 + 2e-4 * a2 + 1e-8 * a1 * a2 - 1 + b2$u))
d <- data.frame(y, a1, a2, g = b2$g)
run("interact-pois", bf(y ~ a1 * a2 + (1 | g)), poisson(), d, cbind(1, a1, a2, a1 * a2), b2$g,
    pdens(y))
# intercept-free
y <- rpois(n, exp(0.5 + 0.3 * b2$z + b2$u))
d <- data.frame(y, xs = 2000 + 5 * b2$z, g = b2$g)
run("nointercept-year", bf(y ~ 0 + xs + (1 | g)), poisson(), d, cbind(d$xs), b2$g, pdens(y))
set.seed(seed + 400)
f3 <- factor(sample(letters[1:3], n, TRUE))
d <- data.frame(y, f3, g = b2$g)
run("nointercept-factor", bf(y ~ 0 + f3 + (1 | g)), poisson(), d, model.matrix(~ 0 + f3), b2$g,
    pdens(y))
# an offset that is large
ex <- 1e6 * exp(rnorm(n, 0, 0.5))
y <- rpois(n, ex * 1e-6 * exp(0.5 + 0.3 * b2$z + b2$u))
d <- data.frame(y, z = b2$z, ex, g = b2$g)
run("offset_1e6-pois", bf(y ~ z + offset(log(ex)) + (1 | g)), poisson(), d, cbind(1, b2$z), b2$g,
    pdens(y), off = log(ex))
# gaussian and student in large units
yg <- 0.3 + 0.5 * b2$z + b2$u + rnorm(n, 0, 0.5)
for (s in c(1e3, 1e6, 1e9)) {
  d <- data.frame(y = yg * s, z = b2$z, g = b2$g)
  gd <- function(s) function(eta, ex, ii) dnorm(d$y[ii], eta, exp(ex), log = TRUE)
}
st <- function(lab, s, nu = 4) {
  yt <- (0.3 + 0.5 * b2$z + b2$u + 0.5 * rt(n, 5)) * s
  d <- data.frame(y = yt, z = b2$z, g = b2$g)
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frm(bf(y ~ z + (1 | g), nu = 4), family = student(), data = d, importance = 300L),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 60))
      invokeRestart("muffleWarning")
    })
  fu <- suppressWarnings(frm(bf(y ~ z + (1 | g), nu = 4), family = student(),
                             data = transform(d, y = y / s), importance = 300L))
  lap <- suppressWarnings(frm(bf(y ~ z + (1 | g), nu = 4), family = student(), data = d))
  if (inherits(fit, "error")) {
    cat(sprintf("%-20s STOP (Laplace ok, ll %.3f): %s\n", lab, as.numeric(logLik(lap)),
                substr(conditionMessage(fit), 1, 210)))
    return(invisible())
  }
  shift <- -n * log(s)
  cat(sprintf("%-20s rounds %d capped %s ll-unit %+.5f mcse %.4f theta diff %.4f warns %d %s\n", lab,
              fit$importance$rounds, fit$importance$capped,
              as.numeric(logLik(fit)) - (as.numeric(logLik(fu)) + shift), fit$importance$mcse,
              fit$opt$par[[length(fit$opt$par)]] - (fu$opt$par[[length(fu$opt$par)]] + log(s)),
              length(w), paste(w, collapse = ";")))
}
for (s in c(10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e9)) st(sprintf("student_y*%g", s), s)
for (s in c(1e3, 1e6, 1e9)) {
  d <- data.frame(y = yg * s, z = b2$z, g = b2$g)
  fg <- suppressWarnings(tryCatch(frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d, importance = 64L),
                                  error = function(e) e))
  lg <- suppressWarnings(frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d))
  if (inherits(fg, "error")) {
    cat(sprintf("gaussian_y*%g STOP: %s\n", s, substr(conditionMessage(fg), 1, 200)))
  } else {
    cat(sprintf("gaussian_y*%g rounds %d ll-lap %+.3g maxdpar/SE %.4f SE ratio %.4f\n", s, fg$importance$rounds,
                as.numeric(logLik(fg)) - as.numeric(logLik(lg)),
                max(abs((fg$opt$par - lg$opt$par) / sqrt(diag(frmtmb:::sdr_of(lg)$cov.fixed)))),
                max(sqrt(diag(frmtmb:::sdr_of(fg)$cov.fixed)) / sqrt(diag(frmtmb:::sdr_of(lg)$cov.fixed)))))
  }
}
# a dispersion regression on a scaled column, and on a large offset
for (cs in list(c(1, 2000), c(1e4, 0), c(100, 1e4), c(1, 1e5))) {
  w0 <- rnorm(n)
  d <- data.frame(y = 0.3 + 0.5 * b2$z + b2$u + rnorm(n, 0, exp(-0.5 + 0.4 * w0)), z = b2$z,
                  w = cs[1] * w0 + cs[2], wc = w0, g = b2$g)
  fr <- suppressWarnings(tryCatch(frm(bf(y ~ z + (1 | g), sigma ~ w), family = gaussian(), data = d,
                                      importance = 64L), error = function(e) e))
  lr <- suppressWarnings(frm(bf(y ~ z + (1 | g), sigma ~ w), family = gaussian(), data = d))
  # Gaussian mu with a dispersion column: Laplace is exact, so it is the reference
  if (inherits(fr, "error")) {
    cat(sprintf("sigma~w %g,%g STOP: %s\n", cs[1], cs[2], substr(conditionMessage(fr), 1, 170)))
  } else {
    cat(sprintf("sigma~w %g,%g rounds %d capped %s ll-lap %+.3g grad %.3g conv %d SE ratio %s\n", cs[1], cs[2],
                fr$importance$rounds, fr$importance$capped,
                as.numeric(logLik(fr)) - as.numeric(logLik(lr)), fr$importance$grad, fr$opt$convergence,
                paste(format(sqrt(diag(frmtmb:::sdr_of(fr)$cov.fixed)) / sqrt(diag(frmtmb:::sdr_of(lr)$cov.fixed)), digits = 4), collapse = " ")))
  }
  y2 <- rgamma(n, shape = exp(0.5 + 0.4 * w0), rate = exp(0.5 + 0.4 * w0) / exp(0.3 + b2$u))
  d2 <- data.frame(y = y2, w = d$w, wc = w0, g = b2$g)
  ga <- suppressWarnings(tryCatch(frm(bf(y ~ 1 + (1 | g), shape ~ w), family = Gamma(link = "log"), data = d2,
                                      importance = 200L), error = function(e) e))
  gc0 <- suppressWarnings(frm(bf(y ~ 1 + (1 | g), shape ~ wc), family = Gamma(link = "log"), data = d2,
                              importance = 200L))
  if (inherits(ga, "error")) {
    cat(sprintf("Gamma shape~w %g,%g STOP: %s\n", cs[1], cs[2], substr(conditionMessage(ga), 1, 170)))
  } else {
    sa <- sqrt(diag(frmtmb:::sdr_of(ga)$cov.fixed))
    sc <- sqrt(diag(frmtmb:::sdr_of(gc0)$cov.fixed))
    cat(sprintf("Gamma shape~w %g,%g rounds %d capped %s ll-centered %+.5f mcse %.4f grad %.3g slopeSE %.4f thetaSE %.4f\n",
                cs[1], cs[2], ga$importance$rounds, ga$importance$capped,
                as.numeric(logLik(ga)) - as.numeric(logLik(gc0)), ga$importance$mcse, ga$importance$grad,
                sa[[3]] * cs[1] / sc[[3]], sa[[4]] / sc[[4]]))
  }
}
offb <- 1000
d <- data.frame(y = yg, z = b2$z, g = b2$g, o = rep(14, n))
fo <- suppressWarnings(frm(bf(y ~ z + (1 | g), sigma ~ 1 + offset(o)), family = gaussian(), data = d, importance = 64L))
lo <- suppressWarnings(frm(bf(y ~ z + (1 | g), sigma ~ 1 + offset(o)), family = gaussian(), data = d))
cat(sprintf("sigma offset 14: rounds %d ll-lap %+.3g betad %.3f SE ratio %s\n", fo$importance$rounds,
            as.numeric(logLik(fo)) - as.numeric(logLik(lo)), fo$opt$par[["betad"]],
            paste(format(sqrt(diag(frmtmb:::sdr_of(fo)$cov.fixed)) / sqrt(diag(frmtmb:::sdr_of(lo)$cov.fixed)), digits = 4), collapse = " ")))
# a family's own parameter with no column on a large scale: an ordinal threshold,
# which absorbs a shifted column: y ~ year + (1 | g) has thresholds near 0.1 * 2000
for (cs in list(c(1, 0), c(5, 2000), c(5, 1e5), c(1000, 1e5))) {
  set.seed(seed + 500)
  lat <- 0.7 * b2$z + b2$u + rlogis(n)
  yo <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
  d <- data.frame(y = yo, w = cs[1] * b2$z + cs[2], wc = b2$z, g = b2$g)
  w <- character(0)
  fo <- withCallingHandlers(
    tryCatch(frm(bf(y ~ w + (1 | g)), family = cumulative(), data = d, importance = 200L), error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 70))
      invokeRestart("muffleWarning")
    })
  fc <- suppressWarnings(frm(bf(y ~ wc + (1 | g)), family = cumulative(), data = d, importance = 200L))
  if (inherits(fo, "error")) {
    cat(sprintf("ordinal %g,%g STOP: %s\n", cs[1], cs[2], substr(conditionMessage(fo), 1, 200)))
  } else {
    sa <- sqrt(diag(frmtmb:::sdr_of(fo)$cov.fixed))
    sc <- sqrt(diag(frmtmb:::sdr_of(fc)$cov.fixed))
    nm <- names(fo$opt$par)
    cat(sprintf("ordinal %g,%g par %s rounds %d capped %s conv %d ll-centered %+.5f mcse %.4f grad %.3g slopeSE %.4f thetaSE %.4f warns %d %s\n",
                cs[1], cs[2], paste(format(fo$opt$par[nm != "beta"], digits = 5), collapse = ","),
                fo$importance$rounds, fo$importance$capped, fo$opt$convergence,
                as.numeric(logLik(fo)) - as.numeric(logLik(fc)), fo$importance$mcse, fo$importance$grad,
                sa[which(nm == "beta")[1]] * cs[1] / sc[which(nm == "beta")[1]],
                sa[length(sa)] / sc[length(sc)], length(w), paste(w, collapse = ";")))
  }
}
