# Release reviewer, 0.73.0, check after the fixes: a hunt for a false
# positive of the new check_weights_unbounded() (a model with a maximum
# that it refuses), and, for refusals near the 1/2-per-censored-row
# threshold, the objective along the spike path (dry_run = "objective",
# which skips the check): a refusal is right only if -fn grows without
# bound as the spike dpar goes to its limit.
#   Rscript dev/rel073-rev-fp.R <lib> <part: fit|path|ext>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
part <- a[2]
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib: frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "| part", part, "\n")
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(sprintf("%-62s REFUSED/ERROR: %s\n", label,
                substr(gsub("[[:space:]]+", " ", r), 1, 110)))
  } else {
    cat(sprintf("%-62s fit: logLik %.4f conv %d betad %s | warnings %d\n",
                label, as.numeric(logLik(r)), r$opt$convergence,
                paste(format(r$estimates$betad, digits = 4), collapse = " "),
                length(w)))
  }
  invisible(r)
}
# -fn along the spike path: the spike dpar's intercept moved by `steps`
# from a start, everything else held
path <- function(label, f, data, family, steps, start_betad, ...) {
  obj <- tryCatch(suppressWarnings(suppressMessages(
    frm(f, family = family, data = data, dry_run = "objective", ...))),
    error = function(e) conditionMessage(e))
  if (is.character(obj)) {
    cat(label, ": objective ERROR", substr(obj, 1, 120), "\n"); return()
  }
  o <- if (!is.null(obj$obj)) obj$obj else obj
  p <- o$par
  i <- which(names(p) == "betad")[1]
  v <- vapply(steps, function(s) {
    q <- p; q[i] <- start_betad + s
    tryCatch(-o$fn(q), error = function(e) NA_real_)
  }, 0)
  cat(sprintf("%-62s -fn at betad %s: %s\n", label,
              paste(start_betad + steps, collapse = ","),
              paste(format(v, digits = 6), collapse = "  ")))
}
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g4 = factor(rep(1:20, each = 4)), g6 = factor(rep(1:16, each = 5)))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$g4] + rnorm(n)
q8 <- quantile(d$y, 0.8)
d$c20 <- ifelse(d$y > q8, "right", "none"); d$y20 <- pmin(d$y, q8)
q1 <- quantile(d$y, 0.125)
d$c87 <- ifelse(d$y > q1, "right", "none"); d$y87 <- pmin(d$y, q1)  # 70 censored
d$gy <- rgamma(n, 3, 3 / exp(0.3 + 0.2 * d$x))
qg <- quantile(d$gy, 0.8)
d$cg <- ifelse(d$gy > qg, "right", "none"); d$gyc <- pmin(d$gy, qg)
d$yb <- plogis(0.2 + 0.4 * d$x); d$yb <- rbeta(n, d$yb * 6, (1 - d$yb) * 6)
qb <- quantile(d$yb, 0.8)
d$cb <- ifelse(d$yb > qb, "right", "none"); d$ybc <- pmin(d$yb, qb)
d$yz <- ifelse(runif(n) < 0.2, 0, d$y)
d$yzz <- ifelse(runif(n) < 0.9, 0, d$y)  # mostly zeros
cen <- d$c20 == "right"
if (part == "fit") {
  cat("== absent cases: a maximum exists, the check must not refuse\n")
  d$w15 <- 1.5; d$w09 <- 0.9; d$w1 <- 1
  run("A gaussian cens 20% w 1.5 (1 | g4)",
      frm(y20 | weights(w15) + cens(c20) ~ x + (1 | g4), data = d))
  run("A2 Gamma cens 20% w 1.5 (1 | g4)",
      frm(gyc | weights(w15) + cens(cg) ~ x + (1 | g4),
          family = Gamma(link = "log"), data = d))
  run("A3 asym_laplace cens 20% (1 | g4)",
      frm(y20 | cens(c20) ~ x + (1 | g4), family = asym_laplace(), data = d))
  run("B gaussian 70 of 80 censored, w 0.9 (1 | obs)",
      frm(y87 | weights(w09) + cens(c87) ~ x + (1 | obs), data = d))
  run("B2 gaussian 70 of 80 censored, w 1 (1 | obs)",
      frm(y87 | cens(c87) ~ x + (1 | obs), data = d))
  run("C quadrature cens 20% w 1 (1 | obs)",
      frm(y20 | cens(c20) ~ x + (1 | obs), data = d, quadrature = TRUE))
  run("C2 quadrature cens 20% w 0.9 (1 | obs)",
      frm(y20 | weights(w09) + cens(c20) ~ x + (1 | obs), data = d,
          quadrature = TRUE))
  d$wmix <- ifelse(cen, 2, 0.9)  # sum 99.2 > 80, uncensored excess -6.4
  cat("   wmix sum", sum(d$wmix), "| uncensored excess", sum(d$wmix[!cen] - 1),
      "\n")
  run("C3 quadrature censored w 2, uncensored w 0.9 (1 | obs)",
      frm(y20 | weights(wmix) + cens(c20) ~ x + (1 | obs), data = d,
          quadrature = TRUE))
  run("D Laplace censored w 3, uncensored w 1 (1 | obs)",
      frm(y20 | weights(wcen3) + cens(c20) ~ x + (1 | obs),
          data = transform(d, wcen3 = ifelse(cen, 3, 1))))
  run("D2 Beta cens 20% w 1 (1 | obs)",
      frm(ybc | cens(cb) ~ x + (1 | obs), family = Beta(), data = d))
  run("D3 Gamma cens 20% w 1 (1 | obs)",
      frm(gyc | cens(cg) ~ x + (1 | obs), family = Gamma(link = "log"),
          data = d))
  run("E zi_asym_laplace (1 | g4)",
      frm(yz ~ x + (1 | g4), family = zero_inflated_asym_laplace(), data = d))
  run("E2 zi_asym_laplace 90% zeros (1 | g4)",
      frm(yzz ~ x + (1 | g4), family = zero_inflated_asym_laplace(), data = d))
  run("E3 zi_asym_laplace 90% zeros (1 | obs) Laplace",
      frm(yzz ~ x + (1 | obs), family = zero_inflated_asym_laplace(), data = d))
  run("E4 zi_asym_laplace 90% zeros (1 | obs) quadrature w 1",
      frm(yzz ~ x + (1 | obs), family = zero_inflated_asym_laplace(), data = d,
          quadrature = TRUE))
  run("E5 zi_asym_laplace 90% zeros (1 | g4) quadrature",
      frm(yzz ~ x + (1 | g4), family = zero_inflated_asym_laplace(), data = d,
          quadrature = TRUE))
  d$wz <- ifelse(d$yzz == 0, 3, 0.9)
  run("E6 zi_asym_laplace quadrature, zeros w 3, nonzero w 0.9 (1 | obs)",
      frm(yzz | weights(wz) ~ x + (1 | obs),
          family = zero_inflated_asym_laplace(), data = d, quadrature = TRUE))
  run("F gaussian interval-censored 20%, w 1.5 on censored only (1 | obs)",
      frm(y20 | weights(wc) + cens(ci, y2i) ~ x + (1 | obs),
          data = transform(d, wc = ifelse(cen, 1.5, 1),
                           ci = ifelse(cen, "interval", "none"),
                           y2i = ifelse(cen, y20 + 1, y20))))
}
if (part == "path") {
  cat("== refusals near the threshold: is -fn unbounded along the spike?\n")
  # gaussian: 16 censored at w 1, uncensored excess 9 (> 16 / 2): refused
  d$wt9 <- ifelse(cen, 1, 1 + 9 / sum(!cen))
  run("G gaussian 16 cens, E_unc 9 (1 | obs) Laplace",
      frm(y20 | weights(wt9) + cens(c20) ~ x + (1 | obs), data = d))
  path("G path", bf(y20 | weights(wt9) + cens(c20) ~ x + (1 | obs)), d,
       gaussian(), c(0, -4, -8, -16, -32, -64), 0)
  d$wt85 <- ifelse(cen, 1, 1 + 8.5 / sum(!cen))
  run("G2 gaussian 16 cens, E_unc 8.5 (1 | obs) Laplace",
      frm(y20 | weights(wt85) + cens(c20) ~ x + (1 | obs), data = d))
  path("G2 path", bf(y20 | weights(wt85) + cens(c20) ~ x + (1 | obs)), d,
       gaussian(), c(0, -4, -8, -16, -32, -64), 0)
  # Gamma, shape to infinity
  run("H Gamma 16 cens, E_unc 9 (1 | obs) Laplace",
      frm(gyc | weights(wt9) + cens(cg) ~ x + (1 | obs),
          family = Gamma(link = "log"), data = d))
  path("H path", bf(gyc | weights(wt9) + cens(cg) ~ x + (1 | obs)), d,
       Gamma(link = "log"), c(0, 4, 8, 16, 32, 64), 1)
  # Beta, phi to infinity
  run("I Beta 16 cens, E_unc 9 (1 | obs) Laplace",
      frm(ybc | weights(wt9) + cens(cb) ~ x + (1 | obs), family = Beta(),
          data = d))
  path("I path", bf(ybc | weights(wt9) + cens(cb) ~ x + (1 | obs)), d,
       Beta(), c(0, 4, 8, 16, 32), 1)
  # asym_laplace, unweighted: 30 censored, 50 uncensored: 50 - 15 > 0
  q6 <- quantile(d$y, 0.625)
  d$c30 <- ifelse(d$y > q6, "right", "none"); d$y30 <- pmin(d$y, q6)
  run("J asym_laplace 30 cens (1 | obs) Laplace",
      frm(y30 | cens(c30) ~ x + (1 | obs), family = asym_laplace(), data = d))
  path("J path", bf(y30 | cens(c30) ~ x + (1 | obs)), d, asym_laplace(),
       c(0, -4, -8, -16, -32), 0)
  # all rows censored but a few: asym_laplace 70 censored, 10 uncensored
  # (10 - 35 < 0: not refused, the filed miss); is it unbounded?
  run("K asym_laplace 70 cens (1 | obs) Laplace",
      frm(y87 | cens(c87) ~ x + (1 | obs), family = asym_laplace(), data = d))
  path("K path", bf(y87 | cens(c87) ~ x + (1 | obs)), d, asym_laplace(),
       c(0, -4, -8, -16, -32), 0)
  # zi asym_laplace, 90% zeros, (1 | obs), Laplace: refused; unbounded?
  path("E3 path", bf(yzz ~ x + (1 | obs)), d, zero_inflated_asym_laplace(),
       c(0, -4, -8, -16, -32), 0)
}
if (part == "bypass") {
  # The refusals near the threshold, fitted with the check switched off
  # in this process only, from the default start and from a start with
  # the spike dpar's intercept far toward its limit. A refusal is a false
  # positive if both fits stop at the same interior maximum.
  utils::assignInNamespace("check_weights_unbounded",
                           function(...) invisible(), "frmtmb")
  d$wt9 <- ifelse(cen, 1, 1 + 9 / sum(!cen))
  d$wt85 <- ifelse(cen, 1, 1 + 8.5 / sum(!cen))
  d$wt05 <- ifelse(cen, 1, 1 + 0.5 / sum(!cen))
  q6 <- quantile(d$y, 0.625)
  d$c30 <- ifelse(d$y > q6, "right", "none"); d$y30 <- pmin(d$y, q6)
  d$w15 <- 1.5
  cases <- list(
    list("G gaussian 16 cens, E_unc 9", bf(y20 | weights(wt9) + cens(c20) ~
           x + (1 | obs)), gaussian(), -12, FALSE),
    list("G2 gaussian 16 cens, E_unc 8.5", bf(y20 | weights(wt85) +
           cens(c20) ~ x + (1 | obs)), gaussian(), -12, FALSE),
    list("G3 gaussian 16 cens, E_unc 0.5, quadrature",
         bf(y20 | weights(wt05) + cens(c20) ~ x + (1 | obs)), gaussian(),
         -12, TRUE),
    list("J asym_laplace 30 cens, unweighted", bf(y30 | cens(c30) ~ x +
           (1 | obs)), asym_laplace(), -12, FALSE),
    list("E3 zi_asym_laplace 90% zeros (1 | obs)", bf(yzz ~ x + (1 | obs)),
         zero_inflated_asym_laplace(), -12, FALSE),
    list("N Gamma cens(cg) 20%, w 1.5", bf(gyc | weights(w15) + cens(cg) ~
           x + (1 | obs)), Gamma(link = "log"), 12, FALSE),
    list("O Beta cens(cb) 20%, w 1.5", bf(ybc | weights(w15) + cens(cb) ~
           x + (1 | obs)), Beta(), 12, FALSE))
  for (cs in cases) {
    run(paste(cs[[1]], "| default start"),
        frm(cs[[2]], family = cs[[3]], data = d, quadrature = cs[[5]]))
    tpl <- tryCatch(par_template(cs[[2]], data = d, family = cs[[3]]),
                    error = function(e) NULL)
    if (!is.null(tpl)) {
      tpl$betad[1] <- cs[[4]]
      run(paste(cs[[1]], "| start betad", cs[[4]]),
          frm(cs[[2]], family = cs[[3]], data = d, quadrature = cs[[5]],
              start = tpl))
    }
  }
}
if (part == "exact") {
  # G3 by hand: a gaussian with (1 | obs) has a closed-form marginal.
  # An uncensored row with weight w integrates to
  # (2 pi s^2)^(-(w-1)/2) w^(-1/2) N(y; mu, sqrt(s^2 / w + tau^2)), and a
  # right-censored one to the integral of Phi((mu + u - c) / s)^w over
  # N(0, tau), which lies in (0, 1]. So the exact log-likelihood grows
  # like -(sum over uncensored rows of (w - 1)) log s as s goes to 0.
  utils::assignInNamespace("check_weights_unbounded",
                           function(...) invisible(), "frmtmb")
  d$wt05 <- ifelse(cen, 1, 1 + 0.5 / sum(!cen))
  f <- suppressWarnings(frm(y20 | weights(wt05) + cens(c20) ~ x + (1 | obs),
                            data = d, quadrature = TRUE))
  b <- f$estimates$beta; tau <- exp(f$estimates$theta[1])
  mu <- b[1] + b[2] * d$x; w <- d$wt05
  exact <- function(ls) {
    s <- exp(ls)
    u <- !cen
    lu <- -(w[u] - 1) / 2 * log(2 * pi * s^2) - 0.5 * log(w[u]) +
      dnorm(d$y20[u], mu[u], sqrt(s^2 / w[u] + tau^2), log = TRUE)
    lc <- vapply(which(cen), function(i) {
      g <- function(z) exp(w[i] * pnorm((mu[i] + z - d$y20[i]) / s,
                                        log.p = TRUE)) * dnorm(z, 0, tau)
      log(integrate(g, -Inf, Inf, rel.tol = 1e-10,
                    subdivisions = 2000L)$value)
    }, 0)
    sum(lu) + sum(lc)
  }
  cat(sprintf("G3 quadrature fit: logLik %.6f at log sigma %.4f (conv %d)\n",
              as.numeric(logLik(f)), f$estimates$betad[1], f$opt$convergence))
  for (ls in c(f$estimates$betad[1], -6, -10, -20, -40, -80, -160)) {
    cat(sprintf("  exact log-likelihood at log sigma %8.3f: %.6f\n", ls,
                exact(ls)))
  }
}
if (part == "profile") {
  # The Laplace objective profiled over everything but log sigma, at
  # uncensored excesses just above and below 16 / 2 (16 censored rows).
  # The check refuses above 8, so a refusal is wrong if the profile is
  # bounded there. The derivation says the slope in -log sigma is
  # E_unc - 8, plus terms that grow more slowly than a power.
  prof <- function(label, E) {
    d$wE <- ifelse(cen, 1, 1 + E / sum(!cen))
    o <- suppressWarnings(suppressMessages(frm(
      y20 | weights(wE) + cens(c20) ~ x + (1 | obs), data = d,
      dry_run = "objective")))
    o <- if (!is.null(o$obj)) o$obj else o
    p <- o$par; i <- which(names(p) == "betad")[1]
    rest <- p[-i]
    out <- c()
    for (ls in c(-2, -4, -6, -8, -10, -12, -14, -16)) {
      fr <- function(r) { q <- p; q[-i] <- r; q[i] <- ls
        v <- tryCatch(o$fn(q), error = function(e) NA_real_)
        if (!is.finite(v)) 1e10 else v }
      op <- nlminb(rest, fr)
      if (op$objective < 1e9) rest <- op$par
      out <- c(out, -op$objective)
    }
    cat(sprintf("%-28s profile -fn at log sigma -2..-16 by 2: %s\n", label,
                paste(format(out, digits = 7), collapse = " ")))
    cat(sprintf("%-28s slopes per unit of -log sigma: %s\n", "",
                paste(format(diff(out) / 2, digits = 3), collapse = " ")))
  }
  prof("E_unc 7.5 (not refused)", 7.5)
  prof("E_unc 8.1 (refused)", 8.1)
  prof("E_unc 8.5 (refused)", 8.5)
  prof("E_unc 9 (refused)", 9)
}
if (part == "ext") {
  suppressMessages(library(frmtmb.eam))
  cat("frmtmb.eam", format(packageVersion("frmtmb.eam")), "\n")
  set.seed(29)
  dd <- ddm_simulate(120, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  dd$code <- 0L
  dd$code[dd$rt > 1.1] <- 1L
  dd$y <- pmin(dd$rt, 1.1)
  dd$obs <- factor(seq_len(nrow(dd))); dd$w <- 1.5
  cat("wiener censored rows", sum(dd$code), "\n")
  run("L wiener cens (1 | obs) w 1.5",
      frm(bf(y | dec(upper) + weights(w) + cens(code) ~ 1 + (1 | obs),
             bias = 0.5), family = wiener(), data = dd))
  set.seed(7)
  dr <- rdm_simulate(150L, v = c(3.0, 2.0), A = 0.5, k = 0.5, ndt = 0.2)
  cut <- stats::quantile(dr$rt, 0.75)
  dr$cens <- as.integer(dr$rt > cut)
  dr$ob <- pmin(dr$rt, cut)
  dr$choice[dr$cens == 1] <- 1L
  dr$obs <- factor(seq_len(nrow(dr))); dr$w <- 1.5
  run("M rdm cens (1 | obs) w 1.5",
      frm(bf(ob | vint(choice) + weights(w) + cens(cens) ~ 1 + (1 | obs)),
          family = rdm(2), data = dr))
}
cat("END\n")
