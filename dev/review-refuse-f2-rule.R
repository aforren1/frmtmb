# Fix round 2, N5: rules for the read-row warning against shifts in exact SE.
#   R0  the lane's: |total by integrate(), scaled from 30 rows| > 1
#   R1  the coordinator's: R0, then per parameter of the read response one more
#       integrate() sweep on the same rows at +1 SE; warn on the largest move
#   R1s the same with ONE sweep, every dpar of the read response at +1 SE
#   R2  no integrate(): Tierney-Kadane's first term per value from the taped
#       inner Hessian at u and u +- h (3 Hessian sweeps, every row), central
#       differences in the outer parameters, predicted shift V grad T in SE.
# A move of 1 SE changes the objective's slope by 1, so a total that moves by d
# over 1 SE predicts a shift of d SE: the threshold on R1 and R2 is 0.5.
# Usage: Rscript dev/review-refuse-f2-rule.R table1 <family> | table2 <slope> |
#        vm | cost. Seeds: table1 11200 + family, table2 11501 + slope index,
#        vm 11600, cost 11502.
args <- commandArgs(trailingOnly = TRUE)
source("dev/review-refuse-lib.R")
fams <- c("Gamma", "exponential", "beta", "exgaussian", "student", "lognormal", "von_mises", "gaussian")
dexg <- function(y, mu, sig, bet) {
  m <- mu - bet
  -log(bet) + (m - y) / bet + sig^2 / (2 * bet^2) + pnorm((y - m) / sig - sig / bet, log.p = TRUE)
}
speco <- function(fam, gridn = 1) {
  sp <- switch(fam,
    Gamma = list(np = 1, g = seq(-14, 7, length.out = 801), inv = exp, lj = identity,
      ld = function(y, eta, p) dgamma(y, exp(p), exp(p) / exp(eta), log = TRUE),
      rgen = function(eta) rgamma(length(eta), 1.3, 1.3 / exp(eta))),
    exponential = list(np = 0, g = seq(-14, 7, length.out = 801), inv = exp, lj = identity,
      ld = function(y, eta, p) dexp(y, exp(-eta), log = TRUE),
      rgen = function(eta) rexp(length(eta), exp(-eta))),
    lognormal = list(np = 1, g = seq(-14, 12, length.out = 801), inv = exp, lj = identity,
      ld = function(y, eta, p) dlnorm(y, eta, exp(p), log = TRUE),
      rgen = function(eta) rlnorm(length(eta), eta, 0.9)),
    beta = list(np = 1, g = seq(-30, 30, length.out = 1601), inv = plogis, lj = function(u) -abs(u) - 2 * log1p(exp(-abs(u))),
      ld = function(y, eta, p) { m <- plogis(eta); dbeta(y, m * exp(p), (1 - m) * exp(p), log = TRUE) },
      rgen = function(eta) rbeta(length(eta), plogis(eta) * 3, (1 - plogis(eta)) * 3)),
    exgaussian = list(np = 2, g = seq(-15, 40, length.out = 1401), inv = identity, lj = function(u) 0 * u,
      ld = function(y, eta, p) dexg(y, eta, exp(p[1]), exp(p[2])),
      rgen = function(eta) eta - 2 + rnorm(length(eta), 0, 0.5) + rexp(length(eta), 1 / 2)),
    student = list(np = 2, g = seq(-80, 80, length.out = 3201), inv = identity, lj = function(u) 0 * u,
      ld = function(y, eta, p) dt((y - eta) / exp(p[1]), 1 + exp(p[2]), log = TRUE) - p[1],
      rgen = function(eta) eta + rt(length(eta), 3)),
    gaussian = list(np = 1, g = seq(-20, 20, length.out = 1201), inv = identity, lj = function(u) 0 * u,
      ld = function(y, eta, p) dnorm(y, eta, exp(p), log = TRUE),
      rgen = function(eta) rnorm(length(eta), eta, 1)),
    von_mises = list(np = 1, g = seq(-pi, pi, length.out = 1001), inv = identity, lj = function(u) 0 * u,
      ld = function(y, eta, p) exp(p) * cos(y - 2 * atan(eta)) - log(2 * pi) - log(besselI(exp(p), 0)),
      rgen = function(eta) {
        m <- 2 * atan(eta); k <- 3
        vapply(m, function(mi) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + mi + pi) %% (2 * pi)) - pi }, 0)
      }))
  sp
}
exact <- function(spec, d, start) {
  rows <- which(is.na(d$y)); ob <- setdiff(seq_len(nrow(d)), rows)
  vals <- spec$inv(spec$g); lj <- spec$lj(spec$g); dg <- diff(spec$g)[1]
  nll <- function(p) {
    by <- p[1:2]; dp <- if (spec$np) p[2 + seq_len(spec$np)]
    bz <- p[(3 + spec$np):(5 + spec$np)]; sz <- exp(p[6 + spec$np])
    ey <- by[1] + by[2] * d$x
    s <- sum(spec$ld(d$y[ob], ey[ob], dp)) +
      sum(dnorm(d$z[ob], bz[1] + bz[2] * d$y[ob] + bz[3] * d$x[ob], sz, log = TRUE))
    for (i in rows) {
      li <- spec$ld(vals, ey[i], dp) + lj + dnorm(d$z[i], bz[1] + bz[2] * vals + bz[3] * d$x[i], sz, log = TRUE)
      li[!is.finite(li)] <- -Inf
      mx <- max(li); s <- s + mx + log(sum(exp(li - mx)) * dg)
    }
    if (!is.finite(s)) return(1e10)
    -s
  }
  o <- nlminb(start, nll, control = list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-11))
  H <- stats::optimHess(o$par, nll)
  list(par = o$par, se = sqrt(diag(solve(H))), ll = -o$objective, fn = nll)
}
# the row errors at any outer point, as the lane's, on the same rows
int_T <- function(fit, par, at, n_read) {
  env <- fit$obj$env
  fit$obj$fn(par)
  lp <- env$last.par
  mp <- which(names(lp) == "miss")
  hk <- Matrix::diag(env$spHess(lp, random = TRUE))[match(mp, env$random)]
  f0 <- env$f(lp)
  e <- vapply(at, function(k) {
    s <- 1 / sqrt(hk[k]); u0 <- lp[mp[k]]
    if (!is.finite(s)) return(NA_real_)
    h <- function(z) vapply(z, function(zi) {
      p <- lp; p[mp[k]] <- u0 + s * zi
      v <- exp(f0 - env$f(p)); if (is.finite(v)) v else 0
    }, 0)
    i <- stats::integrate(h, -Inf, Inf, subdivisions = 100L, rel.tol = 0.1, stop.on.error = FALSE)$value
    if (!is.finite(i) || i <= 0) NA_real_ else 0.5 * log(2 * pi) - log(i)
  }, 0)
  mean(e, na.rm = TRUE) * n_read
}
# Tierney-Kadane's first term, every row, from three taped Hessians
tk_T <- function(fit, par, h = 0.05) {
  env <- fit$obj$env
  fit$obj$fn(par)
  lp <- env$last.par
  mp <- which(names(lp) == "miss")
  ix <- match(mp, env$random)
  d2 <- function(p) Matrix::diag(env$spHess(p, random = TRUE))[ix]
  f2 <- d2(lp)
  # the step in units of each value's own width, so that a well-pinned value
  # is read on its own scale
  st <- h * 1 / sqrt(pmax(f2, 1e-8)) * 4
  up <- lp; up[mp] <- up[mp] + st
  dn <- lp; dn[mp] <- dn[mp] - st
  fp <- d2(up); fm <- d2(dn)
  f3 <- (fp - fm) / (2 * st)
  f4 <- (fp - 2 * f2 + fm) / st^2
  e <- f4 / (8 * f2^2) - 5 * f3^2 / (24 * f2^3)
  e[!is.finite(e)] <- 0
  sum(e)
}
rules <- function(a) {
  p <- a$opt$par
  t0 <- proc.time()[["elapsed"]]
  H <- stats::optimHess(p, a$obj$fn, a$obj$gr)
  V <- solve(H)
  se <- sqrt(diag(V))
  th <- proc.time()[["elapsed"]] - t0
  env <- a$obj$env
  ml <- a$mi_laplace
  n_read <- length(a$estimates$miss)
  at <- if (n_read > 30L) unique(round(seq(1, n_read, length.out = 30L))) else seq_len(n_read)
  T0 <- ml$total %||% NA_real_
  pn <- names(p)
  # the read response's own parameters, by the estimates' names
  bn <- names(a$estimates$beta); dn <- names(a$estimates$betad)
  ibeta <- which(pn == "beta"); idp <- which(pn == "betad")
  own <- c(ibeta[grep("^y_", bn)], idp[grep("^y_", dn)])
  ydp <- idp[grep("^y_", dn)]
  if (!length(ydp)) ydp <- ibeta[grep("^y_.Intercept", bn)]
  t0 <- proc.time()[["elapsed"]]
  r1 <- NA_real_; r1s <- NA_real_; t1 <- NA; t1s <- NA
  if (is.finite(T0) && abs(T0) > 1) {
    Tb <- int_T(a, p, at, n_read)
    t00 <- proc.time()[["elapsed"]]
    p1 <- p; p1[ydp] <- p1[ydp] + se[ydp]
    r1s <- abs(int_T(a, p1, at, n_read) - Tb)
    t1s <- proc.time()[["elapsed"]] - t00
    t00 <- proc.time()[["elapsed"]]
    r1 <- max(vapply(own, function(j) {
      q <- p; q[j] <- q[j] + se[j]
      abs(int_T(a, q, at, n_read) - Tb)
    }, 0))
    t1 <- proc.time()[["elapsed"]] - t00
  }
  t0 <- proc.time()[["elapsed"]]
  T2 <- tk_T(a, p)
  g <- vapply(seq_along(p), function(j) {
    d <- 0.5 * se[j]
    q1 <- p; q1[j] <- q1[j] + d; q2 <- p; q2[j] <- q2[j] - d
    (tk_T(a, q1) - tk_T(a, q2)) / (2 * d)
  }, 0)
  pred <- as.numeric(V %*% g) / se
  t2 <- proc.time()[["elapsed"]] - t0
  a$obj$fn(p)
  list(T0 = T0, r1 = r1, r1s = r1s, T2 = T2, pred = pred, se = se, th = th, t1 = t1, t1s = t1s, t2 = t2, names = pn)
}
line <- function(tag, a, ex, p0) {
  r <- rules(a)
  sh <- (p0 - ex$par) / ex$se
  mx <- max(abs(sh)); pmx <- max(abs(r$pred))
  # p0 is ordered as exact()'s vector: y's, then z's. Find each one's order for
  # the prediction by name is not needed for the maxima
  gapll <- as.numeric(logLik(a)) - ex$ll
  pr("RULE %s truth_max %.3f gap_logLik %.3f T0 %.3f R1 %s R1s %s T2 %.3f R2 %.3f th %.2f t1 %s t1s %s t2 %.2f",
     tag, mx, gapll, r$T0, format(r$r1, digits = 3), format(r$r1s, digits = 3), r$T2, pmx, r$th, format(r$t1, digits = 3), format(r$t1s, digits = 3), r$t2)
}
pars0 <- function(a) {
  est <- a$estimates; nm <- names(est$beta); dn <- names(est$betad)
  by0 <- est$beta[grep("^y_", nm)]; bz0 <- est$beta[grep("^z_", nm)]
  dy <- est$betad[grep("^y_", dn)]; sz0 <- est$betad[grep("^z_", dn)]
  c(by0, dy, bz0[c(1, grep("miy", names(bz0)), grep("_x$", names(bz0)))], sz0)
}
one_fit <- function(fam, d, tag) {
  ff <- switch(fam, beta = Beta(), Gamma = stats::Gamma(link = "log"), gaussian = gaussian(), get(fam, asNamespace("frmtmb"))())
  a <- tr(frm(bf(y | mi() ~ x) + ff + bf(z ~ mi(y) + x) + gaussian(), data = d))
  if (is_err(a)) { pr("RULE %s frm error %s", tag, msg(a, 120)); return(invisible()) }
  p0 <- pars0(a)
  sp <- speco(fam)
  ex <- tr(exact(sp, d, p0))
  if (is_err(ex)) { pr("RULE %s exact error %s", tag, msg(ex, 120)); return(invisible()) }
  rr <- tr({ line(tag, a, ex, p0); TRUE })
  if (is_err(rr)) pr("RULE %s rule error %s", tag, msg(rr, 160))
}
mode <- args[1]
if (mode == "table1") {
  fam <- args[2]; fi <- match(fam, fams); sp <- speco(fam)
  set.seed(11200 + fi)
  for (cfg in list(c(200, 60, 0.7, 0.7), c(200, 60, 0.05, 1), c(200, 100, 0.7, 0.7), c(200, 15, 0.7, 0.7))) {
    n <- cfg[1]; m <- cfg[2]; bz <- cfg[3]; sz <- cfg[4]
    x <- rnorm(n)
    eta <- if (fam == "beta") 0.3 + 0.7 * x else 0.3 + 0.5 * x
    y <- sp$rgen(eta)
    d <- data.frame(x = x, y = y)
    d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, sz)
    d$y[sample(n, m)] <- NA
    one_fit(fam, d, sprintf("t1 %s read %d slope %.2f", fam, m, bz))
  }
} else if (mode == "table2") {
  sl <- as.numeric(args[2])
  set.seed(11501 + round(sl * 100))
  for (mm in c(20, 60)) for (k in 1:8) {
    n <- 200; x <- rnorm(n)
    y <- 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2
    d <- data.frame(x = x, y = y); d$z <- 0.5 + sl * y + 0.3 * x + rnorm(n, 0, 1)
    d$y[sample(n, mm)] <- NA
    one_fit("exgaussian", d, sprintf("t2 exgaussian read %d slope %.2f rep %d", mm, sl, k))
  }
} else if (mode == "vm") {
  # item 2: the window's end, with a reader, on my integrate
  set.seed(11600)
  rvm <- function(mu, k) vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
  for (cfg in list(c(0, 3), c(3.05, 2), c(3.05, 3), c(3.14, 20))) for (mm in c(15, 60)) for (bz in c(0.2, 0.7)) for (k in 1:3) {
    n <- 200; x <- rnorm(n)
    y <- rvm(cfg[1] + 0.05 * x, cfg[2])
    d <- data.frame(x = x, y = y); d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, 0.7)
    rows <- sample(n, mm); d$y[rows] <- NA
    a <- tr(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(), data = d))
    if (is_err(a)) { pr("VM %g %g read %d slope %.1f: %s", cfg[1], cfg[2], mm, bz, msg(a, 100)); next }
    warned <- "frmtmb_mi_laplace_error" %in% "x"
    w <- attr(a, "w")
    sp <- speco("von_mises")
    p0 <- pars0(a)
    ex <- tr(exact(sp, d, p0))
    if (is_err(ex)) { pr("VM exact error %s", msg(ex, 100)); next }
    env <- a$obj$env; pv <- env$last.par.best; ii <- which(names(pv) == "miss"); b0 <- env$f(pv)
    H <- Matrix::diag(env$spHess(pv, random = TRUE))[match(ii, env$random)]
    err <- vapply(seq_along(ii), function(q) {
      h <- function(v) vapply(v, function(s) exp(b0 - env$f(replace(pv, ii[q], s))), 0)
      I <- integrate(h, -40, 40, subdivisions = 4000L, stop.on.error = FALSE)$value
      0.5 * log(2 * pi) - 0.5 * log(H[q]) - log(I)
    }, 0)
    r <- tr(rules(a))
    sh <- max(abs((p0 - ex$par) / ex$se))
    pr("VM mean %.2f kappa %g read %d slope %.1f rep %d: per value my integrate %+.3f, lane's total %.2f; logLik - exact %+.2f; truth_max %.2f; warned %s; R2 %s; T2 %s; R1 %s; conv %d; exact at frm-objective gain %s",
       cfg[1], cfg[2], mm, bz, k, mean(err), a$mi_laplace$total %||% NA, as.numeric(logLik(a)) - ex$ll, sh, any(grepl("Laplace approximation of the", w)),
       if (is_err(r)) "-" else format(max(abs(r$pred)), digits = 3), if (is_err(r)) "-" else format(r$T2, digits = 3), if (is_err(r)) "-" else format(r$r1, digits = 3), a$opt$convergence, format(a$opt$objective - a$obj$fn(replace(a$opt$par, seq_along(ex$par), ex$par[match(seq_along(ex$par), seq_along(ex$par))])), digits = 4))
  }
} else if (mode == "cost") {
  orig <- get("mi_laplace_note", asNamespace("frmtmb"))
  set.seed(11502)
  for (cfg in list(c(200, 60), c(200, 120), c(1000, 300), c(5000, 1500))) {
    n <- cfg[1]; m <- cfg[2]; x <- rnorm(n)
    y <- rgamma(n, 1.3, 1.3 / exp(0.3 + 0.5 * x))
    d <- data.frame(x = x, y = y); d$z <- 0.5 + 0.7 * y + 0.3 * x + rnorm(n, 0, 0.7)
    d$y[sample(n, m)] <- NA
    fm <- function() frm(bf(y | mi() ~ x) + stats::Gamma(link = "log") + bf(z ~ mi(y) + x) + gaussian(), data = d)
    assignInNamespace("mi_laplace_note", function(fit) fit, "frmtmb")
    bare <- median(vapply(1:3, function(i) system.time(a0 <- suppressWarnings(fm()))[["elapsed"]], 0))
    assignInNamespace("mi_laplace_note", orig, "frmtmb")
    a <- suppressWarnings(fm())
    t0 <- median(vapply(1:3, function(i) system.time(suppressWarnings(orig(a)))[["elapsed"]], 0))
    r <- rules(a)
    # the fit already holds a Hessian for its SEs: th is an upper bound, and R2
    # without the outer part is the one-point total
    tt <- median(vapply(1:3, function(i) system.time(tk_T(a, a$opt$par))[["elapsed"]], 0))
    pr("COST n=%d read=%d: bare fit %.2fs; R0 %.2fs (%.0f%%); R1 extra %.2fs (%.0f%%, only past the gate); R1s extra %.2fs (%.0f%%); R2 %.2fs (%.0f%%), of which one total %.3fs; Hessian for SEs %.2fs; p=%d", n, m, bare, t0, 100 * t0 / (bare + t0),
       r$t1, 100 * r$t1 / bare, r$t1s, 100 * r$t1s / bare, r$t2, 100 * r$t2 / bare, tt, r$th, length(a$opt$par))
  }
}
