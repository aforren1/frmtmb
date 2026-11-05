# Fix round 2, N5: the candidate ruling, on the same data as review-refuse-f2-rule.R
# (same seeds, same order, no exact fit here: its truth is joined by tag).
#   R3  gate: |Tierney-Kadane total| > 1, three Hessian sweeps, every row.
#       Past the gate: 10 evenly spaced rows, integrate() at the estimate and at
#       +1 SE on each distributional parameter of the read response (its
#       intercept if it has none); warn on the largest move, scaled, >= 0.3.
# Usage: as the rule script, table1 <family> | table2 <slope> | vm | cost
args0 <- commandArgs(trailingOnly = TRUE)
source_rule <- function() {
  txt <- readLines("dev/review-refuse-f2-rule.R")
  i <- which(startsWith(txt, "mode <- args[1]"))
  eval(parse(text = sub('args <- commandArgs(trailingOnly = TRUE)', 'args <- "none"', txt[seq_len(i - 1L)], fixed = TRUE)), envir = globalenv())
}
source_rule()
r3 <- function(a, rows = 10L, tau = 0.3) {
  p <- a$opt$par
  t0 <- proc.time()[["elapsed"]]
  T2 <- tk_T(a, p)
  tg <- proc.time()[["elapsed"]] - t0
  if (!is.finite(T2)) T2 <- Inf
  mv <- NA_real_; ts <- NA_real_; nsw <- 0L
  if (abs(T2) > 1) {
    t0 <- proc.time()[["elapsed"]]
    H <- stats::optimHess(p, a$obj$fn, a$obj$gr)
    se <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(p)))
    se[!is.finite(se)] <- 0.1
    th <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    n_read <- length(a$estimates$miss)
    at <- if (n_read > rows) unique(round(seq(1, n_read, length.out = rows))) else seq_len(n_read)
    pn <- names(p); bn <- names(a$estimates$beta); dn <- names(a$estimates$betad)
    ibeta <- which(pn == "beta"); idp <- which(pn == "betad")
    own <- idp[grep("^y_", dn)]
    if (!length(own)) own <- ibeta[grep("^y_.Intercept", bn)]
    Tb <- int_T(a, p, at, n_read)
    mv <- max(vapply(own, function(j) {
      q <- p; q[j] <- q[j] + se[j]
      abs(int_T(a, q, at, n_read) - Tb)
    }, 0))
    nsw <- length(own)
    ts <- proc.time()[["elapsed"]] - t0
    a$obj$fn(p)
  }
  list(T2 = T2, mv = mv, warn = isTRUE(abs(T2) > 1 && mv >= tau), tg = tg, ts = ts, nsw = nsw)
}
out <- function(tag, a) {
  r <- r3(a)
  pr("R3 %s T2 %.3f mv %s warn %s tg %.3f ts %s T0 %.3f", tag, r$T2, format(r$mv, digits = 3), r$warn, r$tg, format(r$ts, digits = 3), a$mi_laplace$total %||% NA)
}
m <- args0[1]
if (m == "table1") {
  fam <- args0[2]; fi <- match(fam, fams); sp <- speco(fam)
  ff <- switch(fam, beta = Beta(), Gamma = stats::Gamma(link = "log"), gaussian = gaussian(), get(fam, asNamespace("frmtmb"))())
  set.seed(11200 + fi)
  for (cfg in list(c(200, 60, 0.7, 0.7), c(200, 60, 0.05, 1), c(200, 100, 0.7, 0.7), c(200, 15, 0.7, 0.7))) {
    n <- cfg[1]; mm <- cfg[2]; bz <- cfg[3]; sz <- cfg[4]
    x <- rnorm(n); eta <- if (fam == "beta") 0.3 + 0.7 * x else 0.3 + 0.5 * x
    y <- sp$rgen(eta)
    d <- data.frame(x = x, y = y); d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, sz)
    d$y[sample(n, mm)] <- NA
    a <- suppressWarnings(frm(bf(y | mi() ~ x) + ff + bf(z ~ mi(y) + x) + gaussian(), data = d))
    out(sprintf("t1 %s read %d slope %.2f", fam, mm, bz), a)
  }
} else if (m == "table2") {
  sl <- as.numeric(args0[2])
  set.seed(11501 + round(sl * 100))
  for (mm in c(20, 60)) for (k in 1:8) {
    n <- 200; x <- rnorm(n)
    y <- 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2
    d <- data.frame(x = x, y = y); d$z <- 0.5 + sl * y + 0.3 * x + rnorm(n, 0, 1)
    d$y[sample(n, mm)] <- NA
    a <- suppressWarnings(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y) + x) + gaussian(), data = d))
    out(sprintf("t2 exgaussian read %d slope %.2f rep %d", mm, sl, k), a)
  }
} else if (m == "vm") {
  set.seed(11600)
  rvm <- function(mu, k) vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
  sp <- speco("von_mises")
  for (cfg in list(c(0, 3), c(3.05, 2), c(3.05, 3), c(3.14, 20))) for (mm in c(15, 60)) for (bz in c(0.2, 0.7)) for (k in 1:3) {
    n <- 200; x <- rnorm(n)
    y <- rvm(cfg[1] + 0.05 * x, cfg[2])
    d <- data.frame(x = x, y = y); d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, 0.7)
    rows <- sample(n, mm); d$y[rows] <- NA
    a <- suppressWarnings(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(), data = d))
    out(sprintf("vm mean %.2f kappa %g read %d slope %.1f rep %d", cfg[1], cfg[2], mm, bz, k), a)
    # is frm's optimum the Laplace objective's, or did it stop short? exact's
    # point on frm's own objective, parameter by parameter
    if (cfg[1] > 3) {
      p0 <- pars0(a); ex <- tr(exact(sp, d, p0))
      if (!is_err(ex)) {
        q <- a$opt$par; nm <- names(q)
        bn <- names(a$estimates$beta); dn <- names(a$estimates$betad)
        q[nm == "beta"][grep("^y_", bn)] <- ex$par[1:2]
        q[nm == "betad"][grep("^y_", dn)] <- ex$par[3]
        zb <- grep("^z_", bn); o <- c(zb[1], grep("miy", bn), grep("z_x", bn))
        q[nm == "beta"][o] <- ex$par[4:6]
        q[nm == "betad"][grep("^z_", dn)] <- ex$par[7]
        lap <- tr(a$obj$fn(q))
        pr("   frm objective %.2f at its optimum, %.2f at the exact point (lower is better for frm); exact logLik %.2f; frm logLik %.2f; y0 frm %.3f exact %.3f; zy frm %.3f exact %.3f; kappa frm %.3f exact %.3f",
           a$opt$objective, as.numeric(lap), ex$ll, as.numeric(logLik(a)), p0[1], ex$par[1], p0[5], ex$par[5], exp(p0[3]), exp(ex$par[3]))
        # the same objective from the exact point and from the observed rows
        o1 <- tr(nlminb(q, a$obj$fn, a$obj$gr, control = list(iter.max = 500, eval.max = 800)))
        cc <- tr(suppressWarnings(frm(bf(y ~ x), family = von_mises(), data = d[!is.na(d$y), ])))
        o2 <- if (is_err(cc)) cc else {
          q2 <- a$opt$par
          q2[nm == "beta"][grep("^y_", bn)] <- cc$estimates$beta
          q2[nm == "betad"][grep("^y_", dn)] <- cc$estimates$betad
          tr(nlminb(q2, a$obj$fn, a$obj$gr, control = list(iter.max = 500, eval.max = 800)))
        }
        pr("   local optimum: frm %.2f; nlminb from the exact point %s; from the observed-rows fit %s", a$opt$objective,
           if (is_err(o1)) msg(o1, 40) else format(o1$objective, digits = 6), if (is_err(o2)) msg(o2, 40) else format(o2$objective, digits = 6))
        if (!is_err(o2)) {
          z <- (c(o2$par) - 0)
          pos <- c(which(nm == "beta")[grep("^y_", bn)], which(nm == "betad")[grep("^y_", dn)], which(nm == "beta")[o], which(nm == "betad")[grep("^z_", dn)])
          pr("   from the observed-rows fit: largest shift from exact %.2f SE; objective - (-exact logLik) %.2f", max(abs((z[pos] - ex$par) / ex$se)), o2$objective + ex$ll)
        }
        a$obj$fn(a$opt$par)
      }
    }
  }
} else if (m == "cost") {
  orig <- get("mi_laplace_note", asNamespace("frmtmb"))
  set.seed(11502)
  for (cfg in list(c(200, 60), c(200, 120), c(1000, 300), c(5000, 1500), c(20000, 6000))) {
    n <- cfg[1]; mm <- cfg[2]; x <- rnorm(n)
    y <- rgamma(n, 1.3, 1.3 / exp(0.3 + 0.5 * x))
    d <- data.frame(x = x, y = y); d$z <- 0.5 + 0.7 * y + 0.3 * x + rnorm(n, 0, 0.7)
    d$y[sample(n, mm)] <- NA
    fm <- function() frm(bf(y | mi() ~ x) + stats::Gamma(link = "log") + bf(z ~ mi(y) + x) + gaussian(), data = d)
    assignInNamespace("mi_laplace_note", function(fit) fit, "frmtmb")
    tb <- vapply(1:3, function(i) { t0 <- proc.time()[["elapsed"]]; a0 <<- suppressWarnings(fm()); proc.time()[["elapsed"]] - t0 }, 0)
    assignInNamespace("mi_laplace_note", orig, "frmtmb")
    bare <- median(tb)
    t0s <- median(vapply(1:3, function(i) { t0 <- proc.time()[["elapsed"]]; suppressWarnings(orig(a0)); proc.time()[["elapsed"]] - t0 }, 0))
    tgs <- median(vapply(1:7, function(i) { t0 <- proc.time()[["elapsed"]]; for (k in 1:5) tk_T(a0, a0$opt$par); (proc.time()[["elapsed"]] - t0) / 5 }, 0))
    ss <- median(vapply(1:3, function(i) { r <- r3(a0); if (is.na(r$ts)) 0 else r$ts }, 0))
    pr("COST3 n=%d read=%d: bare %.2fs; current rule %.2fs = %.0f%% of bare; R3 gate %.4fs = %.2f%%; R3 second stage %.2fs = %.0f%% (only past the gate)", n, mm, bare, t0s, 100 * t0s / bare, tgs, 100 * tgs / bare, ss, 100 * ss / bare)
  }
}
