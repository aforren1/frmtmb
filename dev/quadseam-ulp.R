# The three new assertions of test-quadrature-defects.R over 20 data
# sets moved by a few ulps (v * (1 + 2e-15 * k), k = 1..20): the
# continuous response (skew normal, mvbf) or the covariate (bernoulli).
#   Rscript dev/quadseam-ulp.R <lib|base> <escape|runaway|boundary|atmax>
# log dev/quadseam-log/ulp-<tag>-<case>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
case <- a[2]
quiet <- function(expr, cls = character()) {
  n <- stats::setNames(integer(length(cls)), cls)
  v <- withCallingHandlers(expr, warning = function(w) {
    for (k in cls) if (inherits(w, k)) n[[k]] <<- n[[k]] + 1L
    invokeRestart("muffleWarning")
  }, message = function(m) {
    for (k in cls) if (inherits(m, k)) n[[k]] <<- n[[k]] + 1L
    invokeRestart("muffleMessage")
  })
  attr(v, "n") <- n
  v
}
ok <- 0L
t0 <- proc.time()[["elapsed"]]
for (k in 1:20) {
  e <- 1 + 2e-15 * k
  if (case == "escape") {
    set.seed(1)
    n <- 200
    xs <- -abs(rnorm(n)) * 3
    y <- (xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)) * e
    dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
    m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
    bad <- list(betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3))
    fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
    f0 <- quiet(frm(fo, family = skew_normal(), data = dd, quadrature = TRUE))
    f1 <- quiet(frm(fo, family = skew_normal(), data = dd, quadrature = TRUE,
                    start = bad))
    gain <- f1$opt[["stationary_escape"]][["gain"]]
    r <- abs(as.numeric(logLik(f1)) - as.numeric(logLik(f0))) / gain
    pass <- isTRUE(gain > 1) && r < 1e-6
    cat(sprintf("K %2d gain %.4f ratio %.3e %s\n", k, gain, r,
                if (pass) "PASS" else "FAIL"))
  } else if (case == "runaway") {
    # fix round: the mi() runaway of test-quadrature-defects.R (the
    # two-response case of the first round is refused by lane wcheck)
    set.seed(20261111)
    n <- 80
    d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
    d$y <- (1 + 0.5 * d$x + rnorm(n)) * e
    d$y[c(5, 17, 33)] <- NA
    d$w <- 1.5
    f <- quiet(frm(bf(y | weights(w) + mi() ~ x + (1 | obs)),
                   family = gaussian(), data = d, quadrature = TRUE),
               "frmtmb_quad_unsettled")
    u <- attr(f, "n")[["frmtmb_quad_unsettled"]]
    pass <- identical(u, 1L)
    cat(sprintf("K %2d unsettled %d log sigma %.3f %s\n", k, u,
                f$estimates$betad[1], if (pass) "PASS" else "FAIL"))
  } else if (case == "atmax") {
    # fix round: the three fits at their maximum, no warning of any class
    set.seed(1)
    d1 <- data.frame(g = factor(rep(1:10, each = 5)), x = rnorm(50))
    u <- rnorm(10, 0, 0.4)
    d1$k <- rbinom(50, 500, plogis(0.3 + 0.5 * d1$x + u[d1$g]))
    d1$nk <- 500 - d1$k
    d1$x <- d1$x * e
    set.seed(4)
    d2 <- data.frame(g = factor(rep(1:6, each = 4)), x = rnorm(24))
    u <- rnorm(6, 0, 0.5)
    d2$z <- c(-3, -2, -1, 1, 2, 3)[d2$g]
    d2$k <- rbinom(24, 100, plogis(1.5 * d2$z + 0.2 * d2$x + u[d2$g]))
    d2$nk <- 100 - d2$k
    d2$x <- d2$x * e
    set.seed(11)
    d3 <- data.frame(g = factor(rep(1:40, each = 5)), x = rnorm(200))
    d3$y <- rbinom(200, 1, plogis(-0.3 + 0.7 * d3$x + rnorm(40)[d3$g]))
    d3$x <- d3$x * e
    nw <- 0L
    cnt <- function(expr) withCallingHandlers(expr, warning = function(w) {
      nw <<- nw + 1L
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
    ctl <- frmtmb_control(grad_tol = 1e-5)
    cnt(frm(cbind(k, nk) ~ x + (1 | g), family = binomial(), data = d1,
            quadrature = TRUE, control = ctl))
    cnt(frm(cbind(k, nk) ~ z + x + (1 | g), family = binomial(), data = d2,
            quadrature = TRUE, control = ctl))
    cnt(frm(y ~ x + (1 | g), family = bernoulli(), data = d3,
            quadrature = TRUE))
    pass <- nw == 0L
    cat(sprintf("K %2d warnings %d %s\n", k, nw, if (pass) "PASS" else "FAIL"))
  } else {
    set.seed(3)
    d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
    d$y <- rbinom(120, 1, plogis(0.2 + 0.5 * d$x))
    d$x <- d$x * e
    f <- quiet(frm(y ~ x + (1 | g), family = bernoulli(), data = d,
                   quadrature = TRUE), "frmtmb_boundary_fit")
    sdr <- frmtmb:::sdr_of(f)
    v <- diag(sdr$cov.fixed)
    pass <- attr(f, "n")[["frmtmb_boundary_fit"]] >= 1L &&
      identical(unname(sdr$se_lost), "boundary") && is.nan(v[[3L]]) &&
      all(is.finite(v[1:2]))
    cat(sprintf("K %2d theta %.3f %s\n", k, f$estimates$theta[1],
                if (pass) "PASS" else "FAIL"))
  }
  ok <- ok + pass
}
cat(sprintf("ULP %s pass %d of 20 secs %.1f\n", case, ok,
            proc.time()[["elapsed"]] - t0))
