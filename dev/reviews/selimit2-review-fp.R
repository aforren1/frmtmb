# Review of lane selimit2, item 4800: does the joint limit verdict fire on
# fits whose maximum likelihood estimate is finite? Each case prints the
# verdict (fit$cache$se_limit), the warnings, and a profile along every
# joint direction found: the objective at the stop moved t natural units
# along it (the rest held), out to 1e6, and the same with the other
# parameters re-optimized at each t. A finite maximum shows as a rise far
# out.
#   Rscript selimit2-review-fp.R <base|lane> [case ...]
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
which_cases <- a[-1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")

fitw <- function(...) {
  ws <- character(0)
  ms <- character(0)
  fit <- withCallingHandlers(frm(...), warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(fit = fit, ws = ws, ms = ms)
}

report <- function(label, r, profile = TRUE) {
  fit <- r$fit
  p <- fit$opt$par
  cat("\n==", label, "| code", fit$opt$convergence, "| obj",
      format(fit$opt$objective, digits = 12), "\n")
  cat("   par:", paste(names(p), format(p, digits = 6), collapse = "; "), "\n")
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(fit)))),
                 error = function(e) NA)
  cat("   se:", paste(format(se, digits = 4), collapse = " "), "\n")
  for (w in r$ws) cat("   W:", substr(w, 1, 260), "\n")
  for (m in r$ms) cat("   M:", substr(m, 1, 200), "\n")
  li <- fit$cache$se_limit
  if (!is.null(li)) {
    cat("   limit kinds:", paste(names(li$kind), li$kind, collapse = "; "),
        "| joint:", paste(li$joint, collapse = ","), "| sep:",
        paste(li$sep, collapse = ","), "| sep1:",
        paste(li$sep1, collapse = ","), "\n")
  }
  if (!profile) return(invisible(fit))
  jw <- li$jw
  if (is.null(jw)) return(invisible(fit))
  obj <- fit$obj
  f0 <- obj$fn(p)
  u <- frmtmb:::se_limit_units(fit)
  for (k in seq_len(NCOL(jw))) {
    w <- jw[, k]
    w <- w / sqrt(sum((w / u)^2))
    cat("   joint dir", k, ":", paste(names(p)[w != 0],
                                      format(w[w != 0], digits = 4),
                                      collapse = " "), "\n")
    ts <- c(-1e3, -100, -10, -2, 2, 12, 100, 1e3, 1e4, 1e5, 1e6)
    held <- vapply(ts, function(t) obj$fn(p + t * w) - f0, 0)
    cat("   held   t:", paste(sprintf("%g:%.3g", ts, held), collapse = " "),
        "\n")
    # re-optimized: fix the joint coefficients at p + t w, free the rest
    jj <- which(w != 0)
    rest <- setdiff(seq_along(p), jj)
    prof <- vapply(ts, function(t) {
      q <- p + t * w
      if (!length(rest)) return(obj$fn(q) - f0)
      o <- tryCatch(stats::nlminb(q[rest], function(x) {
        q2 <- q
        q2[rest] <- x
        v <- obj$fn(q2)
        if (!is.finite(v)) 1e300 else v
      }, function(x) {
        q2 <- q
        q2[rest] <- x
        drop(obj$gr(q2))[rest]
      }, control = list(iter.max = 300, eval.max = 600)),
      error = function(e) list(objective = NA))
      o$objective - f0
    }, 0)
    cat("   reopt  t:", paste(sprintf("%g:%.3g", ts, prof), collapse = " "),
        "\n")
  }
  invisible(fit)
}

cases <- list()

# B1: bernoulli separated at x = 10 but for one reversed pair a distance
# 2 * delta apart: a finite MLE, flat along (intercept, slope) keeping the
# boundary at 10
cases$bern <- function() {
  for (delta in c(1e-4, 1e-5, 1e-6)) {
    x <- c(5:9, 11:15, 10 - delta, 10 + delta)
    y <- c(rep(0, 5), rep(1, 5), 1, 0)
    d <- data.frame(x = rep(x, 3), y = rep(y, 3))
    report(sprintf("bernoulli reversed pair delta %g (finite MLE)", delta),
           fitw(y ~ x, data = d, family = bernoulli()))
  }
  # the same with the pair in order: quasi-separation, MLE infinite
  x <- c(5:9, 11:15, 10 - 1e-5, 10 + 1e-5)
  y <- c(rep(0, 5), rep(1, 5), 0, 1)
  d <- data.frame(x = rep(x, 3), y = rep(y, 3))
  report("bernoulli pair in order (MLE infinite)",
         fitw(y ~ x, data = d, family = bernoulli()))
}

# P1: poisson y ~ x, a count of 1 at the largest x and one at x - delta:
# finite MLE (a single count at the largest x is infinite)
cases$pois <- function() {
  for (delta in c(1e-3, 1e-4, 1e-5)) {
    x <- c(seq(0, 10, length.out = 30), 10 + 5 - delta, 15)
    y <- c(rep(0, 30), 1, 1)
    d <- data.frame(x = x, y = y)
    report(sprintf("poisson top pair delta %g (finite MLE)", delta),
           fitw(y ~ x, data = d, family = poisson()))
  }
  x <- c(seq(0, 10, length.out = 30), 15)
  d <- data.frame(x = x, y = c(rep(0, 30), 1))
  report("poisson one count at the top (MLE infinite)",
         fitw(y ~ x, data = d, family = poisson()))
}

# Z1: zero-inflated poisson with zi ~ f, level a with one nonzero count
# among 40: zi_a near 1 but finite; and all zero (infinite)
cases$zi <- function() {
  set.seed(3)
  for (k in c(0L, 1L, 2L)) {
    ya <- c(rep(0, 40 - k), rep(3, k))
    d <- data.frame(f = factor(rep(c("a", "b"), each = 40)),
                    y = c(ya, stats::rpois(40, 3)))
    report(sprintf("zip zi ~ f, level a %d nonzero of 40", k),
           fitw(bf(y ~ 1, zi ~ f), data = d, family = zero_inflated_poisson()))
    report(sprintf("zip zi ~ 0 + f, level a %d nonzero of 40", k),
           fitw(bf(y ~ 1, zi ~ 0 + f), data = d,
                family = zero_inflated_poisson()))
    report(sprintf("hurdle hu ~ f, level a %d nonzero of 40", k),
           fitw(bf(y ~ 1, hu ~ f), data = d, family = hurdle_poisson()))
  }
  # hu toward 0 on a level: no zeros in level b
  d <- data.frame(f = factor(rep(c("a", "b"), each = 40)),
                  y = c(stats::rpois(40, 1), 1 + stats::rpois(40, 3)))
  report("hurdle hu ~ f, level b has no zeros (hu_b to 0)",
         fitw(bf(y ~ 1, hu ~ f), data = d, family = hurdle_poisson()))
  # zi with a continuous predictor, near-one zi at one end
  set.seed(5)
  x <- seq(-3, 3, length.out = 120)
  zi <- stats::plogis(4 * x)
  y <- ifelse(stats::runif(120) < zi, 0, stats::rpois(120, 4))
  d <- data.frame(x = x, y = y)
  report("zip zi ~ x, zi = plogis(4 x) (finite)",
         fitw(bf(y ~ 1, zi ~ x), data = d, family = zero_inflated_poisson()))
}

# M1: two-gaussian mixture with a small second component, theta2 ~ z:
# a finite but weakly determined mixing logit
cases$mix <- function() {
  set.seed(11)
  n <- 400
  z <- stats::rnorm(n)
  cls <- stats::runif(n) < stats::plogis(-5 + 0.5 * z)
  y <- ifelse(cls, stats::rnorm(n, 6, 0.5), stats::rnorm(n, 0, 1))
  d <- data.frame(y = y, z = z)
  cat("small-component rows:", sum(cls), "\n")
  report("gaussian mixture, small component, theta2 ~ z",
         fitw(bf(y ~ 1, theta2 ~ z), data = d,
              family = mixture(gaussian(), gaussian())))
}

# O1: test-ordinal-mixture.R #28 and #29 (newly NaN): theta1 ~ z
cases$omx <- function() {
  src <- readLines(file.path(
    "C:/Users/adf44/source/r/frmtmb-wt-selimit2/tests/testthat",
    "test-ordinal-mixture.R"))
  i <- grep("^omx_fa_data <- function", src)
  j <- i + which(src[(i + 1):length(src)] == "}")[1]
  eval(parse(text = src[i:j]))
  d <- omx_fa_data(3)
  for (fam in list(mixture(cumulative(), sratio(), order = "mu"),
                   mixture(cumulative(), cumulative("probit"),
                           order = "mu"))) {
    report(paste("ordinal mixture", fam$family, "ysh ~ 1, theta1 ~ z"),
           fitw(bf(ysh ~ 1, theta1 ~ z), family = fam, data = d))
  }
}

# O2: is the ordinal mixtures' end the supremum? Refit from the generating
# mixing logit (0.3, 1.5) and compare objectives; then a profile over a
# scale s of the theta coefficients, s * (a, b), the thresholds re-optimized
cases$omx2 <- function() {
  src <- readLines(file.path(
    "C:/Users/adf44/source/r/frmtmb-wt-selimit2/tests/testthat",
    "test-ordinal-mixture.R"))
  i <- grep("^omx_fa_data <- function", src)
  j <- i + which(src[(i + 1):length(src)] == "}")[1]
  eval(parse(text = src[i:j]))
  d <- omx_fa_data(3)
  fam <- mixture(cumulative(), sratio(), order = "mu")
  r0 <- fitw(bf(ysh ~ 1, theta1 ~ z), family = fam, data = d)
  fit <- r0$fit
  p <- fit$opt$par
  obj <- fit$obj
  cat("default fit obj", format(fit$opt$objective, digits = 12), "\n")
  it <- which(names(p) == "betad")
  rest <- setdiff(seq_along(p), it)
  for (ab in list(c(0.3, 1.5), c(1, 1), c(5, 5), c(20, 18), c(100, 90),
                  c(1000, 920), c(1e4, 9200), c(1e5, 9.2e4))) {
    q <- p
    q[it] <- ab
    o <- stats::nlminb(q[rest], function(x) {
      q2 <- q
      q2[rest] <- x
      v <- obj$fn(q2)
      if (!is.finite(v)) 1e300 else v
    })
    cat(sprintf("theta (%g, %g), thresholds re-optimized: obj %.10f (vs end %.3g)\n",
                ab[1], ab[2], o$objective, o$objective - fit$opt$objective))
  }
  # free optimization from the generating values
  q <- p
  q[it] <- c(0.3, 1.5)
  o <- stats::nlminb(q, obj$fn, obj$gr,
                     control = list(iter.max = 2000, eval.max = 4000))
  cat("free nlminb from (0.3, 1.5): obj", format(o$objective, digits = 12),
      "theta", format(o$par[it], digits = 6), "code", o$convergence, "\n")
}

# B2: the reversed pair, delta 1e-8 (MLE at slope about 34.6), started
# beyond the stop the default tolerance reaches
cases$bernstart <- function() {
  delta <- 1e-8
  x <- c(5:9, 11:15, 10 - delta, 10 + delta)
  y <- c(rep(0, 5), rep(1, 5), 1, 0)
  d <- data.frame(x = rep(x, 3), y = rep(y, 3))
  for (b in c(25, 30)) {
    report(sprintf("bernoulli reversed pair delta 1e-8, start slope %g", b),
           fitw(y ~ x, data = d, family = bernoulli(),
                start = list(beta = c(-10 * b, b))))
  }
  # the finite maximum: the objective along the slope with the boundary at
  # 10
  f <- function(b) {
    eta <- b * (d$x - 10)
    -sum(stats::dbinom(d$y, 1, stats::plogis(eta), log = TRUE))
  }
  bs <- c(20, 25, 30, 33, 34, 35, 36, 40, 50, 1e3, 1e5, 1e7, 1e8)
  cat("   nll(b), boundary 10:", paste(sprintf("%g:%.12f", bs,
                                               vapply(bs, f, 0)),
                                       collapse = " "), "\n")
}

for (nm in if (length(which_cases)) which_cases else names(cases)) {
  cases[[nm]]()
}
