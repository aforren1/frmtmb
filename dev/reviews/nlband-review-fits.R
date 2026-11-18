# Review of lane nlband, checks 1 to 3 and 5 (fits). Each ROW: the default
# fit (logLik, code, estimates, warning and message classes), the runs
# made inside nl_pole_escape() (ESC), the objective and gradient
# evaluations of every optimize_obj() run (EV), and the global maximum
# (GLOB): an exact profile over b for gaussian bodies with one pole (a
# grid of step 2e-4 over b, then optimize() between the grid's
# neighbors), else the best of local fits (nl_starts = FALSE) from a set
# of starts. NLBAND_ARM=before|after Rscript nlband-review-fits.R <sections>
args <- commandArgs(trailingOnly = TRUE)
secs <- if (length(args)) strsplit(args[1], ",")[[1]] else
  c("seeds", "edge", "inside", "dense", "two", "visits", "jitter",
    "student", "poisson", "bstep", "bband")
source("dev/nlband-common.R")
n_ev <- 0; n_esc <- 0
in_esc <- function() {
  any(vapply(sys.calls(), function(cl) {
    is.call(cl) && identical(cl[[1L]], as.name("nl_pole_escape"))
  }, NA))
}
suppressMessages(trace("optimize_obj", exit = quote({
  n_ev <<- n_ev + (returnValue()$evals %||% 0)
  if (in_esc()) n_esc <<- n_esc + 1
}), print = FALSE, where = ns))
`%||%` <- function(a, b) if (is.null(a)) b else a
fit_row <- function(tag, expr, glob = NULL) {
  n_ev <<- 0; n_esc <<- 0
  t0 <- proc.time()[["elapsed"]]
  r <- tt(expr)
  el <- proc.time()[["elapsed"]] - t0
  if (iserr(r)) {
    cat("ROW", tag, "| ERROR", substr(r$msg, 1, 160), "\n")
    return(invisible(NULL))
  }
  ll <- llv(r)
  g <- if (is.null(glob)) c(NA, NA) else glob
  cat(sprintf(paste0("ROW %s | ll %.6f code %d beta %s | warn %s msg %s | ",
                     "ESC %d EV %d t %.2f | GLOB %.6f at %.4f short %.3f\n"),
              tag, ll, r$opt$convergence,
              paste(sprintf("%.5g", r$estimates$beta), collapse = " "),
              wsum(r), msum(r), n_esc, n_ev, el, g[1], g[2], g[1] - ll))
  for (m in attr(r, "m")) if (grepl("pole", m)) cat("   MSG", m, "\n")
  for (w in attr(r, "w")) if (grepl("pole", w)) cat("   WARN", w, "\n")
  invisible(r)
}
# exact profile log-likelihood of y ~ a / (x + b) + N(0, s), a and s out
prof1 <- function(x, y, b) {
  vapply(b, function(bb) {
    z <- 1 / (x + bb)
    a <- sum(y * z) / sum(z * z)
    rss <- sum((y - a * z)^2)
    n <- length(y)
    -n / 2 * (log(2 * pi * rss / n) + 1)
  }, 1)
}
glob1 <- function(x, y) {
  # 15 points inside every gap between adjacent walls -x_i, and 60
  # log-spaced points beyond each end; then optimize() between the best
  # point's neighbors
  w <- sort(unique(-x))
  r <- diff(range(w))
  u <- seq(0, 1, length.out = 17)[2:16]
  b <- c(as.vector(outer(u, diff(w)) + rep(w[-length(w)], each = 15)),
         w[1] - r * 10^seq(-4, 2, length.out = 60),
         w[length(w)] + r * 10^seq(-4, 2, length.out = 60))
  b <- sort(b)
  p <- prof1(x, y, b)
  p[!is.finite(p)] <- -Inf
  k <- which.max(p)
  o <- optimize(function(bb) prof1(x, y, bb), c(b[max(1, k - 1)],
                b[min(length(b), k + 1)]), maximum = TRUE, tol = 1e-12)
  c(max(o$objective, p[k]), if (o$objective >= p[k]) o$maximum else b[k])
}
# the best local fit over starts (nl_starts = FALSE)
globms <- function(form, d, fam, starts) {
  best <- c(-Inf, NA)
  for (s in starts) {
    r <- tryCatch(suppressWarnings(suppressMessages(
      frm(form, data = d, family = fam, start = list(beta = s),
          control = frmtmb_control(nl_starts = FALSE)))),
      error = function(e) NULL)
    if (is.null(r)) next
    v <- llv(r)
    if (is.finite(v) && v > best[1]) best <- c(v, r$estimates$beta[2])
  }
  best
}
band_x <- function(seed) band_data(seed)$x
hole_starts <- function(a = c(-1, 1)) {
  st <- list()
  for (aa in a) for (b0 in c(-2.25, -2.1, -2, -1.9, -1.75, 0.5, -4.5)) {
    st[[length(st) + 1L]] <- c(aa, b0)
  }
  st
}

if ("seeds" %in% secs) {
  # check 3: the band data, seeds 1 to 3 and 5 fresh seeds 201 to 205
  for (s in c(1:3, 201:205)) {
    d <- band_data(s)
    fit_row(paste("seeds", s), frm(band_form, data = d), glob1(d$x, d$y))
  }
}
if ("edge" %in% secs) {
  # a hole near the edge: k rows below it, pole 0.25, no rows in 0.1-0.4
  for (k in c(6, 12)) for (s in 1:3) {
    set.seed(s)
    x <- c(runif(k, 0, 0.1), runif(400 - k, 0.4, 4))
    d <- data.frame(x = x, y = 1 / (x - 0.25) + rnorm(400, 0, 0.2))
    fit_row(sprintf("edge k%d %d", k, s), frm(band_form, data = d),
            glob1(d$x, d$y))
  }
}
if ("inside" %in% secs) {
  # the band data with k rows put back inside the hole (pole at 2)
  for (k in c(1, 3, 5)) for (s in 1:3) {
    d <- band_data(s)
    set.seed(1000 + s)
    xi <- runif(k, 1.7, 2.3)
    d <- rbind(d, data.frame(x = xi, y = 1 / (xi - 2) + rnorm(k, 0, 0.2)))
    fit_row(sprintf("inside k%d %d", k, s), frm(band_form, data = d),
            glob1(d$x, d$y))
  }
}
if ("dense" %in% secs) {
  # the pole inside dense data: no hole
  for (s in 1:3) {
    set.seed(s)
    x <- runif(400, 0, 4)
    d <- data.frame(x = x, y = 1 / (x - 2) + rnorm(400, 0, 0.2))
    fit_row(paste("dense", s), frm(band_form, data = d), glob1(d$x, d$y))
  }
}
if ("two" %in% secs) {
  # two poles: a / ((x + b) * (x + c)), holes at 1 and 3
  prof2 <- function(x, y, b, c) {
    z <- 1 / ((x + b) * (x + c))
    a <- sum(y * z) / sum(z * z)
    rss <- sum((y - a * z)^2)
    -length(y) / 2 * (log(2 * pi * rss / length(y)) + 1)
  }
  for (s in 1:3) {
    set.seed(s)
    x <- runif(400, 0, 4)
    x <- x[abs(x - 1) > 0.25 & abs(x - 3) > 0.25]
    d <- data.frame(x = x, y = 1 / ((x - 1) * (x - 3)) +
                      rnorm(length(x), 0, 0.2))
    o <- optim(c(-1, -3), function(p) -prof2(d$x, d$y, p[1], p[2]),
               control = list(reltol = 1e-14))
    fit_row(paste("two", s),
            frm(bf(y ~ a / ((x + b) * (x + c)), a + b + c ~ 1, nl = TRUE),
                data = d), c(-o$value, o$par[1]))
  }
}
if ("visits" %in% secs) {
  # 10 visits with jitter, the pole between visits 5 and 6
  for (s in 1:5) {
    set.seed(s)
    x <- rep(1:10, each = 40) + runif(400, -0.2, 0.2)
    d <- data.frame(x = x, y = 1 / (x - 5.5) + rnorm(400, 0, 0.2))
    fit_row(paste("visits", s), frm(band_form, data = d), glob1(d$x, d$y))
  }
}
if ("jitter" %in% secs) {
  # replicates at 60 levels measured with a small error; no pole in the
  # data (truth 5 / (x + 0.5)): every hole start is a false start
  for (s in 1:3) {
    set.seed(s)
    x <- rep(1:60, each = 3) + runif(180, -0.01, 0.01)
    d <- data.frame(x = x, y = 5 / (x + 0.5) + rnorm(180, 0, 0.2))
    fit_row(paste("jitter", s), frm(band_form, data = d), glob1(d$x, d$y))
  }
  # a 30-level two-fold dilution series, 4 replicates, truth b 0.5
  for (s in 1:3) {
    set.seed(s)
    x <- rep(0.001 * 2^(0:29), each = 4)
    d <- data.frame(x = x, y = 5 / (x + 0.5) + rnorm(120, 0, 0.2))
    fit_row(paste("dil2fold", s), frm(band_form, data = d), glob1(d$x, d$y))
  }
}
if ("student" %in% secs) {
  for (s in 1:3) {
    set.seed(s)
    x <- band_x(s)
    d <- data.frame(x = x, y = 1 / (x - 2) + 0.2 * rt(length(x), 3))
    fit_row(paste("student", s), frm(band_form, data = d, family = student()),
            globms(band_form, d, student(), hole_starts()))
  }
}
if ("poisson" %in% secs) {
  for (s in 1:3) {
    set.seed(s)
    x <- band_x(s)
    d <- data.frame(x = x, y = rpois(length(x), exp(1 / (x - 2))))
    fit_row(paste("poisson", s), frm(band_form, data = d, family = poisson()),
            globms(band_form, d, poisson(), hole_starts()))
  }
}
if ("bstep" %in% secs) {
  # bernoulli with a step at 2 and no pole: logit p = 4 (x - 2)
  for (s in 1:3) {
    set.seed(s)
    x <- band_x(s)
    d <- data.frame(x = x, y = rbinom(length(x), 1, plogis(4 * (x - 2))))
    r <- fit_row(paste("bstep", s),
                 frm(band_form, data = d, family = bernoulli()),
                 globms(band_form, d, bernoulli(), hole_starts()))
    # the profile over b inside the hole, a out (bernoulli, logit)
    pb <- function(bb) {
      z <- 1 / (d$x + bb)
      optimize(function(a) sum(dbinom(d$y, 1, plogis(a * z), log = TRUE)),
               c(-200, 200), maximum = TRUE)$objective
    }
    lo <- -min(d$x[d$x > 2]); hi <- -max(d$x[d$x < 2])
    bg <- seq(lo, hi, length.out = 401)[-c(1, 401)]
    pv <- vapply(bg, pb, 1)
    k <- which.max(pv)
    cat(sprintf(paste0("   PROF hole %.4f to %.4f: max %.6f at b %.4f ",
                       "(%d of 399), ends %.6f %.6f\n"),
                lo, hi, pv[k], bg[k], k, pv[1], pv[399]))
  }
}
if ("bband" %in% secs) {
  # the lane's bernoulli band rows (logit p = 2 / (x - 2)): is the kept
  # fit an interior maximum of the profile over b in the hole?
  for (s in 1:3) {
    d <- band_data(s)
    d$nt <- 10L
    d$yb <- rbinom(nrow(d), 10, plogis(1 / (d$x - 2)))
    d$y1 <- rbinom(nrow(d), 1, plogis(2 / (d$x - 2)))
    d$y <- d$y1
    fit_row(paste("bband", s), frm(band_form, data = d, family = bernoulli()),
            globms(band_form, d, bernoulli(), hole_starts()))
    pb <- function(bb) {
      z <- 1 / (d$x + bb)
      optimize(function(a) sum(dbinom(d$y, 1, plogis(a * z), log = TRUE)),
               c(-200, 200), maximum = TRUE)$objective
    }
    lo <- -min(d$x[d$x > 2]); hi <- -max(d$x[d$x < 2])
    bg <- seq(lo, hi, length.out = 401)[-c(1, 401)]
    pv <- vapply(bg, pb, 1)
    k <- which.max(pv)
    cat(sprintf(paste0("   PROF hole %.4f to %.4f: max %.6f at b %.4f ",
                       "(%d of 399), ends %.6f %.6f\n"),
                lo, hi, pv[k], bg[k], k, pv[1], pv[399]))
  }
}
cat("DONE\n")
