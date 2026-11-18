# Review of lane nlband: checks 4 (bitwise, guard absent), 5 (cost at n
# 1e4), 6 (5099: other spellings of a pole) and the code-1 kept fit.
# NLBAND_ARM=before|after Rscript nlband-review-extra.R <sections>
args <- commandArgs(trailingOnly = TRUE)
secs <- if (length(args)) strsplit(args[1], ",")[[1]] else
  c("bitwise", "double", "code1", "cost")
source("dev/nlband-common.R")
full <- function(tag, r) {
  if (iserr(r)) {
    cat("BIT", tag, "| ERROR", substr(r$msg, 1, 160), "\n")
    return(invisible())
  }
  se <- tryCatch(sqrt(diag(suppressWarnings(vcov(r)))),
                 error = function(e) NA)
  cat("BIT", tag, "| ll", format(llv(r), digits = 17), "code",
      r$opt$convergence, "par",
      paste(format(r$opt$par, digits = 17), collapse = " "), "se",
      paste(format(se, digits = 17), collapse = " "), "warn", wsum(r),
      "msg", msum(r), "\n")
}
if ("bitwise" %in% secs) {
  # no hole: x ~ U(0, 4), every row kept, truth b 0.5; seeds 11 to 20;
  # default control and nl_starts = FALSE; and a gamma body
  for (s in 11:20) {
    set.seed(s)
    d <- data.frame(x = runif(400, 0, 4))
    d$y <- 1 / (d$x + 0.5) + rnorm(400, 0, 0.2)
    full(paste("nohole default", s), tt(frm(band_form, data = d)))
    full(paste("nohole off", s), tt(frm(band_form, data = d,
      control = frmtmb_control(nl_starts = FALSE))))
    d$yp <- rpois(400, exp(1 / (d$x + 0.5)))
    full(paste("nohole poisson", s), tt(frm(
      bf(yp ~ a / (x + b), a + b ~ 1, nl = TRUE), data = d,
      family = poisson())))
  }
  # the band data with nl_starts = FALSE (a hole, switch off)
  for (s in 1:3) {
    full(paste("band off", s), tt(frm(band_form, data = band_data(s),
      control = frmtmb_control(nl_starts = FALSE))))
  }
}
if ("double" %in% secs) {
  # a double pole, y = 1 / (x - 2)^2 + noise, the band data's x: the same
  # model spelled a / (x + b)^2 and a * (x + b)^-2
  for (s in 1:3) {
    d <- band_data(s)
    d$y <- 1 / (d$x - 2)^2 + rnorm(nrow(d), 0, 0.2)
    for (f in list(bf(y ~ a / (x + b)^2, a + b ~ 1, nl = TRUE),
                   bf(y ~ a * (x + b)^-2, a + b ~ 1, nl = TRUE))) {
      r <- tt(frm(f, data = d))
      g <- tt(frm(f, data = d, start = list(beta = c(1, -1.9)),
                  control = frmtmb_control(nl_starts = FALSE)))
      cat(sprintf("DBL %s seed %d | default %.6f b %.4f msg %s warn %s | from hole %.6f b %.4f\n",
                  deparse(f$formula[[3]]), s, llv(r),
                  if (iserr(r)) NA else r$estimates$beta[2], msum(r),
                  wsum(r), llv(g), if (iserr(g)) NA else g$estimates$beta[2]))
    }
  }
  # EC50^3 + x^3 with x of both signs: walls at EC50 = -x_i, a cubic
  set.seed(7)
  x <- runif(300, -4, 4)
  x <- x[abs(x - 2) > 0.3]
  d <- data.frame(x = x, y = 2 * x^3 / (x^3 - 8) + rnorm(length(x), 0, 0.2))
  f <- bf(y ~ Emax * x^3 / (EC50^3 + x^3), Emax + EC50 ~ 1, nl = TRUE)
  r <- tt(frm(f, data = d))
  g <- tt(frm(f, data = d, start = list(beta = c(2, -2)),
              control = frmtmb_control(nl_starts = FALSE)))
  cat(sprintf("CUBIC default %.6f est %s msg %s warn %s | from EC50 -2 %.6f est %s\n",
              llv(r), if (iserr(r)) r$msg else
                paste(signif(r$estimates$beta, 5), collapse = " "),
              msum(r), wsum(r), llv(g), if (iserr(g)) g$msg else
                paste(signif(g$estimates$beta, 5), collapse = " ")))
}
if ("code1" %in% secs) {
  # the kept fit of "inside k5 3" (nlband-review-fits.R) ended at code 1
  d <- band_data(3)
  set.seed(1003)
  xi <- runif(5, 1.7, 2.3)
  d <- rbind(d, data.frame(x = xi, y = 1 / (xi - 2) + rnorm(5, 0, 0.2)))
  r <- tt(frm(band_form, data = d))
  cat("CODE1 default ll", format(llv(r), digits = 12), "code",
      r$opt$convergence, "message", r$opt$message, "| warn", wsum(r),
      "msg", msum(r), "\n")
  for (m in attr(r, "m")) cat("   MSG", m, "\n")
  for (w in attr(r, "w")) cat("   WARN", w, "\n")
  g <- tt(frm(band_form, data = d, start = list(beta = c(1, -2.0001))))
  cat("CODE1 given start ll", format(llv(g), digits = 12), "code",
      g$opt$convergence, "message", g$opt$message, "| warn", wsum(g),
      "msg", msum(g), "\n")
  for (m in attr(g, "m")) cat("   MSG", m, "\n")
  for (w in attr(g, "w")) cat("   WARN", w, "\n")
  cat("CODE1 xi", format(sort(xi), digits = 6), "\n")
}
if ("cost" %in% secs && arm == "after") {
  holes_fn <- ns$nl_pole_holes
  set_holes <- function(on) {
    f <- if (on) holes_fn else function(w) NULL
    unlockBinding("nl_pole_holes", ns)
    assign("nl_pole_holes", f, envir = ns)
    lockBinding("nl_pole_holes", ns)
  }
  n_ev <- 0
  suppressMessages(trace("optimize_obj",
    exit = quote(n_ev <<- n_ev + returnValue()$evals),
    print = FALSE, where = ns))
  mk <- list(
    "n1e4 no hole, truth b 0.5" = function(s) {
      set.seed(s); x <- runif(1e4, 0, 4)
      data.frame(x = x, y = 1 / (x + 0.5) + rnorm(1e4, 0, 0.2))
    },
    "n1e4 hole, truth b 0.5" = function(s) {
      set.seed(s); x <- runif(1e4, 0, 4); x <- x[abs(x - 2) > 0.3]
      data.frame(x = x, y = 1 / (x + 0.5) + rnorm(length(x), 0, 0.2))
    },
    "n1e4 hole, pole at 2" = function(s) {
      set.seed(s); x <- runif(1e4, 0, 4); x <- x[abs(x - 2) > 0.3]
      data.frame(x = x, y = 1 / (x - 2) + rnorm(length(x), 0, 0.2))
    })
  for (nm in names(mk)) for (s in 1:2) {
    d <- mk[[nm]](s)
    out <- c()
    for (a in c("on", "off", "on")) {
      set_holes(a == "on")
      n_ev <- 0
      t0 <- proc.time()[["elapsed"]]
      r <- tt(frm(band_form, data = d))
      out <- c(out, sprintf("%s ev %d t %.2f ll %.4f", a, n_ev,
                            proc.time()[["elapsed"]] - t0, llv(r)))
    }
    set_holes(TRUE)
    cat("COST", nm, "seed", s, "|", paste(out, collapse = " | "), "\n")
  }
}
prof1 <- function(x, y, b) {
  vapply(b, function(bb) {
    z <- 1 / (x + bb); a <- sum(y * z) / sum(z * z)
    rss <- sum((y - a * z)^2); n <- length(y)
    -n / 2 * (log(2 * pi * rss / n) + 1)
  }, 1)
}
glob1 <- function(x, y) {
  w <- sort(unique(-x)); r <- diff(range(w))
  u <- seq(0, 1, length.out = 17)[2:16]
  b <- sort(c(as.vector(outer(u, diff(w)) + rep(w[-length(w)], each = 15)),
              w[1] - r * 10^seq(-4, 2, length.out = 60),
              w[length(w)] + r * 10^seq(-4, 2, length.out = 60)))
  p <- prof1(x, y, b); p[!is.finite(p)] <- -Inf
  k <- which.max(p)
  o <- optimize(function(bb) prof1(x, y, bb), c(b[max(1, k - 1)],
                b[min(length(b), k + 1)]), maximum = TRUE, tol = 1e-12)
  c(max(o$objective, p[k]), if (o$objective >= p[k]) o$maximum else b[k])
}
if ("nearedge" %in% secs) {
  # the pole inside the hole but near its edge (hole |x - 2| < 0.3, pole
  # at p): does the keep rule (an end no nearer a wall than the median
  # gap) reject the maximum?
  for (p in c(2.25, 2.29, 2.297)) for (s in 1:2) {
    d <- band_data(s)
    set.seed(s + 500)
    d$y <- 1 / (d$x - p) + rnorm(nrow(d), 0, 0.2)
    r <- tt(frm(band_form, data = d))
    g <- glob1(d$x, d$y)
    gap <- median(diff(sort(unique(d$x))))
    near <- min(abs(d$x - p))
    cat(sprintf(paste0("EDGE pole %.3f seed %d | default %.6f b %.4f msg %s ",
                       "warn %s | GLOB %.6f at %.4f short %.3f | nearest row ",
                       "%.4f median gap %.4f\n"), p, s, llv(r),
                if (iserr(r)) NA else r$estimates$beta[2], msum(r), wsum(r),
                g[1], g[2], g[1] - llv(r), near, gap))
  }
}
if ("student" %in% secs) {
  for (s in 1:3) {
    set.seed(s)
    x <- band_data(s)$x
    d <- data.frame(x = x, y = 1 / (x - 2) + 0.2 * rt(length(x), 3))
    r <- tt(frm(band_form, data = d, family = student()))
    g <- tt(frm(band_form, data = d, family = student(),
                start = list(beta = c(1, -2)),
                control = frmtmb_control(nl_starts = FALSE)))
    for (v in list(list("default", r), list("from (1,-2)", g))) {
      f <- v[[2]]
      cat(sprintf("STUD seed %d %s | ll %.6f code %d | %s | warn %s msg %s\n",
                  s, v[[1]], llv(f), f$opt$convergence,
                  paste(names(f$opt$par), signif(f$opt$par, 6), collapse = " "),
                  wsum(f), msum(f)))
      for (w in attr(f, "w")) cat("   WARN", substr(w, 1, 400), "\n")
    }
  }
}
if ("bsep" %in% secs) {
  # separated bernoulli data: y = 1 exactly for x > 2 (band x)
  for (s in 1:2) {
    d <- band_data(s)
    d$y <- as.integer(d$x > 2)
    r <- tt(frm(band_form, data = d, family = bernoulli()))
    cat(sprintf("BSEP seed %d | ll %.6f code %s beta %s | warn %s msg %s\n",
                s, llv(r), if (iserr(r)) "ERR" else r$opt$convergence,
                if (iserr(r)) r$msg else
                  paste(signif(r$estimates$beta, 5), collapse = " "),
                wsum(r), msum(r)))
    for (w in attr(r, "w")) cat("   WARN", substr(w, 1, 400), "\n")
    for (m in attr(r, "m")) cat("   MSG", substr(m, 1, 400), "\n")
  }
}
if ("prof" %in% secs) {
  # bernoulli profiles over b inside the hole, a out, on the log scale
  pll <- function(d, bb) {
    z <- 1 / (d$x + bb)
    optimize(function(a) {
      e <- a * z
      sum(d$y * plogis(e, log.p = TRUE) + (1 - d$y) * plogis(-e, log.p = TRUE))
    }, c(-60, 60), maximum = TRUE, tol = 1e-10)$objective
  }
  for (kind in c("bstep", "bband")) for (s in 1:3) {
    if (kind == "bstep") {
      set.seed(s); x <- band_data(s)$x
      d <- data.frame(x = x, y = rbinom(length(x), 1, plogis(4 * (x - 2))))
    } else {
      d <- band_data(s)
      yb <- rbinom(nrow(d), 10, plogis(1 / (d$x - 2)))
      d$y <- rbinom(nrow(d), 1, plogis(2 / (d$x - 2)))
    }
    lo <- -min(d$x[d$x > 2]); hi <- -max(d$x[d$x < 2])
    bg <- seq(lo, hi, length.out = 401)[-c(1, 401)]
    pv <- vapply(bg, function(bb) pll(d, bb), 1)
    k <- which.max(pv)
    cat(sprintf(paste0("PROF %s %d hole %.4f to %.4f: max %.6f at b %.4f ",
                       "(%d of 399), ends %.6f %.6f\n"), kind, s, lo, hi,
                pv[k], bg[k], k, pv[1], pv[399]))
  }
}
if ("dilw" %in% secs) {
  # does nl_pole_walls() read the walls of a / (x + b) on a 30-level
  # two-fold dilution series (x to 5.4e5)?
  for (s in 1:3) {
    set.seed(s)
    x <- rep(0.001 * 2^(0:29), each = 4)
    d <- data.frame(x = x, y = 5 / (x + 0.5) + rnorm(120, 0, 0.2))
    f <- suppressWarnings(suppressMessages(frm(band_form, data = d,
      control = frmtmb_control(nl_starts = FALSE))))
    lp <- Filter(function(q) !is.null(q$nl_body), f$frame$linpreds)[[1L]]
    wl <- ns$nl_pole_walls(f$frame, lp, f$estimates, "beta", 2L,
                           ns$nl_denoms(lp$nl_body))
    cat(sprintf("DILW seed %d b %.17g | walls read %s holes %d\n", s,
                f$estimates$beta[2], !is.null(wl),
                if (is.null(wl)) NA else NROW(wl[[1]]$holes)))
  }
}
cat("DONE\n")
