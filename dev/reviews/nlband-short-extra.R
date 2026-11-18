# Short check of lane nlband, follow-ups. supz: the separated-z body on
# both arms, and from inside the hole with nl_starts = FALSE. sweep: a
# wider search for a correct hole fit that nl_hole_located() refuses
# (n 60, narrowest flagged hole, amp 0.15 to 0.4, sd 1 and 2, seeds 1 to
# 15). studg: student() on gaussian band data (nu at its limit at the
# maximum): which boundary message the kept hole fit carries.
# NLBAND_ARM=before|after Rscript dev/reviews/nlband-short-extra.R <secs>
args <- commandArgs(trailingOnly = TRUE)
secs <- strsplit(args[1], ",")[[1]]
source("dev/nlband-common.R")
cls <- function(r, k) {
  v <- attr(r, k)
  if (length(v)) paste(unique(v), collapse = "+") else "-"
}
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
if ("supz" %in% secs) {
  form2 <- bf(y ~ a / (x + b) + c * z, a + b + c ~ 1, nl = TRUE)
  for (s in 1:5) {
    d <- band_data(s)
    set.seed(900 + s)
    d$z <- rbinom(nrow(d), 1, 0.2)
    d$y <- rbinom(nrow(d), 1, plogis(1 / (d$x - 2)))
    d$y[d$z == 1] <- 1L
    r <- tt(frm(form2, data = d, family = bernoulli()))
    g <- tt(frm(form2, data = d, family = bernoulli(),
                start = list(beta = c(1, -2, 0)),
                control = frmtmb_control(nl_starts = FALSE)))
    for (v in list(list("default", r), list("hole start, nl_starts off", g))) {
      f <- v[[2]]
      se <- tryCatch(sqrt(diag(suppressWarnings(vcov(f)))),
                     error = function(e) NA)
      cat(sprintf(paste0("SUPZ seed %d %s | ll %.6f code %d beta %s se %s | ",
                         "warn %s msg %s\n"), s, v[[1]], llv(f),
                  f$opt$convergence,
                  paste(sprintf("%.5g", f$estimates$beta), collapse = " "),
                  paste(signif(se, 4), collapse = " "),
                  cls(f, "wc"), cls(f, "mc")))
      for (w in attr(f, "w")) cat("   WARN", substr(w, 1, 200), "\n")
    }
  }
}
if ("sweep" %in% secs && arm == "after") {
  rows <- list()
  suppressMessages(trace("nl_hole_located", where = ns, print = FALSE,
    exit = quote({
      .row <- c(se = if (exists("s", inherits = FALSE)) s else NA,
                half = (ck$hi - ck$lo) / 2, ok = returnValue(),
                ll = tryCatch(as.numeric(logLik(fit)), error = function(e) NA))
      rows[[length(rows) + 1L]] <<- .row
    })))
  n <- 60
  narrow <- function(seed) {
    for (h in seq(0.2, 1.6, by = 0.02)) {
      set.seed(seed)
      x <- runif(2000, 0, 4)
      x <- x[abs(x - 2) > h][seq_len(n)]
      hs <- ns$nl_pole_holes(sort(unique(-x)))
      if (!is.null(hs) && any(hs[, "lo"] < -2 & hs[, "hi"] > -2)) {
        return(list(x = x, h = h))
      }
    }
    NULL
  }
  tot <- 0; nloc <- 0; nrej <- 0; nrej_ok <- 0; maxratio <- 0
  for (sdv in c(1, 2)) for (amp in c(0.15, 0.2, 0.25, 0.3, 0.4)) {
    for (s in 1:15) {
      nw <- narrow(s)
      if (is.null(nw)) next
      set.seed(500 + s)
      d <- data.frame(x = nw$x)
      d$y <- amp / (d$x - 2) + rnorm(n, 0, sdv)
      g <- glob1(d$x, d$y)
      rows <- list()
      r <- tt(frm(band_form, data = d))
      tot <- tot + 1
      ingap <- g[2] > -2 - nw$h && g[2] < -2 + nw$h
      for (lr in rows) {
        nloc <- nloc + 1
        maxratio <- max(maxratio, lr[["se"]] / lr[["half"]], na.rm = TRUE)
        if (!isTRUE(as.logical(lr[["ok"]]))) {
          nrej <- nrej + 1
          right <- abs(lr[["ll"]] - g[1]) < 1e-6 * max(1, abs(g[1]))
          if (right) nrej_ok <- nrej_ok + 1
          cat(sprintf(paste0("REJ amp %.2f sd %.0f seed %d h %.2f | se %.4g half ",
                             "%.4g | refused refit ll %.6f | returned ll %.6f ",
                             "code %d b %.4f warn %s msg %s | GLOB %.6f at %.4f",
                             " in hole %s | refused was the max %s\n"),
                      amp, sdv, s, nw$h, lr[["se"]], lr[["half"]], lr[["ll"]],
                      llv(r), r$opt$convergence, r$estimates$beta[2],
                      cls(r, "wc"), cls(r, "mc"), g[1], g[2], ingap, right))
        }
      }
    }
  }
  cat(sprintf(paste0("SWEEP fits %d | located calls %d | refused %d | ",
                     "refused that were the maximum %d | max se/half %.3f\n"),
              tot, nloc, nrej, nrej_ok, maxratio))
}
if ("studg" %in% secs) {
  for (s in 1:3) {
    d <- band_data(s)
    r <- tt(frm(band_form, data = d, family = student()))
    g <- tt(frm(band_form, data = d, family = student(),
                start = list(beta = c(1, -2)),
                control = frmtmb_control(nl_starts = FALSE)))
    for (v in list(list("default", r), list("from (1,-2)", g))) {
      f <- v[[2]]
      cat(sprintf("STUDG seed %d %s | ll %.6f code %d par %s | warn %s msg %s\n",
                  s, v[[1]], llv(f), f$opt$convergence,
                  paste(sprintf("%.5g", f$opt$par), collapse = " "),
                  cls(f, "wc"), cls(f, "mc")))
    }
  }
}
cat("DONE\n")
