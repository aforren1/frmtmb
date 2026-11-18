# Lane nlband, fix round 2: the short check's sweep
# (dev/reviews/nlband-short-extra.R sweep: n 60, the narrowest hole about
# 2 that nl_pole_holes() flags, y = amp / (x - 2) + N(0, sd), amp 0.15 to
# 0.4, sd 1 and 2, seeds 1 to 15) on the round-2 build. Each call of
# nl_hole_se() is recorded with its standard error, the hole's half width
# and the walls' span, so the refusals of the round-1 rule (se not below
# half the width) and of the round-2 rule (se not below the span) are
# both counted on the same calls. Each fit prints its logLik, the global
# maximum (exact profile over b) and its pole classes.
# Rscript dev/nlband-sweep.R
source("dev/nlband-common.R")
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
rows <- list()
suppressMessages(trace("nl_hole_se", where = ns, print = FALSE,
  exit = quote({
    .nlb_row <- c(se = returnValue(), half = (ck$hi - ck$lo) / 2,
                  span = ck$span,
                  ll = tryCatch(as.numeric(logLik(fit)),
                                error = function(e) NA))
    assign("rows", c(get("rows", envir = globalenv()), list(.nlb_row)),
           envir = globalenv())
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
tot <- 0; ncall <- 0; rej1 <- 0; rej2 <- 0; maxhalf <- 0; maxspan <- 0
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
    for (lr in rows) {
      ncall <- ncall + 1
      maxhalf <- max(maxhalf, lr[["se"]] / lr[["half"]], na.rm = TRUE)
      maxspan <- max(maxspan, lr[["se"]] / lr[["span"]], na.rm = TRUE)
      r1 <- !isTRUE(lr[["se"]] < lr[["half"]])
      r2 <- !isTRUE(lr[["se"]] < lr[["span"]])
      rej1 <- rej1 + r1
      rej2 <- rej2 + r2
      if (r1 || r2) {
        cat(sprintf(paste0("CALL amp %.2f sd %.0f seed %d | se %.4g half %.4g ",
                           "span %.4g | refused round1 %s round2 %s | refit ",
                           "ll %.6f | returned ll %.6f code %d b %.4f warn %s ",
                           "msg %s | GLOB %.6f at %.4f\n"),
                    amp, sdv, s, lr[["se"]], lr[["half"]], lr[["span"]], r1,
                    r2, lr[["ll"]], llv(r), r$opt$convergence,
                    r$estimates$beta[2], wsum(r), msum(r), g[1], g[2]))
      }
    }
  }
}
cat(sprintf(paste0("SWEEP fits %d | calls %d | refused by the round-1 rule ",
                   "%d | by the round-2 rule %d | max se/half %.3f | max ",
                   "se/span %.3f\n"), tot, ncall, rej1, rej2, maxhalf,
            maxspan))
cat("DONE\n")
