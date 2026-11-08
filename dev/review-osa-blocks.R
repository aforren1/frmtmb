# Review lane osa, short check: residuals(type = "osa") in row blocks (the
# fix-round default) against one call over every row
# (options(frmtmb.osa_block = Inf)) and against the exact one-step
# predictive, on (1 | g) and latent ar() designs, rows group-major or
# interleaved (time-major: row order t = 1 of every group, then t = 2...).
#
#   Rscript dev/review-osa-blocks.R
#
# Exact: as dev/review-osa-exact.R, written from stats:: functions; rows of
# a group are conditioned in data order, as TMB conditions them.
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
fams <- list(
  poisson = list(fam = poisson(), inv = exp,
                 p = function(y, m) dpois(y, m), F = function(y, m) ppois(y, m),
                 r = function(m) rpois(length(m), m)),
  bernoulli = list(fam = bernoulli(), inv = plogis,
                   p = function(y, m) dbinom(y, 1, m),
                   F = function(y, m) pbinom(y, 1, m),
                   r = function(m) rbinom(length(m), 1, m)))
make <- function(G, K, design, order, ff, seed) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(G), each = K)), t = rep(1:K, G),
                  x = rnorm(G * K))
  e <- if (design == "re") rnorm(G, 0, 1)[d$g] else
    as.vector(replicate(G, arima.sim(list(ar = 0.5), K, sd = 0.8)))
  d$y <- ff$r(ff$inv(0.2 + 0.5 * d$x + e))
  if (order == "interleaved") d <- d[order(d$t, d$g), ]
  rownames(d) <- NULL
  d
}
exact <- function(fit, d, ff, design, U) {
  eta <- fit$estimates$beta[[1]] + fit$estimates$beta[[2]] * d$x
  th <- fit$estimates$theta
  if (design == "re") {
    sdu <- exp(th[[1]])
    s <- seq(-8, 8, length.out = 1601) * sdu
    w0 <- dnorm(s, 0, sdu)
    scale <- 1
    Kt <- NULL
  } else {
    scale <- exp(th[[1]])
    phi <- frmtmb:::autocor_pacf(th[[2]])
    sds <- 1 / sqrt(1 - phi^2)
    s <- seq(-8, 8, length.out = 801) * sds
    Kt <- outer(s, s, function(a, b) dnorm(b, phi * a, 1))
    w0 <- dnorm(s, 0, sds)
  }
  lo <- hi <- numeric(nrow(d))
  for (gi in split(seq_len(nrow(d)), d$g)) {
    w <- w0 / sum(w0)
    for (i in gi) {
      m <- ff$inv(eta[i] + scale * s)
      p <- ff$p(d$y[i], m)
      Fh <- ff$F(d$y[i], m)
      hi[i] <- sum(w * Fh)
      lo[i] <- sum(w * (Fh - p))
      w <- w * p
      if (!is.null(Kt)) w <- as.vector(w %*% Kt)
      w <- w / sum(w)
    }
  }
  qnorm(hi - U * (hi - lo))
}
osa <- function(fit) suppressWarnings(residuals(fit, type = "osa")[, "Estimate"])
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  list(v = v, s = proc.time()[["elapsed"]] - t0)
}
cat("BEGIN review-osa-blocks\n")
runs <- list(
  list(20, 10, "re", "group-major"), list(20, 10, "re", "interleaved"),
  list(20, 10, "ar", "group-major"), list(20, 10, "ar", "interleaved"))
for (fn in names(fams)) for (rr in runs) {
  ff <- fams[[fn]]
  d <- make(rr[[1]], rr[[2]], rr[[3]], rr[[4]], ff, 701)
  rhs <- if (rr[[3]] == "re") "x + (1 | g)" else "x + ar(t, g, cov = TRUE)"
  fit <- suppressWarnings(frm(as.formula(paste("y ~", rhs)), family = ff$fam,
                              data = d))
  n <- nrow(d)
  U <- local({ set.seed(123); runif(n) })
  b <- tm(osa(fit))
  op <- options(frmtmb.osa_block = Inf)
  u <- tm(osa(fit))
  options(op)
  ez <- exact(fit, d, ff, rr[[3]], U)
  cat(sprintf(paste("%-9s %-3s %-11s n %4d | max|blocked - one call| %.3g",
                    "| max|blocked - exact| %.3g, one call %.3g |",
                    "blocked %.1fs one call %.1fs\n"),
              fn, rr[[3]], rr[[4]], n, max(abs(b$v - u$v)),
              max(abs(b$v - ez)), max(abs(u$v - ez)), b$s, u$s))
}
# larger, poisson only
for (design in c("re", "ar")) {
  ff <- fams$poisson
  d <- make(60, 10, design, "group-major", ff, 702)
  rhs <- if (design == "re") "x + (1 | g)" else "x + ar(t, g, cov = TRUE)"
  fit <- suppressWarnings(frm(as.formula(paste("y ~", rhs)), family = ff$fam,
                              data = d))
  n <- nrow(d)
  U <- local({ set.seed(123); runif(n) })
  b <- tm(osa(fit))
  op <- options(frmtmb.osa_block = Inf)
  u <- tm(osa(fit))
  options(op)
  ez <- exact(fit, d, ff, design, U)
  cat(sprintf(paste("%-9s %-3s %-11s n %4d | max|blocked - one call| %.3g",
                    "| max|blocked - exact| %.3g, one call %.3g |",
                    "blocked %.1fs one call %.1fs\n"),
              "poisson", design, "group-major", n, max(abs(b$v - u$v)),
              max(abs(b$v - ez)), max(abs(u$v - ez)), b$s, u$s))
}
cat("END review-osa-blocks\n")
