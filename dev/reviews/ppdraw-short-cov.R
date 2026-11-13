# Short check of the fix round: an EXACT reference for the side law on
# interior fits. Data are simulated from a known model, so the predictive
# law of a new observation is known, and a predict() interval's coverage
# of it is exact (no Monte Carlo in the reference). Each fit is tagged by
# whether the side law changed its draws (the frmtmb_draw_held message).
#
#   R: gaussian y ~ x + (1 | g), 6 groups x 20, sigma 1, group sd 0.1 or
#      0.3; a new level at x = 0: y_new ~ N(1, sd_g^2 + 1).
#   S: gaussian y ~ x + (1 + x | g), 15 groups x 8, slope sd 0.1 or 0.2,
#      intercept sd 0.8; a new level at x = 2:
#      y_new ~ N(2, 0.8^2 + 4 sd_s^2 + 1).
#   Z: zero_inflated_poisson on poisson data (no zero inflation), 200 rows;
#      predictive mean and P(y = 0) against the true ones.
#
# Usage: Rscript ppdraw-short-cov.R <lib or "base"> <out.tsv> [nrep]
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
out <- args[2]
nrep <- if (length(args) >= 3) as.integer(args[3]) else 100L
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
q <- function(expr) suppressWarnings(expr)
pred <- function(fit, nd, pe = TRUE, ndraws = 1000L) {
  said <- FALSE
  set.seed(99)
  p <- withCallingHandlers(
    q(predict(fit, newdata = nd, ndraws = ndraws, summary = FALSE,
              allow_new_levels = TRUE, propagate_error = pe)),
    frmtmb_draw_held = function(m) {
      said <<- TRUE
      invokeRestart("muffleMessage")
    },
    message = function(m) invokeRestart("muffleMessage"))
  list(p = p, said = said)
}
cover <- function(p, m, s, lev) {
  a <- (1 - lev) / 2
  qq <- stats::quantile(p, c(a, 1 - a), names = FALSE)
  c(stats::pnorm((qq[2] - m) / s) - stats::pnorm((qq[1] - m) / s),
    qq[2] - qq[1])
}
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...)
for (sdg in c(0.1, 0.3)) for (seed in seq_len(nrep)) {
  set.seed(1000 + seed)
  d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, sdg)[d$g] + stats::rnorm(120)
  fit <- tryCatch(q(suppressMessages(frm(bf(y ~ x + (1 | g)), data = d))),
                  error = function(e) NULL)
  if (is.null(fit)) next
  nd <- data.frame(x = 0, g = factor("new"))
  st <- sqrt(sdg^2 + 1)
  for (pe in c(TRUE, FALSE)) {
    r <- pred(fit, nd, pe)
    c90 <- cover(r$p, 1, st, 0.9)
    c98 <- cover(r$p, 1, st, 0.98)
    add(case = "R", par = sdg, seed = seed, pe = pe, changed = r$said,
        cov90 = c90[1], w90 = c90[2], cov98 = c98[1], w98 = c98[2],
        est = NA, truth = NA, p0 = NA, p0true = NA)
  }
}
for (sds in c(0.1, 0.2)) for (seed in seq_len(nrep %/% 2L)) {
  set.seed(2000 + seed)
  G <- 15
  m <- 8
  d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                  x = stats::rnorm(G * m))
  ug <- stats::rnorm(G, 0, 0.8)
  sg <- stats::rnorm(G, 0, sds)
  d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + stats::rnorm(G * m)
  fit <- tryCatch(q(suppressMessages(frm(bf(y ~ x + (1 + x | g)),
                                         data = d))),
                  error = function(e) NULL)
  if (is.null(fit)) next
  nd <- data.frame(x = 2, g = factor("new"))
  st <- sqrt(0.8^2 + 4 * sds^2 + 1)
  for (pe in c(TRUE, FALSE)) {
    r <- pred(fit, nd, pe)
    c90 <- cover(r$p, 2, st, 0.9)
    c98 <- cover(r$p, 2, st, 0.98)
    add(case = "S", par = sds, seed = seed, pe = pe, changed = r$said,
        cov90 = c90[1], w90 = c90[2], cov98 = c98[1], w98 = c98[2],
        est = NA, truth = NA, p0 = NA, p0true = NA)
  }
}
for (seed in 1:30) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(200))
  mu <- exp(1 + 0.3 * d$x)
  d$y <- stats::rpois(200, mu)
  fit <- q(suppressMessages(frm(bf(y ~ x), family = zero_inflated_poisson(),
                                data = d)))
  for (pe in c(TRUE, FALSE)) {
    said <- FALSE
    set.seed(1)
    p <- withCallingHandlers(
      q(predict(fit, ndraws = 2000, summary = FALSE, propagate_error = pe)),
      frmtmb_draw_held = function(m) {
        said <<- TRUE
        invokeRestart("muffleMessage")
      },
      message = function(m) invokeRestart("muffleMessage"))
    add(case = "Z", par = 0, seed = seed, pe = pe, changed = said,
        cov90 = NA, w90 = NA, cov98 = NA, w98 = NA, est = mean(p),
        truth = mean(mu), p0 = mean(p == 0), p0true = mean(exp(-mu)))
  }
}
res <- do.call(rbind, rows)
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("rows", nrow(res), "\n")
