# Fix round 1, item 1: the two things own.R left open.
# (1) an uncentered column on an ordinal: the Laplace fit against the importance
# fit against the centered column, per data seed, with the log-likelihoods.
# (2) x + m on poisson: the Laplace fit's slope against the centered fit's too.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
fitw <- function(...) {
  w <- character(0)
  f <- withCallingHandlers(tryCatch(frm(...), error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, substr(conditionMessage(x), 1, 60))
                             invokeRestart("muffleWarning")
                           })
  list(f = f, w = w)
}
ses <- function(f) tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)), error = function(e) NA)
cat("(1) ordinal\n")
for (m in c(2000, 1e4)) for (s in 1:4) {
  set.seed(s)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10))
  z <- rnorm(150)
  lat <- 0.7 * z + rnorm(G, 0, 0.7)[g] + rlogis(150)
  y <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
  d <- data.frame(y, w = 5 * z + m, wc = 5 * z, g)
  li <- fitw(bf(y ~ w + (1 | g)), family = cumulative(), data = d)
  lc <- fitw(bf(y ~ wc + (1 | g)), family = cumulative(), data = d)
  ii <- fitw(bf(y ~ w + (1 | g)), family = cumulative(), data = d, importance = 300L)
  ic <- fitw(bf(y ~ wc + (1 | g)), family = cumulative(), data = d, importance = 300L)
  st <- function(x) if (inherits(x$f, "error")) "STOP" else sprintf("%.3f", as.numeric(logLik(x$f)))
  th <- function(x) if (inherits(x$f, "error")) NA else x$f$opt$par[[length(x$f$opt$par)]]
  cat(sprintf("  m %-6g s %d ll: lap-w %s lap-c %s imp-w %s imp-c %s | log sd: %.3f %.3f %.3f %.3f | conv lap-w %d imp-w %d | warns imp-w %d %s\n",
              m, s, st(li), st(lc), st(ii), st(ic), th(li), th(lc), th(ii), th(ic),
              li$f$opt$convergence, if (inherits(ii$f, "error")) NA else ii$f$opt$convergence,
              length(ii$w), paste(ii$w, collapse = ";")))
}
cat("(2) poisson x + m: z of the slope against the centered fit\n")
for (m in c(3e4, 1e5, 3e5, 1e6, 1e7)) for (s in 1:3) {
  set.seed(s)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10))
  z <- rnorm(150)
  y <- rpois(150, exp(0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]))
  d <- data.frame(y, w = z + m, wc = z, g)
  ii <- fitw(bf(y ~ w + (1 | g)), family = poisson(), data = d, importance = 300L)
  li <- fitw(bf(y ~ w + (1 | g)), family = poisson(), data = d)
  rr <- fitw(bf(y ~ wc + (1 | g)), family = poisson(), data = d, importance = 300L)
  se <- ses(rr$f)
  if (inherits(ii$f, "error")) {
    cat(sprintf("  m %-6g s %d STOP %s\n", m, s, substr(conditionMessage(ii$f), 1, 120)))
    next
  }
  cat(sprintf("  m %-6g s %d imp slope z %+.3f theta z %+.3f conv %d SE %.4f warns %d | lap slope z %+.3f theta z %+.3f conv %d warns %d | %s\n",
              m, s, (ii$f$opt$par[[2]] - rr$f$opt$par[[2]]) / se[[2]],
              (ii$f$opt$par[[3]] - rr$f$opt$par[[3]]) / se[[3]], ii$f$opt$convergence,
              ses(ii$f)[[2]] / se[[2]], length(ii$w),
              (li$f$opt$par[[2]] - rr$f$opt$par[[2]]) / se[[2]],
              (li$f$opt$par[[3]] - rr$f$opt$par[[3]]) / se[[3]], li$f$opt$convergence, length(li$w),
              paste(c(ii$w)[seq_len(min(1, length(ii$w)))], collapse = ";")))
}
