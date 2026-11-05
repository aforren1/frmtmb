# Fix round 1, item 1: (a) a family's own parameter with no column, on a large
# scale: ordinal thresholds absorbing an uncentered column; (b) the filed x + 1e5
# Poisson, both arms, Laplace too; (c) a hard design (Laplace far from the optimum)
# on large columns: does the fit correct, or stop at the Laplace optimum with a
# false-convergence warning? (d) gaussian with a dispersion regression, no prior,
# against Laplace, the control for the prior table.
# Usage: REV_ARM=after|before|mutant Rscript dev/review-impsamp-fr1-own.R
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "lib", find.package("frmtmb"), "\n")
fitw <- function(...) {
  w <- character(0)
  f <- withCallingHandlers(tryCatch(frm(...), error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, substr(conditionMessage(x), 1, 75))
                             invokeRestart("muffleWarning")
                           })
  list(f = f, w = w)
}
ses <- function(f) tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)), error = function(e) NA)
cat("(a) ordinal: y ~ w + (1 | g), w = 5 z + m, 3 categories, 300 draws\n")
for (m in c(0, 2000, 1e4, 3e4, 1e5)) for (s in 1:3) {
  set.seed(s)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10))
  z <- rnorm(150)
  lat <- 0.7 * z + rnorm(G, 0, 0.7)[g] + rlogis(150)
  y <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
  d <- data.frame(y, w = 5 * z + m, wc = 5 * z, g)
  a <- fitw(bf(y ~ w + (1 | g)), family = cumulative(), data = d, importance = 300L)
  r <- fitw(bf(y ~ wc + (1 | g)), family = cumulative(), data = d, importance = 300L)
  if (inherits(a$f, "error")) {
    cat(sprintf("  m %-6g s %d STOP %s\n", m, s, substr(conditionMessage(a$f), 1, 150)))
    next
  }
  sa <- ses(a$f); sr <- ses(r$f)
  cat(sprintf("  m %-6g s %d rounds %d capped %s conv %d grad %.3g moved %.3g slopeSE %.4f dslope/SE %+.3f ll-c %+.4f mcse %.4f warns %d %s\n",
              m, s, a$f$importance$rounds, a$f$importance$capped, a$f$opt$convergence, a$f$importance$grad,
              a$f$importance$moved, sa[[1]] / sr[[1]], (a$f$opt$par[[1]] - r$f$opt$par[[1]]) / sr[[1]],
              as.numeric(logLik(a$f)) - as.numeric(logLik(r$f)), a$f$importance$mcse, length(a$w),
              paste(a$w, collapse = ";")))
}
cat("\n(b) poisson x + m: importance, and Laplace on the same column\n")
for (m in c(1e3, 1e4, 3e4, 1e5, 1e6)) for (s in 1:2) {
  set.seed(s)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10))
  z <- rnorm(150)
  y <- rpois(150, exp(0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]))
  d <- data.frame(y, w = z + m, wc = z, g)
  a <- fitw(bf(y ~ w + (1 | g)), family = poisson(), data = d, importance = 300L)
  l <- fitw(bf(y ~ w + (1 | g)), family = poisson(), data = d)
  r <- fitw(bf(y ~ wc + (1 | g)), family = poisson(), data = d, importance = 300L)
  if (inherits(a$f, "error")) {
    cat(sprintf("  m %-6g s %d STOP %s\n", m, s, substr(conditionMessage(a$f), 1, 170)))
    next
  }
  cat(sprintf("  m %-6g s %d imp: rounds %d conv %d moved %.3g slopeSE %.4f slope z %+.3f theta z %+.3f warns %d | Laplace: conv %d slopeSE %.4f warns %d %s\n",
              m, s, a$f$importance$rounds, a$f$opt$convergence, a$f$importance$moved, ses(a$f)[[2]] / ses(r$f)[[2]],
              (a$f$opt$par[[2]] - r$f$opt$par[[2]]) / ses(r$f)[[2]],
              (a$f$opt$par[[3]] - r$f$opt$par[[3]]) / ses(r$f)[[3]],
              length(a$w), l$f$opt$convergence, ses(l$f)[[2]] / ses(r$f)[[2]], length(l$w),
              paste(c(a$w, l$w)[1:min(1, length(c(a$w, l$w)))], collapse = ";")))
}
cat("\n(c) hard: bernoulli sd 4, 30 groups of 5, seed 3, on columns: does it correct?\n")
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
de <- design_list$bigsd_bern(3L)
rf <- suppressWarnings(frm(de$formula, family = de$family, data = de$data, importance = 500L))
lf <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
ser <- ses(rf)
cat(sprintf("  standard column: theta imp %.4f Laplace %.4f (SE %.3f), intercept imp %.3f lap %.3f\n",
            rf$opt$par[[3]], lf$opt$par[[3]], ser[[3]], rf$opt$par[[1]], lf$opt$par[[1]]))
for (cs in list(c(1e3, 0), c(1e4, 0), c(1e5, 0), c(1, 2000), c(100, 1e4), c(1, 1e4), c(1, 1e5))) {
  d <- transform(de$data, w = cs[1] * x + cs[2])
  a <- fitw(bf(y ~ w + (1 | g)), family = bernoulli(), data = d, importance = 500L)
  l <- fitw(bf(y ~ w + (1 | g)), family = bernoulli(), data = d)
  if (inherits(a$f, "error")) {
    cat(sprintf("  %g x + %g STOP %s\n", cs[1], cs[2], substr(conditionMessage(a$f), 1, 170)))
    next
  }
  cat(sprintf("  %g x + %g: theta imp %.4f (z from standard %+.3f) lap %.4f rounds %d conv %d moved %.3g grad %.3g slopeSE ratio %.4f warns %s\n",
              cs[1], cs[2], a$f$opt$par[[3]], (a$f$opt$par[[3]] - rf$opt$par[[3]]) / ser[[3]],
              l$f$opt$par[[3]], a$f$importance$rounds, a$f$opt$convergence, a$f$importance$moved,
              a$f$importance$grad, ses(a$f)[[2]] * cs[1] / ser[[2]], paste(a$w, collapse = ";")))
}
cat("\n(d) gaussian, (x | g), sigma ~ z, no prior, against Laplace\n")
set.seed(5)
G <- 25L
g <- factor(rep(seq_len(G), each = 6))
x <- rnorm(150); z <- rnorm(150)
a0 <- rnorm(G, 0, 0.8); b0 <- rnorm(G, 0, 0.4)
y <- 0.3 + 0.5 * x + a0[g] + b0[g] * x + rnorm(150, 0, exp(-0.5 + 0.3 * z))
d <- data.frame(y, x, z, g)
for (N in c(64L, 512L)) {
  fi <- suppressWarnings(frm(bf(y ~ x + (x | g), sigma ~ z), family = gaussian(), data = d, importance = N))
  fl <- suppressWarnings(frm(bf(y ~ x + (x | g), sigma ~ z), family = gaussian(), data = d))
  cat(sprintf("  N %d: obj-lap %+.3g maxpar %.3g mcse %.4f rounds %d\n", N, fi$opt$objective - fl$opt$objective,
              max(abs(fi$opt$par - fl$opt$par)), fi$importance$mcse, fi$importance$rounds))
}
