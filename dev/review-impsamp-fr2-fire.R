# Fix round 2, item 3: the natural firings the search found, at default rounds. Is
# each message true when read at the count it names and at the counts above? Is the
# firing stable under the draw seed, which a real-fit test would need?
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
pois_big <- function(s) {
  set.seed(s); g <- factor(rep(1:10, each = 2L)); x <- rnorm(20)
  list(d = data.frame(y = rpois(20, exp(2 + rnorm(10, 0, 3)[g])), x = x, g = g),
       f = bf(y ~ x + (x | g)), fam = poisson())
}
sep_tight <- function(s) {
  set.seed(s); g <- factor(rep(1:12, each = 3L)); x <- rnorm(36)
  list(d = data.frame(y = as.integer(2 * x + rnorm(12, 0, 3)[g] > 0), x = x, g = g),
       f = bf(y ~ x + (x | g)), fam = bernoulli())
}
tryfit <- function(m, N, ds = 1L) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(frm(m$f, family = m$fam, data = m$d, importance = N,
                 control = frmtmb_control(importance_seed = ds)), error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 36))
      invokeRestart("muffleWarning")
    })
  if (inherits(r, "error")) {
    sprintf("STOP: %s", sub(".*which is [^.]*[.]?[0-9]* of [0-9]+[.] ", "", substr(conditionMessage(r), 1, 520)))
  } else {
    im <- r$importance
    sprintf("fit rounds %d capped %s ll %.3f mcse %.3f ess %.3f warns: %s", im$rounds, im$capped, as.numeric(logLik(r)),
            im$mcse, im$ess_min, paste(w, collapse = ";"))
  }
}
m <- pois_big(14L)
l <- suppressWarnings(frm(m$f, family = m$fam, data = m$d))
cat("pois_big 14: Laplace ll", as.numeric(logLik(l)), "conv", l$opt$convergence, "par", format(l$opt$par, digits = 4), "\n")
for (N in c(8L, 16L, 32L, 64L, 256L, 1024L)) cat(" N", N, ":", tryfit(m, N), "\n")
for (ds in 2:6) cat(" N 8, draw seed", ds, ":", tryfit(m, 8L, ds), "\n")
m <- pois_big(16L)
for (N in c(32L, 64L, 512L, 2048L)) cat("pois_big 16 N", N, ":", tryfit(m, N), "\n")
m <- sep_tight(24L)
for (ds in 1:4) cat("sep_tight 24 N 8, draw seed", ds, ":", tryfit(m, 8L, ds), "\n")
