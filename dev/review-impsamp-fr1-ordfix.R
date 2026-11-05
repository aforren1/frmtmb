source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs(Sys.getenv("REV_ARM", "after"))
suppressPackageStartupMessages(library(frmtmb))
cat(find.package("frmtmb"), "\n")
ses <- function(f) tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)), error = function(e) NA)
for (m in c(2000, 1e4, 3e4, 1e5)) for (s in c(3L, 5L, 6L, 8L)) {
  set.seed(s)
  g <- factor(rep(seq_len(15L), each = 10))
  z <- rnorm(150)
  lat <- 0.7 * z + rnorm(15, 0, 0.7)[g] + rlogis(150)
  y <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
  d <- data.frame(y, w = 5 * z + m, wc = 5 * z, g)
  w <- character(0)
  a <- withCallingHandlers(tryCatch(frm(bf(y ~ w + (1 | g)), family = cumulative(), data = d, importance = 300L), error = function(e) e),
    warning = function(x) { w <<- c(w, substr(conditionMessage(x), 1, 60)); invokeRestart("muffleWarning") })
  r <- suppressWarnings(frm(bf(y ~ wc + (1 | g)), family = cumulative(), data = d, importance = 300L))
  if (inherits(a, "error")) { cat("m", m, "s", s, "STOP", substr(conditionMessage(a), 1, 100), "\n"); next }
  cat(sprintf("m %-6g s %d rounds %d capped %s conv %d log sd %.3f (c %.3f) grad %.3g slopeSE %.4f ll-c %+.4f warns %d %s\n", m, s, a$importance$rounds,
      a$importance$capped, a$opt$convergence, a$opt$par[[2]], r$opt$par[[2]], a$importance$grad, ses(a)[[1]] / ses(r)[[1]],
      as.numeric(logLik(a)) - as.numeric(logLik(r)), length(w), paste(w, collapse = ";")))
}
