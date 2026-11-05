# Fix round 1: moved < imp_move_tol is in the coefficient's own units. A gaussian
# with y in large units: how often does an exact fit warn that it used all its
# rounds? Data seeds 1:8, y x 10^k.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
for (k in c(2, 3, 4, 5, 6, 7)) {
  capn <- 0L
  ex <- ""
  for (s in 1:8) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10))
    z <- rnorm(150)
    d <- data.frame(y = (0.3 + 0.5 * z + rnorm(15, 0, 0.7)[g] + rnorm(150, 0, 0.5)) * 10^k, z, g)
    f <- suppressWarnings(tryCatch(frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d, importance = 64L),
                                   error = function(e) e))
    l <- suppressWarnings(frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d))
    if (inherits(f, "error")) {
      ex <- paste0(ex, " s", s, " STOP")
      next
    }
    if (isTRUE(f$importance$capped)) {
      capn <- capn + 1L
      ex <- sprintf("%s s%d moves %s grad %.3g ll-lap %+.2g;", ex, s,
                    paste(format(f$importance$moves, digits = 2), collapse = ","), f$importance$grad,
                    as.numeric(logLik(f)) - as.numeric(logLik(l)))
    }
  }
  cat(sprintf("y x 1e%d: capped %d of 8.%s\n", k, capn, ex))
}
