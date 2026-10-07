# Lane optima, item 3: ord_log_interior()'s floor before (the minimum
# form on z = b - a) and after (the positive part of the gap a - b), on
# the tape: values for a >= b, on 100000 random pairs as
# dev/ordmix-p1-m4.R drew them (gaps 1e-12 to 50), plus touching and
# crossed pairs, with the gradient's finiteness.
#   Rscript dev/optima-floor-check.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
lil <- frmtmb:::log_inv_logit
old <- function(a, b) {
  z <- b - a
  cap <- -1e-300
  z <- 0.5 * (z + cap - abs(cap - z))
  -b + RTMB::logspace_sub(0 * z, z) + lil(a) + lil(b)
}
new <- function(a, b) frmtmb:::ord_log_interior(a, b, TRUE)
set.seed(20261006)
n <- 100000
b <- stats::runif(n, -40, 40)
gap <- 10^stats::runif(n, -12, log10(50))
a <- b + gap
ev <- function(f, a, b) {
  k <- length(a)
  tp <- RTMB::MakeTape(function(x) f(x[seq_len(k)], x[k + seq_len(k)]),
                       c(a, b))
  list(v = tp(c(a, b)),
       gfin = all(is.finite(RTMB::MakeTape(function(x) {
         sum(f(x[seq_len(k)], x[k + seq_len(k)]))
       }, c(a, b))$jacobian(c(a, b)))))
}
vo <- ev(old, a, b)
vn <- ev(new, a, b)
cat("a >= b, 100000 pairs: values differ on", sum(vo$v != vn$v),
    "; gradient finite old", vo$gfin, "new", vn$gfin, "\n")
for (case in c("touch", "cross 1e-8", "cross 0.5", "cross 5")) {
  bb <- c(-3, 0.2, 7.5)
  aa <- switch(case, touch = bb, "cross 1e-8" = bb - 1e-8,
               "cross 0.5" = bb - 0.5, "cross 5" = bb - 5)
  o <- ev(old, aa, bb)
  w <- ev(new, aa, bb)
  cat(sprintf("%-11s old: %s (gradient finite %s) | new: %s (gradient finite %s)\n",
              case, paste(format(o$v, digits = 5), collapse = " "), o$gfin,
              paste(format(w$v, digits = 5), collapse = " "), w$gfin))
}
