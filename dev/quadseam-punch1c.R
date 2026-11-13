# Item 828: under quadrature = TRUE the stationary-point escape
# re-optimizes a Gauss-Kronrod tape frozen at the Laplace stall. Port of
# dev/skewinit-punch1c.R (its first half) to the round's libraries: the
# data are wt-skewinit's acceptance seed 1 (make_data(1, FALSE) of
# dev/skewinit-lib.R, inlined) plus a 20-level (1 | g).
#   Rscript dev/quadseam-punch1c.R <lib|base>
# log dev/quadseam-log/punch1c-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
op <- options(digits = 12)
skew <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
set.seed(1)
n <- 200
xs <- -abs(rnorm(n)) * 3
y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
dd <- data.frame(y = y, xs = xs)
set.seed(99)
dd$g <- factor(rep(1:20, length.out = nrow(dd)))
m3 <- skew(dd$y)
bad <- list(betad = c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3))
fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
fitq <- function(q, ...) {
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd,
                            quadrature = q, ...))
  attr(f, "secs") <- proc.time()[["elapsed"]] - t0
  f
}
for (q in c(FALSE, TRUE)) {
  d0 <- fitq(q)
  d1 <- fitq(q, start = bad)
  e <- d1$opt[["stationary_escape"]]
  cat(sprintf(paste0("ARM quadrature=%-5s default %16.9f  rescued %16.9f",
                     "  d %11.3e  esc %s  secs %.2f %.2f\n"),
              q, as.numeric(logLik(d0)), as.numeric(logLik(d1)),
              as.numeric(logLik(d1)) - as.numeric(logLik(d0)),
              if (is.null(e)) "-" else format(e[["gain"]], digits = 5),
              attr(d0, "secs"), attr(d1, "secs")))
  cat(sprintf("    alpha default %.6f rescued %.6f\n",
              d0$estimates$betad[2], d1$estimates$betad[2]))
  if (q) {
    # the rescued fit's tape at its own estimate: a fresh rule anchored
    # there, against the tape the fit carries
    g1 <- max(abs(d1$obj$gr(d1$opt$par)))
    cat(sprintf("    rescued: grad on its own tape %.3e  objective %.9f\n",
                g1, d1$opt$objective))
  }
  alt <- vapply(c(1, 3, 6), function(s) {
    as.numeric(logLik(fitq(q, start = list(betad = c(log(stats::sd(dd$y)),
                                                      s)))))
  }, 0)
  cat(sprintf("    best of starts 1,3,6 %16.9f  spread %11.3e\n",
              max(alt), max(alt) - min(alt)))
}
options(op)
cat("END\n")
