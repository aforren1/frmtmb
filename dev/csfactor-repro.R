# Reproduce defect A (cs() on a factor fitted on integer codes) and
# probe defect B (y ~ x + cs(x)). Seed 405, the wt-predfix reviewer's
# construction from dev/predfix-review2/r2-cs2.R.
#
#   Rscript dev/csfactor-repro.R <lib>
# where <lib> is the library holding frmtmb (rellib-r3 for "before",
# the lane library for "after").
args <- commandArgs(TRUE)
LIB <- if (length(args)) args[[1L]] else
  "C:/Users/adf44/source/r/wt-csfactor-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")

set.seed(405)
n <- 500
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(x, fc, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$fb <- as.numeric(d$fc == "b")
d$fcc <- as.numeric(d$fc == "c")

cat("\n== A: cs(fc) on a factor ==\n")
ff <- tryCatch(frm(bf(yo ~ cs(fc)), family = sratio(), data = d),
               error = function(e) e)
if (inherits(ff, "error")) {
  cat("refused:", conditionMessage(ff), "\n")
} else {
  print(round(fixef(ff)[, 1:2], 4))
  cat("variables():\n")
  print(grep("^bcs", variables(ff), value = TRUE))
  nd3 <- data.frame(fc = factor(c("a", "b", "c")))
  cat("fitted at 3 rows (a, b, c):\n")
  print(round(fitted(ff, newdata = nd3)[, "Estimate", ], 4))
  cat("fitted at the single row fc = factor('c'):\n")
  print(round(fitted(ff, newdata = data.frame(fc = factor("c")))[,
                                                                "Estimate", ],
              4))
  cat("empirical shares by level:\n")
  print(round(prop.table(table(d$fc, d$yo), 1), 4))
  fd <- frm(bf(yo ~ cs(fb) + cs(fcc)), family = sratio(), data = d)
  cat(sprintf("logLik cs(fc) %.4f | hand dummies %.4f | diff %.4f\n",
              as.numeric(logLik(ff)), as.numeric(logLik(fd)),
              as.numeric(logLik(ff)) - as.numeric(logLik(fd))))
  cat("hand-dummy fitted at 3 rows (a, b, c):\n")
  nd3d <- data.frame(fb = c(0, 1, 0), fcc = c(0, 0, 1))
  print(round(fitted(fd, newdata = nd3d)[, "Estimate", ], 4))
  cat("hand-dummy variables():\n")
  print(grep("^bcs", variables(fd), value = TRUE))
  cat("default_prior() rows for the factor fit:\n")
  print(default_prior(ff))
}

cat("\n== A2: cs() on a character column ==\n")
dch <- d
dch$fch <- as.character(d$fc)
print(tryCatch({
  fh <- frm(bf(yo ~ cs(fch)), family = sratio(), data = dch)
  c(logLik = as.numeric(logLik(fh)),
    ncoef = length(grep("^bcs", variables(fh))))
}, error = function(e) conditionMessage(e)))

cat("\n== A3: cs() on a numeric factor (levels 1, 2, 3) ==\n")
dnf <- d
dnf$fnum <- factor(as.integer(d$fc))
print(tryCatch({
  fn <- frm(bf(yo ~ cs(fnum)), family = sratio(), data = dnf)
  c(logLik = as.numeric(logLik(fn)),
    ncoef = length(grep("^bcs", variables(fn))))
}, error = function(e) conditionMessage(e)))

cat("\n== B: y ~ x + cs(x) ==\n")
fb <- tryCatch(frm(bf(yo ~ x + cs(x)), family = sratio(), data = d),
               error = function(e) e)
if (inherits(fb, "error")) {
  cat("refused:", conditionMessage(fb), "\n")
} else {
  print(round(fixef(fb), 6))
  cat("logLik", sprintf("%.6f", as.numeric(logLik(fb))), "\n")
  cat("max |se|", sprintf("%.6g", suppressWarnings(max(fixef(fb)[, 2],
                                                       na.rm = TRUE))), "\n")
}
cat("\ndone\n")
