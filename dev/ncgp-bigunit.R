# M2 outside gp: a large-unit non-gp fit, before and after. Standard errors over the
# scale, on the fit's own unit. Seeds in place. In %a, so a diff reads digit for digit.
#   Rscript dev/ncgp-bigunit.R <lib or base>
a <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] != "base") libs <- c(a[1], libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("# lib:", dirname(find.package("frmtmb")), "\n")
set.seed(15)
g <- gl(12, 8)
z <- rnorm(96)
v <- 2 + 0.5 * z + rnorm(12, 0, 0.7)[g] + 0.5 * rnorm(96)
set.seed(16)
u2 <- rnorm(12, 0, 0.3)
s <- 0.5 + u2[g] * z + 0.3 * rnorm(96)
rows <- list(
  ri = function(c0) frm(bf(y ~ z + (1 | g)), data = data.frame(y = v * c0, z = z, g = g)),
  rs = function(c0) frm(bf(y ~ z + (z | g)), data = data.frame(y = s * c0, z = z, g = g)),
  lm = function(c0) frm(bf(y ~ z), data = data.frame(y = v * c0, z = z)),
  poly = function(c0) frm(bf(y ~ poly(z, 3) + (1 | g)), data = data.frame(y = v * c0, z = z, g = g))
)
for (nm in names(rows)) for (c0 in c(1, 1e3, 1e6, 1e9)) {
  w <- character()
  f <- withCallingHandlers(rows[[nm]](c0), warning = function(cn) {
    w <<- c(w, substr(conditionMessage(cn), 1, 50))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  se <- sqrt(diag(as.matrix(vcov(f, full = TRUE))))
  se[seq_along(f$estimates$beta)] <- se[seq_along(f$estimates$beta)] / c0
  cat(sprintf("BU %-5s c %.0e u %.4g conv %d ll+nlogc %.8f se %s | w %d %s\n", nm, c0,
              (f$par_units %||% 1)[1], f$opt$convergence, as.numeric(logLik(f)) + 96 * log(c0),
              paste(sprintf("%.6g", se), collapse = " "), length(w), paste(unique(w), collapse = ";")))
}
