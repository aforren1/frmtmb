.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
y <- 3; x <- (y - 0.2) / 1.3
pieces <- list(
  l = function(p) ns$gev_parts(y, list(mu = 0.2, sigma = 1.3, xi = p))$l,
  h = function(p) ns$gev_parts(y, list(mu = 0.2, sigma = 1.3, xi = p))$h,
  out = function(p) 1e8 * ns$gev_parts(y, list(mu = 0.2, sigma = 1.3, xi = p))$out,
  eh = function(p) ns$exp_neg_capped(ns$gev_parts(y, list(mu = 0.2, sigma = 1.3, xi = p))$h))
want <- c(l = x, h = -x^2 / 2, out = 0, eh = exp(-x) * x^2 / 2)
for (nm in names(pieces)) {
  tp <- RTMB::MakeTape(pieces[[nm]], 0)
  cat(sprintf("%-4s value %.17g grad %.17g want %.17g\n", nm, tp(0), tp$jacobian(0)[1, 1], want[[nm]]))
}
