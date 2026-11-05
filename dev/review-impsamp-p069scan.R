.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib","C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
res <- integer(0)
for (s in 9:120) {
  d <- E$eq_dat("poisson", seed = s)
  lap <- suppressMessages(suppressWarnings(frm(bf(y ~ x + (1 + x | g)), family = poisson(), data = d)))
  v <- frmtmb:::varcorr_matrices(lap)[[1L]]
  r <- v[1, 2] / sqrt(v[1, 1] * v[2, 2])
  if (is.finite(r) && abs(r) > 0.99) {
    cat(s, r, sqrt(diag(v)), "\n")
    res <- c(res, s)
  }
}
cat("BOUNDARY", res, "\n")
