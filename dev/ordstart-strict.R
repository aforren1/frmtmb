# Lane ordstart, short check item 1: would a stricter "proper" test
# (an exact-Hessian eigenvalue floor, or a category-cover floor of 1e-8)
# reject the six flat points the review found without touching the fits
# 0.81.0 passed? On the 40-fit table, at the default start's point
# (starts = 0) and at the point the lane returns: the smallest eigenvalue
# of the exact outer Hessian, the smallest category cover of
# mixture_ord_degeneracy(), and the fit-end degenerate verdict.
#   Rscript dev/ordstart-strict.R out.tsv
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
fitw <- function(d, fam, ...) {
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x + z), family = fam, data = d, ...),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  H <- f$obj$he(f$opt$par)
  ev <- if (all(is.finite(H))) {
    min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  } else NA_real_
  dg <- frmtmb:::mixture_ord_degeneracy(f, "y")
  c(ll = as.numeric(logLik(f)), code = f$opt$convergence, eig = ev,
    cover = min(dg$cover, na.rm = TRUE),
    degen = any(grepl("ended at a degenerate boundary", w, fixed = TRUE)))
}
rows <- list()
for (nm in c("cum_cum", "cum_sratio")) for (s in 1:20) {
  fam <- switch(nm, cum_cum = mixture(cumulative(), cumulative()),
                cum_sratio = mixture(cumulative(), sratio()))
  d <- mk(s, mix = TRUE)
  a <- fitw(d, fam, control = frmtmb_control(starts = 0))
  b <- fitw(d, fam)
  rows[[length(rows) + 1L]] <- data.frame(
    family = nm, seed = s, ll0 = a[["ll"]], code0 = a[["code"]],
    eig0 = a[["eig"]], cover0 = a[["cover"]], degen0 = a[["degen"]],
    ll = b[["ll"]], code = b[["code"]], eig = b[["eig"]],
    cover = b[["cover"]], degen = b[["degen"]])
  cat(nm, s, "\n")
}
X <- do.call(rbind, rows)
utils::write.table(X, args[1], sep = "\t", quote = FALSE, row.names = FALSE)
