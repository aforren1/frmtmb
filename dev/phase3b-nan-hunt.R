# Where does the censored or contaminated Wiener likelihood give a NaN
# gradient? Tape each new piece and scan wild parameter values, the
# kind an optimizer step reaches. Output to stdout.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib2",
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE,
                                   export_all = TRUE))
q <- c(0.3, 0.45, 0.7, 1.0, 1.0, 1.0)
up <- c(0, 1, 1, 0, 0, 1)
pieces <- list(
  lS = function(p) sum(ddm_rt_lcdf2(q - p[3], p[1], exp(p[2]), 0.5)$lS),
  lF = function(p) sum(ddm_rt_lcdf2(q - p[3], p[1], exp(p[2]), 0.5)$lF),
  Fprob = function(p) sum(exp(ddm_rt_lcdf2(q - p[3], p[1], exp(p[2]),
                                           0.5)$lF)))
set.seed(1)
grid <- rbind(
  expand.grid(v = c(-60, -20, -5, -1, 0, 1e-9, 1, 5, 20, 60),
              la = c(-12, -8, -4, -2, -0.5, 0, 0.5, 1.5, 3, 6, 10),
              t0 = c(-5, 0, 0.1, 0.25, 0.2999, 0.29999999, 0.3, 0.5, 2)))
for (nm in names(pieces)) {
  tp <- RTMB::MakeTape(pieces[[nm]], c(0, 0, 0))
  bad <- 0
  worst <- NULL
  for (i in seq_len(nrow(grid))) {
    p <- as.numeric(grid[i, ])
    val <- tp(p)
    g <- tp$jacobian(p)
    if (!all(is.finite(g)) || !is.finite(val)) {
      bad <- bad + 1
      if (is.null(worst)) worst <- c(p, val, g)
    }
  }
  cat(nm, ": non-finite at", bad, "of", nrow(grid), "points")
  if (!is.null(worst)) cat("; first:", format(worst, digits = 4))
  cat("\n")
}
pieces$dens <- function(p) sum(ddm_lpdf_both(q - p[3], p[1], exp(p[2]), 0.5, up))
for (nm in names(pieces)) {
  tp <- RTMB::MakeTape(pieces[[nm]], c(0, 0, 0))
  ok <- vapply(seq_len(nrow(grid)), function(i) {
    p <- as.numeric(grid[i, ])
    all(is.finite(tp$jacobian(p))) && is.finite(tp(p))
  }, TRUE)
  cat(nm, "bad by la:", format(table(grid$la[!ok])), "\n")
  cat(nm, "bad by t0:", format(table(grid$t0[!ok])), "\n")
}
