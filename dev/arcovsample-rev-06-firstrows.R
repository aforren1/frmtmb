# REVIEW script 06, claim 6: is `?sample-log_lik`'s sentence "with each
# group's first max(p, q) rows taking the unshifted mean" true of the
# code?
#
# brms 2.23.0 initializes Err to zero and fills Err[n + 1, i] from
# J_lag[n], so row 2 of a group already has Err[., 1] = err of row 1.
# With max(p, q) = 2 that row is PARTIALLY shifted, not unshifted. This
# measures which rows of each group actually carry no shift.
#
#   Rscript dev/arcovsample-rev-06-firstrows.R

LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))

set.seed(6161L)
ng <- 4L
dd <- do.call(rbind, lapply(seq_len(ng), function(i) {
  data.frame(g = factor(i, levels = seq_len(ng)), t = seq_len(7L))
}))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + as.numeric(arima.sim(list(ar = c(0.5, 0.25)),
                                                nrow(dd)))
dd <- dd[sample(nrow(dd)), ]
rownames(dd) <- NULL

report <- function(lab, form, p, q) {
  f <- frm(form, family = gaussian(), data = dd)
  d0 <- eval_dpars(f)
  d1 <- arma_cond_dpars(f, d0)
  sh <- abs(as.numeric(d1$y$mu) - as.numeric(d0$y$mu))
  ord <- order(dd$g, dd$t)
  pos <- unlist(lapply(split(seq_along(ord), dd$g[ord]), seq_along))
  shifted <- sh[ord] > 1e-12
  cat("\n---- ", lab, "  max(p, q) = ", max(p, q), "\n", sep = "")
  cat("  within-group position : ",
      paste(sprintf("%2d", pos[1:8]), collapse = " "), " ...\n", sep = "")
  cat("  shift nonzero        : ",
      paste(sprintf("%2s", ifelse(shifted[1:8], "Y", ".")),
            collapse = " "), " ...\n", sep = "")
  tab <- tapply(shifted, pos, function(v) sum(v))
  cnt <- tapply(shifted, pos, length)
  cat("  groups with a nonzero shift, by within-group position:\n")
  for (i in names(tab)) {
    cat("    position ", i, ": ", tab[[i]], " of ", cnt[[i]], "\n",
        sep = "")
  }
  cat("  rows per group with NO shift: ",
      paste(sort(unique(as.integer(names(tab))[vapply(tab, function(v)
        v == 0L, TRUE)])), collapse = ","), "\n", sep = "")
}

report("ar(t, g, p = 1)", bf(y ~ x + ar(t, g, p = 1)), 1L, 0L)
report("ar(t, g, p = 2)", bf(y ~ x + ar(t, g, p = 2)), 2L, 0L)
report("ma(t, g, q = 2)", bf(y ~ x + ma(t, g, q = 2)), 0L, 2L)
report("arma(t, g, p = 2, q = 2)", bf(y ~ x + arma(t, g, p = 2, q = 2)),
       2L, 2L)
report("arma(t, g, p = 3, q = 1)", bf(y ~ x + arma(t, g, p = 3, q = 1)),
       3L, 1L)
cat("\nDONE\n")
