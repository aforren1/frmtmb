# Lane wt-arcovsample: is pp_check(type = "loo_pit_overlay") broken for
# every model, or only for the one this lane changed? The before-run hit
# "One of 'lw' and 'psis_object' must be specified" on a cov = FALSE
# ARMA fit, where log_lik() was refused; a plain model tells which.
#   Rscript dev/arcovsample-ppcheck.R
LIB <- if (identical(Sys.getenv("ARCOVSAMPLE_REF"), "true")) {
  "C:/Users/adf44/source/r/rellib-r3"
} else {
  "C:/Users/adf44/source/r/wt-arcovsample-lib"
}
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("library      :", LIB, "\n")
cat("frmtmb       :", format(packageVersion("frmtmb")), "\n")
cat("frmtmb.sample:", format(packageVersion("frmtmb.sample")), "\n")

set.seed(31)
dd <- data.frame(x = stats::rnorm(60))
dd$y <- 1 + 0.5 * dd$x + stats::rnorm(60)
ds <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ x), family = gaussian(), data = dd, chains = 1,
             iter = 400, refresh = 0, seed = 3)))

for (ty in c("dens_overlay", "loo_pit_overlay", "loo_pit", "loo_intervals",
             "loo_ribbon")) {
  v <- tryCatch(suppressWarnings(suppressMessages(pp_check(ds, type = ty))),
                error = identity)
  cat("[", ty, "] ",
      if (inherits(v, "condition")) {
        paste0(class(v)[1L], ": ", conditionMessage(v))
      } else paste(class(v)[1L], "OK"), "\n", sep = "")
}
