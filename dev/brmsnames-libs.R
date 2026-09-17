## Library order for lane brmsnames, one place so no script drifts.
## `base` reads the shared reference build of aa9227e (read-only);
## `lane` puts this lane's private library first. pinlib must beat the
## user library, or a fresh Stan compile dies (dev/lane-rules.md).
brmsnames_libs <- function(arm) {
  pin <- "C:/Users/adf44/source/r/pinlib"
  usr <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
  ref <- "C:/Users/adf44/source/r/rellib-r3"
  own <- "C:/Users/adf44/source/r/brmsnames-lib"
  .libPaths(switch(arm,
                   base = c(ref, pin, usr),
                   lane = c(own, ref, pin, usr),
                   stop("arm must be base or lane")))
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               mustWork = FALSE))
  invisible(.libPaths())
}
