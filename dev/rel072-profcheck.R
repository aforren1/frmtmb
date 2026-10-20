# (Copy of dev/rel071-profcheck.R for 0.72.0.) Do the functions dev/rel070-profile.R wraps keep the signatures it
# assumes on this build? A wrapper with a stale signature would log
# nothing or break the fit.
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
want <- list(se_check = c("fit", "control"), se_deferred_report = c("fit", "sdr"),
             se_report = c("fit", "lost", "act"), separation_check = "fit",
             check_convergence = c("fit", "control"),
             nl_flat_message = c("obj", "opt", "frame", "fit"),
             se_explained_pars = "fit", se_par_grouped = NULL, fit_assembled = NULL)
for (nm in names(want)) {
  f <- get0(nm, envir = ns)
  if (is.null(f)) { cat(nm, ": MISSING\n"); next }
  a <- names(formals(f))
  ok <- is.null(want[[nm]]) || identical(a[seq_along(want[[nm]])], want[[nm]])
  cat(sprintf("%-20s %s  (%s)\n", nm, if (ok) "ok" else "CHANGED", paste(a, collapse = ", ")))
}
