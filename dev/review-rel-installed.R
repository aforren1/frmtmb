# Release review 0.80.0: is the installed code (rellib-r19) the tree's?
# Stronger than dev/rel080-installed.R: functions are compared whole
# (formals with their defaults, and body), assignments by `=` count, and a
# top-level value that evaluates in base R alone (a constant such as
# me_net_tol or opt_basis_kappa) is compared by value.
#
#   Rscript dev/review-rel-installed.R
L <- "C:/Users/adf44/source/r/rellib-r19"
.libPaths(c(L, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
pk <- c(frmtmb = "R", frmtmb.sample = "extensions/frmtmb.sample/R",
        frmtmb.latent = "extensions/frmtmb.latent/R",
        frmtmb.coupling = "extensions/frmtmb.coupling/R",
        frmtmb.eam = "extensions/frmtmb.eam/R",
        frmtmb.learn = "extensions/frmtmb.learn/R",
        frmtmb.ode = "extensions/frmtmb.ode/R",
        frmtmb.spatial = "extensions/frmtmb.spatial/R",
        frmtmb.survival = "extensions/frmtmb.survival/R")
strip <- function(f) {
  f <- utils::removeSource(f)
  environment(f) <- baseenv()
  f
}
for (p in names(pk)) {
  suppressMessages(requireNamespace(p, quietly = TRUE))
  ns <- asNamespace(p)
  nf <- 0L; nv <- 0L
  bad <- character(); absent <- character()
  for (f in list.files(file.path(root, pk[[p]]), pattern = "[.][Rr]$",
                       full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) {
      if (!(is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") &&
              is.name(e[[2L]]))) next
      nm <- as.character(e[[2L]])
      rhs <- e[[3L]]
      if (!exists(nm, envir = ns, inherits = FALSE)) {
        absent <- c(absent, nm)
        next
      }
      if (bindingIsActive(nm, ns)) next
      inst <- get(nm, envir = ns, inherits = FALSE)
      if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) {
        nf <- nf + 1L
        src <- eval(rhs, baseenv())
        if (!is.function(inst) ||
            !identical(deparse(strip(src)), deparse(strip(inst)))) {
          bad <- c(bad, nm)
        }
      } else {
        v <- tryCatch(eval(rhs, new.env(parent = baseenv())),
                      error = function(e) NULL)
        if (is.null(v) || is.function(v)) next
        nv <- nv + 1L
        if (!identical(v, inst)) bad <- c(bad, paste0(nm, " (value)"))
      }
    }
  }
  cat(sprintf("%-16s %s %s  functions %d  values %d  differ %d  absent %d\n",
              p, dirname(find.package(p)), as.character(packageVersion(p)),
              nf, nv, length(bad), length(absent)))
  if (length(bad)) cat("  DIFFER:", paste(bad, collapse = ", "), "\n")
  if (length(absent)) cat("  ABSENT:", paste(unique(absent), collapse = ", "), "\n")
}
cat("INSTALLED CHECK DONE\n")
