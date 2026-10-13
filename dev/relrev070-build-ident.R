# Reviewer of the 0.70.0 consolidation: is rellib-r8 built from this tree?
# Every function a package's R/ files define at top level is compared,
# deparsed without source references, with the installed namespace's.
#   Rscript dev/relrev070-build-ident.R
.libPaths(c(Sys.getenv("RELREV_LIB", "C:/Users/adf44/source/r/rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
T <- "C:/Users/adf44/source/r/frmtmb-wt-release"
pk <- c(frmtmb = ".", frmtmb.sample = "extensions/frmtmb.sample",
        frmtmb.eam = "extensions/frmtmb.eam",
        frmtmb.coupling = "extensions/frmtmb.coupling",
        frmtmb.latent = "extensions/frmtmb.latent",
        frmtmb.learn = "extensions/frmtmb.learn",
        frmtmb.ode = "extensions/frmtmb.ode",
        frmtmb.spline = "extensions/frmtmb.spline")
dp <- function(f) {
  f <- utils::removeSource(f)
  paste(deparse(f, control = c("keepInteger", "keepNA", "niceNames")),
        collapse = "\n")
}
for (p in names(pk)) {
  ns <- asNamespace(p)
  files <- list.files(file.path(T, pk[[p]], "R"), "[.][Rr]$",
                      full.names = TRUE)
  n <- 0L; same <- 0L; diff <- character(0); absent <- character(0)
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (!is.call(e) || !(identical(e[[1]], as.name("<-")) ||
                           identical(e[[1]], as.name("="))) ||
          !is.name(e[[2]])) next
      rhs <- e[[3]]
      if (!is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
      nm <- as.character(e[[2]])
      n <- n + 1L
      if (!exists(nm, envir = ns, inherits = FALSE)) {
        absent <- c(absent, nm); next
      }
      g <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
      if (!is.function(g)) { absent <- c(absent, nm); next }
      fs <- eval(rhs, envir = baseenv())
      if (identical(dp(fs), dp(g))) same <- same + 1L else
        diff <- c(diff, nm)
    }
  }
  cat(sprintf("%-16s %s: top-level functions %d, identical %d, differ %d (%s), not in namespace %d (%s)\n",
              p, as.character(packageVersion(p)), n, same, length(diff),
              paste(head(diff, 8), collapse = " "), length(absent),
              paste(head(absent, 8), collapse = " ")))
}
