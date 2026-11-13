# Short check of fix round 2: is the installed code the tree's code? Parse each package's
# R/ directory, compare every top-level function, without comments, to the namespace of
# rellib-r21. Also lists the source lines that differ between a comment-free parse.
#   Rscript dev/review-predread-r2-installed.R
L <- "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(L, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-rel082"
pk <- list(frmtmb = "R", frmtmb.sample = "extensions/frmtmb.sample/R", frmtmb.latent = "extensions/frmtmb.latent/R",
           frmtmb.coupling = "extensions/frmtmb.coupling/R", frmtmb.eam = "extensions/frmtmb.eam/R",
           frmtmb.learn = "extensions/frmtmb.learn/R", frmtmb.ode = "extensions/frmtmb.ode/R",
           frmtmb.spatial = "extensions/frmtmb.spatial/R", frmtmb.survival = "extensions/frmtmb.survival/R")
for (p in names(pk)) {
  suppressMessages(requireNamespace(p, quietly = TRUE))
  ns <- asNamespace(p)
  fs <- list.files(file.path(root, pk[[p]]), pattern = "[.]R$", full.names = TRUE)
  n <- 0L
  bad <- character()
  missing <- character()
  for (f in fs) {
    for (e in parse(f, keep.source = FALSE)) {
      if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
          is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
        nm <- as.character(e[[2L]])
        if (!exists(nm, envir = ns, inherits = FALSE)) {
          missing <- c(missing, nm)
          next
        }
        inst <- get(nm, envir = ns, inherits = FALSE)
        if (!is.function(inst)) next
        src <- eval(e[[3L]], envir = baseenv())
        a <- deparse(body(src), control = NULL)
        b <- deparse(body(inst), control = NULL)
        if (bindingIsActive(nm, ns)) next
        n <- n + 1L
        if (!identical(a, b) || !identical(names(formals(src)), names(formals(inst)))) bad <- c(bad, nm)
      }
    }
  }
  cat(sprintf("%-16s %s  functions compared %d  differ %d  absent %d  installed %s  version %s\n", p,
              dirname(find.package(p)), n, length(bad), length(missing),
              format(file.mtime(file.path(find.package(p), "DESCRIPTION")), "%m-%d %H:%M:%S"), packageVersion(p)))
  if (length(bad)) cat("   differ:", paste(utils::head(bad, 20), collapse = " "), "\n")
  if (length(missing)) cat("   absent:", paste(utils::head(missing, 20), collapse = " "), "\n")
}
