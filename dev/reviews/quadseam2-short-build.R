# Short check of lane quadseam2's fix round: is the lane library's
# installed code the fix-round source? Each function of the named source
# files is compared, deparsed without srcrefs, with the installed
# namespace's copy.
#   Rscript quadseam2-short-build.R
lib <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-quadseam2"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "version",
    format(packageVersion("frmtmb")), "\n")
cat("latent:", dirname(find.package("frmtmb.latent")), "\n")
ns <- asNamespace("frmtmb")
for (f in c("R/fit.R", "R/autocor-latent.R", "R/se-check.R")) {
  ex <- parse(file.path(wt, f), keep.source = FALSE)
  n <- 0L; bad <- character()
  for (e in ex) {
    if (!is.call(e) || !identical(e[[1]], as.name("<-")) ||
          !is.name(e[[2]])) next
    if (!is.call(e[[3]]) || !identical(e[[3]][[1]], as.name("function"))) next
    nm <- as.character(e[[2]])
    src <- eval(e[[3]], baseenv())
    if (!exists(nm, ns, inherits = FALSE)) { bad <- c(bad, nm); next }
    inst <- utils::removeSource(get(nm, ns))
    n <- n + 1L
    if (!identical(deparse(formals(src)), deparse(formals(inst))) ||
          !identical(deparse(body(src)), deparse(body(inst)))) {
      bad <- c(bad, nm)
    }
  }
  cat(sprintf("BUILD %s functions %d differ %d %s\n", f, n, length(bad),
              paste(bad, collapse = ",")))
}
cat("END\n")
