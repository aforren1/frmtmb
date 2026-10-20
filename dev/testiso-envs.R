# Every top-level environment in the eight namespaces: the candidates for
# session state a test can leave changed. Lane testiso's audit list.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("frmtmb", "frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
            "frmtmb.learn", "frmtmb.ode", "frmtmb.sample", "frmtmb.spline")) {
  ns <- asNamespace(p)
  for (nm in ls(ns, all.names = TRUE)) {
    if (bindingIsActive(nm, ns)) next
    x <- get(nm, envir = ns)
    if (is.environment(x) && !isNamespace(x)) {
      el <- ls(x, all.names = TRUE)
      cat(p, nm, length(el), paste(head(el), collapse = ","), "\n")
    }
  }
}
