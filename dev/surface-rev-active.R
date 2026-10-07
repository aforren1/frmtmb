# Reviewer of lane surface: count the active bindings among each
# namespace's exports (the lane-rules.md count), per arm.
#   Rscript dev/surface-rev-active.R lane|base
source("dev/surface-rev-env.R")
rev_env(commandArgs(TRUE)[1])
for (p in c("frmtmb", "frmtmb.sample")) {
  ns <- asNamespace(p)
  ex <- getNamespaceExports(p)
  act <- ex[vapply(ex, function(n) {
    exists(n, envir = ns, inherits = FALSE) && bindingIsActive(n, ns)
  }, NA)]
  # a re-export lives in the imports env, not the namespace itself
  imp <- parent.env(ns)
  act_imp <- ex[vapply(ex, function(n) {
    !exists(n, envir = ns, inherits = FALSE) &&
      exists(n, envir = imp, inherits = FALSE) && bindingIsActive(n, imp)
  }, NA)]
  cat(p, format(packageVersion(p)), "| exports", length(ex),
      "| active in namespace", length(act),
      "| active re-exports via imports", length(act_imp), "\n")
  cat("  own:", sort(act), "\n")
}
cat("frm_generic_owners", length(frmtmb:::frm_generic_owners), "\n")
cat("sample_generic_owners",
    length(frmtmb.sample:::sample_generic_owners), "\n")
