# (copy of dev/rel072-active.R for 0.73.0, on the merged build)
# Active bindings among the exports of frmtmb and frmtmb.sample on the
# release build (dev/lane-rules.md, the local_mocked_bindings() rule),
# without and with brms attached.
#   Rscript dev/rel073-active.R [brms]
.libPaths(c("C:/Users/adf44/source/r/wt-release-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (identical(commandArgs(TRUE)[1], "brms")) suppressMessages(library(brms))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
for (p in c("frmtmb", "frmtmb.sample")) {
  ns <- asNamespace(p)
  ex <- getNamespaceExports(p)
  n <- sum(vapply(ex, function(x) {
    exists(x, envir = ns, inherits = FALSE) && bindingIsActive(x, ns)
  }, NA))
  cat(p, "active", n, "of", length(ex), "exports\n")
}
