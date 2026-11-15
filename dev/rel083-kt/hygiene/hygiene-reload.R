# Item 3154: does a reload of an extension duplicate its frame checks
# and prior-default providers? Loads frmtmb.eam, frmtmb.ode and
# frmtmb.sample from the base library, then pkgload::load_all()s the
# worktree copy of frmtmb.ode twice and calls frmtmb.sample's .onLoad()
# again, counting the two list registries after each step.
#
# Usage: Rscript dev/hygiene-reload.R <mine|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "mine") "C:/Users/adf44/source/r/rellib-r22" else
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
say <- function(step) {
  cat(sprintf("%-28s checks=%d providers=%d\n", step,
              length(ns$frmtmb_frame_checks$fns),
              length(ns$frmtmb_prior_defaults$providers)))
}
cat("arm", arm, "frmtmb from", dirname(find.package("frmtmb")), "\n")
say("core alone")
for (p in c("frmtmb.eam", "frmtmb.ode", "frmtmb.sample")) {
  suppressMessages(loadNamespace(p))
}
say("eam, ode, sample loaded")
ode_dir <- "dev/rel083-log/kt/hygiene/scratch/regkey-src-mine/frmtmb.ode"
suppressMessages(pkgload::load_all(ode_dir, quiet = TRUE,
                                   export_all = FALSE, helpers = FALSE))
say("load_all(frmtmb.ode) 1")
suppressMessages(pkgload::load_all(ode_dir, quiet = TRUE,
                                   export_all = FALSE, helpers = FALSE))
say("load_all(frmtmb.ode) 2")
sns <- asNamespace("frmtmb.sample")
sns$.onLoad(dirname(getNamespaceInfo("frmtmb.sample", "path")),
            "frmtmb.sample")
say("frmtmb.sample .onLoad() again")
