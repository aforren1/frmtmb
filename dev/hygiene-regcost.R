# F1 of the review: the cost of one registration under the new key, as
# 1000 re-registrations of one closure from a testthat test environment
# (frmtmb's namespace is its topenv()), with frmtmb.eam, frmtmb.ode and
# frmtmb.sample loaded so the list holds their checks too.
#   Rscript dev/hygiene-regcost.R <mine|base>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "mine") "C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (p in c("frmtmb.eam", "frmtmb.ode", "frmtmb.sample")) {
  suppressMessages(loadNamespace(p))
}
ns <- asNamespace("frmtmb")
te <- testthat::test_env("frmtmb")
f <- eval(quote(function(spec, frame) invisible(NULL)), te)
n0 <- length(ns$frmtmb_frame_checks$fns)
gc(); m0 <- sum(gc()[, 2])
t0 <- proc.time()[["elapsed"]]
for (i in 1:1000) frmtmb_register_frame_check(f)
t1 <- proc.time()[["elapsed"]]
gc(); m1 <- sum(gc()[, 2])
cat(sprintf("arm %s: 1000 registrations %.2f s (%.2f ms each), +%.1f MB, checks %d -> %d\n",
            arm, t1 - t0, t1 - t0, m1 - m0, n0,
            length(ns$frmtmb_frame_checks$fns)))
