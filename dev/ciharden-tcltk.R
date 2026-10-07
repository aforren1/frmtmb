# Which call of frmtmb.sample's test-compat-preflight.R:95 loads tcltk,
# whose load on a headless Linux runner warns "no DISPLAY variable so
# Tk is not available"? Trace the namespace loads around the call.
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.sample))
before <- loadedNamespaces()
setHook(packageEvent("tcltk", "onLoad"), function(...) {
  cat("tcltk loading; call stack:\n")
  print(sys.calls())
})
trace(loadNamespace, quote(if (identical(package, "tcltk") ||
                              identical(package, "frmtmb.ode"))
  cat("loadNamespace(", package, ")\n")), print = FALSE)
pf <- frmtmb.sample:::sample_preflight
rows <- data.frame(feature_a = "frm_sample", kind_a = "method",
                   feature_b = "frm_ode()", kind_b = "special",
                   status = "refused", note = "x", stringsAsFactors = FALSE)
form <- bf(y ~ frm_ode(dyn, init = list(d), times = t), a ~ 1, nl = TRUE)
r <- try(pf(form, rows = rows))
cat("new namespaces:", setdiff(loadedNamespaces(), before), "\n")
