# Lane arspat: the frm_compat() rows of the pairs this lane touches.
# ARSPAT_LIB="" reads the base build alone. Output on stdout.
lib <- Sys.getenv("ARSPAT_LIB", "C:/Users/adf44/source/r/wt-arspat-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("frmtmb from", dirname(dirname(system.file(package = "frmtmb"))), "\n")
pairs <- list(c("ar()", "se()"), c("ma()", "se()"), c("arma()", "se()"),
              c("cosy()", "se()"), c("unstr()", "se()"), c("fcor()", "se()"),
              c("sar()", "se()"), c("ar()", "mi()"), c("cosy()", "mi()"),
              c("unstr()", "mi()"), c("fcor()", "mi()"), c("sar()", "mi()"),
              c("ar()", "weights()"), c("cosy()", "weights()"),
              c("sar()", "ar()"), c("ar()", "marginaleffects"),
              c("sar()", "marginaleffects"))
for (p in pairs) {
  r <- tryCatch(frm_compat(p[1], p[2]), error = function(e) {
    cat(p[1], "x", p[2], "ERROR", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(r)) next
  st <- if (is.data.frame(r)) r$status[1] else as.character(r)[1]
  nt <- if (is.data.frame(r) && !is.null(r$note)) substr(r$note[1], 1, 160) else ""
  cat(sprintf("%-8s x %-16s %-12s %s\n", p[1], p[2], st, nt))
}
