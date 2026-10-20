# Reviewer scratch (fix round 1): what the spline rename looks like to a
# user of frm_compat(), and how core names its own argument-setting rows.
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
loadNamespace("frmtmb.spline")
cat("spline from:", dirname(find.package("frmtmb.spline")), "\n")
ft <- frm_compat_features()
print(ft[ft$kind == "method", c("name", "key")], row.names = FALSE)
try1 <- function(expr) tryCatch(expr, error = function(e)
  paste("ERROR:", substr(conditionMessage(e), 1, 160)))
r <- try1(frm_compat("frm_curve(contrast =)"))
cat("new name rows:", if (is.data.frame(r)) nrow(r) else r, "\n")
r <- try1(frm_compat("frm_curve(contrast =)", "s()"))
cat("new name x s():", if (is.data.frame(r)) r$status else r, "\n")
r <- try1(frm_compat("frm_curve_contrast"))
cat("old name:", if (is.data.frame(r)) paste(nrow(r), "rows") else r, "\n")
r <- try1(frm_compat("frm_curve(contrast=)"))
cat("no-space spelling:", if (is.data.frame(r)) paste(nrow(r), "rows") else r,
    "\n")
