# Lane eamhier: what happens when a caller asks ndt_time() for a
# standard error.
#
# This row needed a population non-decision time in SECONDS with an
# interval on it, and the harness had to build one by hand out of
# predict(dpar = "ndt", se.fit = TRUE) and the bound. That is worth a
# measurement rather than an opinion: either the package already
# supports it, or it refuses by name, or it does something worse.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-ndttime-se.R
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

d <- eamhier_data(20260910L, "A", ns = 6L, nt = 80L)
key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]

for (grp in c(TRUE, FALSE)) {
  form <- eamhier_form("A", if (grp) "pg" else "gl")
  fit <- frm(form, family = wiener(), data = d, se = TRUE)
  tag <- if (grp) "ndt_group(s)" else "global bound"
  cat("\n== ", tag, " ==\n", sep = "")
  cat("ndt_time(newdata) : ",
      paste(formatC(as.numeric(ndt_time(fit, newdata = key)),
                    digits = 5, format = "f"), collapse = " "), "\n")
  r <- tryCatch(ndt_time(fit, newdata = key, se.fit = TRUE),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat("ndt_time(se.fit = TRUE) ERROR: ", conditionMessage(r), "\n")
  } else {
    cat("ndt_time(se.fit = TRUE) returned a ", class(r)[1L],
        " of length ", length(r), ":\n", sep = "")
    print(utils::str(r))
  }
  # what the harness does instead, and what the delta method says it
  # should be: the bound is DATA, so the standard error of a time is
  # the standard error of the fraction times that bound exactly
  p <- suppressWarnings(
    stats::predict(fit, newdata = key, dpar = "ndt", type = "response",
                   se.fit = TRUE))
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  fl <- if (is.null(bd[["floors"]])) bd[["ub"]] else bd[["floors"]]
  cat("by hand, time  : ",
      paste(formatC(as.numeric(p$fit) * fl, digits = 5, format = "f"),
            collapse = " "), "\n")
  cat("by hand, se    : ",
      paste(formatC(as.numeric(p$se.fit) * fl, digits = 5,
                    format = "f"), collapse = " "), "\n")
}
