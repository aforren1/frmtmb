# Review of lane parity2, item 3922: conditional_effects() default effects
# on formulas the lane did not run, both arms. Each default display is
# then drawn; an error on a display is reported.
#   Rscript dev/reviews/parity2-review-ce.R base|lane
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/rellib-r22", libs[-1])
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, " lib:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/parity2/parity2-review-data.R")
try_ <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) cat(label, ": ERROR ", conditionMessage(r), "\n",
                                sep = "")
  for (x in unique(w)) cat("   warning: ", substr(x, 1, 160), "\n", sep = "")
  invisible(r)
}
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("\n## ", nm, "\n", sep = "")
  fit <- try_(paste(nm, "fit"), frm(s$f, family = s$fam, data = rd))
  if (inherits(fit, "error")) next
  ce <- try_(paste(nm, "ce"), conditional_effects(fit))
  if (!inherits(ce, "error")) {
    cat(nm, " ce: ", paste(names(ce), collapse = " "), "\n", sep = "")
  }
  if (!is.null(s$dpar)) {
    ce2 <- try_(paste(nm, "ce dpar"), conditional_effects(fit, dpar = s$dpar))
    if (!inherits(ce2, "error")) {
      cat(nm, " ce dpar=", s$dpar, ": ", paste(names(ce2), collapse = " "),
          "\n", sep = "")
    }
  }
}
