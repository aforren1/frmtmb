# A mixture of a core and an extension component, saved and read where
# the extension is absent. Usage: <save|absent>
arm <- commandArgs(trailingOnly = TRUE)[1]
first <- if (arm == "absent") "C:/Users/adf44/source/r/rev-splitsurv-lib" else
  "C:/Users/adf44/source/r/wt-splitsurv-lib"
.libPaths(c(first, "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
O <- "C:/Users/adf44/source/r/frmtmb-wt-splitsurv/dev/review-splitsurv-out/mix.rds"
show <- function(lab, expr) {
  r <- tryCatch({ v <- suppressWarnings(suppressMessages(expr)); "OK" },
    error = function(e) paste0(if (inherits(e, "frmtmb_extension_missing"))
      "REFUSED[ext] " else "ERROR ", substr(conditionMessage(e), 1, 140)))
  cat(sprintf("  %-12s %s\n", lab, r))
}
if (arm == "save") {
  suppressMessages(library(frmtmb.survival))
  set.seed(4)
  dd <- data.frame(x = rnorm(300))
  dd$time <- ifelse(runif(300) < 0.5, rlnorm(300, 0.5 * dd$x, 0.4),
                    rweibull(300, 1.5, 3))
  bk <- range(log(dd$time))
  fm <- tryCatch(suppressWarnings(frm(bf(time ~ x),
    family = mixture(lognormal(), royston_parmar(df = 1, knots = numeric(0),
                                                 bknots = bk)), data = dd)),
    error = function(e) conditionMessage(e))
  if (is.character(fm)) { cat("mixture fit failed:", fm, "\n"); quit() }
  cat("mixture fit logLik", format(as.numeric(logLik(fm)), digits = 10),
      "; family package:", frm_family_package(family(fm)),
      "; component packages:", attr(family(fm), "frmtmb_component_packages"),
      "\n")
  saveRDS(fm, O)
}
fm <- readRDS(O)
cat("arm", arm, "\n")
show("print", capture.output(print(fm)))
show("summary", summary(fm))
show("fitted", fitted(fm))
show("update", update(fm, . ~ 1))
show("logLik", logLik(fm))
