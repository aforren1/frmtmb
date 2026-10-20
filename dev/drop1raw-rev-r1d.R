# Reviewer, re-check of fix round 1, section D alone (section D of
# dev/drop1raw-rev-r1.R stopped on the reviewer's own control typo), and
# what getCall() of an update records on both builds.
# Arg: lane | base. Log: dev/drop1raw-rev-out/r1d-<arm>.txt (seed 9)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
ll <- function(f) as.numeric(logLik(f))
tr <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-50s ERROR %s\n", lab,
                substr(gsub("\n", " ", conditionMessage(r)), 1, 170)))
  } else if (inherits(r, "frmtmb_fit")) {
    cat(sprintf("%-50s nobs %d logLik %.10f\n", lab, nobs(r), ll(r)))
  } else {
    cat(sprintf("%-50s %s\n", lab, paste(format(r), collapse = " ")))
  }
  invisible(r)
}
set.seed(9)
m <- 240
e <- data.frame(x = rnorm(m), x3 = rnorm(m),
                g = factor(sample(sprintf("g%02d", 1:12), m, TRUE)))
e$z <- 1 + e$x + rnorm(12)[e$g] + rnorm(m)
f0 <- q(frm(bf(z ~ x + x3 + (1 | g)), data = e,
            control = frmtmb_control(check_se = "ignore")))
if (arm == "lane") {
  tr("update(fit): check_se of the refit", update(f0)$control$check_se)
  tr("update(fit, control = <stop>): check_se",
     update(f0, control = frmtmb_control(check_se = "stop"))$control$check_se)
  fb <- q(frm(bf(z ~ x), data = e))
  st <- fb$estimates
  fs <- local({
    st_l <- st
    q(frm(bf(z ~ x), data = e, start = st_l))
  })
  tr("start fit: update(fs)", update(fs))
  tr("start fit: update(fs, ~ . + x3)", update(fs, ~ . + x3))
  tr("reference frm(z ~ x + x3)", frm(bf(z ~ x + x3), data = e))
  tr("start fit: update(fs, family = student())",
     update(fs, family = student()))
  tr("reference frm(z ~ x, student())", frm(bf(z ~ x), family = student(),
                                           data = e))
  tr("start fit: update(fs, ~ . + x3, start = NULL)",
     update(fs, ~ . + x3, start = NULL))
}
u <- q(update(f0, ~ . - x3))
cat("getCall(update(fit, ~ . - x3)):", substr(deparse1(getCall(u)), 1, 300),
    "\n")
cat("nchar:", nchar(deparse1(getCall(u))), "\n")
u0 <- q(update(f0, family = student()))
cat("getCall(update(fit, family = student())):",
    substr(deparse1(getCall(u0)), 1, 300), "\n")
