# Release review of 0.74.0: a user who updates frmtmb to 0.74.0 and
# keeps an extension built for 0.73.0 (frmtmb.spline 0.10.2, which
# 0.74.0 retires; frmtmb.sample 0.21.0, whose floor is now 0.74.0).
# Library: rev-rel074-lib/mixed (frmtmb 0.74.0, survival 0.1.0 from
# rellib-r12; spline 0.10.2 and sample 0.21.0 from rellib-r11).
#   Rscript dev/review-rel074-stale.R <case>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/rev-rel074-lib/mixed",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("case", a[1], "| frmtmb", format(packageVersion("frmtmb")), "\n")
try_ <- function(label, expr) {
  r <- tryCatch({ v <- withCallingHandlers(expr, warning = function(w) {
    cat("  WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
    "ok" }, error = function(e) paste("ERROR", gsub("\\s+", " ", conditionMessage(e))))
  cat(sprintf("%-40s %s\n", label, substr(r, 1, 300)))
}
feat <- function() {
  f <- frm_compat_features()
  f[f$feature %in% c("frm_curve", "frm_curve_contrast", "royston_parmar",
                     "rp_floored", "cox"), , drop = FALSE]
}
if (a[1] == "spline") {
  before <- feat()
  try_("library(frmtmb.spline) 0.10.2", library(frmtmb.spline))
  cat("search():", paste(grep("frmtmb", search(), value = TRUE), collapse = " "), "\n")
  cat("frm_curve found in:", environmentName(environment(get("frm_curve"))), "\n")
  after <- feat()
  cat("compat features before:\n"); print(before)
  cat("compat features after:\n"); print(after)
  try_("then library(frmtmb.survival)", library(frmtmb.survival))
  cat("royston_parmar found in:", environmentName(environment(get("royston_parmar"))), "\n")
  set.seed(1)
  dd <- data.frame(x = runif(200, 0, 10)); dd$y <- sin(dd$x) + rnorm(200, 0, 0.3)
  fit <- frm(bf(y ~ s(x)), family = gaussian(), data = dd)
  try_("frm_curve() (whichever is first)", frm_curve(fit, newdata = data.frame(x = 1:5), seed = 1))
  try_("frmtmb::frm_curve()", frmtmb::frm_curve(fit, newdata = data.frame(x = 1:5), seed = 1))
}
if (a[1] == "sample") {
  try_("library(frmtmb.sample) 0.21.0", library(frmtmb.sample))
  cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "\n")
}
if (a[1] == "spline-survival") {
  try_("library(frmtmb.survival)", library(frmtmb.survival))
  try_("then library(frmtmb.spline) 0.10.2", library(frmtmb.spline))
  cat("royston_parmar found in:", environmentName(environment(get("royston_parmar"))), "\n")
}
