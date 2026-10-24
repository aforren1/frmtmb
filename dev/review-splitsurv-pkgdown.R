.libPaths(c("C:/Users/adf44/source/r/wt-splitsurv-lib","C:/Users/adf44/source/r/rellib-r11","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("pkgdown", format(packageVersion("pkgdown")), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-splitsurv"
for (p in c(wt, file.path(wt, "extensions/frmtmb.survival"))) {
  cat("== check_pkgdown", basename(p), "\n")
  r <- tryCatch({pkgdown::check_pkgdown(p); "ok"}, error = function(e) conditionMessage(e))
  print(r)
}
dst <- file.path(wt, "dev/review-splitsurv-out/redir")
unlink(dst, recursive = TRUE)
pkg <- pkgdown::as_pkgdown(wt, override = list(destination = dst))
r <- tryCatch({pkgdown:::build_redirects(pkg); "ok"}, error = function(e) conditionMessage(e))
print(r)
f <- list.files(file.path(dst, "frmtmb.spline"), recursive = TRUE)
print(f)
cat(readLines(file.path(dst, "frmtmb.spline/index.html")), sep = "\n")
cat(readLines(file.path(dst, "reference/cox_baseline.html")), sep = "\n")
