# Release review of 0.74.0, item 3: frm_curve() on a gam-like fit, the
# user's path. Arm "old": rellib-r11 with library(frmtmb.spline) 0.10.2;
# arm "new": a library holding frmtmb 0.74.0 alone (no extension). Each
# arm saves its outputs; "cmp" compares them with identical().
#   Rscript dev/review-rel074-curve.R old|new <out.rds> [libs]
#   Rscript dev/review-rel074-curve.R cmp <old.rds> <new.rds>
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "cmp") {
  o <- readRDS(a[2]); n <- readRDS(a[3])
  for (k in union(names(o), names(n))) {
    x <- o[[k]]; y <- n[[k]]
    if (!identical(unclass(x), unclass(y)) && is.list(x)) for (nm in union(names(x), names(y))) if (!identical(x[[nm]], y[[nm]])) cat("   differs in", nm, class(x[[nm]])[1], if (is.numeric(x[[nm]])) max(abs(unlist(x[[nm]]) - unlist(y[[nm]])), na.rm = TRUE) else "", "
")
    cat(sprintf("%-22s identical %s | all.equal %s\n", k,
                identical(unclass(x), unclass(y)),
                isTRUE(all.equal(unclass(x), unclass(y), check.attributes = FALSE))))
  }
  quit(save = "no")
}
if (a[1] == "old") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r11", user))
  suppressMessages({library(frmtmb); library(frmtmb.spline)})
} else {
  .libPaths(c(a[3], user))
  suppressMessages(library(frmtmb))
}
cat("frmtmb", as.character(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "| loaded:", paste(grep("^frmtmb", loadedNamespaces(), value = TRUE),
                       collapse = " "), "\n")
set.seed(20261023)
n <- 300
dd <- data.frame(x = runif(n, 0, 10), g = factor(sample(letters[1:2], n, TRUE)))
dd$y <- sin(dd$x) + 0.3 * (dd$g == "b") + rnorm(n, 0, 0.4)
fit <- frm(bf(y ~ g + s(x, k = 10)), family = gaussian(), data = dd)
nd <- data.frame(x = seq(0.5, 9.5, length.out = 40), g = "a")
out <- list()
run <- function(k, expr) {
  out[[k]] <<- tryCatch(expr, error = function(e) paste("ERROR", conditionMessage(e)))
  r <- out[[k]]
  cat(k, ":", if (is.character(r) && length(r) == 1) r else class(r)[1], "\n")
}
run("logLik", as.numeric(logLik(fit)))
run("curve", frm_curve(fit, newdata = nd, seed = 1))
run("curve_trans", frm_curve(fit, newdata = nd, transform = TRUE, seed = 1))
run("deriv1", frm_curve_deriv(fit, var = "x", newdata = nd, seed = 1))
run("deriv2", frm_curve_deriv(fit, var = "x", order = 2L, newdata = nd, seed = 1))
run("feat_max", frm_curve_feature(fit, var = "x", type = "maximum",
                                  newdata = data.frame(x = seq(0.5, 3.5, length.out = 30), g = "a")))
run("feat_cross", frm_curve_feature(fit, var = "x", type = "crossing", at = 0,
                                    newdata = data.frame(x = seq(2, 4.5, length.out = 30), g = "a")))
run("contrast", frm_curve(fit, newdata = nd,
                          contrast = transform(nd, g = "b"), seed = 1))
print(out$feat_max)
saveRDS(out, a[2])
