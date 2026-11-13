# Lane arspat, backlog 4320: what reaches frmtmb when marginaleffects'
# predictions(newdata =) runs on an ar(cov = FALSE) fit.
# Output: dev/arspat-out/4320.txt
lib <- Sys.getenv("ARSPAT_LIB", "C:/Users/adf44/source/r/wt-arspat-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(marginaleffects)
})
cat("marginaleffects", format(packageVersion("marginaleffects")), "\n")
set.seed(2)
dd <- data.frame(t = rep(1:20, 3), h = factor(rep(1:3, each = 20)),
                 x = rnorm(60))
dd$y <- dd$x + rnorm(60)
f <- frm(y ~ x + ar(time = t, gr = h), data = dd)
nd <- dd[c(1, 2, 21), ]
seen <- list()
trace(frmtmb:::predict.frmtmb_fit, quote({
  nd_ <- get("newdata", envir = parent.frame())
  if (!is.null(nd_)) {
    assign("seen", c(get("seen", envir = globalenv()),
                     list(nd_[, intersect(c("t", "h", "x", "y"), names(nd_)),
                              drop = FALSE])), envir = globalenv())
  }
}), print = FALSE)
r <- tryCatch(predictions(f, newdata = nd), error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
untrace(frmtmb:::predict.frmtmb_fit)
cat("predict() calls with newdata:", length(seen), "\n")
for (s in seen) print(utils::head(s, 8))
cat("\nfrmtmb's own predict on nd:\n")
print(tryCatch(predict(f, newdata = nd), error = function(e) conditionMessage(e)))
cat("\nfitted on nd:\n")
print(tryCatch(fitted(f, newdata = nd), error = function(e) conditionMessage(e)))
