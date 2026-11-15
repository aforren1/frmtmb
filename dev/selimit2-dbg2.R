# Lane selimit2: why the asym_laplace fitted() warning keeps its old words.
#   Rscript dev/selimit2-dbg2.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(5)
g <- factor(rep(1:10, each = 10))
x <- stats::rnorm(100)
dd <- data.frame(g = g, x = x, y = 1 + 0.5 * x + stats::rnorm(100))
f <- suppressWarnings(frm(y ~ x, data = dd, family = asym_laplace()))
jc <- ns$get_joint_cov(f)
str(jc[c("lost_pos", "lost_why", "null", "limit_range")])
print(diag(jc$V))
trace("se_pred_warn", quote({str(bad); print(vapply(sys.calls(), function(cl) paste(deparse(cl[[1L]]), collapse = ""), ""))}), print = FALSE, where = ns)
r <- withCallingHandlers(fitted(f), warning = function(w) {
  cat("W:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
})
