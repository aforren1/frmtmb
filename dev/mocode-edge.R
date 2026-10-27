## mo() edge cases in frmtmb: ordered factors with unused levels, gaps,
## two values, one value, negative integers, newdata outside the range.
##   Rscript dev/mocode-edge.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("== frmtmb", format(packageVersion("frmtmb")), arm, "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
info <- function(lab, x, nd = NULL) {
  set.seed(3)
  d <- data.frame(x = x, y = rnorm(length(x)))
  if (is.numeric(x)) d$y <- d$y + x
  fit <- tryCatch(q(frm(y ~ mo(x), data = d)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-26s ERROR: %s\n", lab, conditionMessage(fit)))
    return(invisible())
  }
  mi <- fit$frame$linpreds[[1]]$mo[[1]]
  cat(sprintf("%-26s D = %d; codes = %s; min = %s; levels = %s\n", lab,
              mi$D, paste(mi$codes, collapse = " "), format(mi$min),
              paste(mi$levels, collapse = ",")))
  if (!is.null(nd)) {
    r <- tryCatch(q(fitted(fit, newdata = data.frame(x = nd))[, 1]),
                  error = function(e) e)
    cat(sprintf("%-26s newdata %s: %s\n", "", paste(nd, collapse = ","),
                if (inherits(r, "error")) paste("ERROR:", conditionMessage(r))
                else paste(signif(r, 6), collapse = " ")))
  }
}
xs <- rep(1:4, 10)
info("1..4", xs, nd = c(1, 4))
info("1..4, newdata 0", xs, nd = 0)
info("1..4, newdata 5", xs, nd = 5)
info("1..4, newdata 2.4", xs, nd = 2.4)
info("-2..1", xs - 3L, nd = c(-2, 1))
info("gaps 1,2,5", rep(c(1, 2, 5), 12), nd = c(1, 5))
info("two values 1,2", rep(1:2, 15))
info("two values 0,1", rep(0:1, 15))
info("one value 3", rep(3, 20))
info("ordered, unused level e",
     factor(rep(c("a", "b", "c", "d"), 10), levels = c("a", "b", "c", "d", "e"),
            ordered = TRUE))
info("ordered, unused first",
     factor(rep(c("b", "c", "d"), 10), levels = c("a", "b", "c", "d"),
            ordered = TRUE))
info("ordered, unused middle",
     factor(rep(c("a", "c", "d"), 10), levels = c("a", "b", "c", "d"),
            ordered = TRUE))
