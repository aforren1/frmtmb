# Release review of 0.82.0, probe 4: ranef() and coef() of two terms on
# one factor, by block label. usage: probe4.R <lib>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
set.seed(11)
d <- data.frame(x = rnorm(120), z = rnorm(120),
                g = factor(rep(1:12, each = 10)),
                h = factor(rep(1:6, 20)))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] + 0.3 * d$x * rnorm(12)[d$g] +
  rnorm(6)[d$h] + rnorm(120)
show <- function(label, expr) {
  cat("\n##", label, "\n")
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cat("  W:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
    }, message = function(m) {
      cat("  M:", substr(conditionMessage(m), 1, 100), "\n")
      invokeRestart("muffleMessage")
    }),
    error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      cat("  call:", deparse1(conditionCall(e)), "\n")
      NULL
    })
  if (!is.null(v)) str(v, max.level = 2, give.attr = FALSE)
  invisible(v)
}
f2 <- frm(y ~ x + (1 | g) + (0 + x | g), data = d)
r <- show("ranef(f2)", ranef(f2))
show("names(ranef(f2))", names(r))
show("labels of blocks", vapply(f2$frame$re_blocks, function(b) {
  paste(b$label %||% NA, collapse = "/")
}, ""))
show("ranef(f2)$g", r$g)
show("ranef(f2)[['0 + x | g']]", r[["0 + x | g"]])
show("ranef(f2)[['1 | g']]", r[["1 | g"]])
show("coef(f2)", coef(f2))
show("coef(f2)$g", coef(f2)$g)
f3 <- frm(y ~ x + (1 | g) + (0 + x | g) + (1 | h), data = d)
show("ranef(f3) names", names(ranef(f3)))
show("coef(f3)$h", coef(f3)$h)
f4 <- frm(bf(y ~ x + (1 | g), sigma ~ (1 | g)), data = d)
show("ranef(f4)", ranef(f4))
show("ranef(f4)[['sigma: 1 | g']]", ranef(f4)[["sigma: 1 | g"]])
show("coef(f4)", coef(f4))
