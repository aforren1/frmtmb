## Reviewer probe: the value emmeans and marginaleffects give a mo()
## variable they do not vary, on the seed-5 data of y ~ mo(k) + z.
##   Rscript dev/review-mocode-interop.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), arm, "\n")
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
fz <- suppressWarnings(suppressMessages(frm(y ~ mo(k) + z, data = d)))
show <- function(lab, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e)
  cat(sprintf("%-30s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR:", substr(conditionMessage(r), 1, 120)) else r))
}
show("emmeans ref_grid k", {
  rg <- emmeans::ref_grid(fz)
  paste(format(rg@grid$k), collapse = " ")
})
show("emmeans(~ z) emmean", {
  s <- summary(emmeans::emmeans(fz, ~ z))
  paste(format(s$emmean), "at z =", format(s$z))
})
show("emmeans(~ k) emmean", {
  s <- summary(emmeans::emmeans(fz, ~ k))
  paste(format(s$k), "->", format(signif(s$emmean, 6)), collapse = "; ")
})
show("marginaleffects datagrid k", {
  g <- marginaleffects::datagrid(model = fz)
  format(g$k)
})
show("predictions(datagrid())", {
  p <- marginaleffects::predictions(fz, newdata = "mean")
  paste("k =", format(p$k), "estimate =", format(p$estimate))
})
show("avg_slopes(variables = 'k')", {
  s <- marginaleffects::avg_slopes(fz, variables = "k")
  format(s$estimate)
})
show("avg_comparisons(k +1)", {
  s <- marginaleffects::avg_comparisons(fz, variables = "k")
  format(s$estimate)
})
