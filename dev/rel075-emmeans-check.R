# Release review finding 1: emmeans on a fit with a part this session
# cannot use, and on an intact fit, on one build. The broken fit is
# built as test-extension-home-mixed.R builds it (its family said to be
# owned by a package no session has). The intact outputs are saved, so
# two builds can be compared with identical().
#   Rscript dev/rel075-emmeans-check.R <lib> <out.rds>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans)})
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
set.seed(20261024)
d <- data.frame(g = factor(rep(1:12, each = 8)),
                f = factor(rep(c("a", "b", "c"), 32)), x = rnorm(96))
d$y <- 1 + 0.5 * d$x + c(0, 0.4, -0.3)[d$f] + rnorm(12, 0, 0.6)[d$g] +
  rnorm(96, 0, 0.5)
fit <- frm(y ~ f + x + (1 | g), family = gaussian(), data = d)
intact <- list(
  ref_grid = summary(ref_grid(fit)),
  emmeans = summary(emmeans(fit, ~ f)),
  pairs = summary(pairs(emmeans(fit, ~ f))),
  slope = summary(emtrends(fit, ~ 1, var = "x")))
print(intact$emmeans)
broken <- fit
fam <- broken$spec$responses[[1L]]$family
attr(fam, "frmtmb_package") <- "frmtmbabsenttestpkga"
broken$spec$responses[[1L]]$family <- fam
broken$frame$spec <- broken$spec
broken$cache <- new.env()
for (nm in c("ref_grid", "emmeans")) {
  e <- tryCatch({
    if (nm == "ref_grid") ref_grid(broken) else emmeans(broken, ~ f)
    NULL
  }, error = identity)
  cat("BROKEN", nm, "::", if (is.null(e)) "OK" else
    paste(class(e)[1], "::", substr(conditionMessage(e), 1, 160)), "\n")
}
saveRDS(intact, a[2])
