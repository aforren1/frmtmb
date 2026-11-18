# Review of lane squant: residuals(type = "osa") on a right-censored fit
# errors for quantile_elf() (squant-review-integ.txt). Is it the family or
# every censored fit? The same data with gaussian() and asym_laplace().
# Run (worktree root), lane build and base:
#   Rscript dev/reviews/squant-review-osa.R lane > dev/reviews/squant-review-osa-lane.txt 2>&1
#   Rscript dev/reviews/squant-review-osa.R base > dev/reviews/squant-review-osa-base.txt 2>&1
arm <- commandArgs(TRUE)[1]
lp <- c("C:/Users/adf44/source/r/rellib-r22",
        "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") lp <- c("C:/Users/adf44/source/r/wt-squant-lib", lp)
.libPaths(lp)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(601)
n <- 300
d <- data.frame(x = runif(n))
d$y <- 1 + 2 * d$x + (0.5 + d$x) * rnorm(n)
d$ci <- 0
d$ci[d$y > 3.5] <- 1
d$ub <- pmin(d$y, 3.5)
fams <- list(gaussian = quote(gaussian()), asym_laplace = quote(asym_laplace()),
             student = quote(student()))
if (arm == "lane") fams$quantile_elf <- quote(quantile_elf())
for (nm in names(fams)) {
  fo <- if (nm %in% c("asym_laplace", "quantile_elf")) {
    bf(ub | cens(ci) ~ x, quantile = 0.5)
  } else bf(ub | cens(ci) ~ x)
  r <- tryCatch({
    f <- suppressWarnings(frm(fo, family = eval(fams[[nm]]), data = d))
    o <- residuals(f, type = "osa")
    sprintf("ok, %d rows, %d NA", nrow(o), sum(is.na(o[, 1])))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 200)))
  cat(nm, ":", r, "\n")
}
cat("DONE\n")
