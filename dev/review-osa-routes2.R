# Review lane osa: student() and gaussian() under residual correlation
# terms, where the likelihood is not the rowwise loop that holds the "cdf"
# hook.
#   Rscript dev/review-osa-routes2.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-routes2", arm, "\n")
set.seed(502)
d <- expand.grid(t = 1:10, g = factor(1:15))
d <- d[order(d$g, d$t), ]
d$x <- rnorm(nrow(d))
e <- as.vector(replicate(15, arima.sim(list(ar = 0.5), 10, sd = 0.8)))
d$y <- 1 + 0.5 * d$x + e
d$yt <- 1 + 0.5 * d$x + as.vector(replicate(15, arima.sim(list(ar = 0.5), 10,
                                                          sd = 0.8,
                                                          rand.gen = function(n, ...) rt(n, 5) * 0.8)))
d$cnt <- rpois(nrow(d), exp(0.3 + 0.4 * d$x + e))
cases <- list(
  "gaussian ar(cov = TRUE)" = list(y ~ x + ar(t, g, cov = TRUE), gaussian()),
  "student ar(cov = TRUE)" = list(yt ~ x + ar(t, g, cov = TRUE), student()),
  "student ar()" = list(yt ~ x + ar(t, g), student()),
  "gaussian trunc ar(cov = TRUE)" = list(y | trunc(lb = -20) ~ x + ar(t, g, cov = TRUE), gaussian()),
  "poisson ar() (cov = FALSE)" = list(cnt ~ x + ar(t, g), poisson()))
for (nm in names(cases)) {
  f <- tryCatch(suppressWarnings(frm(cases[[nm]][[1]], family = cases[[nm]][[2]],
                                     data = d)), error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-32s FIT ERROR %s\n", nm, substr(conditionMessage(f), 1, 150)))
    next
  }
  r <- tryCatch(suppressWarnings(residuals(f, type = "osa")[, "Estimate"]),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-32s ERROR %s\n", nm, substr(gsub("\n", " ", conditionMessage(r)), 1, 220)))
  } else {
    cat(sprintf("%-32s nonfinite %d mean %.3f sd %.3f ks.p %.3g\n", nm,
                sum(!is.finite(r)), mean(r, na.rm = TRUE), sd(r, na.rm = TRUE),
                suppressWarnings(ks.test(r[is.finite(r)], "pnorm")$p.value)))
  }
}
cat("END review-osa-routes2\n")
