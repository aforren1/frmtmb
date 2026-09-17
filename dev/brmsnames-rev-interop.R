## Reviewer, claim 7: third-party readers on a fit, base against lane.
##   Rscript dev/brmsnames-rev-interop.R base|lane|compare   data seed 3
arm <- commandArgs(trailingOnly = TRUE)[1L]
lib <- switch(arm, lane = c("C:/Users/adf44/source/r/brmsnames-lib",
                            "C:/Users/adf44/source/r/rellib-r3"),
              "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
out <- sprintf("dev/stan-cache/brmsnames-rev-interop-%s.rds", arm)
q <- function(e) suppressWarnings(suppressMessages(e))
if (arm != "compare") {
  q(library(frmtmb))
  set.seed(3)
  d <- data.frame(x = rnorm(200), f = factor(sample(letters[1:3], 200, TRUE)),
                  g = factor(rep(1:20, 10)))
  d$y <- 1 + d$x + rnorm(20)[d$g] + rnorm(200)
  fit <- q(frm(bf(y ~ x + f + (1 + x | g)), family = gaussian(), data = d))
  t1 <- function(e) tryCatch(q(e), error = function(err) paste("ERROR:", conditionMessage(err)))
  r <- list(
    find_parameters = t1(insight::find_parameters(fit)),
    get_parameters = t1(insight::get_parameters(fit)),
    get_variance = t1(insight::get_variance(fit)),
    get_statistic = t1(insight::get_statistic(fit)),
    find_random = t1(insight::find_random(fit)),
    emmeans = t1(as.data.frame(emmeans::emmeans(fit, ~ f))),
    avg_slopes = t1(as.data.frame(marginaleffects::avg_slopes(fit, variables = "x"))),
    summary_print = t1(capture.output(print(summary(fit)))),
    print = t1(capture.output(print(fit)))
  )
  saveRDS(r, out)
  for (k in names(r)) cat(k, ":", substr(paste(format(r[[k]]), collapse = " "), 1, 200), "\n")
} else {
  b <- readRDS(sprintf("dev/stan-cache/brmsnames-rev-interop-%s.rds", "base"))
  l <- readRDS(sprintf("dev/stan-cache/brmsnames-rev-interop-%s.rds", "lane"))
  for (k in names(b)) {
    cat(sprintf("%-16s identical %s\n", k, identical(b[[k]], l[[k]])))
    if (!identical(b[[k]], l[[k]])) print(utils::head(all.equal(b[[k]], l[[k]]), 5))
  }
}
