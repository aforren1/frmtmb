## Reviewer probe: frm_sample() draws of y ~ mo(k), k in 1..4, saved by
## frmtmb 0.75.0 (rellib-r13), read with the lane library.
##   Rscript dev/review-mocode-draws.R base   # writes the rds
##   Rscript dev/review-mocode-draws.R lane   # reads it
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), arm, "lib:",
    dirname(dirname(system.file(package = "frmtmb.sample"))), "\n")
rds <- "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/review-mocode-log/old-draws.rds"
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
if (identical(arm, "base")) {
  ds <- suppressWarnings(suppressMessages(frm_sample(
    y ~ mo(k), data = d, chains = 1, iter = 400, seed = 1, refresh = 0)))
  saveRDS(ds, rds)
  cat("saved\n")
  quit("no")
}
ds <- readRDS(rds)
run <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  cat(sprintf("%-28s %s\n", lab, if (inherits(r, "error")) {
    paste("REFUSED/ERROR:", substr(gsub("\n", " ", conditionMessage(r)),
                                    1, 110))
  } else {
    v <- suppressWarnings(as.numeric(unlist(r)))
    paste("RETURNED:", paste(signif(head(v[is.finite(v)], 4), 6),
                             collapse = " "))
  }))
}
run("fixef()", fixef(ds))
run("summary()", summary(ds)$fixed)
run("variables()", length(variables(ds)))
run("posterior_epred()", posterior_epred(ds))
run("posterior_epred(newdata)",
    posterior_epred(ds, newdata = data.frame(k = 1:4)))
run("posterior_predict()", posterior_predict(ds))
run("log_lik()", log_lik(ds))
run("loo()", loo(ds)$estimates)
run("fitted()", fitted(ds))
run("conditional_effects()", conditional_effects(ds)[[1]]$estimate__)
run("bayes_R2()", bayes_R2(ds))
