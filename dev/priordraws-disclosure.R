# Does frm_sample()'s default-prior disclosure say "b (flat)" when the
# call gives class "b" a prior? On rellib-r14 (the base) only.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + 0.5 * d$x + rnorm(40)
withCallingHandlers(
  invisible(suppressWarnings(frm_sample(
    bf(y ~ x), data = d, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "b"), chains = 1,
    iter = 100, seed = 1, refresh = 0))),
  message = function(m) {
    cat("MESSAGE:", conditionMessage(m))
    invokeRestart("muffleMessage")
  })
