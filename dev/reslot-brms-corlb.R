# Lane reslot: is brms 2.23.0's program for lb/ub on class cor valid Stan?
# rstan::stanc() parses and type-checks only; nothing is compiled.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
dd <- data.frame(g = factor(rep(1:8, each = 10)), x = rnorm(80))
dd$y <- rnorm(80)
for (pr in list(set_prior("lkj(2)", class = "cor", lb = 0),
                set_prior("", class = "cor", group = "g", lb = -0.5,
                          ub = 0.5),
                set_prior("lkj(2)", class = "cor"))) {
  code <- stancode(bf(y ~ x + (1 + x | g)), data = dd, prior = pr)
  r <- tryCatch({rstan::stanc(model_code = code); "parses"},
                error = function(e) conditionMessage(e))
  cat(format(pr$prior), format(pr$lb), format(pr$ub), ":",
      substr(gsub("\n", " ", r), 1, 300), "\n")
}
