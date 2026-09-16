## The measured defect, against the base install (frmtmb 0.57.0 in the
## shared reference library). Run:
##   Rscript dev/argspell-defect.R
## Seed 2505 for the design, 25051 for the response draw; the same
## fixture as tests/testthat/test-arg-refusal.R.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(utils::packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n\n")

set.seed(2505)
dd <- data.frame(x = stats::rnorm(120),
                 g = factor(rep(1:8, each = 15)), y = 0)
dd$y <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                     newparams = list(Intercept = 1, x = 0.5, sigma = 0.4,
                                      sd_g__Intercept = 1.2),
                     nsim = 1, seed = 25051)[[1]]
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

cond <- fitted(fit)
swallowed <- suppressWarnings(fitted(fit, re_formula = NA))
pop <- suppressWarnings(predict(fit, re.form = NA, type = "response"))

cat(sprintf("fitted(fit)[1]                   = %.6f\n", cond[1]))
cat(sprintf("fitted(fit, re_formula = NA)[1]  = %.6f\n", swallowed[1]))
cat(sprintf("predict(re.form = NA)[1]         = %.6f\n", pop[1]))
cat(sprintf("identical(fitted, swallowed)     = %s\n",
            identical(unname(cond), unname(swallowed))))
cat(sprintf("max |fitted - population|        = %.6f\n",
            max(abs(cond - pop))))
cat(sprintf("misspelling swallowed silently   = %s\n",
            identical(unname(cond),
                      unname(suppressWarnings(fitted(fit, re_frmula = NA))))))
cat("\nranef structure:\n")
str(ranef(fit), max.level = 2)
cat("\nresponse names:", names(fit$spec$responses), "\n")
cat("dpars:", names(fit$spec$responses[[1]]$dpars), "\n")
