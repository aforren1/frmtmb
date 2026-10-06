# Reviewer of lane fixes, claim 4b: duplicate coef rows, a coef row with
# resp in a multivariate model, and frm_sample()'s use of a coef row.
#   Rscript dev/fixes-rev-prior2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(9090)
n <- 250
d <- data.frame(x = rnorm(n, 1, 1))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -0.2) + (u > 0.8) + (u > 1.8) + (u > 2.8)
d$w <- 1L + (u + rnorm(n) > 0.5) + (u + rnorm(n) > 1.5)
tr <- function(lab, expr) {
  r <- tryCatch({force(expr); "OK"}, error = function(e)
    paste("ERROR:", substr(conditionMessage(e), 1, 200)))
  cat(sprintf("%-40s %s\n", lab, r))
}
p2 <- set_prior("normal(0.7, 0.4)", class = "Intercept", coef = "2")
tr("duplicate coef = 2 rows", frm(y ~ x, family = cumulative(), data = d,
                                  prior = p2 + set_prior("normal(0, 1)",
                                                         class = "Intercept",
                                                         coef = "2")))
tr("mv, coef = 2 with resp = w",
   frm(bf(y ~ x) + bf(w ~ x), family = cumulative(), data = d,
       prior = set_prior("normal(0, 1)", class = "Intercept", coef = "2",
                         resp = "w")))
tr("mv, coef = 3 with resp = w (w has 2)",
   frm(bf(y ~ x) + bf(w ~ x), family = cumulative(), data = d,
       prior = set_prior("normal(0, 1)", class = "Intercept", coef = "3",
                         resp = "w")))
fit <- frm(y ~ x, family = cumulative(), data = d, prior = p2)
tr("frm_sample with a coef row", {
  ds <- suppressMessages(suppressWarnings(frm_sample(fit, chains = 1,
                                                     iter = 200, refresh = 0,
                                                     seed = 3)))
  print(round(posterior_summary(ds)[1:6, 1:2], 3))
  ds
})
