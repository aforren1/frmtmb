## Fix 3, run: predictive_interval() answers under either spelling and
## the two agree; pp_check() still refuses the alias.
.libPaths(c("C:/Users/adf44/source/r/asrev-lib2",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
set.seed(3)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- 0.3 + 0.5 * dd$x + rnorm(8, 0, 1.2)[dd$g] + rnorm(80, 0, 0.4)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- suppressMessages(frm_sample(fit, chains = 1, iter = 400,
                                  refresh = 0, seed = 7))
p <- function(lab, e) {
  r <- tryCatch(force(e), error = function(x) x)
  cat(sprintf("%-48s %s\n", lab,
              if (inherits(r, "error"))
                paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
              else paste("ok", paste(class(r), collapse = "/"))))
  invisible(r)
}
cat("---- ours: predictive_interval ----\n")
set.seed(9)
a <- p("predictive_interval(ds, re_formula = NA)",
       predictive_interval(ds, re_formula = NA, ndraws = 12))
set.seed(9)
b <- p("predictive_interval(ds, re.form = NA)",
       predictive_interval(ds, re.form = NA, ndraws = 12))
set.seed(9)
d <- p("predictive_interval(ds)  [default]",
       predictive_interval(ds, ndraws = 12))
cat("  identical(alias, primary)  : ", identical(a, b), "\n")
cat("  alias differs from default : ",
    !isTRUE(all.equal(unname(a), unname(d))), "\n")
p("predictive_interval(re_formula = NA, re.form = NA)",
  predictive_interval(ds, re_formula = NA, re.form = NA))
cat("\n---- ours: pp_check ----\n")
p("pp_check(ds, re.form = NA)", pp_check(ds, re.form = NA, ndraws = 5))
p("pp_check(fit, re.form = NA)", pp_check(fit, re.form = NA, ndraws = 5))
p("pp_check(ds, nosucharg = 1)", pp_check(ds, nosucharg = 1, ndraws = 5))
p("pp_check(ds, re_formula = NA)",
  pp_check(ds, re_formula = NA, ndraws = 5))
cat("\nDONE\n")
