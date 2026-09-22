# Reviewer, priority 3: interop entry points the lane's own script does
# NOT run. Each is a route by which another package can reach
# predict(), fitted() or vcov() on a frmtmb fit without any call site
# in this repository.
#
#   Rscript dev/shapes-rev-interop2.R lane|base

arg <- commandArgs(trailingOnly = TRUE)
which_lib <- if (identical(arg[1], "base")) "base" else "lane"
lib <- if (which_lib == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
dd <- rev_data()
cat("lib:", lib, "\n")
fg <- frm(bf(y ~ x + f) + gaussian(), data = dd)
fm <- frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd)
fo <- frm(bf(ord ~ x) + cumulative(), data = dd)

p <- function(lab, e) {
  r <- tryCatch(withCallingHandlers(suppressMessages(e),
                                    warning = function(w) {
                                      cat("   [warn] ",
                                          substr(conditionMessage(w), 1, 90),
                                          "\n", sep = "")
                                      invokeRestart("muffleWarning") }),
                error = function(c) paste0("ERROR: ",
                                           substr(conditionMessage(c), 1, 90)))
  if (is.character(r) && length(r) == 1 && startsWith(r, "ERROR")) {
    cat(sprintf("  %-46s %s\n", lab, r)); return(invisible(NULL))
  }
  v <- if (is.data.frame(r)) {
    num <- vapply(r, is.numeric, TRUE)
    paste0("df ", nrow(r), "x", ncol(r), " sum=",
           format(sum(unlist(r[num]), na.rm = TRUE), digits = 10))
  } else if (is.numeric(r)) {
    paste0("num ", length(r), " sum=",
           format(sum(r, na.rm = TRUE), digits = 10))
  } else paste(class(r)[1], length(r))
  cat(sprintf("  %-46s %s\n", lab, v))
}

cat("\n== marginaleffects ==\n")
p("avg_predictions(gaussian)",
  as.data.frame(marginaleffects::avg_predictions(fg)))
p("avg_comparisons(gaussian)",
  as.data.frame(marginaleffects::avg_comparisons(fg)))
p("comparisons(gaussian, variables='x')",
  as.data.frame(marginaleffects::comparisons(fg, variables = "x")))
p("hypotheses(gaussian, 'b2 = 0')",
  as.data.frame(marginaleffects::hypotheses(fg, "b2 = 0")))
p("predictions(gaussian, type='link')",
  as.data.frame(marginaleffects::predictions(fg, type = "link")))
p("slopes(mixed)", as.data.frame(marginaleffects::slopes(fm)))
p("avg_slopes(ordinal)", as.data.frame(marginaleffects::avg_slopes(fo)))

cat("\n== insight ==\n")
p("get_predicted(gaussian)",
  as.numeric(insight::get_predicted(fg)))
p("get_predicted(gaussian, 'link')",
  as.numeric(insight::get_predicted(fg, predict = "link")))
p("get_predicted(gaussian, 'prediction')",
  as.numeric(insight::get_predicted(fg, predict = "prediction")))
p("get_predicted(ordinal)", as.numeric(insight::get_predicted(fo)))
p("get_fitted(gaussian)", as.numeric(insight::get_fitted(fg)))
p("get_residuals(gaussian)", as.numeric(insight::get_residuals(fg)))
p("get_sigma(gaussian)", as.numeric(insight::get_sigma(fg)))
p("get_df(gaussian)", as.numeric(insight::get_df(fg)))
p("get_variance(mixed)", unlist(insight::get_variance(fm)))
p("model_info(gaussian)", unlist(insight::model_info(fg)))
p("find_parameters(ordinal)", unlist(insight::find_parameters(fo)))

cat("\n== emmeans ==\n")
p("emmeans(gaussian, ~f)",
  as.data.frame(summary(emmeans::emmeans(fg, ~ f))))
p("emtrends(gaussian, ~f, 'x')",
  as.data.frame(summary(emmeans::emtrends(fg, ~ f, var = "x"))))
p("emmeans(mixed, ~x, at)", as.data.frame(summary(
  emmeans::emmeans(fm, "x", at = list(x = c(-1, 1))))))
p("emmeans(ordinal, ~x)", as.data.frame(summary(
  emmeans::emmeans(fo, "x", at = list(x = 0)))))

cat("\n== base R and broom-style seams ==\n")
p("as.data.frame(fit)", as.data.frame(fg))
p("coef(fit)", unlist(coef(fg)))
p("model.matrix(fit)", as.numeric(model.matrix(fg)))
p("simulate(fit, nsim = 2)", as.data.frame(simulate(fg, nsim = 2)))
p("stats::confint(fit)", as.data.frame(confint(fg)))
p("stats::nobs(fit)", as.numeric(stats::nobs(fg)))
p("sandwich::sandwich-style vcov(cluster)",
  as.numeric(vcov(fg, cluster = ~ g)))
