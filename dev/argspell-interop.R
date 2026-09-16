## The reviewer's interop harness (dev/asrev-interop.R), rerun against
## the LANE library after five previously exempt interop methods were
## given the dots guard. Only the library path differs.
## REVIEW 2.5e, risk 1: does the dots refusal over-refuse?
## Exercise the real third-party packages against a real fit.
LIB <- "C:/Users/adf44/source/r/argspell-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
cat("version    :", as.character(utils::packageVersion("frmtmb")), "\n\n")

set.seed(2505)
dd <- data.frame(x = stats::rnorm(120),
                 f = factor(rep(c("a", "b", "c"), length.out = 120)),
                 g = factor(rep(1:8, each = 15)), y = 0)
dd$y <- frm_simulate(bf(y ~ x + f + (1 | g)) + gaussian(), dd,
                     newparams = list(Intercept = 1, x = 0.5,
                                      fb = 0.3, fc = -0.2,
                                      sigma = 0.4, sd_g__Intercept = 1.2),
                     nsim = 1, seed = 25051)[[1]]
fit <- frm(bf(y ~ x + f + (1 | g)) + gaussian(), data = dd)
cat("fit ok\n\n")

try_it <- function(label, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("FAIL  %-46s %s\n", label,
                gsub("\n", " ", conditionMessage(r))))
  } else {
    cat(sprintf("ok    %-46s %s\n", label,
                paste(class(r), collapse = "/")))
  }
  invisible(r)
}

cat("---- emmeans ----\n")
if (requireNamespace("emmeans", quietly = TRUE)) {
  cat("emmeans", as.character(packageVersion("emmeans")), "\n")
  try_it("emmeans(fit, ~ f)", emmeans::emmeans(fit, ~ f))
  try_it("emmeans(fit, pairwise ~ f)", emmeans::emmeans(fit, pairwise ~ f))
  try_it("ref_grid(fit)", emmeans::ref_grid(fit))
  try_it("emtrends(fit, ~ f, var='x')",
         emmeans::emtrends(fit, ~ f, var = "x"))
  try_it("recover_data via ref_grid(data=)",
         emmeans::ref_grid(fit, data = dd))
  try_it("emmeans(fit, ~f, weights='cells')",
         emmeans::emmeans(fit, ~ f, weights = "cells"))
} else cat("emmeans NOT INSTALLED\n")

cat("\n---- insight ----\n")
if (requireNamespace("insight", quietly = TRUE)) {
  cat("insight", as.character(packageVersion("insight")), "\n")
  try_it("get_parameters(fit)", insight::get_parameters(fit))
  try_it("get_parameters(fit, effects='all')",
         insight::get_parameters(fit, effects = "all"))
  try_it("n_parameters(fit)", insight::n_parameters(fit))
  try_it("n_parameters(fit, effects='fixed')",
         insight::n_parameters(fit, effects = "fixed"))
  try_it("find_formula(fit)", insight::find_formula(fit))
  try_it("find_random(fit)", insight::find_random(fit))
  try_it("find_statistic(fit)", insight::find_statistic(fit))
  try_it("link_function(fit)", insight::link_function(fit))
  try_it("link_inverse(fit)", insight::link_inverse(fit))
  try_it("get_varcov(fit)", insight::get_varcov(fit))
  try_it("get_predicted(fit)", insight::get_predicted(fit))
  try_it("model_info(fit)", insight::model_info(fit))
  try_it("get_data(fit)", insight::get_data(fit))
  try_it("get_variance(fit)", insight::get_variance(fit))
} else cat("insight NOT INSTALLED\n")

cat("\n---- parameters / performance ----\n")
if (requireNamespace("parameters", quietly = TRUE)) {
  cat("parameters", as.character(packageVersion("parameters")), "\n")
  try_it("parameters::model_parameters(fit)",
         parameters::model_parameters(fit))
} else cat("parameters NOT INSTALLED\n")
if (requireNamespace("performance", quietly = TRUE)) {
  cat("performance", as.character(packageVersion("performance")), "\n")
  try_it("performance::r2(fit)", performance::r2(fit))
} else cat("performance NOT INSTALLED\n")

cat("\n---- marginaleffects ----\n")
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  cat("marginaleffects", as.character(packageVersion("marginaleffects")), "\n")
  try_it("avg_slopes(fit)", marginaleffects::avg_slopes(fit))
  try_it("slopes(fit)", marginaleffects::slopes(fit))
  try_it("predictions(fit)", marginaleffects::predictions(fit))
  try_it("avg_predictions(fit, by='f')",
         marginaleffects::avg_predictions(fit, by = "f"))
  try_it("avg_comparisons(fit)", marginaleffects::avg_comparisons(fit))
  try_it("hypotheses(avg_slopes)",
         marginaleffects::hypotheses(marginaleffects::avg_slopes(fit),
                                     "b1 = 0"))
} else cat("marginaleffects NOT INSTALLED\n")

cat("\n---- other dots consumers on the fit surface ----\n")
try_it("broom.mixed::tidy", {
  if (!requireNamespace("broom.mixed", quietly = TRUE)) stop("not installed")
  broom.mixed::tidy(fit)
})
try_it("modelsummary", {
  if (!requireNamespace("modelsummary", quietly = TRUE)) stop("not installed")
  modelsummary::modelsummary(fit, output = "data.frame")
})
try_it("stats::confint(fit)", stats::confint(fit))
try_it("stats::update(fit)", stats::update(fit))
try_it("anova(fit)", stats::anova(fit))
try_it("car::Anova(fit)", {
  if (!requireNamespace("car", quietly = TRUE)) stop("not installed")
  car::Anova(fit)
})
cat("\nDONE\n")
