# Reviewer, lane wt-conditions (round 2): interop entry points on a
# gaussian, a binomial and a mixed fit, run on one library arm, and the
# error paths where the caller inspects the class of a caught condition.
# Run once per arm and diff the two outputs.
#   Rscript dev/conditions-rev-interop.R base|lane
arm <- commandArgs(TRUE)[1L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib"
            else "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(width = 120, digits = 6)
set.seed(20260917)
d <- data.frame(y = rnorm(80), x = rnorm(80), f = gl(4, 20),
                g = factor(rep(1:10, 8)))
d$n <- 10L; d$k <- rbinom(80, 10, plogis(0.3 * d$x))
fits <- list(
  gaussian = frm(y ~ x + f, data = d, family = gaussian()),
  binomial = frm(k | trials(n) ~ x + f, data = d, family = binomial()),
  mixed = frm(y ~ x + f + (1 | g), data = d, family = gaussian()))
show <- function(label, expr) {
  cat("\n##", label, "\n")
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat("ERROR [", paste(class(e), collapse = "/"), "]",
          conditionMessage(e), "\n")
      NULL
    }),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    },
    message = function(c) {
      w <<- c(w, paste("MSG", conditionMessage(c)))
      invokeRestart("muffleMessage")
    })
  if (length(w)) cat("CONDITIONS:", unique(w), sep = "\n  ")
  if (!is.null(r)) {
    out <- tryCatch(utils::capture.output(print(r)), error = function(e)
      paste("print failed:", conditionMessage(e)))
    cat(head(out, 12), sep = "\n")
  }
}
nd_one <- data.frame(x = c(0, 1), f = factor(c("1", "1")),
                     g = factor(c("1", "1")), n = 10L)
nd_bad <- data.frame(x = c(0, 1), f = factor(c("zz", "zz")),
                     g = factor(c("1", "1")), n = 10L)
for (nm in names(fits)) {
  fit <- fits[[nm]]
  cat("\n==========", nm, "\n")
  show("emmeans f", as.data.frame(emmeans::emmeans(fit, ~ f)))
  show("emmeans contrast", as.data.frame(emmeans::contrast(
    emmeans::emmeans(fit, ~ f), "pairwise")))
  show("emmeans bad spec", emmeans::emmeans(fit, ~ nosuch))
  show("mfx avg_predictions", marginaleffects::avg_predictions(fit))
  show("mfx avg_slopes", marginaleffects::avg_slopes(fit, variables = "x"))
  show("mfx predictions bad newdata",
       marginaleffects::predictions(fit, newdata = nd_bad))
  show("mfx predictions bad type",
       marginaleffects::predictions(fit, type = "bogus"))
  show("insight get_predicted", head(as.data.frame(
    insight::get_predicted(fit, ci = 0.95))))
  show("insight get_predicted single-level factor data",
       insight::get_predicted(fit, data = nd_one, ci = 0.95))
  show("insight get_predicted bad level data",
       insight::get_predicted(fit, data = nd_bad, ci = 0.95))
  show("insight .get_predicted_ci_modelmatrix bad data",
       insight:::.get_predicted_ci_modelmatrix(fit, data = nd_bad))
  show("insight get_modelmatrix bad data",
       insight::get_modelmatrix(fit, data = nd_bad))
  show("insight find_formula", insight::find_formula(fit))
  show("insight get_varcov", insight::get_varcov(fit))
  show("insight get_parameters", insight::get_parameters(fit))
  show("model.matrix bad data", model.matrix(fit, data = nd_bad))
  show("broom tidy", broom::tidy(fit))
}
