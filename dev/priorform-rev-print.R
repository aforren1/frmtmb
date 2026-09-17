# Reviewer, lane wt-priorform: set_prior() recycling and the print layout
# against brms 2.23.0, character for character.
#   Rscript dev/priorform-rev-print.R brms|lane
# Writes dev/priorform-rev-print-<mode>.txt
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (mode == "lane") "C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (mode == "brms") suppressMessages(library(brms)) else
  suppressMessages(library(frmtmb))
calls <- c(
  'set_prior("normal(0,1)")',
  'set_prior("normal(0,1)", coef = "x")',
  'set_prior("cauchy(0,1)", class = "sd", group = "x")',
  'set_prior("cauchy(0, 1)", class = "sd", group = "g", coef = "Intercept")',
  'set_prior("normal(0, 1)", class = "b", dpar = "sigma")',
  'set_prior("normal(0, 1)", class = "b", nlpar = "a", coef = "x")',
  'set_prior("normal(0, 1)", class = "b", resp = "y1", coef = "x")',
  'set_prior("normal(0, 1)", class = "sd", group = "g", dpar = "phi")',
  'set_prior("normal(0, 1)", class = "sd", group = "g", resp = "y1", dpar = "sigma", nlpar = "a")',
  'set_prior("normal(0, 1)", lb = 0)',
  'set_prior("normal(0, 1)", ub = 1)',
  'set_prior("normal(0, 1)", lb = 0, ub = 1)',
  'set_prior("", class = "b", coef = "x", lb = 0)',
  'set_prior("student_t(3, 0, 2.5)", class = "sigma")',
  'set_prior("student_t(3, 0, 2.5)", class = "Intercept")',
  'set_prior("lkj(2)", class = "cor")',
  'set_prior("normal(0, 2)", class = c("b", "sd"))',
  'set_prior(c("normal(0, 1)", "cauchy(0, 2)"), class = c("b", "sd"))',
  'set_prior("normal(0, 1)", coef = c("x", "z", "w", "v"), class = c("b", "b"))',
  'set_prior("normal(0, 1)", coef = c("x", "z", "w"), class = c("b", "b"))',
  'set_prior("normal(0, 1)", coef = character(0))',
  'set_prior(character(0))',
  'set_prior("normal(0, 1)", coef = NA)',
  'set_prior("normal(0, 1)", lb = c(0, 1), coef = c("x", "z"))',
  'set_prior("normal(0, 1)", lb = c(0, NA), coef = c("x", "z"))',
  'prior(normal(0, 1), class = sd) + prior(cauchy(0, 2), class = b)',
  'c(set_prior("normal(0, 1)"), set_prior("normal(0, 2)", coef = "x"))',
  'set_prior("normal(0,1)")$prior',
  'set_prior("normal(0, 2)", class = c("b", "sd"))$class',
  'set_prior("normal(0, 1)", lb = 0)$lb',
  'set_prior("normal(0, 1)")$ub',
  'set_prior("normal(0, 1)")$source',
  'set_prior("normal(0, 1)")$coef',
  'empty_prior()',
  'nrow(as.data.frame(empty_prior()))',
  'empty_prior() + set_prior("normal(0, 1)")',
  'as.brmsprior(data.frame(prior = "normal(0,1)", x = "test", coef = c("a", "b")))',
  'as.brmsprior(data.frame(prior = c("normal(0,1)", ""), lb = c(NA, "0")))',
  'names(as.data.frame(set_prior("normal(0, 1)")))'
)
lines <- character(0)
for (cl in calls) {
  o <- tryCatch(utils::capture.output(print(eval(parse(text = cl)))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  lines <- c(lines, paste0("### ", cl), o)
}
writeLines(lines, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-priorform/dev/priorform-rev-print-%s.txt", mode))
