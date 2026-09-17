## Reviewer check for lane wt-famlink, priority 2: a fit saved by the
## BASE build (dev/famlink-rev-savedfit-base.rds, written by
## dev/famlink-rev-interop.R base, seed 20260916) read by either arm.
## Usage: Rscript dev/famlink-rev-savedfit.R <lane|base>
ARM <- commandArgs(trailingOnly = TRUE)[1]
source("dev/famlink-rev-common.R")
suppressMessages({library(emmeans); library(insight); library(marginaleffects)})
fits <- readRDS("dev/famlink-rev-savedfit-base.rds")
try1 <- function(label, expr) {
  v <- tryCatch(withCallingHandlers({force(expr); "ok"},
                  warning = function(w) invokeRestart("muffleWarning"),
                  message = function(m) invokeRestart("muffleMessage")),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-22s %s\n", label, substr(v, 1, 150)))
}
for (nm in names(fits)) {
  fit <- fits[[nm]]
  cat(nm, "\n")
  try1("family(fit)$link", { v <- family(fit)$link; cat("    value class:", class(v), "\n") })
  try1("family(fit)$links", family(fit)$links)
  try1("print(family)", capture.output(print(family(fit))))
  try1("print(fit)", capture.output(print(fit)))
  try1("summary(fit)", capture.output(summary(fit)))
  try1("predict(fit)", predict(fit))
  try1("predict(newdata)", predict(fit, newdata = fit$data[1:5, ]))
  try1("fitted(fit)", fitted(fit))
  try1("residuals(fit)", residuals(fit))
  try1("simulate(fit)", simulate(fit, nsim = 1, seed = 1))
  try1("confint(fit)", confint(fit))
  try1("update(fit)", update(fit))
  try1("model_info", insight::model_info(fit))
  try1("emmeans", emmeans::emmeans(fit, "f"))
  try1("avg_slopes", marginaleffects::avg_slopes(fit))
  try1("conditional_effects", conditional_effects(fit))
  try1("loo", loo(fit))
  try1("get_prior", get_prior(fit))
}
