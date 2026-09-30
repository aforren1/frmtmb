# Reviewer, lane ordinal: reach of the constant disc coefficient that
# every ordinal family now carries. Usage: Rscript ... <base|lane>
# Data seed 20261003. Output: dev/ordinal-rev-log-reach-<arm>.txt
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(emmeans)})
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261003)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
d$yg <- rnorm(n) + d$x
try_ <- function(lab, expr) {
  cat("\n##", lab, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  print(r)
}
fits <- list(cumulative = frm(y ~ x, family = cumulative(), data = d),
             hurdle = frm(yh ~ x, family = hurdle_cumulative(), data = d))
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n==========", nm, "\n")
  try_("fixef(flatten = TRUE)", fixef(f, flatten = TRUE))
  try_("documented invariant: flatten names in rownames(confint())",
       all(names(fixef(f, flatten = TRUE)) %in% rownames(confint(f))))
  try_("confint rownames", rownames(confint(f)))
  bs <- frm_bootstrap(f, nsim = 20, seed = 4)
  try_("bootstrap t colnames", colnames(bs$t))
  try_("print(bootstrap)", utils::capture.output(print(bs)))
  try_("confint(bootstrap)", confint(bs))
  try_("confint(bootstrap, type = 'bca')", confint(bs, type = "bca"))
  try_("summary(bootstrap)", summary(bs))
  try_("emmeans", as.data.frame(emmeans(f, ~ x)))
  try_("summary()$fixed_dpars", summary(f)$fixed_dpars)
  try_("tidy-like: fixef()", fixef(f))
  try_("frm_linpred(dpar = 'disc')", head(frm_linpred(f, dpar = "disc")))
  try_("variables", variables(f))
}
mv <- frm(bf(y ~ x) + bf(yg ~ x), family = list(cumulative(), gaussian()),
          data = d)
try_("mv coef names", names(coef(mv)))
try_("mv coef y_disc", coef(mv)[["y_disc"]])
try_("mv emmeans resp = y", as.data.frame(emmeans(mv, ~ x, resp = "y")))
try_("mv fixef flatten", fixef(mv, flatten = TRUE))
