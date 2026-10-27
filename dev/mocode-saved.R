## A fit saved by frmtmb 0.75.0 (rellib-r13), read by the lane's build.
##   Rscript dev/mocode-saved.R base   # writes dev/mocode-out/old-*.rds
##   Rscript dev/mocode-saved.R lane   # reads them
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
out <- "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/mocode-out"
cat("== frmtmb", format(packageVersion("frmtmb")), arm, "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$k0 <- d$k - 1L
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
if (identical(arm, "base")) {
  saveRDS(q(frm(y ~ mo(k), data = d)), file.path(out, "old-k.rds"))
  saveRDS(q(frm(y ~ mo(k0), data = d)), file.path(out, "old-k0.rds"))
  cat("saved\n")
  quit("no")
}
try_ <- function(lab, e) {
  r <- tryCatch(q(e), error = function(e) e)
  cat(sprintf("%-34s %s\n", lab, if (inherits(r, "error")) {
    paste("ERROR:", substr(conditionMessage(r), 1, 150))
  } else paste(signif(head(as.numeric(r), 4), 7), collapse = " ")))
}
for (v in c("k", "k0")) {
  fit <- readRDS(file.path(out, paste0("old-", v, ".rds")))
  nd <- data.frame(x = if (v == "k") 1:4 else 0:3)
  names(nd) <- v
  cat("-- 0.75.0 fit of y ~ mo(", v, ")\n", sep = "")
  try_("fitted()", fitted(fit)[, 1])
  try_("fitted(newdata)", fitted(fit, newdata = nd)[, 1])
  try_("predict(newdata)", predict(fit, newdata = nd)[, 1])
  try_("simulate()", simulate(fit, nsim = 1, seed = 1)[[1]])
  try_("residuals()", residuals(fit))
  try_("conditional_effects()",
       conditional_effects(fit)[[1]][["estimate__"]])
  try_("logLik()", logLik(fit))
  try_("fixef()", fixef(fit)[, 1])
  try_("coef()", unlist(coef(fit)))
  try_("AIC()", AIC(fit))
  try_("BIC()", BIC(fit))
  try_("vcov()", vcov(fit))
  try_("confint()", confint(fit))
  try_("summary()", summary(fit)$fixed[, 1])
  try_("drop1()", drop1(fit)[, 2])
  try_("anova(fit, update(fit))", anova(fit, update(fit))[, 1])
  try_("frm_sample()", fixef(frm_sample(fit, chains = 1, iter = 200,
                                        seed = 1, refresh = 0))[, 1])
  try_("update() (refit)", fixef(update(fit))[, 1])
}
