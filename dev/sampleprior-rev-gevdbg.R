source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
ns <- asNamespace("frmtmb.sample")
f <- eval(CASES$gev_xi1$frm)
od <- cap_obj(f[[1]], f[[2]], f[[3]])
of <- cap_obj(f[[1]], f[[2]], f[[3]], prior = "flat")
fit <- suppressMessages(suppressWarnings(ns$sample_assemble(
      f[[1]], f[[2]], f[[3]], data2 = list(), start = NULL,
      fit_control = frmtmb_control(), na.action = stats::na.omit, REML = FALSE)))
ob <- ns$prior_augmented_obj(fit, list())
set.seed(7)
p0 <- od$env$par
for (k in 1:6) {
  q <- p0 + rnorm(length(p0), 0, 0.4)
  cat(k, "q:", round(q, 4), " od", od$env$f(q), " of", of$env$f(q), " bare-retape", ob$env$f(q),
      " lp(od-bare)", -(od$env$f(q) - ob$env$f(q)), "\n")
}
