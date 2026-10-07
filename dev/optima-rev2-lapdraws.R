# Reviewer of lane optima, re-check (c, d): frm_sample(laplace = TRUE) on
# a mo() model with random effects, where the added simo_[D] column
# enters draws_is_laplace()'s column count, and the remaining readers
# (print, summary, ranef, mcmc_plot, pairs) on the draws.
#   Rscript dev/optima-rev2-lapdraws.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
try1 <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  cat(sprintf("  %-28s %s\n", label, if (inherits(r, "error")) {
    paste("ERROR", substr(gsub("\n", " ", conditionMessage(r)), 1, 160))
  } else "ok"))
  invisible(r)
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(5)
n <- 300
d <- data.frame(x1 = sample(0:3, n, TRUE), x2 = sample(0:4, n, TRUE),
                g = factor(sample(1:15, n, TRUE)))
d$y <- c(0, 1, 1, 2)[d$x1 + 1] + c(0, 0.3, 0.6, 0.6, 1)[d$x2 + 1] +
  rnorm(15, sd = 0.5)[d$g] + rnorm(n)
fit <- frm(bf(y ~ mo(x1) + mo(x2) + (1 | g)), data = d, family = gaussian())
for (lap in c(FALSE, TRUE)) {
  cat("== laplace =", lap, "==\n")
  s <- q(frm_sample(fit, chains = 2, iter = 600, warmup = 300, seed = 4,
                    cores = 1, refresh = 0, laplace = lap))
  cat("  draws_is_laplace:", frmtmb.sample:::draws_is_laplace(s),
      "| ncol", ncol(s$draws), "\n")
  pe <- try1("posterior_epred", posterior_epred(s))
  if (is.matrix(pe)) cat("  epred finite:", all(is.finite(pe)), "dim",
                         dim(pe), "\n")
  try1("print", capture.output(print(s)))
  try1("summary", summary(s))
  try1("ranef", ranef(s))
  try1("log_lik", log_lik(s))
  try1("loo", q(loo(s)))
  try1("mcmc_plot", mcmc_plot(s))
  try1("pairs", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    pairs(s)
  })
  try1("conditional_effects", q(conditional_effects(s, "x2")))
  try1("predict(newdata)", q(predict(s, newdata = d[1:5, ])))
}
