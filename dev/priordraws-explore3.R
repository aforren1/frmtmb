# The spec fields the prior-draw plan reads, on three models
# (rellib-r14).
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
show <- function(fit) {
  for (r in names(fit$spec$responses)) {
    rs <- fit$spec$responses[[r]]
    cat("resp", r, "names:", paste(names(rs), collapse = " "), "\n")
    cat("  family dpars:", paste(rs$family$dpars, collapse = " "), "\n")
    cat("  primary_dpars:", paste(rs$primary_dpars, collapse = " "), "\n")
    cat("  nlpars:", paste(rs$nlpars, collapse = " "), "\n")
  }
  cat("linpreds:", paste(vapply(fit$frame$linpreds, function(lp)
    paste0(lp$resp, ":", lp$dpar), ""), collapse = " "), "\n")
  print(frmtmb::brms_coef_table(fit))
  tb <- validate_prior(NULL, fit,
                       route = "sample")
  print(as.data.frame(tb))
}
show(frm(bf(y ~ x + mo(m) + (1 + x | g), sigma ~ z), data = d,
         family = student()))
show(frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d,
         family = gaussian(), start = list(beta = c(1, 0.1))))
show(frm(bf(y ~ x) + bf(z ~ x) + set_rescor(TRUE), data = d))
show(frm(bf(y ~ x), data = d, family = mixture(gaussian(), gaussian())))
