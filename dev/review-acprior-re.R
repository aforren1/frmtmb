# Reviewer: the random-effect analogue of the lane's double slot. A
# 3 x 3 block (1 + x + z | g): does one class theta prior on ONE
# correlation retire the joint LKJ on all three, on frm() and on
# frm_sample()'s defaults?
# Usage: Rscript dev/review-acprior-re.R before|after
# Log: dev/review-acprior-log/re-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(8)
d <- data.frame(g = rep(1:15, each = 8), x = rnorm(120), z = rnorm(120))
d$y <- 1 + d$x + rnorm(15)[d$g] + rnorm(120)
fo <- y ~ x + z + (1 + x + z | g)
fit <- frm(fo, data = d, dry_run = "objective")
cat("theta names:", names(unlist(fit$frame$par_template$theta)), "\n")
ent <- function(ri) {
  for (e in ri$entries) {
    if (!identical(e$comp, "theta")) next
    cat("  entry theta", paste(e$idx, collapse = ","), e$dist$kind,
        if (identical(e$dist$kind, "trans")) e$dist$inner$kind, "\n")
  }
}
pr <- c(set_prior("lkj(2)", class = "cor"),
        set_prior("normal(0, 1)", class = "theta", coef = "theta_4"))
cat("--- frm: cor lkj(2) then theta_4 normal\n")
ent(frmtmb:::resolve_prior_input(fit, pr))
cat("--- frm_sample: defaults alone\n")
ent(suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))$ri)
cat("--- frm_sample: defaults + theta_4 normal\n")
ent(suppressMessages(frmtmb.sample:::sample_resolve_priors(fit,
  set_prior("normal(0, 1)", class = "theta", coef = "theta_4")))$ri)
cat("--- frm_sample: defaults + sd normal(0, 1), coef x + theta_2 normal\n")
ent(suppressMessages(frmtmb.sample:::sample_resolve_priors(fit,
  c(set_prior("normal(0, 1)", class = "sd", group = "g", coef = "x"),
    set_prior("normal(0, 5)", class = "theta", coef = "theta_2"))))$ri)
