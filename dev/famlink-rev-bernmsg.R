## Reviewer check for lane wt-famlink, priority 6: how many times does
## the bernoulli message fire across refit entry points, on a model
## where brms messages once per brm() call? Seed 20260916, n = 80.
## Usage: Rscript dev/famlink-rev-bernmsg.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
set.seed(20260916)
d <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
d$y <- rbinom(80, 1, plogis(0.3 + 0.6 * d$x))
d$one <- 1
count <- function(label, expr) {
  n <- 0L
  v <- tryCatch(withCallingHandlers(expr,
    message = function(m) {
      if (grepl("bernoulli", conditionMessage(m))) n <<- n + 1L
      invokeRestart("muffleMessage")
    }, warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-34s messages %3d %s\n", label, n,
              if (is.character(v) && length(v) == 1) substr(v, 1, 90) else ""))
  invisible(v)
}
fit <- count("frm binomial trials(1)",
             frm(y | trials(one) ~ x, data = d, family = binomial()))
count("frm bernoulli (control)", frm(y ~ x, data = d, family = bernoulli()))
count("update()", update(fit))
count("update(formula.)", update(fit, . ~ . + g))
count("refit()", lme4::refit(fit, d$y))
count("frm_bootstrap(nsim = 5)", frm_bootstrap(fit, nsim = 5, seed = 1))
count("frm_allfit()", frm_allfit(fit))
count("profile(parm = 1)", profile(fit, parm = 1))
count("confint(method=profile)", confint(fit, parm = 1, method = "profile"))
count("simulate()", simulate(fit, nsim = 2, seed = 1))
count("predict(newdata)", predict(fit, newdata = d[1:5, ]))
count("loo()", loo(fit))
count("get_prior()", get_prior(y | trials(one) ~ x, data = d, family = binomial()))
count("frm(dry_run = 'frame')", frm(y | trials(one) ~ x, data = d,
                                     family = binomial(), dry_run = "frame"))
count("frm_multiple(list of 3)", frm_multiple(y | trials(one) ~ x,
                                               data = list(d, d, d),
                                               family = binomial()))
count("mvbf two binomial responses", {
  d$y2 <- d$y
  frm(mvbf(bf(y | trials(one) ~ x), bf(y2 | trials(one) ~ x)) + binomial(),
      data = d)
})
count("mixture(binomial, binomial) trials(1)",
      frm(y | trials(one) ~ 1, data = d, family = mixture(binomial, binomial)))
