# Reviewer of lane fixes, claim 2: what brms 2.23.0 does on the update
# of ledger row brmsfit-methods:955. testmode returns fit2 unchanged
# (update.R: the recompile branch skips brm()), so the model the update
# describes is built here by hand and sampled with fit2's priors and
# with brms's default (flat) priors on the coefficients.
#   Rscript dev/fixes-rev-u955-brms.R
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
fit2 <- brms:::brmsfit_example2
up <- update(fit2, formula. = bf(count ~ a + b, nl = TRUE), testmode = TRUE)
cat("testmode returns fit2 unchanged:",
    identical(deparse1(formula(up)$formula), deparse1(formula(fit2)$formula)),
    "\n")
dat <- fit2$data
fo <- bf(count | weights(AgeSD) ~ a + b, a ~ Age + (1 | ID1 | patient),
         b ~ Age + (1 | ID1 | patient), nl = TRUE)
pr <- c(prior(normal(2, 2), nlpar = "a"), prior(normal(0, 3), nlpar = "b"))
summ <- function(fit, lab) {
  s <- summary(fit)$fixed
  cat("==", lab, "\n")
  print(round(s[, c("Estimate", "Est.Error", "Rhat", "Bulk_ESS")], 3))
  np <- nuts_params(fit)
  cat("divergences:", sum(np$Value[np$Parameter == "divergent__"]), "\n")
  dr <- as_draws_df(fit)
  cat("posterior cor(b_a_Intercept, b_b_Intercept):",
      round(cor(dr$b_a_Intercept, dr$b_b_Intercept), 3),
      " sd of the sum:", round(sd(dr$b_a_Intercept + dr$b_b_Intercept), 3),
      "\n")
}
for (arm in c("fit2 priors", "flat priors")) {
  t0 <- Sys.time()
  f <- tryCatch(brm(fo, data = dat, family = Gamma("identity"),
                    prior = if (arm == "fit2 priors") pr else NULL,
                    chains = 2, iter = 1000, cores = 1, refresh = 0,
                    seed = 955),
                error = function(e) e)
  if (inherits(f, "error")) cat(arm, "ERROR:", conditionMessage(f), "\n") else
    summ(f, arm)
  cat("elapsed", format(Sys.time() - t0), "\n")
}
