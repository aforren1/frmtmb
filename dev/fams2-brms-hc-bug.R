# brms 2.23.0's hurdle_cumulative_logit_lpmf() (the generic logit density,
# used when disc is modeled or cs() is present) tests y == nthres + 2 for
# the top category, which never occurs, so the top category y = nthres + 1
# falls to the interior branch and reads thres[nthres + 1]. Evaluate the
# compiled program's log density at a point to see what Stan does there.
# Seed 18; the same data as row 16g of test-brms-likelihood.R.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(brms); library(rstan)})
set.seed(18)
n <- 300
d <- data.frame(x = rnorm(n))
u <- stats::rlogis(n) + 0.8 * d$x
yc <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
d$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * d$x), 0L, yc)
cat("rows in the top category:", sum(d$y == max(d$y)), "\n")
run <- function(link, form) {
  bform <- bf(form, disc ~ 0 + x)
  pr <- get_prior(bform, data = d, family = hurdle_cumulative(link))
  pr$prior <- ""
  code <- stancode(bform, data = d, family = hurdle_cumulative(link),
                   prior = pr)
  sdat <- standata(bform, data = d, family = hurdle_cumulative(link),
                   prior = pr)
  cat("\n== link", link, "==\n")
  cat("lpmf used:", regmatches(code, regexpr("target \\+= [a-z_]+", code)),
      "\n")
  mod <- stan_model(model_code = code)
  sf <- suppressMessages(sampling(mod, data = sdat, chains = 0))
  pars <- list(b = array(0.8, 1), Intercept = c(-1, 0.3, 1.5), hu = 0.35,
               b_disc = array(0.1, 1))
  up <- unconstrain_pars(sf, pars)
  r <- tryCatch(log_prob(sf, up, adjust_transform = FALSE),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("log_prob at a legal point:", format(r), "\n")
}
run("logit", y ~ x)
run("probit", y ~ x)
