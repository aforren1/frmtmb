## brms 2.23.0: can a mo() term on a single value (Jmo = 0) be sampled?
## The ported row brmsfit-methods:924 updates fixture 1 with Exp = 4 on
## every row under testmode = TRUE, which never samples.
##   Rscript dev/mocode-brms-const.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
set.seed(1)
d <- data.frame(x = rep(4, 40), z = rnorm(40))
d$y <- d$z + rnorm(40)
r <- tryCatch(suppressMessages(brm(y ~ z + mo(x), data = d, chains = 1,
                                   iter = 400, seed = 1, refresh = 0)),
              error = function(e) e)
if (inherits(r, "error")) {
  cat("brm ERROR:", conditionMessage(r), "\n")
} else {
  cat("brm returned; draws:", posterior::ndraws(r), "\n")
  print(summary(r)$fixed)
  print(tryCatch(posterior::summarise_draws(posterior::subset_draws(
    posterior::as_draws(r), variable = "bsp_mox")), error = function(e) e))
}
