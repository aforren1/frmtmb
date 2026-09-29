# Does a fresh Stan program compile here at all? Isolates the
# "invalid connection" failure of dev/csfactor-lp.R from the cs() work.
.libPaths(c("C:/Users/adf44/source/r/wt-csfactor-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
cat("R_MAKEVARS_USER =", Sys.getenv("R_MAKEVARS_USER"), "\n")
cat("StanHeaders", as.character(packageVersion("StanHeaders")),
    "rstan", as.character(packageVersion("rstan")), "\n")
code <- "
data { int<lower=0> N; vector[N] y; }
parameters { real mu; real<lower=0> s; }
model { target += normal_lpdf(y | mu, s); }
"
r <- tryCatch({
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  sf <- suppressMessages(rstan::sampling(m, data = list(N = 3,
                                                        y = c(1, 2, 3)),
                                         chains = 0))
  paste("OK lp =", rstan::log_prob(sf, rstan::unconstrain_pars(
    sf, list(mu = 0, s = 1)), adjust_transform = FALSE))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, "\n")
