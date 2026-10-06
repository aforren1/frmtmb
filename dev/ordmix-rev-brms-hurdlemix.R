# Reviewer of lane ordmix: brms 2.23.0's mixture(cumulative(),
# hurdle_cumulative()), which the lane refuses. Does brms's compiled
# program evaluate on a response with zeros (the hurdle's), and on one
# without? Seed 1, n = 200. Stan cache: dev/ordmix-rev-stan-cache.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(rstan)
})
cache_dir <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-rev-stan-cache"
stan_mod <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- sample(1:4, n, TRUE)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
fam <- mixture(cumulative(), hurdle_cumulative())
for (resp in if (length(commandArgs(TRUE))) commandArgs(TRUE) else c("yh", "y")) {
  f <- bf(as.formula(paste(resp, "~ x")))
  code <- suppressMessages(stancode(f, data = d, family = fam))
  sdat <- suppressMessages(standata(f, data = d, family = fam))
  mod <- stan_mod(code)
  sf <- suppressMessages(sampling(mod, data = sdat, chains = 0))
  np <- rstan::get_num_upars(sf)
  set.seed(2)
  r <- tryCatch(rstan::log_prob(sf, rnorm(np, 0, 0.3)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("response %s (zeros: %d): nthres %d, log_prob %s\n", resp,
              sum(d[[resp]] == 0), sdat$nthres, substr(as.character(r), 1, 300)))
}
