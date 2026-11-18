# Smoke test: both tools load and sample a tiny model.
lib <- "C:/Users/adf44/source/r/wt-samplers-lib"
Sys.setenv(R_USER_CACHE_DIR = file.path(lib, "_cache"))
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(stanli); library(SparseNUTS); library(RTMB)
})
set.seed(1)
d <- data.frame(x = rnorm(50))
d$y <- 1 + 0.5 * d$x + rnorm(50)
t0 <- proc.time()[[3]]
code <- brms::make_stancode(y ~ x, data = d)
sdat <- brms::make_standata(y ~ x, data = d)
m <- stanli_model(code = code, data = unclass(sdat))
t1 <- proc.time()[[3]]
fit <- sample_model(m, chains = 4, seed = 1, refresh = 0,
                    parallel_chains = 1)
t2 <- proc.time()[[3]]
print(summary(fit))
sf <- as_stanfit(fit)
bf <- brms::brm(y ~ x, data = d, backend = "mock", mock_fit = sf,
                rename = TRUE)
print(brms::fixef(bf))
cat("stanli setup", t1 - t0, "sample", t2 - t1, "\n")

dat <- list(y = c(28, 8, -3, 7, -1, 1, 18, 12),
            sigma = c(15, 10, 16, 11, 9, 11, 10, 18))
f <- function(pars) {
  getAll(dat, pars)
  theta <- mu + exp(logtau) * eta
  -(sum(dnorm(eta, 0, 1, log = TRUE)) +
      sum(dnorm(y, theta, sigma, log = TRUE)) + logtau)
}
obj <- MakeADFun(f, list(mu = 0, logtau = 0, eta = rep(1, 8)),
                 random = "eta", silent = TRUE)
t3 <- proc.time()[[3]]
sf2 <- sample_snuts(obj, refresh = 0, seed = 1, cores = 1, chains = 2,
                    model_name = "schools")
cat("snuts", proc.time()[[3]] - t3, "\n")
print(sf2)
str(sf2$monitor)
