# Probe 2: the three value checks that did not hold on the first record
# (log_lik:271 discrete_weibull, log_lik:462 cumulative, log_lik:478
# hurdle_cumulative), on brms's own preps rebuilt here with the blocks'
# seeds, value by value. Log: dev/bin2port-out/probe2.txt.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/brmsport-bin2-adapters.R")
tab <- function(a, b) {
  cls <- function(x) ifelse(is.nan(x), "NaN", ifelse(is.na(x), "NA",
                       ifelse(is.infinite(x), "-Inf", "finite")))
  print(table(frmtmb = cls(a), brms = cls(b)))
}

## log_lik:271, the block at tests.log_lik.R:233, seed 233 (the
## generated block's withr::local_seed(<block line>))
blk <- function(line) {
  # the generator seeds each block with its first line
  set.seed(line)
}
src <- "dev/brms-suite/brms/tests/testthat/tests.log_lik.R"
cat("block line of log_lik:271:",
    grep("count and survival", readLines(src)), "\n")
gen <- readLines("dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R")
i0 <- grep("count and survival models works correctly", gen)
cat(gen[i0 + 1], "\n")
seed <- as.integer(sub(".*local_seed[(]([0-9]+)L[)].*", "\\1", gen[i0 + 1]))
set.seed(seed)
ns <- 25
nobs <- 10
trials <- sample(10:30, nobs, replace = TRUE)
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(
  eta = matrix(rnorm(ns * nobs), ncol = nobs),
  shape = rgamma(ns, 4),
  xi = runif(ns, -1, 0.5),
  phi = rgamma(ns, 1)
)
prep$dpars$sigma <- 1 / prep$dpars$shape
prep$dpars$nu <- prep$dpars$shape + 1
prep$data <- list(
  Y = rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1)),
  trials = trials
)
i <- sample(nobs, 1)
prep$dpars$mu <- 1 / (1 + exp(-prep$dpars$eta))
a <- brms_adapt_log_lik("discrete_weibull", i, prep)
b <- brms:::log_lik_discrete_weibull(i, prep = prep)
cat("y =", prep$data$Y[i], "\n")
tab(a, b)
d <- data.frame(mu = prep$dpars$mu[, i], shape = prep$dpars$shape,
                frmtmb = a, brms = b)
print(head(d[order(d$mu), ], 25), digits = 6)
# brms's own pmf: P(Y = y) = mu^(y^shape) - mu^((y+1)^shape)
y <- prep$data$Y[i]
mu <- prep$dpars$mu[, i]
sh <- prep$dpars$shape
exact <- log(mu^(y^sh) - mu^((y + 1)^sh))
lexact <- (y^sh) * log(mu) + log(-expm1(((y + 1)^sh - y^sh) * log(mu)))
cat("closed form on the log scale (no underflow) vs frmtmb: max abs diff",
    format(max(abs(lexact - a)), digits = 3), "\n")
print(head(cbind(mu, sh, naive = exact, logscale = lexact, frmtmb = a,
                 brms = b)[order(mu), ], 8), digits = 6)

## log_lik:462 and :478, the ordinal block (tests.log_lik.R:446): all
## thresholds 0, so the two middle categories have probability 0
gen_line <- grep("log_lik for ordinal models", gen)
seed <- as.integer(sub(".*local_seed[(]([0-9]+)L[)].*", "\\1",
                       gen[gen_line + 1]))
set.seed(seed)
ns <- 50
nobs <- 8
nthres <- 3
ncat <- nthres + 1
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(
  mu = array(rnorm(ns * nobs), dim = c(ns, nobs)),
  disc = rexp(ns),
  hu = rbeta(ns, 1, 1)
)
prep$thres$thres <- array(0, dim = c(ns, nthres))
prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
prep$family$link <- "logit"
a <- suppressWarnings(sapply(1:nobs, brms_adapt_log_lik, fun = "cumulative",
                             prep = prep))
b <- sapply(1:nobs, brms:::log_lik_cumulative, prep = prep)
cat("cumulative, all thresholds 0, Y =", prep$data$Y, "\n")
for (j in 1:4) {
  cat("category", prep$data$Y[j], ": frmtmb", format(a[1:3, j], digits = 4),
      " brms", format(b[1:3, j], digits = 4), "\n")
}
tab(a, b)
a <- suppressWarnings(brms_adapt_log_lik("hurdle_cumulative", 3, prep))
b <- brms:::log_lik_hurdle_cumulative(3, prep = prep)
cat("hurdle_cumulative, y =", prep$data$Y[3], ": frmtmb",
    format(a[1:3], digits = 4), " brms", format(b[1:3], digits = 4), "\n")
tab(a, b)
# the same rows with distinct thresholds, where the category has mass
prep$thres$thres <- matrix(rep(c(-1, 0, 1), each = ns), ns)
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "cumulative", prep = prep)
b <- sapply(1:nobs, brms:::log_lik_cumulative, prep = prep)
cat("cumulative, thresholds -1, 0, 1: max abs diff",
    format(max(abs(a - b)), digits = 3), "\n")
a <- brms_adapt_log_lik("hurdle_cumulative", 3, prep)
b <- brms:::log_lik_hurdle_cumulative(3, prep = prep)
cat("hurdle_cumulative, thresholds -1, 0, 1: max abs diff",
    format(max(abs(a - b)), digits = 3), "\n")
