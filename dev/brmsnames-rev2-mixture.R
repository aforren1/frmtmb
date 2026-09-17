## Reviewer recheck, MAJOR 3: a two-component mixture's `theta1`. In brms
## `theta1` is the mixing PROBABILITY of component 1. What does the lane's
## `theta1` hold, on the fit, on draws and in frm_simulate(newparams =)?
##   Rscript dev/brmsnames-rev2-mixture.R base|lane|brms
## Data seed 64: 90% of rows from component 1, so theta1 is about 0.9,
## whose log odds (2.2) cannot be mistaken for it.
arm <- commandArgs(trailingOnly = TRUE)[1L]
libs <- list(
  base = c("C:/Users/adf44/source/r/rellib-r3"),
  lane = c("C:/Users/adf44/source/r/brmsnames-lib",
           "C:/Users/adf44/source/r/rellib-r3"),
  brms = character(0))
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
set.seed(64)
n <- 400
z <- runif(n) < 0.9
d <- data.frame(x = rnorm(n))
d$y <- ifelse(z, -2 + 0.3 * d$x + rnorm(n, 0, 0.8),
              2.5 + 0.3 * d$x + rnorm(n, 0, 1.2))
cat("empirical share of component 1 (mean -2):", mean(z), "\n")
if (arm == "brms") {
  stopifnot(packageVersion("StanHeaders") == "2.32.10")
  q(library(brms))
  pri <- c(set_prior("normal(-2, 1)", class = "Intercept", dpar = "mu1"),
           set_prior("normal(2.5, 1)", class = "Intercept", dpar = "mu2"))
  b <- brm(bf(y ~ x), data = d, family = mixture(gaussian, gaussian),
           prior = pri, chains = 1, iter = 400, warmup = 200, seed = 1,
           refresh = 0, backend = "rstan")
  saveRDS(b, "dev/stan-cache/brmsnames-rev2-brms-mixture.rds")
  print(variables(b))
  print(posterior_summary(b, variable = c("theta1", "theta2", "sigma1",
                                          "sigma2")))
  print(hypothesis(b, "theta1 = 0.5", class = NULL)$hypothesis)
} else {
  q(library(frmtmb)); q(library(frmtmb.sample))
  fit <- q(frm(bf(y ~ x), family = mixture(gaussian(), gaussian()),
               data = d))
  cat("== estimates ==\n"); print(fit$estimates$betad)
  cat("variables(fit):", try1(variables(fit)), "\n")
  cat("== hypothesis(fit, theta1 = 0.5) ==\n")
  print(try1(if (arm == "lane") {
    hypothesis(fit, "theta1 = 0.5", class = NULL)$hypothesis
  } else hypothesis(fit, "theta1 = 0.5")))
  set.seed(5)
  ds <- q(frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 3))
  cat("== draws column means ==\n"); print(colMeans(ds$draws))
  cat("== posterior mixing probability of component 1, pp_mixture ==\n")
  pm <- try1(q(pp_mixture(ds)))
  if (is.array(pm)) print(colMeans(pm[, 1, ])) else print(pm)
  if (arm == "lane") {
    cat("== hypothesis(ds, theta1 = 0.5, class = NULL) ==\n")
    print(try1(q(hypothesis(ds, "theta1 = 0.5", class = NULL)$hypothesis)))
    cat("== posterior_summary(ds, variable = theta1) ==\n")
    print(try1(q(posterior_summary(ds, variable = "theta1"))))
    cat("== frm_simulate(newparams = list(theta1 = 0.7)) share below 0 ==\n")
    np <- as.list(q(sapply(c("b_mu1_Intercept", "b_mu1_x", "b_mu2_Intercept",
                             "b_mu2_x"), function(v) 0)))
    np$b_mu1_Intercept <- -5; np$b_mu2_Intercept <- 5
    np$sigma1 <- 0.5; np$sigma2 <- 0.5; np$theta1 <- 0.7
    sim <- try1(q(frm_simulate(bf(y ~ x), family = mixture(gaussian(),
                                                          gaussian()),
                               data = d, newparams = np, seed = 2)))
    if (is.data.frame(sim)) print(mean(sim[[1]] < 0)) else print(sim)
    cat("== the same with theta1 = qlogis(0.7), read as a log ratio ==\n")
    np$theta1 <- stats::qlogis(0.7)
    sim <- try1(q(frm_simulate(bf(y ~ x), family = mixture(gaussian(),
                                                          gaussian()),
                               data = d, newparams = np, seed = 2)))
    if (is.data.frame(sim)) print(mean(sim[[1]] < 0)) else print(sim)
  }
}
