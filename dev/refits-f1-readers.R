# F1 of the review: on a bernoulli fit coded 1/2 (1 = failure), does every
# reader of a new response take 1 as a failure? Seed 9, n = 40 (the
# reviewer's probe5 data). brms 2.23.0: standata() of an empty fit on new
# data, which codes against the stored resp_levels.
source("dev/refits-lib.R")
set.seed(9)
n <- 40
x <- rnorm(n)
s <- rbinom(n, 1, plogis(0.3 + x))
d12 <- data.frame(x = x, y = s + 1)
d01 <- data.frame(x = x, y = s)
f12 <- frm(y ~ x, data = d12, family = bernoulli())
f01 <- frm(y ~ x, data = d01, family = bernoulli())
nd12 <- transform(d12, y = 1)
nd01 <- transform(d01, y = 0)
v <- rep(1, n); v[1] <- 2
go("refit: c(2, 1, 1, ...) on f12 = c(1, 0, 0, ...) on f01",
   identical(fixef(refit(f12, v)), fixef(refit(f01, v - 1))))
go("refit: rep(1, n) on f12", tryCatch(refit(f12, rep(1, n)),
   error = function(e) "refused (reads two ways)"))
go("residuals(newdata y = 1) on f12 = (y = 0) on f01",
   identical(residuals(f12, newdata = nd12, type = "response"),
             residuals(f01, newdata = nd01, type = "response")))
go("log_lik(newdata y = 1) on f12 = (y = 0) on f01", {
  suppressPackageStartupMessages(library(frmtmb.sample))
  set.seed(1); a <- log_lik(f12, newdata = nd12, ndraws = 50)
  set.seed(1); b <- log_lik(f01, newdata = nd01, ndraws = 50)
  identical(unname(a), unname(b))
})
go("bayes_R2(newdata y = 1) on f12 = (y = 0) on f01", {
  set.seed(2); a <- bayes_R2(f12, newdata = nd12)
  set.seed(2); b <- bayes_R2(f01, newdata = nd01)
  identical(unname(a), unname(b))
})
go("pp_check(newdata y = 1): the y it plots", {
  set.seed(3); p <- pp_check(f12, newdata = nd12, ndraws = 5)
  unique(p$data$value[p$data$is_y])
})
go("brms 2.23.0 standata(empty fit, newdata y = 1)$Y", {
  suppressPackageStartupMessages(library(brms))
  b12 <- brm(y ~ x, data = d12, family = bernoulli(), empty = TRUE)
  unique(standata(b12, newdata = nd12)$Y)
})
## the draws readers, on frmtmb.sample draws of each fit (same seed)
suppressPackageStartupMessages(library(frmtmb.sample))
s12 <- go("frm_sample f12", suppressWarnings(frm_sample(
  f12, chains = 1, iter = 400, warmup = 200, seed = 4, refresh = 0)))
s01 <- go("frm_sample f01", suppressWarnings(frm_sample(
  f01, chains = 1, iter = 400, warmup = 200, seed = 4, refresh = 0)))
go("draws equal", identical(as.matrix(s12), as.matrix(s01)))
go("log_lik(draws, newdata y = 1) on f12 = (y = 0) on f01",
   identical(unname(log_lik(s12, newdata = nd12)),
             unname(log_lik(s01, newdata = nd01))))
go("log_lik(draws, newdata y = 1) on f12 = (y = 1) on f01",
   identical(unname(log_lik(s12, newdata = nd12)),
             unname(log_lik(s01, newdata = transform(d01, y = 1)))))
go("bayes_R2(draws, newdata) on f12 = on f01", {
  nd12b <- d12; nd12b$y[1:5] <- 1
  nd01b <- d01; nd01b$y[1:5] <- 0
  set.seed(2); a <- bayes_R2(s12, newdata = nd12b)
  set.seed(2); b <- bayes_R2(s01, newdata = nd01b)
  identical(unname(a), unname(b))
})
