# Reviewer: frmtmb.sample on addition-term expressions against the
# precomputed columns. Data seed 101, sampler seed 3, 1 chain, 400 iter.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/formrobust-rev-stan-cache")
suppressMessages(library(frmtmb.sample))
same <- function(label, a, b) cat(sprintf("[%s] identical: %s\n", label,
                                          identical(a, b)))
set.seed(101)
n <- 60
d <- data.frame(x = rnorm(n), wt = runif(n, 0.5, 2), time = runif(n, 1, 3))
d$y <- rnorm(n, 0.5 * d$x)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * d$time)
d$w2 <- d$wt * 2; d$t2 <- d$time * 2; d$lbm <- min(d$y) - 1
s <- function(f, fam = gaussian()) suppressWarnings(suppressMessages(
  frm_sample(f, family = fam, data = d, chains = 1, iter = 400,
             refresh = 0, seed = 3)))
pairs <- list(
  weights = list(bf(y | weights(wt * 2) ~ x), bf(y | weights(w2) ~ x),
                 gaussian()),
  rate = list(bf(yc | rate(time * 2) ~ x), bf(yc | rate(t2) ~ x),
              poisson()),
  trunc = list(bf(y | trunc(lb = min(y) - 1) ~ x),
               bf(y | trunc(lb = lbm) ~ x), gaussian()))
nd <- d[1:8, ]
for (nm in names(pairs)) {
  p <- pairs[[nm]]
  a <- s(p[[1]], p[[3]]); b <- s(p[[2]], p[[3]])
  same(paste(nm, "draws"), as.matrix(a$draws), as.matrix(b$draws))
  same(paste(nm, "log_lik"), log_lik(a), log_lik(b))
  set.seed(1); pa <- posterior_predict(a, newdata = nd)
  set.seed(1); pb <- posterior_predict(b, newdata = nd)
  same(paste(nm, "posterior_predict(newdata)"), pa, pb)
  same(paste(nm, "posterior_epred(newdata)"), posterior_epred(a, newdata = nd),
       posterior_epred(b, newdata = nd))
  la <- tryCatch(loo(a)$estimates, error = function(e) conditionMessage(e))
  lb <- tryCatch(loo(b)$estimates, error = function(e) conditionMessage(e))
  same(paste(nm, "loo"), la, lb)
}
