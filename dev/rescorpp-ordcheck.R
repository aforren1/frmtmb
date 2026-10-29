## Lane rescorpp, fix round 1: what each reader does on a multivariate
## model with an ordinal response (data as rc_case("ord") in
## test-rescor-draws.R, seed 16), and loo()/loo_linpred() with several
## responses on the g3-like case.
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-common.R")
set.seed(16); n <- 120
x <- rnorm(n)
Z <- matrix(rnorm(3 * n), n)
dd <- data.frame(x = x, y2 = -1 + 0.5 * x + 1.5 * Z[, 2])
dd$o <- cut(x + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
f <- frm(mvbf(bf(o ~ x, family = cumulative()),
              bf(y2 ~ x, family = gaussian()), rescor = FALSE), data = dd)
ds <- suppressWarnings(suppressMessages(frm_sample(f, chains = 1,
  iter = 600, refresh = 0, seed = 4)))
show("posterior_linpred(ds) dim", dim(posterior_linpred(ds, ndraws = 5)))
show("posterior_epred(ds) dim", dim(posterior_epred(ds, ndraws = 5)))
show("fitted(ds) dim", dim(fitted(ds, ndraws = 5)))
show("predict(ds)", dim(predict(ds, ndraws = 5)))
show("predict(ds, summary = FALSE) dim",
     dim(predict(ds, ndraws = 5, summary = FALSE)))
show("residuals(ds)", dim(residuals(ds, ndraws = 5)))
show("loo(ds, resp = c('o', 'y2')) elpd",
     loo(ds, resp = c("o", "y2"))$estimates["elpd_loo", ])
show("loo_linpred(ds, resp = c('y2', 'o')) dim",
     suppressMessages(dim(loo_linpred(ds, resp = c("y2", "o")))))
