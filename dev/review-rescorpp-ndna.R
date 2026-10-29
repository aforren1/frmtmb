## Review rescorpp: newdata with a missing predictor that only one
## response reads, under rescor: do the joint layers stay row-aligned?
##   Rscript dev/review-rescorpp-ndna.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
suppressMessages(loadNamespace("brms"))
set.seed(81); n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n))
E <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .7, .7, 1), 2))
d$y1 <- 1 + 2 * d$x + E[, 1]
d$y2 <- -1 + 3 * d$z + 0.3 * E[, 2]
fit <- quiet(frm(mvbf(bf(y1 ~ x), bf(y2 ~ z)) + set_rescor(TRUE), data = d))
ds <- quiet(frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 2))
nd <- data.frame(x = c(-2, 0, 2, 4), z = c(1, NA, -1, 2))
show("epred y1 resp", colMeans(posterior_epred(ds, newdata = nd, resp = "y1")))
show("epred y2 resp", colMeans(posterior_epred(ds, newdata = nd, resp = "y2")))
show("epred joint", apply(posterior_epred(ds, newdata = nd), c(2, 3), mean))
set.seed(1)
pj <- tryCatch(posterior_predict(ds, newdata = nd), error = function(e) e)
show("pp joint col means", if (inherits(pj, "error")) pj else apply(pj, c(2, 3), mean))
show("pp joint dim", if (inherits(pj, "error")) conditionMessage(pj) else dim(pj))
show("pp y2 alone col means", colMeans(posterior_predict(ds, newdata = nd, resp = "y2")))
show("core predict", predict(fit, newdata = nd))
