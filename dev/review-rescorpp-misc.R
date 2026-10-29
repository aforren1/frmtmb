## Review rescorpp: name spellings, duplicated resp, a rescor model with
## NA rows, and core predict() of a rescor fit with sigma ~ x.
##   Rscript dev/review-rescorpp-misc.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
set.seed(61); n <- 100
x <- rnorm(n)
E <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .7, .7, 1), 2))
d <- data.frame(x = x, y_a = 1 + x + E[, 1], y_b = 2 * E[, 2])
d$x[c(3, 7)] <- NA
fit <- quiet(frm(mvbf(bf(y_a ~ x), bf(y_b ~ x)) + set_rescor(TRUE), data = d,
                 na.action = na.exclude))
ds <- quiet(frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 5))
cat("responses:", names(fit$spec$responses), " nobs", nobs(fit), "\n")
show("pp dimnames", dimnames(posterior_predict(ds, ndraws = 5)))
show("pp dim", dim(posterior_predict(ds, ndraws = 5)))
show("epred dim", dim(posterior_epred(ds, ndraws = 5)))
show("pp resp = 'ya' (brms spelling)", dim(posterior_predict(ds, resp = "ya", ndraws = 5)))
show("log_lik resp = 'ya' (brms spelling)", dim(log_lik(ds, resp = "ya", ndraws = 5)))
show("loo resp = 'ya'", loo(ds, resp = "ya")$estimates[1, ])
show("log_lik resp = 'y_a'", dim(log_lik(ds, resp = "y_a", ndraws = 5)))
show("pp resp = c('y_a','y_a')", dim(posterior_predict(ds, resp = c("y_a", "y_a"), ndraws = 5)))
show("pp resp = c('ya','y_a')", dim(posterior_predict(ds, resp = c("ya", "y_a"), ndraws = 5)))
show("log_lik resp = c('y_a','y_a') head", log_lik(ds, resp = c("y_a", "y_a"), ndraws = 2)[, 1:3])
show("log_lik resp = character(0)", dim(log_lik(ds, resp = character(0), ndraws = 2)))
show("pp resp = character(0)", dim(posterior_predict(ds, resp = character(0), ndraws = 2)))
show("fitted dimnames", dimnames(fitted(ds, ndraws = 5))[2:3])
show("predict resp ya dim", dim(predict(ds, resp = "ya", ndraws = 5)))
show("residuals no resp dim", dim(residuals(ds, ndraws = 5)))
## core predict on the ML fit of a sigma ~ x rescor model, resp subset
set.seed(62)
d2 <- data.frame(x = rnorm(300))
E2 <- MASS::mvrnorm(300, c(0, 0, 0), matrix(c(1, .8, .4, .8, 1, .3, .4, .3, 1), 3))
d2$y1 <- d2$x + exp(d2$x) * E2[, 1]
d2$y2 <- 3 * E2[, 2]
d2$y3 <- 0.1 * E2[, 3]
f2 <- quiet(frm(mvbf(bf(y1 ~ x, sigma ~ x), bf(y2 ~ x), bf(y3 ~ x)) +
                  set_rescor(TRUE), data = d2))
R2 <- rescor_matrix(f2)
set.seed(63)
p <- quiet(predict(f2, summary = FALSE, ndraws = 4000, propagate_error = FALSE))
cat("core predict summary=FALSE class", class(p), "names", names(p) %||% "-", "dim",
    if (is.list(p)) paste(lapply(p, dim), collapse = " ") else dim(p), "\n")
## residual correlation per row at the estimate
ep <- fitted(f2)
getE <- function(r) {
  pr <- if (is.list(p)) p[[r]] else p[, , r]
  sweep(pr, 2, as.numeric(if (length(dim(ep)) == 3) ep[, "Estimate", r] else ep[[r]][, "Estimate"]))
}
e1 <- getE("y1"); e2 <- getE("y2"); e3 <- getE("y3")
rc <- function(a, b) mean(vapply(seq_len(ncol(a)), function(i) cor(a[, i], b[, i]), 1))
cat(sprintf("core predict mean row cor: 12 %.4f (R %.4f)  13 %.4f (R %.4f)  23 %.4f (R %.4f)\n",
            rc(e1, e2), R2[1, 2], rc(e1, e3), R2[1, 3], rc(e2, e3), R2[2, 3]))
sdr <- apply(e1, 2, sd) / eval_dpars(f2)$y1$sigma
cat("core predict sd(e1)/sigma_i: range", round(range(sdr), 3), "mean", round(mean(sdr), 4), "\n")
