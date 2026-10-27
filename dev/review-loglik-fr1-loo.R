# Reviewer of lane loglik, fix round 1: loo(), waic() and loo_epred() of
# frmtmb.sample on an unstr() block against the loo package on the same
# log_lik() matrix. Seed 111. Output: dev/review-loglik-log/fr1-loo.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(111)
d <- data.frame(g = rep(1:15, each = 4), t = rep(1:4, 15))
d <- d[!(d$g == 3 & d$t == 2), ]
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d))
d <- d[sample(nrow(d)), ]
fit <- frm(y ~ x + unstr(t, g), data = d, dry_run = "objective")
lab <- c(brms_par_labels(fit), "lp__")
K <- 300L
m <- matrix(0, K, length(lab), dimnames = list(NULL, lab))
m[, "b_Intercept"] <- rnorm(K, 1, 0.1)
m[, "b_x"] <- rnorm(K, 0.5, 0.1)
m[, "sigma"] <- exp(rnorm(K, 0, 0.05))
for (j in grep("^thetaac", lab)) m[, j] <- rnorm(K, 0.2, 0.05)
ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                class = "frmtmb_draws")
ll <- log_lik(ds)
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
lo <- suppressWarnings(loo(ds))
lr <- suppressWarnings(loo::loo(ll, r_eff = loo::relative_eff(
  exp(ll), chain_id = rep(1L, K))))
cat(sprintf("loo estimates rel %.3e, pointwise elpd rel %.3e\n",
            rel(lo$estimates, lr$estimates),
            rel(lo$pointwise[, "elpd_loo"], lr$pointwise[, "elpd_loo"])))
wa <- suppressWarnings(waic(ds))
wr <- suppressWarnings(loo::waic(ll))
cat(sprintf("waic estimates rel %.3e\n", rel(wa$estimates, wr$estimates)))
ep <- posterior_epred(ds)
ps <- suppressWarnings(loo::psis(-ll, r_eff = loo::relative_eff(
  exp(ll), chain_id = rep(1L, K))))
le <- suppressWarnings(suppressMessages(loo_epred(ds)))
ler <- loo::E_loo(ep, ps, type = "mean", log_ratios = -ll)$value
cat(sprintf("loo_epred vs E_loo rel %.3e\n", rel(le, ler)))
