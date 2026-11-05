# Short check, item 3: predict() with theta held. Est.Error and intervals against
# ar1() and brms, on fitted rows and on new rows. Data: review's seeds 1, 3, 7, 19,
# plus a stronger latent with the same n, plus a larger n. Also interval coverage
# of a fresh latent path on new rows, by simulation from the truth.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
sim <- function(seed, s, ar = .4, G = 10, K = 6, b0 = .5) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- rnorm(n)
  d$ft <- factor(d$t)
  e <- as.vector(replicate(G, arima.sim(list(ar = ar), K, sd = s)))
  d$y <- rpois(n, exp(b0 + .3 * d$x + e))
  d
}
rat <- function(a, b) c(med = median(a / b), lo = unname(quantile(a / b, .1)), hi = unname(quantile(a / b, .9)))
cases <- list(c(1, .25), c(3, .25), c(7, .25), c(19, .25), c(11, .6), c(12, .6), c(21, 1), c(22, 1))
res <- list()
for (cc in cases) {
  d <- sim(cc[1], cc[2])
  fl <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  fr <- suppressWarnings(frm(y ~ x + ar1(ft + 0 | g), family = poisson(), data = d))
  set.seed(1)
  pl <- suppressWarnings(predict(fl, ndraws = 4000))
  pr <- suppressWarnings(predict(fr, ndraws = 4000))
  fb <- brms::bf(y ~ x + ar(t, g, cov = TRUE))
  code <- brms::stancode(fb, d, family = brms::brmsfamily("poisson"))
  sdat <- brms::standata(fb, d, family = brms::brmsfamily("poisson"))
  mod <- hel$brms_stan_model(code)
  sf <- suppressWarnings(rstan::sampling(mod, data = sdat, chains = 4, iter = 1500, warmup = 750, seed = 11, refresh = 0))
  scaf <- brms::brm(fb, d, family = brms::brmsfamily("poisson"), empty = TRUE)
  scaf$fit <- sf
  bfit <- brms::rename_pars(scaf)
  pb <- brms::posterior_predict(bfit)
  b_sd <- apply(pb, 2, sd)
  b_lo <- apply(pb, 2, quantile, .025)
  b_hi <- apply(pb, 2, quantile, .975)
  # brms re-sorts by g and t: put frm's rows into its order
  o <- order(d$g, d$t)
  el <- pl[o, "Est.Error"]; er <- pr[o, "Est.Error"]
  wl <- (pl[o, "Q97.5"] - pl[o, "Q2.5"]); wr <- (pr[o, "Q97.5"] - pr[o, "Q2.5"]); wb <- b_hi - b_lo
  # new rows, a new group
  nd <- data.frame(t = 1:6, g = factor(rep("new", 6)), x = 0, ft = factor(1:6, levels = levels(d$ft)))
  set.seed(2)
  nl <- suppressWarnings(predict(fl, newdata = nd, ndraws = 4000, allow_new_levels = TRUE))
  nr <- suppressWarnings(predict(fr, newdata = nd, ndraws = 4000, allow_new_levels = TRUE, sample_new_levels = "gaussian"))
  pbn <- brms::posterior_predict(bfit, newdata = nd, allow_new_levels = TRUE)
  th <- fl$estimates$theta
  vl <- sqrt(diag(as.matrix(vcov(fl, full = TRUE))))
  res[[length(res) + 1]] <- data.frame(seed = cc[1], s_true = cc[2], sderr = exp(th[1]), ar = th[2] / sqrt(1 + th[2]^2),
    se_log_s = vl["theta_1"], brms_sderr = median(rstan::extract(sf, "sderr")[[1]]),
    ee_ratio_ar1 = median(el / er), ee_ratio_brms = median(el / b_sd), ar1_ratio_brms = median(er / b_sd),
    width_ratio_ar1 = median(wl / wr), width_ratio_brms = median(wl / wb), ar1width_brms = median(wr / wb),
    new_ee_frm = mean(nl[, "Est.Error"]), new_ee_ar1 = mean(nr[, "Est.Error"]), new_ee_brms = mean(apply(pbn, 2, sd)),
    new_q975_frm = mean(nl[, "Q97.5"]), new_q975_ar1 = mean(nr[, "Q97.5"]), new_q975_brms = mean(apply(pbn, 2, quantile, .975)),
    new_mean_frm = mean(nl[, "Estimate"]), new_mean_brms = mean(pbn))
  print(res[[length(res)]], digits = 3, row.names = FALSE)
}
o <- do.call(rbind, res)
print(o, digits = 3, row.names = FALSE)
cat("median over cases: fitted rows, Est.Error latent/ar1", median(o$ee_ratio_ar1), " latent/brms", median(o$ee_ratio_brms), " ar1/brms", median(o$ar1_ratio_brms), "\n")
cat("interval width latent/ar1", median(o$width_ratio_ar1), " latent/brms", median(o$width_ratio_brms), " ar1/brms", median(o$ar1width_brms), "\n")
cat("new rows Est.Error frm/brms", median(o$new_ee_frm / o$new_ee_brms), " ar1/brms", median(o$new_ee_ar1 / o$new_ee_brms), " Q97.5 frm/brms", median(o$new_q975_frm / o$new_q975_brms), "\n")

# Coverage: 300 datasets, a fresh series from the truth, predict() on the new
# rows, how often does the row's 95% interval hold the new count?
covr <- function(seed, s) {
  d <- sim(seed, s, b0 = 1)
  f <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  nd <- data.frame(t = 1:6, g = factor(rep("new", 6)), x = 0)
  newy <- rpois(6, exp(1 + as.numeric(arima.sim(list(ar = .4), 6, sd = s))))
  p <- suppressWarnings(predict(f, newdata = nd, ndraws = 400, allow_new_levels = TRUE))
  # the same with the parameter draw of a free theta, capped at the 1e300 issue
  c(hit = mean(newy >= p[, "Q2.5"] & newy <= p[, "Q97.5"]), se = sqrt(diag(as.matrix(vcov(f, full = TRUE))))["theta_1"], mid = mean(p[, "Estimate"]))
}
for (s in c(.25, .6)) {
  v <- sapply(5001:5150, covr, s = s)
  cat("coverage on a new series, s =", s, ": mean hit", mean(v["hit", ]), " median se log s", median(v[2, ], na.rm = TRUE), "\n")
}
