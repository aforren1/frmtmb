## Review rescorpp, short check of fix round 1.
##   Rscript dev/review-rescorpp-short.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
source(file.path(wt, "dev/rescorpp-data.R"))
source(file.path(wt, "dev/review-rescorpp-data.R"))
suppressMessages(loadNamespace("brms"))
tq <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)),
                              error = function(e) paste("ERROR:", conditionMessage(e)))
dsc <- function(x) if (is.character(x)) substr(x[1], 1, 160) else
  paste0("dim=", paste(dim(x) %||% length(x), collapse = "x"),
         " dn3=", paste(if (length(dim(x)) >= 3) dimnames(x)[[3]] else "-", collapse = ","))
pair <- function(tag, fb, ff) cat(sprintf("%-34s\n   brms  : %s\n   frmtmb: %s\n", tag,
                                          dsc(fb), dsc(ff)))
samp <- function(f) quiet(frm_sample(f, chains = 2, iter = 1000, refresh = 0, seed = 21))

cat("\n#### 1. ordinal + gaussian predict()\n")
dord <- rv_data("ord")
bo <- readRDS(file.path(rv_out, "brms-ord.rds"))
dso <- samp(frm(rv_frm_form("ord"), data = dord))
for (cl in list(list(), list(resp = c("o1", "y2")), list(resp = c("y2", "o1")))) {
  tg <- paste(cl$resp %||% "none", collapse = ",")
  pair(paste("predict summary resp", tg),
       tq(do.call(predict, c(list(bo, ndraws = 30), cl))),
       tq(do.call(predict, c(list(dso, ndraws = 30), cl))))
  set.seed(1)
  pair(paste("predict summary=FALSE resp", tg),
       tq(do.call(predict, c(list(bo, ndraws = 30, summary = FALSE), cl))),
       tq(do.call(predict, c(list(dso, ndraws = 30, summary = FALSE), cl))))
}
set.seed(2); a <- predict(dso, ndraws = 30, summary = FALSE)
set.seed(2); b <- posterior_predict(dso, ndraws = 30)
cat("predict(summary = FALSE) identical to posterior_predict():", identical(a, b),
    "; codes in o1 layer:", paste(sort(unique(as.vector(a[, , "o1"]))), collapse = ","), "\n")
pair("predict resp o1", tq(predict(bo, resp = "o1", ndraws = 30)),
     tq(predict(dso, resp = "o1", ndraws = 30)))

cat("\n#### 3. resp spellings, duplicates, empty, against brms\n")
dus <- rv_data("us")
bu <- readRDS(file.path(rv_out, "brms-us.rds"))
dsu <- samp(frm(rv_frm_form("us"), data = dus))
g2 <- rescorpp_data("g3")
b2f <- readRDS(file.path(wt, "dev/rescorpp-brmsfits/g2f.rds"))
ds2f <- samp(frm(eval(rescorpp_frm_forms[["g2f"]]), family = "gaussian", data = g2))
b3 <- readRDS(file.path(wt, "dev/rescorpp-brmsfits/g3.rds"))
ds3 <- samp(frm(eval(rescorpp_frm_forms[["g3"]]), family = "gaussian", data = g2))
for (rs in list("ya", "y_a", c("ya", "ya"), c("ya", "y_a"), character(0))) {
  tg <- paste0("c(", paste0("'", rs, "'", collapse = ","), ")")
  pair(paste("us log_lik", tg), tq(brms::log_lik(bu, resp = rs, ndraws = 20)),
       tq(log_lik(dsu, resp = rs, ndraws = 20)))
  pair(paste("us posterior_predict", tg),
       tq(brms::posterior_predict(bu, resp = rs, ndraws = 20)),
       tq(posterior_predict(dsu, resp = rs, ndraws = 20)))
}
for (rs in list(c("y1", "y1"), character(0))) {
  tg <- paste0("c(", paste0("'", rs, "'", collapse = ","), ")")
  for (m in c("log_lik", "posterior_predict", "posterior_epred", "fitted")) {
    pair(paste("g2f (no rescor)", m, tg),
         tq(do.call(m, list(b2f, resp = rs, ndraws = 20))),
         tq(do.call(m, list(ds2f, resp = rs, ndraws = 20))))
  }
}
## legitimate calls that pass through draws_resp_set()
leg <- list(
  pp_check_y2 = function() class(pp_check(ds3, resp = "y2", ndraws = 3))[1],
  bayes_R2_two = function() rownames(bayes_R2(ds3, resp = c("y1", "y2"))),
  loo_R2 = function() rownames(loo_R2(ds3)),
  residuals_pearson = function() dim(residuals(ds3, type = "pearson", ndraws = 10)),
  perr_epred_two = function() dim(predictive_error(ds3, method = "posterior_epred",
                                                   resp = c("y3", "y2"), ndraws = 10)),
  pint_y3 = function() dim(predictive_interval(ds3, resp = "y3", ndraws = 10)),
  ce_y2 = function() names(conditional_effects(ds3, resp = "y2", ndraws = 10)),
  loo_none = function() round(loo(ds3)$estimates[1, 1], 3),
  waic_two = function() round(waic(ds3, resp = c("y1", "y2"))$estimates[1, 1], 3),
  loo_pred_none = function() dim(loo_predict(ds3)),
  fitted_sigma = function() dim(fitted(ds3, dpar = "sigma", ndraws = 10)),
  linpred_all = function() dim(posterior_linpred(ds3, resp = c("y1", "y2", "y3"),
                                                 ndraws = 10)))
for (k in names(leg)) cat(sprintf("legit %-18s %s\n", k, paste(tq(leg[[k]]()), collapse = " ")))
## a univariate model with a resp argument (brms ignores it)
d1 <- data.frame(x = g2$x, y = g2$y1)
ds1 <- quiet(frm_sample(frm(bf(y ~ x), data = d1), chains = 1, iter = 400,
                        refresh = 0, seed = 3))
for (rs in list("y", "foo", c("y", "y"), character(0))) {
  cat(sprintf("univariate pp resp=%-10s %s | log_lik %s\n", paste(rs, collapse = ","),
              dsc(tq(posterior_predict(ds1, resp = rs, ndraws = 5))),
              dsc(tq(log_lik(ds1, resp = rs, ndraws = 5)))))
}

cat("\n#### 4. predict_rescor_draw() by name\n")
fit3 <- quiet(frm(eval(rescorpp_frm_forms[["g3"]]), family = "gaussian", data = g2))
dp <- eval_dpars(fit3)
ok <- rep(TRUE, nrow(g2))
set.seed(1)
a <- predict_rescor_draw(fit3, fit3$spec$responses[c("y3", "y1")], dp, ok)
cat("subset call: names", names(a), " mean y3", round(mean(a$y3), 3), "(mu_y3",
    round(mean(dp$y3$mu), 3), ") sd(y3 - mu_y3)", round(sd(a$y3 - dp$y3$mu), 3),
    "(sigma_y3", round(dp$y3$sigma[1], 3), ")\n")
set.seed(1)
b <- predict_rescor_draw(fit3, fit3$spec$responses[c("y3", "y1")], dp[c("y3", "y1")], ok)
cat("identical to the call with dps already in that order:", identical(a, b), "\n")
cat("missing name:", tq(predict_rescor_draw(fit3, fit3$spec$responses, dp[c("y1", "y2")], ok)), "\n")
cat("unnamed, wrong length:", tq(predict_rescor_draw(fit3, fit3$spec$responses,
                                                     unname(dp[1:2]), ok)), "\n")

cat("\n#### 5. loo_* with several responses: layers and weights by hand\n")
hand <- function(obj, rs, pred, chain) {
  ll <- log_lik(obj, resp = rs)
  r_eff <- loo::relative_eff(exp(ll), chain_id = chain)
  ps <- suppressWarnings(loo::psis(-ll, r_eff = r_eff))
  vapply(rs, function(r) suppressWarnings(loo::E_loo(pred(obj, r), ps, type = "mean",
                                                     log_ratios = -ll)$value),
         numeric(ncol(ll)))
}
for (rs in list(c("y3", "y1"), c("y1", "y3"), c("y1", "y2", "y3"))) {
  ff <- suppressWarnings(suppressMessages(loo_epred(ds3, resp = rs)))
  hf <- hand(ds3, rs, function(o, r) posterior_epred(o, resp = r), rep(1:2, each = 500))
  # brms 2.23.0 stops on a permuted rescor resp; compare it in model order
  rb <- sort(rs)
  fb <- suppressWarnings(suppressMessages(brms::loo_epred(b3, resp = rb)))
  hb <- hand(b3, rb, function(o, r) brms::posterior_epred(o, resp = r),
             rep(1:2, each = 1000))
  cat(sprintf("loo_epred resp=%-8s frmtmb dim %s dn %s, max|frmtmb - hand| %.3g; brms dim %s dn %s, max|brms - hand| %.3g\n",
              paste(rs, collapse = ","), paste(dim(ff), collapse = "x"),
              paste(dimnames(ff)[[3]], collapse = ","),
              max(abs(ff[, 1, ] - hf)), paste(dim(fb), collapse = "x"),
              paste(dimnames(fb)[[3]], collapse = ","), max(abs(fb[, 1, ] - hb))))
}
## the joint weights differ from per-response weights, so the check bites
ff <- suppressWarnings(suppressMessages(loo_epred(ds3, resp = c("y3", "y1"))))
h1 <- hand(ds3, "y1", function(o, r) posterior_epred(o, resp = r), rep(1:2, each = 500))
cat("y1 layer, joint weights vs y1-only weights: max diff",
    signif(max(abs(ff[, 1, 2] - h1)), 3), "\n")
lp <- suppressWarnings(suppressMessages(loo_linpred(ds3, resp = c("y2", "y3"))))
lpi <- suppressWarnings(suppressMessages(loo_predictive_interval(ds3, resp = c("y2", "y3"))))
lpr <- suppressWarnings(suppressMessages(loo_predict(ds3, resp = c("y2", "y3"))))
cat("loo_linpred", dim(lp), dimnames(lp)[[3]], "| loo_predictive_interval", dim(lpi),
    dimnames(lpi)[[2]], dimnames(lpi)[[3]], "| loo_predict", dim(lpr), "\n")
bl <- suppressWarnings(suppressMessages(brms::loo_predictive_interval(b3, resp = c("y2", "y3"))))
cat("brms loo_predictive_interval", dim(bl), dimnames(bl)[[2]], dimnames(bl)[[3]], "\n")

cat("\n#### 2. tail-share statistic across predictive seeds (st case of the test)\n")
set.seed(13); n <- 120
x <- rnorm(n)
Z <- matrix(rnorm(3 * n), n) %*% chol(matrix(c(1, .8, .5, .8, 1, .3, .5, .3, 1), 3))
dd <- data.frame(x = x)
Z <- Z / sqrt(rchisq(n, 6) / 6)
dd$y1 <- 1 + x + Z[, 1]; dd$y2 <- -1 + 0.5 * x + 1.5 * Z[, 2]
fst <- quiet(frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), family = student(),
                 data = dd))
zs <- NULL
for (ss in 4:5) {
  dss <- quiet(frm_sample(fst, chains = 1, iter = 1200, refresh = 0, seed = ss))
  dr <- as.data.frame(posterior::as_draws_df(dss))
  ep <- posterior_epred(dss)
  for (s in 1:25) {
    set.seed(s)
    e <- posterior_predict(dss) - ep
    for (j in 1:2) {
      u <- e[, , j] / dr[[paste0("sigma_y", j)]]
      hit <- abs(u) > stats::qt(0.995, dr$nu)
      zs <- rbind(zs, data.frame(samp = ss, seed = s, j = j,
                                 z = (mean(hit) - 0.01) / sqrt(0.01 * 0.99 / length(hit))))
    }
  }
}
cat("tail z over", nrow(zs), "runs (2 sampler seeds x 25 predictive seeds x 2 margins):",
    "mean", round(mean(zs$z), 3), "sd", round(sd(zs$z), 3), "max |z|",
    round(max(abs(zs$z)), 3), "; share |z| > 3:", mean(abs(zs$z) > 3), "\n")
