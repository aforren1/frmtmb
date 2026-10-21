# Reviewer of lane loodraws: each refusal with its condition ABSENT, and
# the brms 2.23.0 call it stands for. brms side on the fits of
# dev/loodraws-rev-brmsfit.R; frmtmb side on the same draws
# (dev/loodraws-rev-same.R's construction).
#   Rscript dev/loodraws-rev-guards.R  -> dev/loodraws-rev-out/guards.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
options(mc.cores = 1)
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb.sample")), "\n")
B <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/loodraws-rev-out/brmsfits.rds")
dd <- B$dd
q <- function(e) suppressWarnings(suppressMessages(e))
tc <- function(e) tryCatch({
  v <- q(e)
  paste("ok", if (!is.null(dim(v))) paste(dim(v), collapse = "x") else
    paste(format(v, digits = 6), collapse = " "))
}, error = function(err) paste("ERROR:", substr(gsub("\n", " ",
  conditionMessage(err)), 1, 160)))
mk_ds <- function(bfit, fit) {
  labs <- frmtmb::brms_par_labels(fit)
  M <- as.matrix(bfit)
  if ("theta_1" %in% labs) M <- cbind(M, theta_1 = M[, "sd_g__Intercept"])
  structure(list(stanfit = bfit$fit,
                 draws = cbind(M[, labs, drop = FALSE], lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
g <- B$fits$gauss
bern <- B$fits$bern
ri <- B$fits$ri
dg <- mk_ds(g, q(frm(bf(y ~ x), family = gaussian(), data = dd)))
dbern <- mk_ds(bern, q(frm(bf(y01 ~ x), family = bernoulli(), data = dd)))
dri <- mk_ds(ri, q(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)))
row <- function(label, b, f) {
  cat(sprintf("%-46s brms: %-34s frmtmb: %s\n", label, substr(b, 1, 34), f))
}

cat("\n== brms's default values of the refused arguments\n")
row("loo_epred(allow_new_levels = FALSE)",
    tc(loo_epred(ri, allow_new_levels = FALSE)),
    tc(loo_epred(dri, allow_new_levels = FALSE)))
row("loo_epred(sample_new_levels = 'uncertainty')",
    tc(loo_epred(ri, sample_new_levels = "uncertainty")),
    tc(loo_epred(dri, sample_new_levels = "uncertainty")))
row("loo_epred(re_formula = NULL)", tc(loo_epred(ri, re_formula = NULL)),
    tc(loo_epred(dri, re_formula = NULL)))
row("loo_epred(re_formula = NA)", tc(loo_epred(ri, re_formula = NA)),
    tc(loo_epred(dri, re_formula = NA)))
row("loo(allow_new_levels = FALSE) [pre-existing]",
    tc(loo(ri, allow_new_levels = FALSE)$estimates[1, 1]),
    tc(loo(dri, allow_new_levels = FALSE)$estimates[1, 1]))
row("loo_R2(allow_new_levels = FALSE)",
    tc(loo_R2(ri, allow_new_levels = FALSE, seed = 1)),
    tc(loo_R2(dri, allow_new_levels = FALSE, seed = 1)))
row("loo_model_weights(newdata = NULL)",
    tc(loo_model_weights(g, ri, newdata = NULL)),
    tc(loo_model_weights(dg, dri, newdata = NULL)))
row("loo_model_weights(allow_new_levels = FALSE)",
    tc(loo_model_weights(g, ri, allow_new_levels = FALSE)),
    tc(loo_model_weights(dg, dri, allow_new_levels = FALSE)))
ids <- 1:1000
row("loo_model_weights(draw_ids = 1:1000)",
    tc(loo_model_weights(g, ri, draw_ids = ids)),
    tc(loo_model_weights(dg, dri, draw_ids = ids)))
row("loo_epred(draw_ids = 1:1000)", tc(loo_epred(g, draw_ids = ids)),
    tc(loo_epred(dg, draw_ids = ids)))
row("loo_model_weights(ndraws = 1000)",
    tc(loo_model_weights(g, ri, ndraws = 1000)),
    tc(loo_model_weights(dg, dri, ndraws = 1000)))
row("loo_model_weights(r_eff_list = )",
    tc(loo_model_weights(g, ri, r_eff_list = list(1, 1))),
    tc(loo_model_weights(dg, dri, r_eff_list = list(1, 1))))
row("loo_model_weights(method='pseudobma', BB=F)",
    tc(loo_model_weights(g, ri, method = "pseudobma", BB = FALSE)),
    tc(loo_model_weights(dg, dri, method = "pseudobma", BB = FALSE)))
row("loo_epred(save_psis = TRUE)", tc(loo_epred(g, save_psis = TRUE)),
    tc(loo_epred(dg, save_psis = TRUE)))
row("loo_epred(pointwise = TRUE)", tc(loo_epred(g, pointwise = TRUE)),
    tc(loo_epred(dg, pointwise = TRUE)))
row("loo_epred(cores = 1)", tc(loo_epred(g, cores = 1)),
    tc(loo_epred(dg, cores = 1)))
row("loo_epred(type = 'quantile', probs = 0)",
    tc(loo_epred(g, type = "quantile", probs = 0)),
    tc(loo_epred(dg, type = "quantile", probs = 0)))
row("loo_epred(type = 'mean', probs = 2) [ignored]",
    tc(loo_epred(g, probs = 2)), tc(loo_epred(dg, probs = 2)))
row("loo_predictive_interval(prob = 0.5)",
    tc(loo_predictive_interval(g, prob = 0.5)),
    tc(loo_predictive_interval(dg, prob = 0.5)))
row("loo_R2(seed = 1, junk = 1)", tc(loo_R2(g, seed = 1, junk = 1)),
    tc(loo_R2(dg, seed = 1, junk = 1)))

cat("\n== scale on a non-identity link (bernoulli)\n")
s_resp <- q(loo_linpred(bern, scale = "response"))
s_def <- q(loo_linpred(bern))
s_tr <- q(loo_linpred(bern, transform = TRUE))
cat("brms loo_linpred(scale = 'response') identical to default:",
    identical(s_resp, s_def), "; to transform = TRUE:",
    identical(s_resp, s_tr), "; range", signif(range(s_resp), 3), "\n")
cat("frmtmb loo_linpred(scale = 'response'):",
    tc(loo_linpred(dbern, scale = "response")), "\n")

cat("\n== responses with an underscore: frmtmb's name and brms's\n")
dd2 <- dd
dd2$y_a <- dd2$y
fm <- q(frm(mvbf(bf(y_a ~ x), bf(y2 ~ x), rescor = FALSE), data = dd2))
bt <- brms::brmsterms(bf(y_a ~ x) + bf(y2 ~ x) + set_rescor(FALSE))
cat("brms response names:", bt$responses, "\n")
set.seed(5)
tpl <- fm$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fm$estimates[[cp]]))
M <- matrix(rep(est, each = 200) + rnorm(200 * length(est), 0, 0.03), 200,
            dimnames = list(NULL, frmtmb::brms_par_labels(fm)))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fm), lp__ = 0)
dmv <- structure(list(stanfit = NULL, draws = M, fit = fm),
                 class = "frmtmb_draws")
cat("loo_epred(resp = 'ya'):", tc(loo_epred(dmv, resp = "ya")), "\n")
cat("loo_epred(resp = 'y_a'):", tc(loo_epred(dmv, resp = "y_a")), "\n")
cat("loo_R2(resp = 'ya'):", tc(loo_R2(dmv, resp = "ya", seed = 1)), "\n")
cat("loo_R2(resp = 'y_a'):", tc(loo_R2(dmv, resp = "y_a", seed = 1)), "\n")
cat("loo_epred() dimnames[[3]]:",
    dimnames(q(loo_epred(dmv)))[[3]], "\n")
cat("colnames(loo_R2()):", rownames(q(loo_R2(dmv, seed = 1))), "\n")
cat("posterior_epred(resp = 'ya'):",
    tc(posterior_epred(dmv, resp = "ya")), "\n")
cat("bayes_R2(resp = 'y_a'):", tc(bayes_R2(dmv, resp = "y_a")), "\n")

cat("\n== add_criterion with two criteria and a loo_R2 argument\n")
cat("brms add_criterion(c('loo', 'loo_R2'), seed = 1):",
    tc(names(add_criterion(g, c("loo", "loo_R2"), seed = 1)$criteria)),
    "\n")
cat("frmtmb add_criterion(c('loo', 'loo_R2'), seed = 1):",
    tc(names(add_criterion(dg, c("loo", "loo_R2"), seed = 1)$criteria)),
    "\n")

cat("\n== brms posterior_average/model_weights with a stored loo\n")
n_ll <- 0
suppressMessages(trace(brms:::log_lik.brmsfit, print = FALSE,
                       tracer = quote(n_ll <<- n_ll + 1)))
g2 <- q(add_criterion(g, "loo", ndraws = 500))
r2 <- q(add_criterion(ri, "loo", ndraws = 500))
n_ll <- 0
w_st <- q(model_weights(g2, r2, weights = "stacking"))
cat("log_lik calls in model_weights(stacking) with stored loo:", n_ll, "\n")
w_fresh <- q(model_weights(g, ri, weights = "stacking"))
cat("weights with stored loo", format(w_st, digits = 8), "| without",
    format(w_fresh, digits = 8), "| identical", identical(unname(w_st),
                                                          unname(w_fresh)),
    "\n")
suppressMessages(untrace(brms:::log_lik.brmsfit))

cat("\n== frm_sample(thin = 2): the chain structure of real thinned draws\n")
st <- tryCatch({
  d <- q(frm_sample(bf(y ~ x), family = gaussian(), data = dd, chains = 2,
                    iter = 1200, warmup = 400, thin = 2, seed = 1,
                    refresh = 0))
  ll <- log_lik(d)
  cid <- attr(ll, "chain_id")
  ps <- q(psis(d))
  rh <- loo::relative_eff(exp(ll), chain_id = rep(1:2, each = nrow(ll) / 2))
  sprintf("draws %d, nchains %d, chain_id %s, psis r_eff == two-chain relative_eff %s",
          nrow(d$draws), frmtmb.sample:::draws_nchains(d),
          if (is.null(cid)) "NULL" else paste(table(cid), collapse = "/"),
          isTRUE(all.equal(attr(ps, "r_eff"), rh)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(st, "\n")
cat("done\n")
