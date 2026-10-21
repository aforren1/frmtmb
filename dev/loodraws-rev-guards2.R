# Reviewer of lane loodraws: the second half of dev/loodraws-rev-guards.R
# (response spellings, add_criterion with two criteria, brms's
# model_weights with a stored loo, real thinned draws).
#   Rscript dev/loodraws-rev-guards2.R -> dev/loodraws-rev-out/guards2.txt
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
ri <- B$fits$ri
dg <- mk_ds(g, q(frm(bf(y ~ x), family = gaussian(), data = dd)))

cat("\n== responses with an underscore: frmtmb's name and brms's\n")
dd2 <- dd
dd2$y_a <- dd2$y
fm <- q(frm(mvbf(bf(y_a ~ x), bf(y2 ~ x), rescor = FALSE), data = dd2))
bmv <- brms::mvbf(brms::bf(y_a ~ x), brms::bf(y2 ~ x), rescor = FALSE)
cat("brms response names:", bmv$responses, "\n")
cat("brms stancode names the response:",
    any(grepl("Y_ya", brms::stancode(bmv, data = dd2))), "\n")
set.seed(5)
tpl <- fm$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fm$estimates[[cp]]))
M <- matrix(rep(est, each = 200) + rnorm(200 * length(est), 0, 0.03), 200,
            dimnames = list(NULL, frmtmb::brms_par_labels(fm)))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fm), lp__ = 0)
dmv <- structure(list(stanfit = NULL, draws = M, fit = fm),
                 class = "frmtmb_draws")
cat("frmtmb draw labels:", colnames(M), "\n")
cat("loo_epred(resp = 'ya'):", tc(loo_epred(dmv, resp = "ya")), "\n")
cat("loo_epred(resp = 'y_a'):", tc(loo_epred(dmv, resp = "y_a")), "\n")
cat("loo_R2(resp = 'ya'):", tc(loo_R2(dmv, resp = "ya", seed = 1)), "\n")
cat("loo_R2(resp = 'y_a'):", tc(loo_R2(dmv, resp = "y_a", seed = 1)), "\n")
cat("loo_epred() dimnames[[3]]:", dimnames(q(loo_epred(dmv)))[[3]], "\n")
cat("rownames(loo_R2()):", rownames(q(loo_R2(dmv, seed = 1))), "\n")
cat("posterior_epred(resp = 'ya'):", tc(posterior_epred(dmv, resp = "ya")),
    "\n")
cat("posterior_epred(resp = 'y_a'):",
    tc(posterior_epred(dmv, resp = "y_a")), "\n")
cat("log_lik(resp = 'ya'):", tc(log_lik(dmv, resp = "ya")), "\n")
cat("loo(resp = 'ya'):", tc(loo(dmv, resp = "ya")$estimates[1, 1]), "\n")
cat("bayes_R2(resp = 'ya'):", tc(bayes_R2(dmv, resp = "ya")), "\n")
cat("bayes_R2(resp = 'y_a'):", tc(bayes_R2(dmv, resp = "y_a")), "\n")

cat("\n== add_criterion with two criteria and a loo_R2 argument\n")
cat("brms add_criterion(c('loo', 'loo_R2'), seed = 1):",
    tc(names(add_criterion(g, c("loo", "loo_R2"), seed = 1)$criteria)),
    "\n")
cat("frmtmb add_criterion(c('loo', 'loo_R2'), seed = 1):",
    tc(names(add_criterion(dg, c("loo", "loo_R2"), seed = 1)$criteria)),
    "\n")
cat("frmtmb add_criterion(c('loo', 'loo_R2')):",
    tc(names(add_criterion(dg, c("loo", "loo_R2"))$criteria)), "\n")

cat("\n== brms model_weights with a stored loo\n")
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
    format(w_fresh, digits = 8), "| identical",
    identical(unname(w_st), unname(w_fresh)), "\n")
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
  sprintf(paste("draws %d, nchains %d, chain_id %s, psis r_eff equals",
                "two-chain relative_eff %s"),
          nrow(d$draws), frmtmb.sample:::draws_nchains(d),
          if (is.null(cid)) "NULL" else paste(table(cid), collapse = "/"),
          isTRUE(all.equal(attr(ps, "r_eff"), rh)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(st, "\n")
cat("done\n")
