# Lane sampfix, script 03: every exported draws method on (a) draws
# from frm_sample(laplace = TRUE) and (b) ordinary draws with the
# stanfit removed (stanfit = NULL, how a test supplies draws without a
# sampler). Records, per call, ERROR with the message, or the result's
# class and its count of non-finite numbers.
#
#   Rscript dev/sampfix-03-sweep.R <lane|ref>

arm <- commandArgs(trailingOnly = TRUE)[1L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("ARM ", arm, "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a)) b else a

nonfinite <- function(r) {
  if (is.list(r) && !is.data.frame(r)) {
    return(sum(vapply(r, nonfinite, numeric(1))))
  }
  if (is.data.frame(r)) r <- unlist(r[vapply(r, is.numeric, NA)])
  if (!is.numeric(r)) return(0)
  sum(!is.finite(as.numeric(r)))
}

run <- function(lab, expr) {
  w <- character(0)
  pdf(NULL)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) invokeRestart("muffleMessage"))
  if (inherits(r, "ggplot")) try(print(r), silent = TRUE)
  dev.off()
  out <- if (inherits(r, "error")) {
    paste0("ERROR: ", gsub("[\r\n]+", " ", substr(conditionMessage(r), 1, 110)))
  } else {
    nf <- tryCatch(nonfinite(r), error = function(e) NA)
    paste0("OK ", class(r)[1L], if (!is.na(nf) && nf > 0) {
      paste0(" NONFINITE=", nf)
    })
  }
  if (length(w)) out <- paste0(out, "  [", length(w), " warn: ",
                               substr(unique(w)[1L], 1, 50), "]")
  cat(sprintf("  %-34s %s\n", lab, out))
}

set.seed(1212L)
dd <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(6L, 0, 0.5)[as.integer(dd$g)] +
  rnorm(nrow(dd), 0, 0.7)
nd <- data.frame(x = c(-1, 1), g = factor(c(1, 2), levels = 1:6))

q <- function(e) e
calls <- function(ds) list(
  print = quote(capture.output(print(ds))),
  summary = quote(summary(ds)),
  fixef = quote(fixef(ds)),
  ranef = quote(ranef(ds)),
  coef = quote(coef(ds)),
  VarCorr = quote(VarCorr(ds)),
  "print(VarCorr)" = quote(capture.output(print(VarCorr(ds)))),
  prior_summary = quote(capture.output(prior_summary(ds))),
  "hypothesis b" = quote(hypothesis(ds, "x > 0")),
  "hypothesis sd" = quote(hypothesis(ds, "sd_g__Intercept > 0",
                                     class = NULL)),
  "hypothesis ranef" = quote(hypothesis(ds, "Intercept > 0",
                                        scope = "ranef", group = "g")),
  posterior_epred = quote(posterior_epred(ds, ndraws = 20)),
  "posterior_epred NA" = quote(posterior_epred(ds, ndraws = 20,
                                               re_formula = NA)),
  "posterior_epred newdata" = quote(posterior_epred(ds, newdata = nd,
                                                    ndraws = 20)),
  "posterior_epred newdata NA" = quote(posterior_epred(
    ds, newdata = nd, ndraws = 20, re_formula = NA)),
  posterior_linpred = quote(posterior_linpred(ds, ndraws = 20)),
  "posterior_linpred NA" = quote(posterior_linpred(ds, ndraws = 20,
                                                   re_formula = NA)),
  posterior_predict = quote(posterior_predict(ds, ndraws = 20)),
  "posterior_predict NA" = quote(posterior_predict(ds, ndraws = 20,
                                                   re_formula = NA)),
  predictive_error = quote(predictive_error(ds, ndraws = 20)),
  predictive_interval = quote(predictive_interval(ds, ndraws = 20)),
  posterior_interval = quote(posterior_interval(ds)),
  posterior_summary = quote(posterior_summary(ds)),
  fitted = quote(fitted(ds, ndraws = 20)),
  "fitted NA" = quote(fitted(ds, ndraws = 20, re_formula = NA)),
  predict = quote(predict(ds, ndraws = 20)),
  residuals = quote(residuals(ds, ndraws = 20)),
  log_lik = quote(log_lik(ds)),
  loo = quote(loo(ds)),
  waic = quote(waic(ds)),
  psis = quote(psis(ds)),
  bayes_R2 = quote(bayes_R2(ds)),
  "pp_check dens" = quote(pp_check(ds, ndraws = 5)),
  "pp_check loo_pit_overlay" = quote(pp_check(ds, type = "loo_pit_overlay")),
  "pp_check loo_intervals" = quote(pp_check(ds, type = "loo_intervals")),
  pp_mixture = quote(pp_mixture(ds)),
  conditional_effects = quote(conditional_effects(ds)),
  mcmc_plot = quote(mcmc_plot(ds)),
  pairs = quote(pairs(ds)),
  plot = quote(plot(ds)),
  rhat = quote(rhat(ds)),
  neff_ratio = quote(neff_ratio(ds)),
  nuts_params = quote(nuts_params(ds)),
  log_posterior = quote(log_posterior(ds)),
  as_draws = quote(as_draws(ds)),
  as_draws_matrix = quote(as_draws_matrix(ds)),
  as_draws_array = quote(as_draws_array(ds)),
  as_draws_df = quote(as_draws_df(ds)),
  as_draws_list = quote(as_draws_list(ds)),
  as_draws_rvars = quote(as_draws_rvars(ds)),
  as.matrix = quote(as.matrix(ds)),
  as.array = quote(as.array(ds)),
  as.data.frame = quote(as.data.frame(ds)),
  as.mcmc = quote(as.mcmc(ds)),
  ndraws = quote(ndraws(ds)),
  nchains = quote(nchains(ds)),
  niterations = quote(niterations(ds)),
  nvariables = quote(nvariables(ds)),
  variables = quote(variables(ds)),
  nobs = quote(nobs(ds)),
  formula = quote(formula(ds)),
  family = quote(family(ds)),
  getCall = quote(getCall(ds)),
  ngrps = quote(ngrps(ds)),
  parnames = quote(parnames(ds)),
  posterior_samples = quote(posterior_samples(ds)),
  nsamples = quote(nsamples(ds)),
  stancode = quote(stancode(ds)),
  standata = quote(standata(ds)),
  restructure = quote(restructure(ds)),
  update = quote(update(ds)),
  loo_compare = quote(loo_compare(ds, ds)),
  kfold = quote(kfold(ds)),
  reloo = quote(reloo(ds)),
  bridge_sampler = quote(bridge_sampler(ds))
)

fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
lap <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 3,
             laplace = TRUE)))
full <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 3)))
nul <- full
nul$stanfit <- NULL

for (nm in c("lap", "nul")) {
  ds <- get(nm)
  cat("\n==== ", nm, ": ", paste(colnames(ds$draws)[1:5], collapse = " "),
      " ... (", ncol(ds$draws), " columns)\n", sep = "")
  cl <- calls(ds)
  for (k in names(cl)) {
    set.seed(9)
    run(k, eval(cl[[k]]))
  }
}
cat("\nDONE\n")
