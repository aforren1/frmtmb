# Lane loodraws, fix round 1 (copied from the reviewer's dev/loodraws-rev-same.R,
# library switched to the lane's): brms 2.23.0 and frmtmb.sample on the SAME
# posterior draws. Each brms fit's draws (as.matrix(), chain-major) are put
# into a frmtmb_draws over the ML fit of the same model, with the brms
# stanfit as its `stanfit`, so the chain structure is brms's too. Then
# every quantity that is a function of the draws must agree to rounding,
# and loo_R2() at one seed must agree draw by draw. No Monte Carlo.
#   Rscript dev/loodraws-rev-same.R
# Reads dev/loodraws-rev-out/brmsfits.rds (dev/loodraws-rev-brmsfit.R);
# log dev/loodraws-log/same.txt.
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
options(mc.cores = 1)
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/loodraws-rev-out"
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "brms",
    format(packageVersion("brms")), "loo", format(packageVersion("loo")),
    "\n")
cat("loo_epred resolves to",
    environmentName(environment(get("loo_epred"))), "\n")
B <- readRDS(file.path(out_dir, "brmsfits.rds"))
dd <- B$dd
q <- function(e) suppressWarnings(suppressMessages(e))
# relative difference, by the largest magnitude of the reference
rd <- function(a, b) {
  if (!identical(dim(a), dim(b))) return(paste("DIM", paste(dim(a),
    collapse = "x"), "vs", paste(dim(b), collapse = "x")))
  m <- max(abs(b), na.rm = TRUE)
  signif(max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE) /
           if (m > 0) m else 1, 3)
}
same_dn <- function(a, b) identical(dimnames(a), dimnames(b))
frm_of <- list(
  gauss = function() frm(bf(y ~ x), family = gaussian(), data = dd),
  wt = function() frm(bf(y | weights(w) ~ x), family = gaussian(), data = dd),
  bern = function() frm(bf(y01 ~ x), family = bernoulli(), data = dd),
  binom = function() frm(bf(s | trials(nt) ~ x), family = binomial(),
                         data = dd),
  ri = function() frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd),
  mv = function() frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
                      data = dd),
  thin = function() frm(bf(y ~ x), family = gaussian(), data = dd))
mk_ds <- function(bfit, fit) {
  labs <- frmtmb::brms_par_labels(fit)
  M <- as.matrix(bfit)
  # the group SD: the conditional log_lik and epred never read it, so
  # brms's sd stands in for frmtmb's own column
  if ("theta_1" %in% labs && "sd_g__Intercept" %in% colnames(M)) {
    M <- cbind(M, theta_1 = M[, "sd_g__Intercept"])
  }
  miss <- setdiff(labs, colnames(M))
  if (length(miss)) stop("labels missing from brms: ", paste(miss, collapse = " "))
  structure(list(stanfit = bfit$fit,
                 draws = cbind(M[, labs, drop = FALSE], lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
DS <- list()
for (nm in names(B$fits)) {
  bf_ <- B$fits[[nm]]
  fit <- q(frm_of[[nm]]())
  ds <- mk_ds(bf_, fit)
  DS[[nm]] <- ds
  cat("\n==== ", nm, ": ", nrow(ds$draws), " draws, chains ",
      frmtmb.sample:::draws_nchains(ds), "\n", sep = "")
  mv <- nm == "mv"
  rr <- if (mv) list(NULL, "y", "y2") else list(NULL)
  for (r in rr) {
    tag <- if (is.null(r)) "" else paste0(" resp=", r)
    llb <- log_lik(bf_, resp = r)
    llf <- log_lik(ds, resp = r)
    cat("log_lik", tag, " rel diff ", rd(llf, llb), "\n", sep = "")
    if (!mv || !is.null(r)) {
      cat("posterior_epred", tag, " rel diff ",
          rd(posterior_epred(ds, resp = r), posterior_epred(bf_, resp = r)),
          "\n", sep = "")
    }
    for (ty in c("mean", "var", "quantile")) {
      pr <- if (ty == "quantile") c(0.025, 0.5, 0.975) else 0.5
      eb <- q(loo_epred(bf_, type = ty, probs = pr, resp = r))
      ef <- q(loo_epred(ds, type = ty, probs = pr, resp = r))
      cat(sprintf("loo_epred%s %-8s rel diff %-10s dimnames same %s dim %s\n",
                  tag, ty, rd(ef, eb), same_dn(ef, eb),
                  paste(dim(ef), collapse = "x")))
      if (!same_dn(ef, eb)) {
        cat("  frmtmb dimnames:"); str(dimnames(ef))
        cat("  brms dimnames:  "); str(dimnames(eb))
      }
    }
    lb <- q(loo_linpred(bf_, type = "var", resp = r))
    lf <- q(loo_linpred(ds, type = "var", resp = r))
    cat("loo_linpred", tag, " var rel diff ", rd(lf, lb), " dimnames same ",
        same_dn(lf, lb), "\n", sep = "")
    if (nm %in% c("bern", "binom")) {
      lb <- q(loo_linpred(bf_, transform = TRUE, resp = r))
      lf <- q(loo_linpred(ds, transform = TRUE, resp = r))
      cat("loo_linpred(transform = TRUE) rel diff ", rd(lf, lb), "\n", sep = "")
    }
    # psis_object = brms's own: the prediction is then the only input
    psb <- q(loo(bf_, resp = r, save_psis = TRUE))$psis_object
    eb <- q(loo_epred(bf_, psis_object = psb, resp = r))
    ef <- q(loo_epred(ds, psis_object = psb, resp = r))
    cat("loo_epred(psis_object = brms's)", tag, " identical ",
        identical(ef, eb), " rel diff ", rd(ef, eb), "\n", sep = "")
    psf <- q(loo(ds, resp = r, save_psis = TRUE))$psis_object
    cat("psis r_eff rel diff ", rd(attr(psf, "r_eff"), attr(psb, "r_eff")),
        "  pareto k rel diff ", rd(psf$diagnostics$pareto_k,
                                    psb$diagnostics$pareto_k), "\n", sep = "")
    # loo_predict: posterior_predict is random on each side; the shape,
    # the names and, for a discrete response, how many rows agree
    set.seed(11)
    pb <- q(loo_predict(bf_, type = "quantile", probs = c(0.25, 0.75),
                        resp = r))
    set.seed(11)
    pf <- q(loo_predict(ds, type = "quantile", probs = c(0.25, 0.75),
                        resp = r))
    cat("loo_predict quantile dim ", paste(dim(pf), collapse = "x"), " / ",
        paste(dim(pb), collapse = "x"), " dimnames same ", same_dn(pf, pb),
        " rows equal ", sum(apply(abs(pf - pb) < 1e-12, 1, all)), " of ",
        nrow(pf), " max diff / sd(y) ",
        signif(max(abs(pf - pb)) / stats::sd(as.numeric(
          posterior_predict(bf_, resp = r, ndraws = 200))), 3),
        "\n", sep = "")
    set.seed(12)
    ib <- q(loo_predictive_interval(bf_, prob = 0.5, resp = r))
    set.seed(12)
    iff <- q(loo_predictive_interval(ds, prob = 0.5, resp = r))
    cat("loo_predictive_interval(0.5) colnames ",
        paste(colnames(iff), collapse = " "), " / ",
        paste(colnames(ib), collapse = " "), "\n", sep = "")
    # loo_R2: same draws, same seed, brms's order of rexp()
    if (!(mv && is.null(r))) {
      Rb <- q(loo_R2(bf_, resp = r, seed = 7, summary = FALSE))
      Rf <- q(loo_R2(ds, resp = r, seed = 7, summary = FALSE))
      cat("loo_R2 draws rel diff ", rd(Rf, Rb), " colnames ",
          paste(colnames(Rf), collapse = " "), " / ",
          paste(colnames(Rb), collapse = " "), "\n", sep = "")
      Sb <- q(loo_R2(bf_, resp = r, seed = 7, robust = TRUE,
                     probs = c(0.1, 0.9)))
      Sf <- q(loo_R2(ds, resp = r, seed = 7, robust = TRUE,
                     probs = c(0.1, 0.9)))
      cat("loo_R2 robust summary rel diff ", rd(Sf, Sb), " dimnames same ",
          same_dn(Sf, Sb), "\n", sep = "")
    } else {
      Rb <- q(loo_R2(bf_, seed = 7, summary = FALSE))
      Rf <- q(loo_R2(ds, seed = 7, summary = FALSE))
      cat("loo_R2 (mv, all resp) draws rel diff ", rd(Rf, Rb), " colnames ",
          paste(colnames(Rf), collapse = " "), " / ",
          paste(colnames(Rb), collapse = " "), "\n", sep = "")
    }
  }
}

cat("\n==== loo_model_weights (gauss vs ri) on the same draws\n")
g <- B$fits$gauss; ri <- B$fits$ri
for (cfg in list(list(method = "stacking"),
                 list(method = "pseudobma", BB = FALSE),
                 list(method = "pseudobma", BB = TRUE, BB_n = 500))) {
  set.seed(5)
  wb <- q(do.call(loo_model_weights, c(list(g, ri), cfg)))
  set.seed(5)
  wf <- q(do.call(loo_model_weights, c(list(DS$gauss, DS$ri), cfg)))
  cat(sprintf("%-28s brms %s | frmtmb %s | rel diff %s | class %s / %s\n",
              paste(unlist(cfg), collapse = " "),
              paste(format(as.numeric(wb), digits = 8), collapse = " "),
              paste(format(as.numeric(wf), digits = 8), collapse = " "),
              rd(as.numeric(wf), as.numeric(wb)), class(wf)[1], class(wb)[1]))
}

cat("\n==== r_eff on a subset of draws: brms vs frmtmb, same draw_ids\n")
set.seed(13)
ids <- sort(sample(4000, 1000))
for (nm in c("gauss", "ri")) {
  bf_ <- B$fits[[nm]]; ds <- DS[[nm]]
  pb <- q(loo(bf_, draw_ids = ids, save_psis = TRUE))
  pf <- q(loo(ds, draw_ids = ids, save_psis = TRUE))
  llb <- log_lik(bf_, draw_ids = ids)
  one_chain <- loo::relative_eff(exp(llb), chain_id = rep(1L, length(ids)))
  cat(nm, ": brms r_eff = relative_eff(one chain): ",
      isTRUE(all.equal(attr(pb$psis_object, "r_eff"), one_chain)),
      "; range ", paste(signif(range(attr(pb$psis_object, "r_eff")), 3),
                        collapse = " "),
      "; frmtmb r_eff range ",
      paste(signif(range(attr(pf$psis_object, "r_eff")), 3), collapse = " "),
      "\n", sep = "")
  cat("  elpd_loo brms ", format(pb$estimates["elpd_loo", 1], digits = 10),
      " frmtmb ", format(pf$estimates["elpd_loo", 1], digits = 10),
      " diff ", signif(pf$estimates["elpd_loo", 1] -
                         pb$estimates["elpd_loo", 1], 3),
      " (brms se_elpd ", signif(pb$estimates["elpd_loo", 2], 3), ")\n",
      sep = "")
  eb <- q(loo_epred(bf_, draw_ids = ids))
  ef <- q(loo_epred(ds, draw_ids = ids))
  cat("  loo_epred(draw_ids) rel diff ", rd(ef, eb), "\n", sep = "")
  # pp_check: brms runs loo(draw_ids =) and hands the log weights on
  LW <- NULL
  suppressMessages(trace(bayesplot::ppc_loo_pit_overlay, print = FALSE,
                         tracer = quote(LW <<- lw)))
  invisible(q(pp_check(bf_, type = "loo_pit_overlay", draw_ids = ids)))
  lwb <- LW
  LW <- NULL
  invisible(q(pp_check(ds, type = "loo_pit_overlay", draw_ids = ids)))
  lwf <- LW
  LW <- NULL
  invisible(q(pp_check(bf_, type = "loo_pit_overlay")))
  lwb_all <- LW
  LW <- NULL
  invisible(q(pp_check(ds, type = "loo_pit_overlay")))
  lwf_all <- LW
  suppressMessages(untrace(bayesplot::ppc_loo_pit_overlay))
  lw_loo <- stats::weights(pf$psis_object, log = TRUE)
  cat("  pp_check lw, draw_ids: frmtmb vs brms max abs diff ",
      signif(max(abs(lwf - lwb)), 3), "; frmtmb pp_check vs frmtmb loo() ",
      signif(max(abs(lwf - lw_loo)), 3), "; brms pp_check vs brms loo() ",
      signif(max(abs(lwb - stats::weights(pb$psis_object, log = TRUE))), 3),
      "\n", sep = "")
  cat("  pp_check lw, all draws: frmtmb vs brms max abs diff ",
      signif(max(abs(lwf_all - lwb_all)), 3), "\n", sep = "")
}

cat("\n==== ndraws: frmtmb's evenly spaced rows, brms's random ones\n")
ds <- DS$gauss
rows <- frmtmb.sample:::draws_subsample(ds, 1000, NULL)
ll <- log_lik(ds, draw_ids = rows)
cat("chain_id attr on ndraws = 1000:", is.null(attr(ll, "chain_id")),
    "(TRUE = none, r_eff 1)\n")
r1 <- loo::relative_eff(exp(ll), chain_id = rep(1L, 1000))
cat("one-chain relative_eff on frmtmb's evenly spaced rows: range",
    signif(range(r1), 3), "\n")
ef1 <- q(loo_epred(ds, ndraws = 1000))
ps1 <- q(loo::psis(-ll, r_eff = r1))
ref1 <- loo::E_loo(posterior_epred(ds, draw_ids = rows), ps1)$value
cat("loo_epred(ndraws = 1000) r_eff 1 vs one-chain r_eff: max rel diff",
    rd(ef1[, 1], ref1), "\n")

cat("\n==== brms argument behaviors\n")
g <- B$fits$gauss
tc <- function(e) tryCatch({ q(e); "ok" }, error = function(err)
  paste("ERROR:", conditionMessage(err)))
cat("brms loo_epred(save_psis = TRUE):", tc(loo_epred(g, save_psis = TRUE)),
    "\n")
cat("brms loo_linpred(scale = 'linear') identical to default:",
    identical(q(loo_linpred(g, scale = "linear")), q(loo_linpred(g))), "\n")
cat("brms loo_epred(junk = 1):", tc(loo_epred(g, junk = 1)), "\n")
cat("frmtmb loo_epred(junk = 1):", tc(loo_epred(DS$gauss, junk = 1)), "\n")
cat("brms loo_epred(probs = 0.3) [type mean] identical to default:",
    identical(q(loo_epred(g, probs = 0.3)), q(loo_epred(g))), "\n")
cat("brms loo_predictive_interval(prob = 1.5):",
    tc(loo_predictive_interval(g, prob = 1.5)), "\n")
cat("frmtmb loo_predictive_interval(prob = 1.5):",
    tc(loo_predictive_interval(DS$gauss, prob = 1.5)), "\n")
cat("brms loo_R2(ndraws = 500) seed 7 twice identical:",
    identical(q(loo_R2(g, ndraws = 500, seed = 7)),
              q(loo_R2(g, ndraws = 500, seed = 7))), "\n")
cat("brms loo_model_weights(g, ri, ndraws = 500):",
    tc(loo_model_weights(g, ri, ndraws = 500)), "\n")
cat("frmtmb loo_model_weights(ndraws = 500):",
    tc(loo_model_weights(DS$gauss, DS$ri, ndraws = 500)), "\n")
cat("brms loo_epred(resp = 'y') on a univariate model:",
    tc(loo_epred(g, resp = "y")), "\n")
cat("frmtmb loo_epred(resp = 'y') on a univariate model:",
    tc(loo_epred(DS$gauss, resp = "y")), "\n")
cat("brms loo_epred(ndraws = 500, moment_match = FALSE):",
    tc(loo_epred(g, ndraws = 500, moment_match = FALSE)), "\n")
cat("brms loo_epred(cores = 1):", tc(loo_epred(g, cores = 1)), "\n")
cat("frmtmb loo_epred(cores = 1):", tc(loo_epred(DS$gauss, cores = 1)), "\n")
cat("brms loo_R2(seed = 7) twice identical (stored criterion absent):",
    identical(q(loo_R2(g, seed = 7)), q(loo_R2(g, seed = 7))), "\n")
cat("done\n")
