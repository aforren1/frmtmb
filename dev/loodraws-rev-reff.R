# Reviewer of lane loodraws: r_eff on a subset of draws. brms 2.23.0 and
# frmtmb.sample on the SAME draws (dev/loodraws-rev-same.R's
# construction) and the same draw_ids, for loo(), loo_epred() and
# pp_check(type = "loo_pit_overlay"); and frmtmb's ndraws subsample under
# r_eff 1 and under brms's one-chain relative_eff.
#   Rscript dev/loodraws-rev-reff.R -> dev/loodraws-rev-out/reff.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
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
rd <- function(a, b) signif(max(abs(as.numeric(a) - as.numeric(b))) /
                              max(abs(as.numeric(b))), 3)
tc <- function(e) tryCatch({ q(e); "ok" }, error = function(err)
  paste("ERROR:", substr(conditionMessage(err), 1, 70)))
mk_ds <- function(bfit, fit) {
  labs <- frmtmb::brms_par_labels(fit)
  M <- as.matrix(bfit)
  if ("theta_1" %in% labs) M <- cbind(M, theta_1 = M[, "sd_g__Intercept"])
  structure(list(stanfit = bfit$fit,
                 draws = cbind(M[, labs, drop = FALSE], lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
DS <- list(
  gauss = mk_ds(B$fits$gauss, q(frm(bf(y ~ x), family = gaussian(),
                                    data = dd))),
  ri = mk_ds(B$fits$ri, q(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                              data = dd))))
set.seed(13)
for (k in c(1000L, 3000L)) {
  ids <- sort(sample(4000, k))
  cat("\n==== draw_ids: a sorted random", k, "of 4000 (seed 13)\n")
  for (nm in names(DS)) {
    bf_ <- B$fits[[nm]]
    ds <- DS[[nm]]
    llb <- log_lik(bf_, draw_ids = ids)
    one <- loo::relative_eff(exp(llb), chain_id = rep(1L, k))
    pb <- q(loo(bf_, draw_ids = ids, save_psis = TRUE))
    cat(nm, ": brms loo() r_eff is one-chain relative_eff: ",
        isTRUE(all.equal(attr(pb$psis_object, "r_eff"), one)),
        " (range ", paste(signif(range(one), 3), collapse = " "), ")\n",
        sep = "")
    cat("  frmtmb loo(draw_ids =): ", tc(loo(ds, draw_ids = ids)),
        "; psis(draw_ids =): ", tc(psis(ds, draw_ids = ids)), "\n", sep = "")
    # frmtmb's convention, r_eff 1, built by hand on the same matrix
    ll <- log_lik(ds, draw_ids = ids)
    cat("  frmtmb log_lik chain_id attr NULL: ", is.null(attr(ll, "chain_id")),
        "; log_lik equal to brms's: ", identical(unname(ll), unname(llb)),
        "\n", sep = "")
    p1 <- q(loo::loo(ll, r_eff = 1, save_psis = TRUE))
    cat("  elpd_loo brms (one-chain r_eff) ",
        format(pb$estimates["elpd_loo", 1], digits = 10),
        " | r_eff 1 ", format(p1$estimates["elpd_loo", 1], digits = 10),
        " | diff ", signif(p1$estimates["elpd_loo", 1] -
                             pb$estimates["elpd_loo", 1], 3),
        " | brms SE ", signif(pb$estimates["elpd_loo", 2], 3),
        " | max |pareto k diff| ",
        signif(max(abs(p1$diagnostics$pareto_k -
                         pb$diagnostics$pareto_k)), 3), "\n", sep = "")
    eb <- q(loo_epred(bf_, draw_ids = ids))
    ef <- q(loo_epred(ds, draw_ids = ids))
    cat("  loo_epred(draw_ids =) frmtmb vs brms rel diff ", rd(ef, eb),
        "\n", sep = "")
    LW <- NULL
    suppressMessages(trace(bayesplot::ppc_loo_pit_overlay, print = FALSE,
                           tracer = quote(LW <<- lw)))
    invisible(q(pp_check(bf_, type = "loo_pit_overlay", draw_ids = ids)))
    lwb <- LW
    LW <- NULL
    invisible(q(pp_check(ds, type = "loo_pit_overlay", draw_ids = ids)))
    lwf <- LW
    suppressMessages(untrace(bayesplot::ppc_loo_pit_overlay))
    cat("  pp_check(draw_ids =) lw: frmtmb vs brms max abs diff ",
        if (is.null(lwf)) "NA (no lw captured)" else
          signif(max(abs(lwf - lwb)), 3),
        "; frmtmb pp_check vs r_eff-1 weights ",
        if (is.null(lwf)) "NA" else
          signif(max(abs(lwf - stats::weights(p1$psis_object, log = TRUE))),
                 3), "\n", sep = "")
  }
}

cat("\n==== ndraws: frmtmb's evenly spaced rows\n")
ds <- DS$gauss
for (k in c(1000L, 2000L, 3000L)) {
  rows <- frmtmb.sample:::draws_subsample(ds, k, NULL)
  ll <- log_lik(ds, draw_ids = rows)
  r1 <- loo::relative_eff(exp(ll), chain_id = rep(1L, k))
  ep <- posterior_epred(ds, draw_ids = rows)
  a <- q(loo_epred(ds, ndraws = k))
  b <- loo::E_loo(ep, q(loo::psis(-ll, r_eff = r1)))$value
  l1 <- q(loo(ds, ndraws = k))
  lb <- q(loo::loo(ll, r_eff = r1))
  cat(sprintf(paste("ndraws %d: one-chain r_eff range %s; loo_epred r_eff 1",
                    "vs one-chain max rel diff %s; elpd r_eff 1 %.6f",
                    "one-chain %.6f\n"),
              k, paste(signif(range(r1), 3), collapse = " "), rd(a[, 1], b),
              l1$estimates["elpd_loo", 1], lb$estimates["elpd_loo", 1]))
}
cat("done\n")
