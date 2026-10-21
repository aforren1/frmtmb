# The verdict moves of the 0.73.0 consolidation in bin 2 of the ported
# brms suite. dev/brmsport-bin2-verdicts-2t.R and -2p.R source this
# after the 0.71.0 and 0.72.0 moves (dev/rel073-resolve.R adds the
# calls). Each move is a row the merged build records differently from
# 0.72.0, with the reason and the lane whose change moved it.
#
# Returns function(X, kind): kind "2t" edits the 2T table, "2p" the 2P
# row verdicts, "values" the value-check verdicts.
function(X, kind) {
  drop <- function(X, ids, pkg = NULL) {
    sel <- X$id %in% ids & (is.null(pkg) | X$pkg %in% pkg)
    stopifnot(sum(sel) == length(ids))
    cat("bin-2", kind, "verdicts removed (now hold):", ids, "\n")
    X[!sel, , drop = FALSE]
  }
  set <- function(X, id, pkg, verdict, class, reason) {
    i <- which(X$id == id & X$pkg == pkg)
    stopifnot(length(i) == 1L)
    X$verdict[i] <- verdict
    X$class[i] <- class
    X$reason[i] <- reason
    cat("bin-2", kind, "verdict moved:", id, pkg, "to", verdict, "/",
        class, "\n")
    X
  }
  if (kind == "2t") {
    S <- "frmtmb.sample"
    F <- "frmtmb"
    m <- function(k) paste0("brmsfit-methods:", k)
    # lane loodraws: brms's loo_*() functions on draws. 13 of the 17
    # rows hold
    X <- drop(X, m(c(487, 490, 497, 511, 512, 519, 522, 533, 547, 548,
                     555, 562, 563)), S)
    for (k in c(508, 544)) {
      X <- set(X, m(k), S, "divergence", "policy", paste(
        "the loo_*() functions refuse newdata by name, as loo() does: a",
        "leave-one-out weight is the observation's own column of",
        "log_lik(), so it exists only for a row the model was fitted on.",
        "brms computes the weights on the new rows' own response (this row",
        "with allow_new_levels = TRUE). Lane loodraws's decision 1"))
    }
    X <- set(X, m(526), S, "divergence", "policy", paste(
      "loo_linpred(scale = ) is refused by name: brms 2.23.0 reads it",
      "nowhere, and its result is the default's, on a logit link too",
      "(lane loodraws, decision 3; its review re-ran brms)"))
    X <- set(X, m(523), S, "divergence", "inherited", paste(
      "loo_linpred() of an ordinal fit with a cs() term returns one",
      "matrix where brms returns one layer per threshold (dim(llp)[3] is",
      "3 in brms, NA here): posterior_linpred() of draws has no cs()",
      "threshold layer, and loo_linpred() inherits it",
      "(dev/test-backlog.md, filed at 0.73.0)"))
    # loo() of draws accepts brms's default newdata = NULL since lane
    # loodraws, and a frmtmb_draws has no $data, so the row's
    # newdata = fit6$data is NULL and the call is loo(fit6): it holds
    # hollowly (the review of lane loodraws, "Row 469 is hollow")
    X <- set(X, m(469), S, "divergence", "hollow", paste(
      "policy (newdata): holds hollowly. A frmtmb_draws has no $data",
      "(names stanfit, draws, fit), so fit6$data is NULL and the call is",
      "loo(fit6), which since 0.73.0 accepts brms's default newdata =",
      "NULL; with the data given (newdata = model.frame(fit6)) loo() of",
      "draws refuses newdata by name, as log_lik() of draws does. brms",
      "computes the elpd on newdata"))
    # on the fit (the restatement, not counted): core now has the
    # generics and refuses by name, so the reason's "absent" is stale
    for (k in c(487, 490, 497, 508, 511, 512, 519, 522, 523, 526, 533,
                544, 547, 548, 555, 562, 563)) {
      X <- set(X, m(k), F, "divergence", "no-draws", paste(
        "a leave-one-out quantity over posterior draws, which a",
        "maximum-likelihood fit has none of: since 0.73.0 frmtmb's own",
        "generic refuses it on a fit by name and names frm_sample()"))
    }
    for (k in 478:480) {
      X <- set(X, m(k), F, "divergence", "no-draws", paste(
        "loo_subsample() approximates loo() over posterior draws, which a",
        "maximum-likelihood fit has none of: since 0.73.0 frmtmb's own",
        "generic refuses it on a fit by name"))
    }
  } else if (kind == "2p") {
    # lane nb2dm: dirichlet_multinomial() is a frmtmb family, and both
    # of its rows hold (log_lik:510's value check holds against brms's
    # log_lik_dirichlet_multinomial(), 2.13e-14 on 400 values)
    X <- drop(X, c("log_lik:510", "posterior_predict:355"))
    # lane cdfs: the 15 p-function rows of families that gained a CDF,
    # and the censored von_mises row, hold, each with a value check
    # against brms's own function
    X <- drop(X, c(paste0("distributions:",
                          c(6, 52, 64, 83, 125, 150, 164, 169, 174, 179,
                            185, 195, 200, 205, 210)),
                   "log_lik:392"))
  }
  X
}
