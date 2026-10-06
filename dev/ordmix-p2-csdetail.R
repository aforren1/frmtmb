# Punch round 2, B3: the cs() fits of dev/ordmix-p2-crit-log/cs_*.txt
# where the warning and the review's label disagree: per component, the
# whole-component and per-threshold sharpening, and the estimates with
# their standard errors.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cut4 <- function(lat) 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
gen <- function(cfg, s) {
  set.seed(s)
  n <- if (cfg == "cs_sr_acat") 300 else 400
  x <- rnorm(n)
  z <- rnorm(n)
  sample(c("a", "b"), n, TRUE)
  cls <- if (cfg == "cs_sr_acat") rbinom(n, 1, 0.4) else
    rbinom(n, 1, plogis(-0.4 + 0.5 * z))
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) +
    (if (cfg == "cs_sr_acat") rlogis(n) else rlogis(n) / exp(0.3 * z * cls))
  data.frame(x = x, y = cut4(lat))
}
fitof <- function(cfg, d) {
  if (cfg == "cs_sr_acat") {
    frm(bf(y ~ cs(x)), family = mixture(sratio(), acat()), data = d)
  } else {
    frm(bf(y ~ x, mu2 ~ cs(x)), family = mixture(cumulative(), sratio()),
        data = d)
  }
}
# per-threshold probes, as mixture_ord_degeneracy() makes them
probes <- function(fit) {
  resp <- names(fit$spec$responses)[1]
  rspec <- fit$spec$responses[[resp]]
  fam <- rspec$family
  mx <- fam[["mix"]][["ord"]]
  dp <- ns$with_cs_offsets(fit, rspec, ns$eval_dpars(fit))[[resp]]
  av <- fit$frame[["aterm_values"]][[resp]]
  ex <- ns$fit_extras(fit, resp)
  y <- fit$frame[["y"]][[resp]]
  n <- length(y)
  ll_at <- function(d) {
    lp <- fam[["mix"]][["log_pi"]](d)
    L <- vapply(seq_len(mx$K), function(k) {
      rep(as.numeric(lp[[k]]), length.out = n) +
        fam[["mix"]][["comp_lpdf"]](y, d, av, k, ex)
    }, numeric(n))
    m <- apply(L, 1, max)
    sum(m + log(rowSums(exp(L - m))))
  }
  ll0 <- ll_at(dp)
  for (k in seq_len(mx$K)) {
    cs <- dp[[ns$cs_slot(paste0("mu", k))]]
    if (is.null(cs)) next
    lp <- fit$frame[["linpreds"]][[ns$linpred_key(resp, paste0("mu", k))]]
    ob <- ns$ord_lp_block(fit$frame, fit$spec, lp)
    tau <- as.numeric(ns$ord_threshold_values(ob$fam,
                                              fit$estimates[[ob$comp]]))
    eta <- as.numeric(dp[[paste0("mu", k)]])
    for (j in seq_len(ncol(cs))) {
      d3 <- dp
      cs3 <- cs
      cs3[, j] <- 2 * cs[, j] + eta - tau[j]
      d3[[ns$cs_slot(paste0("mu", k))]] <- cs3
      # the identity probe, factor 1, must give 0
      d1 <- dp
      cs1 <- cs
      cs1[, j] <- 1 * cs[, j] + 0 * (eta - tau[j])
      d1[[ns$cs_slot(paste0("mu", k))]] <- cs1
      # rows whose category the step decides: category >= j in sratio
      cat(sprintf("  component %d threshold %d: sharpen %.4g (identity %.3g); rows at or past step %d: %d\n",
                  k, j, ll_at(d3) - ll0, ll_at(d1) - ll0, j, sum(y >= j)))
    }
  }
}
for (cs in list(c("cs_mu2", 3), c("cs_sr_acat", 12), c("cs_sr_acat", 1),
                c("cs_sr_acat", 9), c("cs_sr_acat", 11))) {
  cfg <- cs[1]
  s <- as.integer(cs[2])
  fit <- suppressWarnings(fitof(cfg, gen(cfg, s)))
  cat("==", cfg, "seed", s, "\n")
  print(frmtmb:::mixture_ord_degeneracy(fit, names(fit$spec$responses)[1]))
  probes(fit)
  fx <- suppressWarnings(fixef(fit))
  print(signif(fx[, c("Estimate", "Est.Error")], 3))
}
