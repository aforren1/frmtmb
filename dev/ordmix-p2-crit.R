# Punch round 2, B3: candidate statistics for the degenerate-component
# warning, measured on three populations of fits:
#   margin_<link>_<c>_<n>: the review's strong-predictor designs
#     (dev/ordmix-rev2-margin.R, seeds 1..20); sound = every standard
#     error finite and both slopes within 50% of the truth;
#   rev_<...>, id_<...>: the 160 + 140 fits of dev/ordmix-p1-degen.R,
#     labeled by the review's rule (|estimate| > 30 or a non-finite
#     standard error).
# Per fit and component: reach (the round-1 statistic), collapsed, and
# the fraction of rows the component predicts as a point mass, the
# largest category probability above 1 - eps, for eps = 1e-6, 1e-9 and
# 1e-12. A step-function component predicts nearly every row with
# certainty; a sound component with a strong predictor does so only on
# the rows far from every threshold, and a heavy-tailed link on almost
# none. Usage: Rscript dev/ordmix-p2-crit.R <config>
# Output: one REP line per fit on stdout (the driver redirects to
# dev/ordmix-p2-crit-log/<config>.txt).
cfg <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cut4 <- function(lat) 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
# the category probabilities of component k, n x ncat
comp_probs <- function(fit, resp) {
  rspec <- fit$spec$responses[[resp]]
  fam <- rspec$family
  mx <- fam[["mix"]][["ord"]]
  dp <- ns$with_cs_offsets(fit, rspec, ns$eval_dpars(fit))[[resp]]
  av <- fit$frame[["aterm_values"]][[resp]]
  ex <- ns$fit_extras(fit, resp)
  n <- length(as.numeric(dp[["mu1"]]))
  codes <- ns$ord_code0(fam) + seq_len(ns$ordinal_ncat(fit, resp)) - 1L
  lapply(seq_len(mx$K), function(k) {
    vapply(codes, function(cd) {
      exp(fam[["mix"]][["comp_lpdf"]](rep.int(cd, n), dp, av, k, ex))
    }, numeric(n))
  })
}
stats_of <- function(fit) {
  resp <- names(fit$spec$responses)[1]
  dg <- ns$mixture_ord_degeneracy(fit, resp)
  P <- comp_probs(fit, resp)
  pm <- function(eps) {
    vapply(P, function(M) {
      mx <- apply(M, 1, max)
      mean(mx > 1 - eps, na.rm = TRUE)
    }, 0)
  }
  sh <- sharpen(fit, resp)
  list(reach = dg$reach, collapsed = dg$collapsed, pm6 = pm(1e-6),
       pm9 = pm(1e-9), pm12 = pm(1e-12), dll2 = sh$dll2,
       ll1err = sh$ll1err, near = sh$near)
}
# The mixture's log-likelihood with component k's discrimination scaled
# by s, which scales every latent distance of that component: dll2 is
# ll(2 disc_k) - ll(disc_k). A component at a step-function boundary
# loses nothing by sharpening (dll2 >= 0); a sound one loses. near: for
# each identified threshold, the latent distance to the nearest row,
# the largest over thresholds: a threshold the data place has rows near
# it. ll1err checks the evaluation against logLik().
sharpen <- function(fit, resp) {
  rspec <- fit$spec$responses[[resp]]
  fam <- rspec$family
  mx <- fam[["mix"]][["ord"]]
  dp <- ns$with_cs_offsets(fit, rspec, ns$eval_dpars(fit))[[resp]]
  av <- fit$frame[["aterm_values"]][[resp]]
  ex <- ns$fit_extras(fit, resp)
  y <- fit$frame[["y"]][[resp]] %||% fit$frame[["y"]]
  n <- length(as.numeric(dp[["mu1"]]))
  ll_of <- function(d) {
    lpi <- fam[["mix"]][["log_pi"]](d)
    L <- vapply(seq_len(mx$K), function(k) {
      rep(as.numeric(lpi[[k]]), length.out = n) +
        fam[["mix"]][["comp_lpdf"]](y, d, av, k, ex)
    }, numeric(n))
    m <- apply(L, 1, max)
    sum(m + log(rowSums(exp(L - m))))
  }
  ll1 <- ll_of(dp)
  dll2 <- vapply(seq_len(mx$K), function(k) {
    d2 <- dp
    nm <- paste0("disc", k)
    d2[[nm]] <- 2 * (if (is.null(dp[[nm]])) 1 else dp[[nm]])
    ll_of(d2) - ll1
  }, 0)
  near <- vapply(seq_len(mx$K), function(k) {
    lp <- fit$frame[["linpreds"]][[ns$linpred_key(resp, paste0("mu", k))]]
    ob <- ns$ord_lp_block(fit$frame, fit$spec, lp)
    tau <- as.numeric(ns$ord_threshold_values(ob$fam,
                                              fit$estimates[[ob$comp]]))
    eta <- as.numeric(dp[[paste0("mu", k)]])
    disc <- rep(as.numeric(if (is.null(dp[[paste0("disc", k)]])) 1 else
      dp[[paste0("disc", k)]]), length.out = n)
    cs <- dp[[ns$cs_slot(paste0("mu", k))]]
    th <- ob$fam[["thres"]]
    lay <- ns$thres_layout(if (is.null(th[["nthres"]])) length(tau) else
      th[["nthres"]])
    gi <- if (isTRUE(th[["grouped"]])) as.integer(av[["thres_gr"]]) else
      rep(1L, n)
    un <- th[["unident"]]
    out <- 0
    for (g in seq_len(lay$G)) {
      jg <- lay$start[g]:lay$end[g]
      r <- which(gi == g)
      if (!length(r)) next
      for (jj in seq_along(jg)) {
        if (jg[jj] %in% un) next
        dd <- disc[r] * (tau[jg[jj]] - eta[r] -
                           if (is.null(cs)) 0 else cs[r, jj])
        out <- max(out, min(abs(dd), na.rm = TRUE))
      }
    }
    out
  }, 0)
  list(dll2 = dll2, ll1err = ll1 - as.numeric(logLik(fit)), near = near)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
fmt <- function(v) paste(signif(v, 4), collapse = ",")
one <- function(s, fit_fun, label_fun) {
  w <- character(0)
  fit <- tryCatch(withCallingHandlers(fit_fun(), warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("REP cfg=%s seed=%d ERROR %s\n", cfg, s,
                substr(conditionMessage(fit), 1, 120)))
    return(invisible())
  }
  st <- stats_of(fit)
  fx <- suppressWarnings(fixef(fit))
  cat(sprintf(paste0("REP cfg=%s seed=%d label=%s se_finite=%s ",
                     "maxabs=%.4g maxse=%.4g reach=%s collapsed=%s pm6=%s pm9=%s ",
                     "pm12=%s dll2=%s near=%s ll1err=%.3g degwarn=%s\n"),
              cfg, s, label_fun(fit, fx), all(is.finite(fx[, "Est.Error"])),
              max(abs(fx[, "Estimate"])), max(fx[, "Est.Error"]), fmt(st$reach),
              paste(st$collapsed, collapse = ","), fmt(st$pm6),
              fmt(st$pm9), fmt(st$pm12), fmt(st$dll2), fmt(st$near),
              st$ll1err, any(grepl("degenerate boundary", w, fixed = TRUE))))
}
rev_label <- function(fit, fx) {
  if (any(abs(fx[, "Estimate"]) > 30) || any(!is.finite(fx[, "Est.Error"])))
    "degenerate" else "sound"
}
if (startsWith(cfg, "margin_")) {
  p <- strsplit(sub("^margin_", "", cfg), "_")[[1]]
  link <- p[1]
  cc <- as.numeric(p[2])
  n <- as.integer(p[3])
  noise <- switch(link,
    logit = function(n) rlogis(n),
    probit = function(n) rnorm(n),
    cloglog = function(n) log(-log(1 - runif(n))),
    cauchit = function(n) rcauchy(n))
  truth <- c(1.5, -0.8) * cc
  cut_at <- c(-1.5, 0, 1.5) * max(1, cc / 2)
  for (s in 1:20) {
    set.seed(s)
    x <- rnorm(n)
    cls <- rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, truth[1] * x + 1, truth[2] * x - 1) + noise(n)
    d <- data.frame(x = x, y = 1L + (lat > cut_at[1]) + (lat > cut_at[2]) +
                      (lat > cut_at[3]))
    one(s, function() frm(bf(y ~ x), family = mixture(cumulative(link),
                                                      cumulative(link)),
                          data = d),
        function(fit, fx) {
          b <- fx[c("mu1_x", "mu2_x"), "Estimate"]
          rel <- min(max(abs(b - truth) / abs(truth)),
                     max(abs(rev(b) - truth) / abs(truth)))
          if (all(is.finite(fx[, "Est.Error"])) && rel <= 0.5) "sound" else
            if (any(abs(fx[, "Estimate"]) > 30 * cc) ||
                any(!is.finite(fx[, "Est.Error"]))) "degenerate" else "poor"
        })
  }
} else if (startsWith(cfg, "rev_")) {
  p <- strsplit(sub("^rev_", "", cfg), "_")[[1]]
  n <- as.integer(p[length(p)])
  fam <- if (identical(p[2], "cum")) mixture(cumulative(), cumulative()) else
    mixture(cumulative(), sratio())
  for (s in 1:40) {
    set.seed(s)
    x <- rnorm(n)
    cls <- rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
    d <- data.frame(x = x, y = cut4(lat))
    one(s, function() frm(bf(y ~ x), family = fam, data = d), rev_label)
  }
} else if (startsWith(cfg, "cs_")) {
  # cs() designs, the two of the suites where round 2's first rule was
  # silent (dev/ordmix-p2-log-suitefires.txt): y ~ cs(x) with
  # mixture(sratio(), acat()) on the gated brms test's process (n =
  # 300), and mu2 ~ cs(x) with mixture(cumulative(), sratio()) on
  # test-ordinal-mixture.R's (n = 400); seeds 1..20
  for (s in 1:20) {
    set.seed(s)
    n <- if (cfg == "cs_sr_acat") 300 else 400
    x <- rnorm(n)
    z <- rnorm(n)
    sample(c("a", "b"), n, TRUE)
    cls <- if (cfg == "cs_sr_acat") rbinom(n, 1, 0.4) else
      rbinom(n, 1, plogis(-0.4 + 0.5 * z))
    lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) +
      (if (cfg == "cs_sr_acat") rlogis(n) else rlogis(n) / exp(0.3 * z * cls))
    d <- data.frame(x = x, y = cut4(lat))
    if (cfg == "cs_sr_acat") {
      one(s, function() frm(bf(y ~ cs(x)), family = mixture(sratio(), acat()),
                            data = d), rev_label)
    } else {
      one(s, function() frm(bf(y ~ x, mu2 ~ cs(x)),
                            family = mixture(cumulative(), sratio()),
                            data = d), rev_label)
    }
  }
} else {
  gen_id <- function(seed, n = 500) {
    set.seed(seed)
    x <- rnorm(n)
    z <- rnorm(n)
    cls <- rbinom(n, 1, 0.4)
    lat_none <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
    lat_mu <- ifelse(cls == 1, 2 * x, -1.5 * x) + rlogis(n)
    lat_pr <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rnorm(n)
    y <- cut4(lat_none)
    hu <- plogis(ifelse(cls == 1, -2, -0.4) + 0.8 * z)
    yh <- ifelse(runif(n) < hu, 0L, y)
    data.frame(x, z, y, yh, ymu = cut4(lat_mu), ypr = cut4(lat_pr),
               gr = factor(sample(letters[1:2], n, TRUE)))
  }
  sp <- switch(sub("^id_", "", cfg),
    cum2 = list(f = bf(y ~ x), fam = mixture(cumulative(), cumulative())),
    probit2 = list(f = bf(ypr ~ x),
                   fam = mixture(cumulative("probit"), cumulative("probit"))),
    cum_sr_clog = list(f = bf(y ~ x),
                       fam = mixture(cumulative(), sratio("cloglog"))),
    mu2 = list(f = bf(ymu ~ x),
               fam = mixture(cumulative(), cumulative(), order = "mu")),
    hurdle_huz = list(f = bf(yh ~ x, hu1 ~ z, hu2 ~ z),
                      fam = mixture(hurdle_cumulative(), hurdle_cumulative())),
    sr_acat_th = list(f = bf(y ~ x, theta1 ~ z),
                      fam = mixture(sratio(), acat())),
    gr_cum2 = list(f = bf(y | thres(gr = gr) ~ x),
                   fam = mixture(cumulative(), cumulative())))
  for (s in 1:20) {
    d <- gen_id(s)
    one(s, function() frm(sp$f, family = sp$fam, data = d), rev_label)
  }
}
