# brms 2.23.0 against frmtmb.sample on the same models and data: the
# loo_*() functions, judged against the Monte Carlo SE each side's
# PSIS weights give. Also records brms's shapes and two brms behaviors
# the findings cite (an ignored `scale`, and the draws ndraws picks).
#   Rscript dev/loodraws-parity.R brms|frmtmb|compare
# Writes dev/loodraws-out/parity-<side>.rds; log
# dev/loodraws-log/parity-<side>.txt. Seeds: data 1, samplers 1.
side <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/loodraws-out"
dir.create(out_dir, showWarnings = FALSE)
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 6)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.7)[dd$g])
dd$y2 <- rnorm(n, -1 + 0.3 * dd$x)
dd$k <- rpois(n, exp(0.2 + 0.4 * dd$x))
q <- function(e) suppressWarnings(suppressMessages(e))
# per-observation MCSE of a self-normalized importance-sampling mean,
# with loo's own relative efficiency (loo::E_loo gives no MCSE)
mcse_mean <- function(x, ps) {
  w <- stats::weights(ps, log = FALSE)
  mu <- colSums(w * x)
  sqrt(colSums(w^2 * sweep(x, 2, mu)^2) / attr(ps, "r_eff"))
}
models <- list(
  gaussian = list(f = y ~ x, fam = "gaussian", resp = NULL),
  poisson = list(f = k ~ x, fam = "poisson", resp = NULL),
  ranef = list(f = y ~ x + (1 | g), fam = "gaussian", resp = NULL),
  mvbf = list(f = NULL, fam = "gaussian", resp = "y2"))
if (side == "brms") {
  suppressMessages(library(brms))
  cat("brms", format(packageVersion("brms")), "\n")
  res <- list()
  for (nm in names(models)) {
    m <- models[[nm]]
    form <- if (nm == "mvbf") {
      bf(y ~ x) + bf(y2 ~ x) + set_rescor(FALSE)
    } else bf(m$f)
    t0 <- proc.time()[["elapsed"]]
    fit <- q(brm(form, data = dd, family = m$fam, chains = 4, iter = 2000,
                 seed = 1, refresh = 0, cores = 1))
    cat(nm, "fit s:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
    r <- m$resp
    ll <- log_lik(fit, resp = r)
    ps <- q(loo::psis(-ll, r_eff = loo::relative_eff(
      exp(ll), chain_id = rep(1:4, each = 1000))))
    ep <- posterior_epred(fit, resp = r)
    lp <- posterior_linpred(fit, resp = r)
    res[[nm]] <- list(
      fit = fit,
      loo_epred = q(loo_epred(fit, resp = r)),
      loo_linpred = q(loo_linpred(fit, resp = r)),
      loo_epred_q = q(loo_epred(fit, resp = r, type = "quantile",
                                probs = c(0.25, 0.75))),
      loo_pi = {set.seed(2); q(loo_predictive_interval(fit, resp = r))},
      loo_R2 = q(loo_R2(fit, resp = r, seed = 7)),
      loo_R2_draws = q(loo_R2(fit, resp = r, seed = 7, summary = FALSE)),
      mcse_epred = mcse_mean(ep, ps), mcse_linpred = mcse_mean(lp, ps),
      ll = ll, ep = ep)
    if (nm == "mvbf") {
      res[[nm]]$mv_epred_dim <- dim(q(loo_epred(fit)))
      res[[nm]]$mv_epred_dimnames <- dimnames(q(loo_epred(fit)))
      res[[nm]]$mv_post_epred_dim <- dim(posterior_epred(fit))
      res[[nm]]$mv_post_predict_dim <- dim(posterior_predict(fit))
      res[[nm]]$mv_loo_R2 <- q(loo_R2(fit, seed = 7))
    }
    if (nm == "gaussian") {
      # `scale` reaches nothing in brms 2.23.0
      res[[nm]]$scale_ignored <- identical(
        q(loo_linpred(fit, scale = "response")), q(loo_linpred(fit)))
      # ndraws: log_lik() and posterior_epred() each draw their own
      # subsample (validate_draw_ids() calls sample()), in that order
      set.seed(42)
      got <- q(loo_epred(fit, ndraws = 100))
      set.seed(42)
      i1 <- sample(4000, 100)
      i2 <- sample(4000, 100)
      ll1 <- log_lik(fit)[i1, ]
      ps1 <- q(loo::psis(-ll1, r_eff = loo::relative_eff(
        exp(ll1), chain_id = rep(1L, 100))))
      ref <- matrix(loo::E_loo(posterior_epred(fit)[i2, ], ps1)$value)
      colnames(ref) <- "mean"
      res[[nm]]$ndraws_other_draws <- identical(got, ref)
      res[[nm]]$ndraws_same_draws <- {
        ps2 <- q(loo::psis(-log_lik(fit)[i2, ], r_eff = 1))
        isTRUE(all.equal(got[, 1],
                         loo::E_loo(posterior_epred(fit)[i2, ], ps2)$value))
      }
    }
    print(utils::str(res[[nm]][c("loo_epred", "loo_R2")]))
  }
  res$lmw <- q(loo_model_weights(res$gaussian$fit, res$ranef$fit))
  set.seed(5)
  res$lmw_pb <- q(loo_model_weights(res$gaussian$fit, res$ranef$fit,
                                    method = "pseudobma"))
  res$lmw_names <- names(res$lmw)
  for (nm in names(models)) res[[nm]]$fit <- NULL
  saveRDS(res, file.path(out_dir, "parity-brms.rds"))
} else if (side == "frmtmb") {
  suppressMessages({library(frmtmb); library(frmtmb.sample)})
  cat("lib:", dirname(system.file(package = "frmtmb.sample")),
      format(packageVersion("frmtmb.sample")), "\n")
  res <- list()
  ds <- list()
  for (nm in names(models)) {
    m <- models[[nm]]
    form <- if (nm == "mvbf") {
      mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE)
    } else bf(m$f)
    fam <- if (m$fam == "poisson") poisson() else gaussian()
    t0 <- proc.time()[["elapsed"]]
    d <- q(frm_sample(form, family = fam, data = dd, chains = 4,
                      iter = 2000, warmup = 1000, seed = 1, refresh = 0))
    cat(nm, "sample s:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
    ds[[nm]] <- d
    r <- m$resp
    ll <- log_lik(d, resp = r)
    ps <- q(psis(d, resp = r))
    ep <- posterior_epred(d, resp = r)
    lp <- posterior_linpred(d, resp = r)
    res[[nm]] <- list(
      loo_epred = q(loo_epred(d, resp = r)),
      loo_linpred = q(loo_linpred(d, resp = r)),
      loo_epred_q = q(loo_epred(d, resp = r, type = "quantile",
                                probs = c(0.25, 0.75))),
      loo_pi = {set.seed(2); q(loo_predictive_interval(d, resp = r))},
      loo_R2 = q(loo_R2(d, resp = r, seed = 7)),
      loo_R2_draws = q(loo_R2(d, resp = r, seed = 7, summary = FALSE)),
      mcse_epred = mcse_mean(ep, ps), mcse_linpred = mcse_mean(lp, ps),
      ll = ll, ep = ep)
    if (nm == "mvbf") {
      res[[nm]]$mv_epred_dim <- dim(q(loo_epred(d)))
      res[[nm]]$mv_epred_dimnames <- dimnames(q(loo_epred(d)))
      res[[nm]]$mv_post_epred_dim <- dim(posterior_epred(d))
      res[[nm]]$mv_post_predict_dim <- dim(posterior_predict(d))
      res[[nm]]$mv_loo_R2 <- q(loo_R2(d, seed = 7))
    }
  }
  res$lmw <- q(loo_model_weights(ds$gaussian, ds$ranef))
  set.seed(5)
  res$lmw_pb <- q(loo_model_weights(ds$gaussian, ds$ranef,
                                    method = "pseudobma"))
  res$lmw_names <- names(res$lmw)
  saveRDS(res, file.path(out_dir, "parity-frmtmb.rds"))
} else {
  b <- readRDS(file.path(out_dir, "parity-brms.rds"))
  f <- readRDS(file.path(out_dir, "parity-frmtmb.rds"))
  zline <- function(lab, fv, bv, fse, bse) {
    z <- (fv - bv) / sqrt(fse^2 + bse^2)
    cat(sprintf("%-26s n=%d mean z=%6.2f sd z=%5.2f max|z|=%5.2f |z|>3: %d\n",
                lab, length(z), mean(z), stats::sd(z), max(abs(z)),
                sum(abs(z) > 3)))
  }
  for (nm in names(models)) {
    zline(paste(nm, "loo_epred mean"), f[[nm]]$loo_epred[, 1],
          b[[nm]]$loo_epred[, 1], f[[nm]]$mcse_epred, b[[nm]]$mcse_epred)
    zline(paste(nm, "loo_linpred mean"), f[[nm]]$loo_linpred[, 1],
          b[[nm]]$loo_linpred[, 1], f[[nm]]$mcse_linpred,
          b[[nm]]$mcse_linpred)
    cat(nm, "shapes  epred", dim(f[[nm]]$loo_epred), "/",
        dim(b[[nm]]$loo_epred), " quantile cols",
        colnames(f[[nm]]$loo_epred_q), "/", colnames(b[[nm]]$loo_epred_q),
        " pi cols", colnames(f[[nm]]$loo_pi), "/", colnames(b[[nm]]$loo_pi),
        "\n")
    cat(nm, "loo_R2 frmtmb:", sprintf("%.4f", f[[nm]]$loo_R2[1, ]),
        " rows", rownames(f[[nm]]$loo_R2), "\n")
    cat(nm, "loo_R2 brms:  ", sprintf("%.4f", b[[nm]]$loo_R2[1, ]),
        " rows", rownames(b[[nm]]$loo_R2), "\n")
    cat(nm, "loo_R2 (f - b) / brms Est.Error:",
        sprintf("%.3f", (f[[nm]]$loo_R2[1, 1] - b[[nm]]$loo_R2[1, 1]) /
                  b[[nm]]$loo_R2[1, 2]),
        " draws dims", dim(f[[nm]]$loo_R2_draws), "/",
        dim(b[[nm]]$loo_R2_draws), " colnames",
        colnames(f[[nm]]$loo_R2_draws), "/", colnames(b[[nm]]$loo_R2_draws),
        "\n")
  }
  cat("mv loo_epred dim frmtmb", f$mvbf$mv_epred_dim, "brms",
      b$mvbf$mv_epred_dim, "\n")
  cat("mv loo_epred dimnames frmtmb:",
      vapply(f$mvbf$mv_epred_dimnames, function(x) paste(x, collapse = ","),
             ""), "\n")
  cat("mv loo_epred dimnames brms:  ",
      vapply(b$mvbf$mv_epred_dimnames, function(x) paste(x, collapse = ","),
             ""), "\n")
  cat("mv posterior_epred(resp = NULL) dim frmtmb", f$mvbf$mv_post_epred_dim,
      "brms", b$mvbf$mv_post_epred_dim, "\n")
  cat("mv posterior_predict(resp = NULL) dim frmtmb",
      f$mvbf$mv_post_predict_dim, "brms", b$mvbf$mv_post_predict_dim, "\n")
  cat("mv loo_R2 rows frmtmb", rownames(f$mvbf$mv_loo_R2), "brms",
      rownames(b$mvbf$mv_loo_R2), "\n")
  print(rbind(frmtmb = f$mvbf$mv_loo_R2[, 1], brms = b$mvbf$mv_loo_R2[, 1],
              brms_err = b$mvbf$mv_loo_R2[, 2]))
  cat("loo_model_weights stacking frmtmb", sprintf("%.4f", f$lmw),
      "brms", sprintf("%.4f", b$lmw), "\n")
  cat("loo_model_weights pseudobma frmtmb", sprintf("%.4f", f$lmw_pb),
      "brms", sprintf("%.4f", b$lmw_pb), "\n")
  cat("classes", class(f$lmw), "/", class(b$lmw), "\n")
  cat("brms: scale ignored by loo_linpred():", b$gaussian$scale_ignored, "\n")
  cat("brms: loo_epred(ndraws = 100) weights other draws:",
      b$gaussian$ndraws_other_draws, " equals the same-draw answer:",
      b$gaussian$ndraws_same_draws, "\n")
  # exactness across implementations: brms's internals on frmtmb's
  # matrices are in dev/loodraws-probe-exact.R
}
