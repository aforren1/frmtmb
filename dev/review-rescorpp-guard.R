## Review rescorpp: guard absent. Univariate (and rescor = FALSE first
## layer) outputs must be bit-identical between rellib-r15 and the lane,
## on the draws route and the fit route.
##   Rscript dev/review-rescorpp-guard.R base make
##   Rscript dev/review-rescorpp-guard.R <base|lane> eval
##   Rscript dev/review-rescorpp-guard.R lane compare
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
mode <- commandArgs(TRUE)[2]
dsfile <- file.path(rv_out, "guard-ds.rds")
if (mode == "make") {
  set.seed(51); n <- 120
  x <- rnorm(n); g <- rep(1:8, length.out = n)
  d <- data.frame(x = x, g = g)
  d$yg <- 1 + x + rnorm(8)[g] + rnorm(n)
  d$yp <- rpois(n, exp(0.3 + 0.5 * x))
  d$ym <- ifelse(runif(n) < 0.4, rnorm(n, -2), rnorm(n, 2)) + 0.3 * x
  d$yc <- factor(sample(c("a", "b", "c"), n, TRUE))
  d$yo <- cut(x + rnorm(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
  yl <- exp(1 + 0.5 * x + rnorm(n, 0, 0.5))
  d$cen <- ifelse(yl > 6, "right", "none")
  d$yl <- pmin(yl, 6)
  d$ys <- 1 + x + exp(0.4 * x) * rt(n, 5)
  d$y2 <- 0.5 * x + rnorm(n)
  specs <- list(
    gauss_re = list(f = bf(yg ~ x + (1 | g)), fam = gaussian()),
    poisson = list(f = bf(yp ~ x), fam = poisson()),
    mixture = list(f = bf(ym ~ x), fam = mixture(gaussian(), gaussian())),
    categorical = list(f = bf(yc ~ x), fam = categorical()),
    censored = list(f = bf(yl | cens(cen) ~ x), fam = lognormal()),
    ordinal = list(f = bf(yo ~ x), fam = cumulative()),
    student_sig = list(f = bf(ys ~ x, sigma ~ x), fam = student()),
    mv_norescor = list(f = mvbf(bf(yg ~ x), bf(y2 ~ x)) + set_rescor(FALSE),
                       fam = gaussian()))
  res <- list()
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    cat("fit", nm, "\n")
    fit <- tryCatch(quiet(frm(sp$f, family = sp$fam, data = d)),
                    error = function(e) {cat("  ERROR", conditionMessage(e), "\n"); NULL})
    if (is.null(fit)) next
    ds <- tryCatch(quiet(frm_sample(fit, chains = 1, iter = 400, refresh = 0,
                                    seed = 3)),
                   error = function(e) {cat("  ERROR", conditionMessage(e), "\n"); NULL})
    res[[nm]] <- list(fit = fit, ds = ds, d = d)
  }
  saveRDS(res, dsfile)
  quit()
}
if (mode == "eval") {
  res <- readRDS(dsfile)
  ev <- list()
  tq <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)),
                                error = function(e) paste("ERROR", conditionMessage(e)))
  for (nm in names(res)) {
    fit <- res[[nm]]$fit; ds <- res[[nm]]$ds; d <- res[[nm]]$d
    nd <- d[1:6, ]
    r1 <- if (nm == "mv_norescor") "yg" else NULL
    m <- list(
      pp = function() posterior_predict(ds),
      pp_resp = function() posterior_predict(ds, resp = r1),
      pp_nd = function() posterior_predict(ds, newdata = nd),
      epred = function() posterior_epred(ds),
      epred_resp = function() posterior_epred(ds, resp = r1),
      linpred = function() posterior_linpred(ds),
      perr = function() predictive_error(ds),
      perr_resp = function() predictive_error(ds, resp = r1),
      fitted = function() fitted(ds),
      fitted_resp = function() fitted(ds, resp = r1),
      predict = function() predict(ds),
      predict_resp = function() predict(ds, resp = r1),
      resid = function() residuals(ds),
      pint = function() predictive_interval(ds),
      pint_resp = function() predictive_interval(ds, resp = r1),
      loglik = function() log_lik(ds),
      loglik_resp = function() log_lik(ds, resp = r1),
      loo = function() loo(ds)$pointwise,
      waic = function() waic(ds)$pointwise,
      loo_pred = function() loo_predict(ds),
      loo_ep = function() loo_epred(ds),
      bR2 = function() bayes_R2(ds),
      looR2 = function() loo_R2(ds),
      ce = function() as.data.frame(conditional_effects(ds, ndraws = 50)[[1]]),
      ce_pred = function() as.data.frame(conditional_effects(
        ds, ndraws = 50, method = "posterior_predict")[[1]]),
      ppc = function() pp_check(ds, ndraws = 5, resp = r1)$data,
      ## fit route
      f_predict = function() predict(fit, nsim = 50),
      f_fitted = function() fitted(fit),
      f_resid = function() residuals(fit),
      f_sim = function() simulate(fit, nsim = 3, seed = 9),
      f_loglik = function() logLik(fit))
    ev[[nm]] <- lapply(m, function(fn) {set.seed(11); tq(fn())})
    cat(nm, "evaluated\n")
  }
  saveRDS(ev, file.path(rv_out, paste0("guard-eval-", arm, ".rds")))
  quit()
}
if (mode == "compare") {
  a <- readRDS(file.path(rv_out, "guard-eval-base.rds"))
  b <- readRDS(file.path(rv_out, "guard-eval-lane.rds"))
  ni <- 0; nn <- 0
  for (nm in names(a)) {
    for (k in names(a[[nm]])) {
      x <- a[[nm]][[k]]; y <- b[[nm]][[k]]
      if (nm == "mv_norescor" && k %in% c("pp", "epred", "linpred", "perr",
                                          "fitted", "predict", "resid")) {
        # the old 2D first response against the new first layer
        yy <- if (is.array(y) && length(dim(y)) == 3L) {
          if (k %in% c("fitted", "predict", "resid")) y[, , 1] else unname(y[, , 1])
        } else y
        if (k %in% c("fitted", "predict", "resid") && is.array(yy)) {
          dimnames(yy) <- dimnames(x)
        }
        same <- identical(x, yy)
        cat(sprintf("%-12s %-12s %s (first layer vs old; new dim %s)\n", nm, k,
                    if (same) "IDENTICAL" else "DIFFERS",
                    paste(dim(y), collapse = "x")))
        if (!same) {cat("   base:", substr(paste(format(x)[1:2], collapse=" "), 1, 150),
                        "\n   lane:", substr(paste(format(yy)[1:2], collapse=" "), 1, 150), "\n")}
        if (same) ni <- ni + 1 else nn <- nn + 1
        next
      }
      same <- identical(x, y)
      if (same) ni <- ni + 1 else nn <- nn + 1
      cat(sprintf("%-12s %-12s %s\n", nm, k, if (same) "IDENTICAL" else "DIFFERS"))
      if (!same) {
        cat("   base:", substr(paste(utils::capture.output(str(x)), collapse = " "), 1, 200), "\n")
        cat("   lane:", substr(paste(utils::capture.output(str(y)), collapse = " "), 1, 200), "\n")
      }
    }
  }
  cat("GUARD identical", ni, "differs", nn, "\n")
}
