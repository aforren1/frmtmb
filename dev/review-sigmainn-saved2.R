# Short check, fix round 1: the saved 0.76.0 fits of
# dev/review-sigmainn-saved.R (saved-r14-fits.rds), every reader run on
# rellib-r14 (arm r14) and on the lane build (arm lane); values saved, so
# "answering" methods can be compared for identity. Part B: fresh fits,
# insight::get_sigma() on both builds.
# Usage: Rscript dev/review-sigmainn-saved2.R r14|lane|compare
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
out <- file.path(wt, "dev", "review-sigmainn-out")
if (arm == "compare") {
  A <- readRDS(file.path(out, "saved2-r14.rds"))
  B <- readRDS(file.path(out, "saved2-lane.rds"))
  for (m in names(B)) {
    cat("\n==", m, "\n")
    for (k in names(B[[m]])) {
      a <- A[[m]][[k]]; b <- B[[m]][[k]]
      st <- if (inherits(b, "refused")) paste("REFUSED on lane:", substr(b, 1, 70)) else
        if (inherits(a, "refused")) paste("refused on r14, lane answers") else
        if (identical(a, b)) "answers, IDENTICAL to r14" else {
          na <- suppressWarnings(as.numeric(unlist(a))); nb <- suppressWarnings(as.numeric(unlist(b)))
          if (length(na) == length(nb) && length(na))
            sprintf("answers, DIFFERS: max rel %.2e", max(abs(na - nb) / pmax(abs(na), 1e-300), na.rm = TRUE))
          else "answers, DIFFERS (shape)"
        }
      cat(sprintf("  %-32s %s\n", k, st))
    }
  }
  quit(save = "no")
}
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
fits <- readRDS(file.path(out, "saved-r14-fits.rds"))
nd <- fits$plain$frame$data_frame
grab <- function(expr) {
  set.seed(5)
  tryCatch(suppressWarnings(suppressMessages({
    v <- force(expr)
    if (inherits(v, "emmGrid")) v <- as.data.frame(summary(v))
    if (is.list(v) && !is.data.frame(v)) v <- lapply(v, function(z) if (is.data.frame(z)) as.data.frame(z) else z)
    v
  })), error = function(e) structure(conditionMessage(e), class = "refused"))
}
res <- list()
for (m in if (nzchar(Sys.getenv("REV_PARTB"))) character(0) else names(fits)) {
  f <- fits[[m]]
  r <- list()
  r$print <- grab(capture.output(print(f)))
  r$summary_fixed <- grab(summary(f)$fixed)
  r$sigma <- grab(sigma(f))
  r$VarCorr <- grab(VarCorr(f))
  r$fixef <- grab(fixef(f))
  r$coef <- grab(coef(f))
  r$vcov <- grab(vcov(f))
  r$vcov_full <- grab(vcov(f, full = TRUE))
  r$confint <- grab(confint(f))
  r$confint_x <- grab(confint(f, parm = "x"))
  r$logLik <- grab(logLik(f))
  r$fitted <- grab(fitted(f))
  r$fitted_mu <- grab(fitted(f, dpar = "mu"))
  r$fitted_sigma <- grab(fitted(f, dpar = "sigma"))
  r$frm_linpred <- grab(frm_linpred(f))
  r$frm_linpred_se <- grab(frm_linpred(f, se.fit = TRUE))
  r$frm_linpred_sigma <- grab(frm_linpred(f, dpar = "sigma"))
  r$predict <- grab(predict(f))
  r$resid_response <- grab(residuals(f, type = "response"))
  r$resid_pearson <- grab(residuals(f, type = "pearson"))
  r$simulate <- grab(simulate(f, nsim = 1, seed = 1))
  r$emmeans_mu <- grab(emmeans::emmeans(f, ~ 1))
  r$emmeans_x <- grab(emmeans::emmeans(f, ~ x, at = list(x = c(-1, 1))))
  r$emmeans_sigma <- grab(emmeans::emmeans(f, ~ 1, dpar = "sigma"))
  r$ce_mu <- grab(conditional_effects(f, effects = "x")[[1]][, c("estimate__", "lower__", "upper__")])
  r$ce_pp <- grab(conditional_effects(f, effects = "x", method = "posterior_predict")[[1]][, c("estimate__", "lower__", "upper__")])
  r$ce_sigma <- grab(conditional_effects(f, effects = "x", dpar = "sigma")[[1]][, c("estimate__", "lower__", "upper__")])
  r$me_predictions <- grab(as.data.frame(marginaleffects::predictions(f))[, c("estimate", "std.error")])
  r$me_avg_slopes <- grab(as.data.frame(marginaleffects::avg_slopes(f))[, c("estimate", "std.error")])
  r$insight_get_sigma <- grab(as.numeric(insight::get_sigma(f)))
  r$insight_get_parameters <- grab(insight::get_parameters(f))
  r$insight_get_varcov <- grab(insight::get_varcov(f))
  r$insight_get_predicted <- grab(as.numeric(insight::get_predicted(f)))
  r$hypothesis_x <- grab(hypothesis(f, "x > 0")$hypothesis)
  r$cooks <- grab(cooks.distance(f))
  r$nobs <- grab(nobs(f))
  r$AIC <- grab(AIC(f))
  r$update_logLik <- grab(as.numeric(logLik(update(f))))
  res[[m]] <- r
}
if (!nzchar(Sys.getenv("REV_PARTB"))) saveRDS(res, file.path(out, paste0("saved2-", arm, ".rds")))
# Part B: insight::get_sigma() on fresh fits
set.seed(3)
d <- data.frame(x = rnorm(120), x2 = rnorm(120), g = factor(rep(1:12, 10)), t = rep(1:10, each = 12))
d$y <- 1 + d$x + rnorm(120, sd = exp(0.2 * d$x2))
d$cnt <- rpois(120, exp(0.3 + 0.2 * d$x))
d$bin <- rbinom(120, 1, plogis(0.3 * d$x))
d$y2 <- d$y + rnorm(120)
cases <- list(
  gaussian = function() frm(y ~ x, data = d),
  sigmax = function() frm(bf(y ~ x, sigma ~ x2), data = d),
  student = function() frm(y ~ x, family = student(), data = d),
  poisson = function() frm(cnt ~ x, family = poisson(), data = d),
  binomial = function() frm(bin ~ x, family = bernoulli(), data = d),
  re = function() frm(y ~ x + (1 | g), data = d),
  arcov = function() frm(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d),
  mv = function() frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = TRUE), data = d))
for (nm in names(cases)) {
  f <- cases[[nm]]()
  v <- tryCatch(withCallingHandlers(insight::get_sigma(f), warning = function(w) {
    cat("   [warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }),
    error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("get_sigma %-9s %s | sigma() %s\n", nm,
              paste(format(unclass(v), digits = 8), collapse = " "),
              paste(format(suppressWarnings(sigma(f)), digits = 8), collapse = " ")))
}
