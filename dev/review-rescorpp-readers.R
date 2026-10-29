## Review rescorpp: other readers of multivariate draws, saved 0.77.0
## draws, and new levels under rescor.
##   Rscript dev/review-rescorpp-readers.R base   (saves 0.77.0 draws)
##   Rscript dev/review-rescorpp-readers.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
source(file.path(wt, "dev/rescorpp-data.R"))
source(file.path(wt, "dev/review-rescorpp-data.R"))
d <- rescorpp_data("g3")
saved <- file.path(rv_out, "ds-g3-saved-base.rds")
if (arm == "base") {
  fit <- quiet(frm(eval(rescorpp_frm_forms[["g3"]]), family = "gaussian",
                   data = d))
  ds <- quiet(frm_sample(fit, chains = 2, iter = 600, refresh = 0, seed = 21))
  ds <- add_criterion(ds, c("loo", "loo_R2"))
  saveRDS(ds, saved)
}
ds <- readRDS(saved)
cat("saved draws from rellib-r15, read on", arm, "\n")
show("pp no resp dim", dim(posterior_predict(ds, ndraws = 20)))
show("epred no resp dim", dim(posterior_epred(ds, ndraws = 20)))
show("fitted no resp dim", dim(fitted(ds, ndraws = 20)))
show("log_lik resp y1 [1:2,1:3]", log_lik(ds, resp = "y1")[1:2, 1:3])
show("loo stored", loo(ds)$estimates)
show("loo resp y1", loo(ds, resp = "y1")$estimates)
show("loo_R2 stored", ds$criteria$loo_R2[1:2, ])
show("loo_R2 now", loo_R2(ds, summary = FALSE)[1:2, ])
show("hypothesis y1_x > 0", hypothesis(ds, "y1_x > 0")$hypothesis)
show("hypothesis rescor", hypothesis(ds, "rescor__y1__y2 > 0.5",
                                     class = NULL)$hypothesis)
show("bayes_R2", bayes_R2(ds))
show("conditional_effects no resp names",
     names(conditional_effects(ds, ndraws = 30)))
show("conditional_effects resp y2 names",
     names(conditional_effects(ds, resp = "y2", ndraws = 30)))
show("conditional_effects resp c(y1,y2)",
     names(conditional_effects(ds, resp = c("y1", "y2"), ndraws = 30)))
show("conditional_effects predict resp c(y1,y2)",
     names(conditional_effects(ds, resp = c("y1", "y2"), ndraws = 30,
                               method = "posterior_predict")))
show("emmeans", as.data.frame(emmeans::emmeans(ds, ~ x, at = list(x = 0))))
show("emmeans resp y2", as.data.frame(emmeans::emmeans(ds, ~ x, resp = "y2",
                                                       at = list(x = 0))))
show("marginaleffects predictions", head(as.data.frame(
  marginaleffects::predictions(ds, newdata = d[1:2, ]))))
show("marginaleffects avg_slopes", as.data.frame(
  marginaleffects::avg_slopes(ds, variables = "x")))
show("pp_check resp y2", class(pp_check(ds, resp = "y2", ndraws = 5)))
show("pp_check no resp", class(pp_check(ds, ndraws = 5)))
show("predictive_interval resp y3", dim(predictive_interval(ds, resp = "y3")))
show("posterior_smooths / summary", class(summary(ds)))
show("update() exists", getS3method("update", "frmtmb_draws", optional = TRUE))

if (arm == "lane") {
  ## new levels under rescor with (1 | p | g): brms draws the new level's
  ## effects jointly across responses
  suppressMessages(loadNamespace("brms"))
  dre <- rv_data("re")
  b <- readRDS(file.path(rv_out, "brms-re.rds"))
  dsr <- readRDS(file.path(rv_out, "ds-re-lane.rds"))
  nd <- data.frame(x = 0, g = 99)
  for (snl in c("gaussian", "uncertainty")) {
    set.seed(3)
    pb <- tryCatch(brms::posterior_predict(b, newdata = nd, allow_new_levels = TRUE,
                                           sample_new_levels = snl),
                   error = function(e) e)
    set.seed(3)
    pf <- tryCatch(quiet(posterior_predict(dsr, newdata = nd, allow_new_levels = TRUE,
                                           sample_new_levels = snl)),
                   error = function(e) e)
    eb <- tryCatch(brms::posterior_epred(b, newdata = nd, allow_new_levels = TRUE,
                                         sample_new_levels = snl),
                   error = function(e) e)
    ef <- tryCatch(quiet(posterior_epred(dsr, newdata = nd, allow_new_levels = TRUE,
                                         sample_new_levels = snl)),
                   error = function(e) e)
    f <- function(p) if (inherits(p, "error")) conditionMessage(p) else
      sprintf("dim %s  cor(y1,y2) %.3f  sd %s", paste(dim(p), collapse = "x"),
              cor(p[, 1, 1], p[, 1, 2]), paste(round(apply(p[, 1, ], 2, sd), 3),
                                               collapse = ","))
    cat("\nnew level, sample_new_levels =", snl, "\n  brms   pp:", f(pb),
        "\n  frmtmb pp:", f(pf), "\n  brms   epred:", f(eb),
        "\n  frmtmb epred:", f(ef), "\n")
  }
  ## existing level with newdata: cor of draws at g = 1
  set.seed(4)
  pf1 <- posterior_predict(dsr, newdata = data.frame(x = 0, g = 1))
  pb1 <- brms::posterior_predict(b, newdata = data.frame(x = 0, g = 1))
  cat("existing level g = 1: cor(y1,y2) frmtmb", round(cor(pf1[, 1, 1], pf1[, 1, 2]), 3),
      " brms", round(cor(pb1[, 1, 1], pb1[, 1, 2]), 3), "\n")
}
