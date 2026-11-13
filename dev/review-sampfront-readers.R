# Review of lane sampfront: which readers take a subset of the fitted rows,
# and what add_criterion() stores from one.
#   Rscript dev/review-sampfront-readers.R
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb.sample", find.package("frmtmb.sample"), "\n")
full_draws <- function(fit, n = 40L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
show <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  cat(sprintf("%-46s ", lab))
  if (inherits(r, "error")) {
    cat("REFUSED:", substr(gsub("\\s+", " ", conditionMessage(r)), 1, 110),
        "\n")
  } else {
    cat("ok:", class(r)[1], if (!is.null(dim(r))) paste(dim(r), collapse = "x")
        else if (inherits(r, "loo")) paste("n =", nrow(r$pointwise))
        else length(r), "\n")
  }
  invisible(r)
}
set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rnorm(48))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] * 0.7 + rnorm(48)
x <- full_draws(q(bf(y ~ x + (1 | g)), data = d))
x0 <- full_draws(q(bf(y ~ 1 + (1 | g)), data = d))
x$data <- d; x0$data <- d
s <- d[1:10, ]
show("log_lik", log_lik(x, newdata = s))
show("logLik", logLik(x, newdata = s))
show("loo", loo(x, newdata = s))
show("waic", waic(x, newdata = s))
show("psis", psis(x, s))
show("loo(x, x0, newdata)", loo(x, x0, newdata = s))
show("loo_model_weights", loo_model_weights(x, x0, newdata = s))
show("model_weights(stacking)", model_weights(x, x0, newdata = s))
show("model_weights(loo)", model_weights(x, x0, weights = "loo", newdata = s))
show("posterior_average(weights = 'loo', ctrl)",
     posterior_average(x, x0, weights = "loo", variable = "b_Intercept",
                       control = list(newdata = s)))
show("pp_check(type = 'loo_pit_overlay')",
     pp_check(x, type = "loo_pit_overlay", newdata = s))
show("loo_R2", loo_R2(x, newdata = s))
show("loo_predict", loo_predict(x, newdata = s))
show("loo_epred", loo_epred(x, newdata = s))
show("loo_linpred", loo_linpred(x, newdata = s))
show("loo_predictive_interval", loo_predictive_interval(x, newdata = s))
show("kfold", kfold(x, newdata = s))
show("loo_compare(x, x0)", loo_compare(x, x0))
xa <- add_criterion(x, "loo", newdata = s)
show("add_criterion(loo, newdata = s) stored", xa$criteria$loo)
show("then loo(xa) (no newdata)", loo(xa))
x0a <- add_criterion(x0, "loo")
show("loo_compare(xa, x0a)", loo_compare(xa, x0a))
# the common ways to take a subset that do not keep row names
sf <- dplyr::filter(d, g == "2")
cat("dplyr::filter(d, g == '2') row names:", rownames(sf), "\n")
show("log_lik(newdata = dplyr::filter(d, g == '2'))", log_lik(x, newdata = sf))
show("log_lik(newdata = dplyr::slice(d, 9:16))",
     log_lik(x, newdata = dplyr::slice(d, 9:16)))
show("log_lik(newdata = subset(d, g == '2'))",
     log_lik(x, newdata = subset(d, g == "2")))
show("log_lik(newdata = d[d$g == '2', ])", log_lik(x, newdata = d[d$g == "2", ]))
cat("DONE\n")
