# Reviewer, claim 6: frmtmb.sample on a multivariate subset model.
# (a) every posterior method without resp: refused naming 'resp', or
# what it returns; (b) with one resp, the draws answers against core's
# ML machinery at the same draw: the objective at the draw's parameter
# vector (an IDENTITY: no prior, no random effect, so -fn = sum of the
# pointwise log-likelihoods) and fitted() of the fit set to that draw.
# Seeds: data 707, sampler 3. Log: dev/aterms2-rev-log-07-sample.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
q <- function(expr) suppressWarnings(suppressMessages(expr))
show <- function(label, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:",
                                                   conditionMessage(e)))
  cat(sprintf("%-32s ", label))
  if (is.character(r) && length(r) == 1L) cat(substr(r, 1, 170), "\n")
  else cat("returns", class(r)[1], paste(if (is.null(dim(r))) length(r)
                                         else dim(r), collapse = "x"), "\n")
  invisible(r)
}
set.seed(707)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                s1 = rep(c(TRUE, FALSE), n / 2),
                s2 = c(rep(TRUE, 40), rep(FALSE, 20)),
                time = runif(n, 0.5, 3))
d$y1 <- 1 + d$x + rnorm(n)
d$y2 <- rpois(n, exp(0.3 - 0.4 * d$z) * d$time)
d$y2[!d$s2] <- NA
d$z[!d$s2] <- NA
fit <- q(frm(bf(y1 | subset(s1) ~ x) + gaussian() +
               bf(y2 | subset(s2) + rate(time) ~ z) + poisson(), data = d))
ds <- q(frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3))
cat("draws", ndraws(ds), "\n")

cat("\n== (a) without resp\n")
show("log_lik", log_lik(ds))
show("loo", loo(ds))
show("waic", waic(ds))
show("psis", psis(ds))
show("loo_subsample", loo_subsample(ds))
show("posterior_epred", posterior_epred(ds))
show("posterior_linpred", posterior_linpred(ds))
show("posterior_predict", posterior_predict(ds))
show("predictive_error", predictive_error(ds))
show("predictive_interval", predictive_interval(ds))
show("fitted", fitted(ds))
show("predict", predict(ds))
show("residuals", residuals(ds))
show("pp_check", pp_check(ds))
show("bayes_R2", bayes_R2(ds))
show("conditional_effects", conditional_effects(ds))
show("kfold(K = 2)", kfold(ds, K = 2))
show("nobs", nobs(ds))
b2 <- q(bayes_R2(ds)); b1 <- q(bayes_R2(ds, resp = "y1"))
cat("  bayes_R2 all rows:", format(b2[, 1]), " resp y1:", format(b1[, 1]), "
")

cat("\n== (a') with resp = 'y1'\n")
show("log_lik", log_lik(ds, resp = "y1"))
show("loo", loo(ds, resp = "y1"))
show("waic", waic(ds, resp = "y1"))
show("fitted", fitted(ds, resp = "y1"))
show("predict", predict(ds, resp = "y1"))
show("residuals", residuals(ds, resp = "y1"))
show("predictive_error", predictive_error(ds, resp = "y1"))
show("bayes_R2", bayes_R2(ds, resp = "y1"))
show("pp_check", pp_check(ds, resp = "y1"))
show("log_lik(resp = c(y1, y2))", log_lik(ds, resp = c("y1", "y2")))

cat("\n== (b) against core at the same draw\n")
idx <- frmtmb.sample:::draws_par_index(ds$fit)
ll1 <- log_lik(ds, resp = "y1")
ll2 <- log_lik(ds, resp = "y2")
ep1 <- posterior_epred(ds, resp = "y1")
ep2 <- posterior_epred(ds, resp = "y2")
nd <- d[c(1:12, 41:46), ]
nd$time <- seq(0.5, 9, length.out = nrow(nd))
ep2n <- posterior_epred(ds, resp = "y2", newdata = nd)
worst <- c(ll = 0, ep1 = 0, ep2 = 0, ep2n = 0)
for (i in c(1L, 17L, 150L, ndraws(ds))) {
  row <- frmtmb.sample:::draws_internal_matrix(ds, i)[1L, ]
  if (i == 1L) cat("  draw columns:", names(row), " opt par", length(fit$opt$par), "
")
  row <- row[setdiff(seq_along(row), grep("^lp__$", names(row)))]
  tot <- -fit$obj$fn(row)
  s <- sum(ll1[i, ]) + sum(ll2[i, ])
  worst["ll"] <- max(worst["ll"], abs(s - tot) / abs(tot))
  sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
  f1 <- fitted(sh, resp = "y1")[, "Estimate"]
  f2 <- fitted(sh, resp = "y2")[, "Estimate"]
  f2n <- fitted(sh, resp = "y2", newdata = nd)[, "Estimate"]
  worst["ep1"] <- max(worst["ep1"], max(abs(ep1[i, ] - f1)) / max(abs(f1)))
  worst["ep2"] <- max(worst["ep2"], max(abs(ep2[i, ] - f2)) / max(abs(f2)))
  worst["ep2n"] <- max(worst["ep2n"],
                       max(abs(ep2n[i, ] - f2n)) / max(abs(f2n)))
}
cat("  columns: ll1", ncol(ll1), "(s1 rows", sum(d$s1), ") ll2", ncol(ll2),
    "(s2 rows", sum(d$s2), ") ep2n", ncol(ep2n), "(newdata s2 rows",
    sum(nd$s2), ")\n")
cat("  max relative residuals over 4 draws (eps =",
    .Machine$double.eps, "):\n")
print(worst)
