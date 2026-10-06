# The post-fit paths of a cs() cumulative() fit on the release build,
# beside brms 2.23.0's own draws where brms has an answer: fitted(),
# predict(), simulate(), conditional_effects(), emmeans(), confint()
# names, update(), and frmtmb.sample's draws (names, posterior_epred()
# against brms's R-side density, posterior_predict()).
#
#   Rscript dev/rel068-cs-probe.R > dev/rel068-log/cs-probe.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-release/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(20261006)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- rlogis(n) + 0.5 * d$x + 0.4 * d$z
d$y <- 1L + (u > -1 + 0.15 * d$x) + (u > 0.3) + (u > 1.5 - 0.15 * d$x)
step <- function(tag, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  cat(sprintf("%-34s %s%s\n", tag,
              if (inherits(r, "err")) paste("ERROR:", r) else "ok",
              if (length(w)) paste0(" [", length(w), " warning(s): ",
                                    substr(w[1], 1, 70), "]") else ""))
  invisible(r)
}
fit <- step("frm", frm(y ~ z + cs(x), family = cumulative(), data = d))
step("summary", capture.output(summary(fit)))
P <- step("fitted", fitted(fit))
cat("  fitted dims", dim(P), "rows summing to 1:",
    sum(abs(rowSums(P[, "Estimate", ]) - 1) < 1e-12), "of", n, "\n")
step("predict", predict(fit))
s <- step("simulate", simulate(fit, nsim = 2, seed = 1))
cat("  simulate NA count", sum(is.na(unlist(s))), "\n")
step("confint names", print(rownames(confint(fit))))
step("conditional_effects", conditional_effects(fit, "x"))
step("emmeans", emmeans::emmeans(fit, ~ z))
step("residuals osa", summary(residuals(fit, type = "osa")[, "Estimate"]))
step("update", update(fit, . ~ . - z))
step("loo-free logLik", logLik(fit))
ds <- step("frm_sample", frm_sample(y ~ z + cs(x), family = cumulative(),
                                    data = d, chains = 1, iter = 300,
                                    refresh = 0, seed = 3,
                                    prior = set_prior("normal(0, 2)",
                                                      class = "b")))
if (!inherits(ds, "err")) {
  cat("  variables:", grep("bcs|Intercept", variables(ds), value = TRUE),
      "\n")
  ep <- step("posterior_epred", posterior_epred(ds, ndraws = 3))
  cat("  epred dims", dim(ep), "\n")
  pp <- step("posterior_predict", posterior_predict(ds, ndraws = 3))
  cat("  predict NA", sum(is.na(pp)), "of", length(pp), "\n")
  dr <- as.matrix(ds$draws)[1, ]
  tau <- dr[paste0("b_Intercept[", 1:3, "]")]
  b <- dr[paste0("bcs_x[", 1:3, "]")]
  thr <- matrix(tau, n, 3, byrow = TRUE) - outer(d$x, b)
  eta <- dr[["b_z"]] * d$z
  ref <- brms:::dcumulative(1:4, eta = eta, thres = thr, disc = 1,
                            link = "logit")
  e1 <- posterior_epred(ds, draw_ids = 1)
  cat("  epred draw 1 vs brms dcumulative: max abs diff",
      format(max(abs(e1[1, , ] - ref)), digits = 3), "\n")
}
cat("PROBE DONE\n")
