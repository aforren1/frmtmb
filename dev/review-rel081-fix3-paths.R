# Release review of 0.81.0, short check of the fixes: which criterion
# paths of frmtmb.sample reach loo_check_finite() (class
# frmtmb_loglik_not_finite) on draws with -Inf log_lik() cells, and which
# stop elsewhere or answer. The draws are the new test's (seed 349).
#   Rscript dev/review-rel081-fix3-paths.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb.sample")), "sample",
    format(packageVersion("frmtmb.sample")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(349)
d <- data.frame(x = stats::rnorm(40))
d$y <- rg(40, 1 + 0.3 * d$x, 1, 0.3)
fit <- suppressWarnings(frm(y ~ x, family = gen_extreme_value(), data = d))
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, seed = 349, refresh = 0)))
i <- which.min(d$y)
dr <- ds$draws
k <- which(dr[, "xi"] > 0.05)[1:3]
for (j in k) {
  mu <- d$y[i] + dr[j, "sigma"] / dr[j, "xi"] + 1
  dr[j, "b_Intercept"] <- mu - dr[j, "b_x"] * d$x[i]
}
mv <- ds
mv$draws <- dr
cat("draws with xi > 0.05:", sum(dr[, "xi"] > 0.05), "of", nrow(dr), "\n")
ll <- log_lik(mv)
cat("-Inf cells:", sum(ll == -Inf), "\n")
try_one <- function(label, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-38s ERROR class %s: %s\n", label,
                if (inherits(r, "frmtmb_loglik_not_finite"))
                  "frmtmb_loglik_not_finite" else class(r)[1],
                substr(gsub("\\s+", " ", conditionMessage(r)), 1, 110)))
  } else {
    cat(sprintf("%-38s ok: %s\n", label, paste(class(r), collapse = "/")))
  }
}
y <- d$y
try_one("loo(mv)", loo(mv))
try_one("waic(mv)", waic(mv))
try_one("psis(mv)", psis(mv))
try_one("psis(mv, r_eff = 1)", psis(mv, r_eff = 1))
try_one("loo(mv, ds) compare", loo(mv, ds))
try_one("loo_compare(mv, ds)", loo_compare(mv, ds))
try_one("loo_model_weights(mv, ds)", loo_model_weights(mv, ds))
try_one("add_criterion(mv, 'loo')", add_criterion(mv, "loo"))
try_one("add_criterion(mv, 'waic')", add_criterion(mv, "waic"))
try_one("loo_predict(mv)", loo_predict(mv))
try_one("loo_linpred(mv)", loo_linpred(mv))
try_one("loo_epred(mv)", loo_epred(mv))
try_one("loo_R2(mv)", loo_R2(mv))
try_one("pp_check(mv, 'loo_pit_overlay')", pp_check(mv, type = "loo_pit_overlay"))
try_one("pp_check(mv, 'loo_intervals')", pp_check(mv, type = "loo_intervals"))
try_one("pp_check(mv, 'loo_pit_qq')", pp_check(mv, type = "loo_pit_qq"))
try_one("loo_moment_match(mv)", loo_moment_match(mv))
try_one("kfold(mv)", kfold(mv))
try_one("loo::waic(log_lik(mv)) (user matrix)", loo::waic(ll))
try_one("loo(ds) control", loo(ds))
try_one("pp_check(ds, 'loo_pit_overlay') ctrl", pp_check(ds, type = "loo_pit_overlay"))
cat("PATHS DONE\n")
