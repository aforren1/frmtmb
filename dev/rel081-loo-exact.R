# Merge check, frmtmb.sample loo.R (hollow x support): log_lik() and
# loo() of draws with newdata equal to the fitted rows must read the
# exact density, as the call without newdata does. A GEV fit is sampled,
# then three draws are moved so that the smallest response lies outside
# their support (t <= 0); there the exact density is -Inf and the fitted
# objective's barrier is finite (about -1e150 or larger). Both calls must
# give -Inf in those cells and agree everywhere.
#
#   Rscript dev/rel081-loo-exact.R [lib]     (default rellib-r20)
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/rellib-r20"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(349)
d <- data.frame(x = stats::rnorm(40))
d$y <- rg(40, 1 + 0.3 * d$x, 1, 0.3)
fit <- suppressWarnings(frm(y ~ x, family = gen_extreme_value(), data = d))
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, seed = 349, refresh = 0)))
dm <- as.data.frame(ds)
cat("draws:", nrow(dm), " columns:", paste(names(dm), collapse = " "), "\n")
cat("draws object elements:", paste(names(ds), collapse = " "), "\n")
# move three draws: mu at the smallest y is put above y + sigma / xi
i <- which.min(d$y)
k <- 1:3
ds2 <- ds
dr <- ds2$draws
cn <- colnames(dr)
cat("internal draws columns:", paste(cn, collapse = " "), "\n")
ll_ref <- log_lik(ds)
cat("unmoved: log_lik finite cells", sum(is.finite(ll_ref)), "of",
    length(ll_ref), "\n")
cat("unmoved: identical(log_lik(ds), log_lik(ds, newdata = ds$data))",
    identical(ll_ref, log_lik(ds, newdata = ds$data)), "\n")
# the moved draws: mu_i = b0 + b1 x_i is set to y_i + sigma / xi + 1,
# so t_i = 1 + xi (y_i - mu_i) / sigma = -xi / sigma < 0
for (j in k) {
  mu_new <- d$y[i] + dr[j, "sigma"] / dr[j, "xi"] + 1
  dr[j, "b_Intercept"] <- mu_new - dr[j, "b_x"] * d$x[i]
}
ds2$draws <- dr
tt <- 1 + dr[k, "xi"] * (d$y[i] - (dr[k, "b_Intercept"] + dr[k, "b_x"] *
                                    d$x[i])) / dr[k, "sigma"]
cat("moved draws", k, ": t at row", i, "=", format(tt, digits = 4), "\n")
ll0 <- log_lik(ds2)
ll1 <- log_lik(ds2, newdata = ds2$data)
ll2 <- log_lik(ds2, newdata = d)
cat("log_lik(draws)[moved, row]:", format(ll0[k, i]), "\n")
cat("log_lik(draws, newdata = draws$data)[moved, row]:", format(ll1[k, i]), "\n")
cat("log_lik(draws, newdata = data)[moved, row]:", format(ll2[k, i]), "\n")
cat("identical(no newdata, newdata = draws$data):", identical(ll0, ll1), "\n")
cat("identical(no newdata, newdata = data):", identical(ll0, ll2), "\n")
cat("cells -Inf:", sum(ll0 == -Inf), " other non-finite:",
    sum(!is.finite(ll0) & ll0 != -Inf), "\n")
cat("unmoved draws equal to the unmoved object:",
    identical(ll0[-k, ], ll_ref[-k, ]), "\n")
# the barrier the objective keeps, for the record: row_lpdf without exact
fam <- gen_extreme_value()
mu_k <- dr[k, "b_Intercept"] + dr[k, "b_x"] * d$x[i]
pk <- list(mu = unname(mu_k), sigma = unname(dr[k, "sigma"]),
           xi = unname(dr[k, "xi"]))
yk <- rep(d$y[i], length(k))
cat("row_lpdf(exact = FALSE) at the moved cells:",
    format(as.numeric(row_lpdf(fam, yk, yk, pk, list(), NULL)), digits = 4), "\n")
cat("row_lpdf(exact = TRUE) at the moved cells:",
    format(as.numeric(row_lpdf(fam, yk, yk, pk, list(), NULL, exact = TRUE))), "\n")
# loo() through both spellings: the pointwise elpd at the row
l0 <- tryCatch(suppressWarnings(loo(ds2)), error = function(e) e)
l1 <- tryCatch(suppressWarnings(loo(ds2, newdata = d)), error = function(e) e)
if (inherits(l0, "error") || inherits(l1, "error")) {
  cat("loo(): no newdata:", if (inherits(l0, "error")) conditionMessage(l0) else "ok",
      "| newdata:", if (inherits(l1, "error")) conditionMessage(l1) else "ok", "\n")
  cat("loo() agree (both error with the same message):",
      inherits(l0, "error") && inherits(l1, "error") &&
        identical(conditionMessage(l0), conditionMessage(l1)), "\n")
} else {
  cat("loo() elpd_loo:", l0$estimates["elpd_loo", "Estimate"],
      l1$estimates["elpd_loo", "Estimate"], " identical estimates:",
      identical(l0$estimates, l1$estimates), " identical pointwise:",
      identical(l0$pointwise, l1$pointwise), "\n")
}
cat("LOO-EXACT DONE\n")
