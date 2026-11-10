# Which value stops loo()/waic()/psis() of draws whose log_lik() holds a
# -Inf cell (a GEV row outside a draw's support, read exactly), and what
# brms 2.23.0 gives for the same draw. The draws are those of
# dev/rel081-loo-exact.R (GEV y ~ x, seed 349, 200 draws, three draws
# moved so the smallest y is below their lower end).
#
#   Rscript dev/rel081-loo-inf.R [lib]     (default rellib-r20)
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/rellib-r20"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "loo",
    format(packageVersion("loo")), "brms", format(packageVersion("brms")),
    "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(349)
d <- data.frame(x = stats::rnorm(40))
d$y <- rg(40, 1 + 0.3 * d$x, 1, 0.3)
fit <- suppressWarnings(frm(y ~ x, family = gen_extreme_value(), data = d))
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, seed = 349, refresh = 0)))
i <- which.min(d$y); k <- 1:3
dr <- ds$draws
for (j in k) {
  mu_new <- d$y[i] + dr[j, "sigma"] / dr[j, "xi"] + 1
  dr[j, "b_Intercept"] <- mu_new - dr[j, "b_x"] * d$x[i]
}
ds2 <- ds; ds2$draws <- dr
ll <- log_lik(ds2)
cat("log_lik cells: -Inf", sum(ll == -Inf), "+Inf", sum(ll == Inf), "NaN",
    sum(is.nan(ll)), "of", length(ll), "\n")
lr <- -ll
cat("-log_lik (the importance log ratios): +Inf", sum(lr == Inf), "NaN",
    sum(is.nan(lr)), "\n")
msg <- function(expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) e)
  if (inherits(r, "error")) paste("ERROR:", conditionMessage(r)) else "ok"
}
cat("loo::validate_ll(log_lik):", msg(loo:::validate_ll(ll)), "\n")
cat("loo::validate_ll(-log_lik):", msg(loo:::validate_ll(lr)), "\n")
cat("frmtmb loo(draws):", msg(loo(ds2)), "\n")
cat("frmtmb waic(draws):", msg(waic(ds2)), "\n")
cat("frmtmb psis(draws):", msg(psis(ds2)), "\n")
re <- rep(1, ncol(ll))
cat("loo::loo.matrix(log_lik):", msg(loo::loo(ll, r_eff = re)), "\n")
cat("loo::psis(-log_lik):", msg(loo::psis(lr, r_eff = re)), "\n")
w <- tryCatch(suppressWarnings(loo::waic(ll)), error = function(e) e)
if (inherits(w, "error")) {
  cat("loo::waic.matrix(log_lik): ERROR:", conditionMessage(w), "\n")
} else {
  cat("loo::waic.matrix(log_lik): ok; elpd_waic", w$estimates["elpd_waic", 1],
      "; pointwise non-finite rows", sum(!is.finite(w$pointwise[, 1])), "\n")
}
# brms 2.23.0's density at the same cells: its R-side GEV density, which
# brms's log_lik_gen_extreme_value() reads
cat("brms:::dgen_extreme_value body:\n")
print(body(brms:::dgen_extreme_value))
mu_k <- dr[k, "b_Intercept"] + dr[k, "b_x"] * d$x[i]
bd <- withCallingHandlers(
  brms:::dgen_extreme_value(rep(d$y[i], 3), mu = unname(mu_k),
                            sigma = unname(dr[k, "sigma"]),
                            xi = unname(dr[k, "xi"]), log = TRUE),
  warning = function(w) {
    cat("brms warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
  })
cat("brms dgen_extreme_value(log = TRUE) at the moved cells:", format(bd), "\n")
llb <- ll; llb[k, i] <- bd
cat("same matrix with brms's values at row", i, "(other off-support cells",
    "set finite): ")
llb[!is.finite(llb) & !is.nan(llb)] <- -50
cat("NaN", sum(is.nan(llb)), "\n")
cat("loo::loo.matrix on it:", msg(loo::loo(llb, r_eff = re)), "\n")
cat("loo::psis(-it):", msg(loo::psis(-llb, r_eff = re)), "\n")
# can brms's own draws leave a fitted row's support? its Stan program
sc <- brms::stancode(y ~ x, family = brms::gen_extreme_value(), data = d)
sl <- strsplit(as.character(sc), "\n")[[1]]
cat("brms stancode lines with xi:\n")
writeLines(paste("  ", grep("xi", sl, value = TRUE)))
cat("LOO-INF DONE\n")
