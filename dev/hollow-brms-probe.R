# Lane hollow: what brms 2.23.0 computes for ledger rows
# brmsfit-methods:394 (Evid.Ratio of a point hypothesis) and :469
# (loo() with newdata). Run with the user library only:
#   Rscript dev/hollow-brms-probe.R > dev/hollow-out/brms-probe.txt 2>&1
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "loo",
    format(packageVersion("loo")), "\n")

fit1 <- brms:::rename_pars(brms:::brmsfit_example1)
cat("\n== Row 394: fixture 1 ==\n")
print(fit1$formula)
cat("sample_prior setting:", deparse(attr(fit1$prior, "sample_prior")),
    "\n")
print(fit1$prior)
pd <- prior_draws(fit1)
cat("prior_draws columns:", paste(names(pd), collapse = " "), "\n")
cat("prior_draws rows:", nrow(pd), "\n")
vars <- variables(fit1)
cat("variables with prior_:", paste(grep("^prior_", vars, value = TRUE),
                                    collapse = " "), "\n")
cat("ndraws:", ndraws(fit1), "\n")

set.seed(394)
hyp <- hypothesis(fit1, "Intercept = 0", class = "sd", group = "visit")
print(hyp$hypothesis, digits = 10)
cat("prior_samples H1 summary:\n")
print(summary(hyp$prior_samples$H1))
cat("posterior samples H1 summary:\n")
print(summary(hyp$samples$H1))

# Rebuild the ratio by hand from the draws brms returned, through
# brms's own density_ratio(), to say WHAT it is.
dr <- brms:::density_ratio(hyp$samples$H1, hyp$prior_samples$H1,
                            point = 0)
cat("density_ratio(post, prior, 0) =", format(dr, digits = 10), "\n")
cat("posterior kde at 0 =",
    format(brms:::density_ratio(hyp$samples$H1, point = 0), digits = 10),
    "\n")
cat("prior kde at 0 =",
    format(brms:::density_ratio(hyp$prior_samples$H1, point = 0),
           digits = 10), "\n")

cat("\n-- the same, directional --\n")
hyp_d <- hypothesis(fit1, "Intercept > 0", class = "sd", group = "visit")
print(hyp_d$hypothesis, digits = 10)
cat("direct posterior odds:",
    format(mean(hyp_d$samples$H1 > 0) / mean(hyp_d$samples$H1 <= 0),
           digits = 10), "\n")

cat("\n-- a coefficient (class b) --\n")
hb <- hypothesis(fit1, c("Trt1 = 0", "Trt1 > 0"))
print(hb$hypothesis, digits = 10)

cat("\n-- a fit WITHOUT prior draws: fixture 2 --\n")
fit2 <- brms:::rename_pars(brms:::brmsfit_example2)
cat("fixture 2 prior_draws ncol:", NCOL(prior_draws(fit2)), "\n")
h2 <- hypothesis(fit2, names(fixef(fit2)[, 1])[2] |>
                   (\(n) c(paste(n, "= 0"), paste(n, "> 0")))())
print(h2$hypothesis, digits = 10)

cat("\n== Row 469: fixture 6 ==\n")
fit6 <- brms:::rename_pars(brms:::brmsfit_example6)
print(fit6$formula)
cat("class(fit6$data):", class(fit6$data), "\n")
cat("dim(fit6$data):", dim(fit6$data), "\n")
cat("names(fit6$data):", names(fit6$data), "\n")
cat("data.name:", fit6$data.name, "\n")
str(fit6$data)
l1 <- suppressWarnings(loo(fit6, cores = 1))
l2 <- suppressWarnings(loo(fit6, cores = 1, newdata = fit6$data))
cat("loo(fit6) estimates:\n"); print(l1$estimates, digits = 12)
cat("loo(fit6, newdata = fit6$data) estimates:\n")
print(l2$estimates, digits = 12)
cat("identical pointwise:", identical(l1$pointwise, l2$pointwise), "\n")
cat("max abs diff pointwise:",
    max(abs(l1$pointwise - l2$pointwise)), "\n")
ll1 <- log_lik(fit6)
ll2 <- log_lik(fit6, newdata = fit6$data)
cat("log_lik identical:", identical(ll1, ll2), " max abs diff:",
    max(abs(ll1 - ll2)), "\n")
print(loo_compare(l1, l2))

cat("\n-- newdata = a subset of rows (first 20) --\n")
nd <- fit6$data[1:20, , drop = FALSE]
l3 <- suppressWarnings(loo(fit6, cores = 1, newdata = nd))
print(l3$estimates, digits = 12)
cat("pointwise rows:", nrow(l3$pointwise), "\n")
cat("pointwise elpd equal to loo(fit6)[1:20]:",
    isTRUE(all.equal(l3$pointwise[, "elpd_loo"],
                     l1$pointwise[1:20, "elpd_loo"])), "\n")
cat("max abs diff vs loo(fit6) rows 1:20:",
    max(abs(l3$pointwise[, "elpd_loo"] - l1$pointwise[1:20, "elpd_loo"])),
    "\n")

cat("\n-- newdata = changed response (y + 1) --\n")
yname <- all.vars(fit6$formula$formula)[1]
cat("response column guess:", yname, "\n")
nd2 <- fit6$data
resp <- names(fit6$data)[1]
cat("first data column:", resp, "\n")
DONE <- TRUE
cat("\nPROBE END\n")
