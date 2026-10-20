# Reviewer: logistic_normal() in mvbf(). (1) The log-likelihood at the
# fit's estimates, written out by hand from brms's density
# (fun_logistic_normal.stan: multi_normal_cholesky of the log ratios
# minus sum(log(y))) plus the gaussian; (2) the code_fix of
# test-refcat-brms.R edits only the generated quantities and a
# perturbed fit still fails; (3) names for a response with an
# underscore and a cbind() response. Lane library.
# Log dev/catmulti-rev-out/mvln.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-catmulti/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(46)
n <- 300
dd <- data.frame(x = stats::rnorm(n))
W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1], -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
Y <- exp(E) / rowSums(exp(E)); colnames(Y) <- c("a", "b", "c")
dd$Y <- Y
dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
ctl <- frmtmb_control(grad_tol = 1e-6, restarts = 3)
fit <- frm(mvbf(bf(Y ~ x, family = logistic_normal()), bf(w ~ x, family = gaussian())),
           data = dd, control = ctl)
cat("logLik(fit):", sprintf("%.10f", as.numeric(logLik(fit))), "\n")
lp <- function(d, r, type = "link") as.numeric(frm_linpred(fit, resp = r, dpar = d, type = type))
mub <- lp("mub", "Y"); muc <- lp("muc", "Y")
sb <- lp("sigmab", "Y", "response"); sc <- lp("sigmac", "Y", "response")
rho <- summary(fit)$spec_pars["lncor_Y(b,c)", "Estimate"]
v <- variables(fit)
cat("variables with lncor:", grep("lncor", v, value = TRUE), "\n")
z1 <- log(Y[, "b"] / Y[, "a"]); z2 <- log(Y[, "c"] / Y[, "a"])
# bivariate normal by hand
u1 <- (z1 - mub) / sb; u2 <- (z2 - muc) / sc
ln <- -log(2 * pi) - log(sb) - log(sc) - 0.5 * log(1 - rho^2) -
  (u1^2 - 2 * rho * u1 * u2 + u2^2) / (2 * (1 - rho^2)) - rowSums(log(Y))
muw <- lp("mu", "w"); sw <- lp("sigma", "w", "response")
lw <- dnorm(dd$w, muw, sw, log = TRUE)
hand <- sum(ln) + sum(lw)
cat(sprintf("hand-written logLik %.10f; difference %.3e (relative %.3e)\n", hand,
            hand - as.numeric(logLik(fit)), (hand - as.numeric(logLik(fit))) / abs(hand)))
# rho from summary is printed to its stored precision; take it from the chart too
th <- fit$estimates[[grep("lncor_raw", names(fit$estimates), value = TRUE)[1]]]
cat("chart name:", grep("lncor_raw", names(fit$estimates), value = TRUE), "theta", th,
    "rho from chart", th / sqrt(1 + th^2), "summary rho", rho, "\n")

## (2) the code_fix
env <- testthat::test_env("frmtmb")
for (f in list.files("tests/testthat", "^helper-.*[.]R$", full.names = TRUE)) sys.source(f, env)
fix <- function(code) {
  out <- gsub("multiply_lower_tri_self_transpose(Llncor);",
              "multiply_lower_tri_self_transpose(Llncor_Y);", code, fixed = TRUE)
  gsub("[nlncor] lncor;", "[nlncor_Y] lncor;", out, fixed = TRUE)
}
bform <- brms::bf(Y ~ x, family = brms::logistic_normal()) +
  brms::bf(w ~ x, family = gaussian()) + brms::set_rescor(FALSE)
code <- brms::make_stancode(bform, data = dd)
a <- strsplit(as.character(code), "\n")[[1]]; b <- strsplit(fix(code), "\n")[[1]]
cat("code lines", length(a), "after fix", length(b), "; changed lines:\n")
ch <- which(a != b)
for (i in ch) cat(sprintf("  %d: %s\n  -> %s\n", i, trimws(a[i]), trimws(b[i])))
gq <- grep("^generated quantities", a)
cat("generated quantities starts at line", gq, "; every change after it:", all(ch > gq), "\n")
r <- tryCatch(env$brms_lp_check(bform, list(brms::logistic_normal(), gaussian()), dd, fit,
                                code_fix = fix), error = function(e) conditionMessage(e))
if (is.character(r)) cat("lp check on the fit: ERROR", r, "\n") else
  cat(sprintf("lp check on the fit: const %.3e grad %.3e\n", r$measured_const, r$max_grad))
# perturb the reported correlation chart by 0.05: brms's density there
# must disagree with frmtmb's logLik at its optimum
fp <- fit
cn <- grep("lncor_raw", names(fp$estimates), value = TRUE)[1]
fp$estimates[[cn]] <- fp$estimates[[cn]] + 0.05
rr <- tryCatch({
  env$brms_lp_check(bform, list(brms::logistic_normal(), gaussian()), dd, fp, code_fix = fix)
  "PASSED (it must not)"
}, error = function(e) paste("FAILED, as it must:", conditionMessage(e)))
cat("perturbed lncor chart:", substr(rr, 1, 300), "\n")
cat("estimate blocks:", names(fit$estimates), "\n")

## (3) names
dd2 <- dd; dd2$Y_a <- dd$Y
f2 <- tryCatch(frm(mvbf(bf(Y_a ~ x, family = logistic_normal()), bf(w ~ x)), data = dd2),
               error = function(e) conditionMessage(e))
if (is.character(f2)) cat("Y_a: ERROR", f2, "\n") else {
  cat("Y_a variables:", grep("lncor|sigma|^b_", variables(f2), value = TRUE), "\n")
  cat("Y_a summary spec rows:", rownames(summary(f2)$spec_pars), "\n")
  dp <- default_prior(f2); print(dp[dp$class %in% c("lncor", "sigma"), c("class", "resp", "dpar")])
}
dd3 <- dd; dd3$a <- Y[, 1]; dd3$b <- Y[, 2]; dd3$c <- Y[, 3]
f3 <- tryCatch(frm(mvbf(bf(cbind(a, b, c) ~ x, family = logistic_normal()), bf(w ~ x)), data = dd3),
               error = function(e) conditionMessage(e))
if (is.character(f3)) cat("cbind: ERROR", f3, "\n") else {
  cat("cbind responses:", names(f3$spec$responses), "\n")
  cat("cbind variables:", grep("lncor|sigma", variables(f3), value = TRUE), "\n")
  dp <- default_prior(f3); print(dp[dp$class %in% c("lncor"), c("class", "resp")])
}
cat("brms names for Y_a (stancode, no compile):\n")
sc <- brms::make_stancode(brms::bf(Y_a ~ x, family = brms::logistic_normal()) +
                            brms::bf(w ~ x) + brms::set_rescor(FALSE), data = dd2)
cat(grep("lncor|sigma", strsplit(sc, "\n")[[1]], value = TRUE)[1:6], sep = "\n")
dp <- brms::default_prior(brms::bf(Y_a ~ x, family = brms::logistic_normal()) +
                            brms::bf(w ~ x) + brms::set_rescor(FALSE), data = dd2)
print(dp[dp$class %in% c("lncor"), c("prior", "class", "resp")])
cat("DONE\n")
