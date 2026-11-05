# Short check of fix round 1, item 2: laplace = TRUE as the low-noise route, and
# check_laplace() beside it. What can a user read from laplace = TRUE at low noise
# (the field, predictions, loo)? Which nugget does check_laplace()'s tape carry,
# and does it report the fit's Laplace correctly there? Non-gaussian routes through
# the factor at 1e-12. Seed 6601, n 40 (the lane's design).
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-lap.R lap.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(6601)
n <- 40
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + 1e-3 * e0)
arg <- list(chains = 2, iter = 600, warmup = 300, seed = 3, refresh = 0, cores = 1)
t0 <- proc.time()[["elapsed"]]
fl <- suppressWarnings(do.call(frm_sample, c(list(bf(y ~ gp(x)), data = d, family = gaussian(),
                                                  laplace = TRUE), arg)))
cat(sprintf("laplace fit, n 40, 2 chains of 600: %.0f s\n", proc.time()[["elapsed"]] - t0))
dd <- posterior::as_draws_df(fl)
cat("columns:", paste(utils::head(posterior::variables(dd), 8), collapse = " "), "\n")
cat(sprintf("sigma median %.3e (target 1.07e-3)\n", median(dd$sigma)))
nd <- data.frame(x = c(0.5, 2.2, 5.1))
try_read <- function(tag, expr) {
  r <- tryCatch({
    v <- suppressWarnings(expr)
    paste(utils::capture.output(str(v, give.attr = FALSE, vec.len = 2))[1], collapse = "")
  }, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 170)))
  cat(sprintf("READ %-34s %s\n", tag, r))
}
try_read("posterior_epred()", brms::posterior_epred(fl))
try_read("posterior_epred(newdata)", brms::posterior_epred(fl, newdata = nd))
try_read("posterior_epred(re_formula = NA)", brms::posterior_epred(fl, re_formula = NA))
try_read("posterior_linpred()", brms::posterior_linpred(fl))
try_read("posterior_predict()", brms::posterior_predict(fl))
try_read("fitted()", fitted(fl))
try_read("predict(newdata)", predict(fl, newdata = nd))
try_read("log_lik()", brms::log_lik(fl))
try_read("loo()", brms::loo(fl))
try_read("frm_curve", frm_curve(fl, "x"))
try_read("conditional_effects()", brms::conditional_effects(fl))
try_read("ranef()", ranef(fl))
try_read("summary()", summary(fl))

# check_laplace(): which nugget does its tape read?
f0 <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
d2 <- data.frame(x = x, y = sin(2 * x) + 0.1 * e0)
f1 <- frm(bf(y ~ gp(x)), family = gaussian(), data = d2)
ns <- asNamespace("frmtmb.sample")
for (lab in c("ordinary noise", "low noise")) {
  fit <- if (lab == "low noise") f0 else f1
  t1 <- proc.time()[["elapsed"]]
  cl <- tryCatch(suppressWarnings(check_laplace(fit, chains = 2, iter = 800, refresh = 0,
                                                seed = 5)),
                 error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("check_laplace, %s, %.0f s:\n", lab, proc.time()[["elapsed"]] - t1))
  if (is.character(cl)) {
    cat(cl, "\n")
  } else {
    print(utils::head(as.data.frame(cl), 6), digits = 3)
  }
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 40, warmup = 20, refresh = 0, seed = 5,
                                    .diagnostic = TRUE))
  bk <- Filter(function(b) identical(b$covstruct, "gp"), ds$fit$frame$re_blocks)[[1]]
  cat("  nugget on the diagnostic run's block:", format(bk$gp_nugget %||% "NULL (frm's 1e-6)"),
      "\n")
  cat(sprintf("  frm() mode: sigma %.3e sdgp %.3f\n", as.numeric(sigma(fit)),
              exp(fit$estimates$theta[[1]])))
}

# laplace = TRUE through the factor, off the gaussian path
d3 <- data.frame(x = x, z = round(runif(n, 0, 6), 3), g = factor(rep(1:5, length.out = n)))
d3$cnt <- rpois(n, exp(2 + sin(2 * x)))
d3$bin <- rbinom(n, 1, plogis(sin(2 * x)))
d3$yy <- sin(2 * x) + 0.2 * cos(d3$z) + rnorm(n, 0, 0.1)
run <- function(tag, expr) {
  t1 <- proc.time()[["elapsed"]]
  r <- tryCatch({
    v <- suppressWarnings(expr)
    dv <- posterior::as_draws_df(v)
    pn <- grep("^sdgp", posterior::variables(dv), value = TRUE)[1]
    m <- posterior::extract_variable_matrix(dv, pn)
    sprintf("%s median %.3f rhat %.3f", pn, median(as.numeric(m)), posterior::rhat(m))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 140)))
  cat(sprintf("LAP %-26s %s, %.0f s\n", tag, r, proc.time()[["elapsed"]] - t1))
}
a2 <- list(chains = 2, iter = 400, warmup = 200, seed = 3, refresh = 0, cores = 1, laplace = TRUE)
run("poisson", do.call(frm_sample, c(list(bf(cnt ~ gp(x)), data = d3, family = poisson()), a2)))
run("bernoulli", do.call(frm_sample, c(list(bf(bin ~ gp(x)), data = d3, family = bernoulli()), a2)))
run("student", do.call(frm_sample, c(list(bf(yy ~ gp(x)), data = d3, family = student()), a2)))
run("two gp terms", do.call(frm_sample, c(list(bf(yy ~ gp(x) + gp(z)), data = d3,
                                                family = gaussian()), a2)))
run("gp + (1|g)", do.call(frm_sample, c(list(bf(yy ~ gp(x) + (1 | g)), data = d3,
                                              family = gaussian()), a2)))
run("hsgp k = 15 (no factor)", do.call(frm_sample, c(list(bf(yy ~ gp(x, k = 15)), data = d3,
                                                           family = gaussian()), a2)))
