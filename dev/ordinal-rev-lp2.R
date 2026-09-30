# Reviewer, lane ordinal: the rows of dev/ordinal-rev-lp.R the helper's
# translator could not take, built by hand: multivariate ordinal models,
# and the hurdle on data with zeros and an unused middle category. Also
# the worker's cratio cloglog NaN on its own data (seed 20260930).
# Data seed 20261001 (as dev/ordinal-rev-lp.R). Output:
# dev/ordinal-rev-log-lp2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "ordrev-stan-cache"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "tests/testthat/helper-brms.R"), envir = env)
ctl <- frmtmb_control(grad_tol = 1e-7, restarts = 3)
set.seed(20261001)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                h = factor(sample(c("p", "q", "r"), n, TRUE)))
ug <- rnorm(8, 0, 0.3)[as.integer(d$g)]
u <- stats::rlogis(n) / exp(0.3 * d$z + ug) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
u2 <- stats::rlogis(n) + 0.5 * d$x
d$y2 <- 1L + (u2 > -1) + (u2 > 0) + (u2 > 1)
d$ym <- d$y
d$ym[d$ym == 3L] <- 2L
d$yhm <- ifelse(runif(n) < 0.25, 0L, d$ym)

brms_prog <- function(bform, family) {
  prior <- env$brms_flat_prior(bform, data = d, family = family)
  code <- brms::make_stancode(bform, data = d, family = family, prior = prior)
  sdat <- env$brms_standata(bform, data = d, family = family, prior = prior)
  sf <- suppressMessages(rstan::sampling(env$brms_stan_model(code),
                                         data = sdat, chains = 0))
  list(code = code, sdat = sdat, sf = sf,
       pars = env$brms_stan_par_names(code))
}
report <- function(lab, b, pars, fit) {
  up <- rstan::unconstrain_pars(b$sf, pars)
  lp <- rstan::log_prob(b$sf, up, adjust_transform = FALSE)
  g <- rstan::grad_log_prob(b$sf, up, adjust_transform = FALSE)
  ours <- as.numeric(logLik(fit))
  cat(sprintf("%s\n  LP const %.4g rel %.3g grad %.3g ours %.10g\n", lab,
              lp - ours, abs(lp - ours) / abs(ours), max(abs(g)), ours))
}
thr <- function(fit, resp) {
  fam <- family(fit)
  if (!is.null(resp)) fam <- fam[[resp]]
  comp <- if (is.null(resp)) "tau_raw" else paste0(resp, "_tau_raw")
  est <- fit$estimates[[comp]] %||% fit$estimates[["tau_raw"]][[resp]]
  frmtmb:::ord_threshold_values(fam, est)
}

cat("== mv 1: cumulative disc ~ 0 + z (y) + sratio sum_to_zero (y2) ==\n")
tryCatch({
  fit <- frm(bf(y ~ x, disc ~ 0 + z) + bf(y2 ~ x), data = d, control = ctl,
             family = list(cumulative(), sratio(threshold = "sum_to_zero")))
  cat("  estimates components:", names(fit$estimates), "\n")
  fe <- fixef(fit, flatten = TRUE); print(fe)
  b <- brms_prog(brms::bf(y ~ x, disc ~ 0 + z) + brms::bf(y2 ~ x) +
                   brms::set_rescor(FALSE),
                 list(brms::cumulative(), brms::sratio(threshold = "sum_to_zero")))
  cat("  brms parameters:", b$pars, "\n")
  t1 <- thr(fit, "y"); t2 <- thr(fit, "y2")
  cat("  sum of y2 thresholds / max:", sum(t2) / max(abs(t2)), "\n")
  pars <- list(b_y = array(fe[["y_x"]], 1),
               Intercept_y = t1 - mean(d$x) * fe[["y_x"]],
               b_disc_y = array(fe[["y_disc_z"]], 1),
               b_y2 = array(fe[["y2_x"]], 1),
               Intercept_y2 = t2)
  stopifnot(setequal(names(pars), b$pars))
  report("  identity", b, pars, fit)
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("\n== mv 2: acat probit (y) + cratio (y2), disc ~ 0 + z on both ==\n")
tryCatch({
  fit <- frm(bf(y ~ x, disc ~ 0 + z) + bf(y2 ~ x, disc ~ 0 + z), data = d,
             control = ctl, family = list(acat("probit"), cratio()))
  fe <- fixef(fit, flatten = TRUE); print(fe)
  b <- brms_prog(brms::bf(y ~ x, disc ~ 0 + z) + brms::bf(y2 ~ x, disc ~ 0 + z) +
                   brms::set_rescor(FALSE),
                 list(brms::acat("probit"), brms::cratio()))
  cat("  brms parameters:", b$pars, "\n")
  pars <- list(b_y = array(fe[["y_x"]], 1),
               Intercept_y = thr(fit, "y") - mean(d$x) * fe[["y_x"]],
               b_disc_y = array(fe[["y_disc_z"]], 1),
               b_y2 = array(fe[["y2_x"]], 1),
               Intercept_y2 = thr(fit, "y2") - mean(d$x) * fe[["y2_x"]],
               b_disc_y2 = array(fe[["y2_disc_z"]], 1))
  stopifnot(setequal(names(pars), b$pars))
  report("  identity", b, pars, fit)
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("\n== hurdle_cumulative equidistant, disc ~ 0 + z, zeros and an unused middle category ==\n")
cat("  table(yhm):", table(d$yhm), "\n")
tryCatch({
  fit <- frm(bf(yhm ~ x, disc ~ 0 + z), data = d, control = ctl,
             family = hurdle_cumulative(threshold = "equidistant"))
  b <- brms_prog(brms::bf(yhm ~ x, disc ~ 0 + z),
                 brms::hurdle_cumulative(threshold = "equidistant"))
  cat("  brms nthres:", b$sdat$nthres, " frmtmb:", family(fit)$thres$nthres, "\n")
  pars <- env$stan_pars_from_fit(fit, b$sdat, b$code)
  report("  identity", b, pars, fit)
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("\n== the worker's cratio cloglog NaN, on its data (seed 20260930) ==\n")
set.seed(20260930)
n <- 300
dw <- data.frame(x = rnorm(n), z = rnorm(n),
                 g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.4 * dw$z) + 0.8 * dw$x
dw$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
fit <- frm(bf(y ~ x, disc ~ 0 + z), family = cratio("cloglog"), data = dw,
           control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
eta <- dw$x * fixef(fit)["x", "Estimate"]
disc <- exp(dw$z * fixef(fit)["disc_z", "Estimate"])
tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
th <- disc * outer(eta, tau, "-")
used <- sweep(matrix(seq_along(tau), n, length(tau), byrow = TRUE), 1,
              pmin(dw$y, length(tau)), "<=")
cat(sprintf("  rows reading a threshold with th > 6.6: %d; max th read %.2f\n",
            sum(rowSums((th > 6.6) & used) > 0), max(th[used])))
