# Reviewer, lane ordinal: the five brms 2.23.0 defects the lane records,
# checked independently, on different data (seed 20261007).
# Output: dev/ordinal-rev-log-defects.txt
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
set.seed(20261007)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 1.2 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$y2 <- 1L + (stats::rlogis(n) + d$x > 0) + (stats::rlogis(n) + d$x > 1)
stanc_ok <- function(code) tryCatch({rstan::stanc(model_code = code); "compiles"},
  error = function(e) paste("stanc ERROR:", substr(gsub("\\s+", " ",
                                                      conditionMessage(e)), 1, 400)))

cat("== 1. cratio cloglog with disc ~ 0 + z: brms gradient ==\n")
fit <- frm(bf(y ~ x, disc ~ 0 + z), family = cratio("cloglog"), data = d,
           control = frmtmb_control(grad_tol = 1e-7, restarts = 3))
bform <- brms::bf(y ~ x, disc ~ 0 + z); bfam <- brms::cratio("cloglog")
prior <- env$brms_flat_prior(bform, data = d, family = bfam)
code <- brms::make_stancode(bform, data = d, family = bfam, prior = prior)
sdat <- env$brms_standata(bform, data = d, family = bfam, prior = prior)
sf <- suppressMessages(rstan::sampling(env$brms_stan_model(code), data = sdat,
                                       chains = 0))
pars <- env$stan_pars_from_fit(fit, sdat, code)
up <- rstan::unconstrain_pars(sf, pars)
lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
cat(sprintf("frmtmb logLik %.10f, brms log_prob %.10f, rel %.3g\n",
            as.numeric(logLik(fit)), lp,
            abs(lp - as.numeric(logLik(fit))) / abs(lp)))
cat("brms gradient:", format(signif(g, 3)), "\n")
cat("frmtmb max|grad|:", max(abs(fit$obj$gr(fit$opt$par))), "\n")
# the mechanism: th = disc (mu - thres_k); exp(-exp(th)) underflows to
# 0 above about 6.6, so q_k = log1m_exp(-exp(th)) is exactly 0 and
# p_k = log1m_exp(q_k) is -Inf, for k < y too, where it is not returned
# but its derivative still enters
eta <- d$x * fixef(fit)["x", "Estimate"]
disc <- exp(d$z * fixef(fit)["disc_z", "Estimate"])
tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
th <- disc * outer(eta, tau, "-")
used <- sweep(matrix(seq_along(tau), n, length(tau), byrow = TRUE), 1,
              pmin(d$y, length(tau)), "<=")
cat(sprintf("rows with th > 6.6 on a threshold the row reads: %d (max th %.2f)\n",
            sum(rowSums((th > 6.6) & used) > 0), max(th[used])))
# the same point without those rows: is the gradient finite?
keep <- rowSums((th > 6.6) & used) == 0
if (any(!keep)) {
  d2 <- d[keep, ]
  sdat2 <- env$brms_standata(bform, data = d2, family = bfam, prior = prior)
  sf2 <- suppressMessages(rstan::sampling(env$brms_stan_model(code),
                                          data = sdat2, chains = 0))
  g2 <- rstan::grad_log_prob(sf2, rstan::unconstrain_pars(sf2, pars),
                             adjust_transform = FALSE)
  cat("brms gradient with those rows dropped:", format(signif(g2, 3)), "\n")
}

cat("\n== 2. cumulative logit sum_to_zero, disc at 1: stanc ==\n")
sc <- brms::make_stancode(brms::bf(y ~ x), data = d,
                          family = brms::cumulative(threshold = "sum_to_zero"))
cat(grep("target \\+=", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
cat(stanc_ok(sc), "\n")
sc <- brms::make_stancode(brms::bf(y ~ 1), data = d,
                          family = brms::cumulative(threshold = "sum_to_zero"))
cat("intercept-only:", stanc_ok(sc), "\n")

cat("\n== 3. acat softit: stanc ==\n")
sc <- brms::make_stancode(brms::bf(y ~ x), data = d, family = brms::acat("softit"))
cat(grep("softit", strsplit(sc, "\n")[[1]], value = TRUE)[1:6], sep = "\n")
cat(stanc_ok(sc), "\n")
sc <- brms::make_stancode(brms::bf(y ~ x), data = d,
                          family = brms::cumulative("softit"))
cat("cumulative softit:", stanc_ok(sc), "\n")

cat("\n== 4. equidistant thresholds in a multivariate model and a mixture ==\n")
sc <- brms::make_stancode(brms::bf(y ~ x) + brms::bf(y2 ~ x) +
                            brms::set_rescor(FALSE), data = d,
                          family = brms::cumulative(threshold = "equidistant"))
L <- strsplit(sc, "\n")[[1]]
cat(grep("delta", L, value = TRUE), sep = "\n")
cat(stanc_ok(sc), "\n")
sc <- brms::make_stancode(brms::bf(y ~ x) + brms::bf(y2 ~ x) +
                            brms::set_rescor(FALSE), data = d,
                          family = list(brms::cumulative(threshold = "equidistant"),
                                        brms::sratio()))
cat("one equidistant response only:", stanc_ok(sc), "\n")
sc <- tryCatch(brms::make_stancode(brms::bf(y ~ x), data = d,
  family = brms::mixture(brms::cumulative(threshold = "equidistant"),
                         brms::cumulative(threshold = "equidistant"))),
  error = function(e) paste("brms ERROR:", conditionMessage(e)))
cat("mixture:", if (startsWith(sc, "brms ERROR")) sc else stanc_ok(sc), "\n")

cat("\n== 5. probit_approx on brms's R side ==\n")
cat("brms:::inv_link(0.7, 'probit_approx') =", brms:::inv_link(0.7, "probit_approx"),
    " pnorm(0.7) =", pnorm(0.7), " Phi_approx(0.7) =",
    plogis(0.07056 * 0.7^3 + 1.5976 * 0.7), "\n")
sc <- brms::make_stancode(brms::bf(y ~ x), data = d,
                          family = brms::acat("probit_approx"))
cat(grep("Phi_approx", strsplit(sc, "\n")[[1]], value = TRUE)[1], "\n")
