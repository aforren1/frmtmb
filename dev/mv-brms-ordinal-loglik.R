# Lane mv, 2026-09-25: an ordinal response in a multivariate model,
# against brms 2.23.0 itself, the same way as
# dev/mv-brms-student-loglik.R: a fixed_param run from an init equal to
# a frmtmb parameter point, then brms::log_lik() summed against
# frmtmb's objective there. brms centers the design, so its threshold
# vector is Intercept_o[k] = tau_k - b * mean(x) (eta enters as
# tau - eta) and the gaussian Intercept_y1 is b0 + b1 * mean(x).
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-mv")
.libPaths(c("/opt/rlib/stan", lib, "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
set.seed(20260925)
n <- 150
x <- rnorm(n)
dd <- data.frame(x = x,
                 o = cut(0.9 * x + rlogis(n), c(-Inf, -1, 0, 1.2, Inf),
                         labels = FALSE),
                 y1 = 1 + 0.5 * x + rnorm(n))
fit <- frm(bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian(), data = dd)
set.seed(3)
p <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.05)
pl <- fit$obj$env$parList(p)
bn <- names(fit$frame$par_template$beta)
b <- function(nm) pl$beta[match(nm, bn)]
raw <- pl[["o_tau_raw"]]
tau <- cumsum(c(raw[1], exp(raw[-1])))
mx <- mean(x)
init <- list(
  b_o = array(b("o_x")), Intercept_o = tau - b("o_x") * mx,
  b_y1 = array(b("y1_x")), Intercept_y1 = b("y1_(Intercept)") +
    b("y1_x") * mx,
  sigma_y1 = exp(pl$betad[1]))
suppressMessages(library(brms))
bfit <- brm(bf(o ~ x, family = cumulative()) +
              bf(y1 ~ x, family = gaussian()) + set_rescor(FALSE),
            data = dd, chains = 1, iter = 1, warmup = 0,
            algorithm = "fixed_param", init = list(init), refresh = 0,
            seed = 1)
ll_brms <- sum(log_lik(bfit))
ll_frm <- -fit$obj$fn(p)
cat(sprintf("brms log_lik summed   %.12f\n", ll_brms))
cat(sprintf("frmtmb objective      %.12f\n", ll_frm))
cat(sprintf("relative difference   %.3e\n",
            abs(ll_brms - ll_frm) / abs(ll_brms)))
dr <- as_draws_df(bfit)
cat("brms b_o_Intercept:", unlist(dr[1, grep("^b_o_Intercept", names(dr))]),
    "\nfrmtmb thresholds: ", tau, "\n")
