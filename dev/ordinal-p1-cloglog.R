# m3: brms's cratio("cloglog") gradient is NaN when a row has
# th = disc * (mu - thres_k) above about 6.6 on a threshold below its
# observed category (exp(-exp(th)) underflows, q_k is 0 and
# log1m_exp(0) is -Inf). Find a data seed of ord_disc_data() without
# such rows, and run the identity check there. Output:
# dev/ordinal-p1-log-cloglog.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordinal/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-brms.R", envir = env)
options(frmtmb.brms_lp_report = TRUE)
ord_disc_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  d$yh <- ifelse(runif(n) < stats::plogis(-1 + 0.3 * d$x), 0L, d$y)
  d
}
ctl <- frmtmb_control(grad_tol = 1e-6, restarts = 3)
for (seed in 20260930 + 0:6) {
  d <- ord_disc_data(seed)
  fit <- frm(bf(y ~ x, disc ~ 0 + z), family = cratio("cloglog"), data = d,
             control = ctl)
  eta <- d$x * fixef(fit)["x", "Estimate"]
  disc <- exp(d$z * fixef(fit)["disc_z", "Estimate"])
  tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
  th <- disc * (outer(eta, tau, "-"))
  below <- outer(d$y, seq_along(tau), ">")
  worst <- max(th[below])
  cat(sprintf("seed %d: largest th below the observed category %.3f\n",
              seed, worst))
  if (worst < 6) {
    r <- env$brms_lp_check(brms::bf(y ~ x, disc ~ 0 + z),
                           brms::cratio("cloglog"), d, fit)
    cat("identity at seed", seed, ": const", r$measured_const,
        " max grad", r$max_grad, "\n")
    break
  }
}
