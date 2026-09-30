# Recovery at a known truth: 95% Wald interval coverage over replicates.
#   Rscript dev/fams2-coverage.R <family> <first seed> <last seed>
# One RDS per replicate in dev/fams2-cov/, so a count comes from the
# results on disk, never from what was launched.
args <- commandArgs(TRUE)
famname <- args[1]
seeds <- as.integer(args[2]):as.integer(args[3])
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-cov"
dir.create(out_dir, showWarnings = FALSE)
n <- 500

# truth on the scale confint() reports: the internal link-scale vector
designs <- list(
  xbeta = list(
    truth = c("(Intercept)" = 0.3, x = 0.6, "phi_(Intercept)" = log(6),
              "kappa_(Intercept)" = -2, kappa_x = 0.5),
    sim = function() {
      d <- data.frame(x = rnorm(n))
      mu <- plogis(0.3 + 0.6 * d$x)
      kap <- exp(-2 + 0.5 * d$x)
      z <- rbeta(n, mu * 6, (1 - mu) * 6)
      d$y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
      d
    },
    fit = function(d) frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d)),
  zibb = list(
    truth = c("(Intercept)" = -0.4, x = 0.5, "phi_(Intercept)" = log(4),
              "zi_(Intercept)" = -1, zi_x = 0.4),
    sim = function() {
      d <- data.frame(x = rnorm(n), tr = sample(4:20, n, TRUE))
      mu <- plogis(-0.4 + 0.5 * d$x)
      yb <- rbinom(n, d$tr, rbeta(n, mu * 4, (1 - mu) * 4))
      d$y <- ifelse(runif(n) < plogis(-1 + 0.4 * d$x), 0L, yb)
      d
    },
    fit = function(d) frm(bf(y | trials(tr) ~ x, zi ~ x),
                          family = zero_inflated_beta_binomial(), data = d)),
  hurdle_cumulative = list(
    truth = c(x = 0.8, "hu_(Intercept)" = -0.6, hu_z = 0.5,
              tau_raw_1 = -1, tau_raw_2 = log(1.3), tau_raw_3 = log(1.2)),
    sim = function() {
      d <- data.frame(x = rnorm(n), z = rnorm(n))
      u <- rlogis(n) + 0.8 * d$x
      d$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * d$z), 0L,
                    1L + (u > -1) + (u > 0.3) + (u > 1.5))
      d
    },
    fit = function(d) frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(),
                          data = d))
)
des <- designs[[famname]]
for (s in seeds) {
  set.seed(s)
  d <- des$sim()
  res <- tryCatch({
    warn <- character()
    f <- withCallingHandlers(des$fit(d), warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    ci <- confint(f)
    ci <- ci[names(des$truth), , drop = FALSE]
    list(seed = s, family = famname, est = ci[, "est"], lwr = ci[, "lwr"],
         upr = ci[, "upr"], truth = des$truth, warn = warn,
         cover = ci[, "lwr"] <= des$truth & des$truth <= ci[, "upr"])
  }, error = function(e) list(seed = s, family = famname,
                              error = conditionMessage(e)))
  saveRDS(res, file.path(out_dir, sprintf("%s-%04d.rds", famname, s)))
}
