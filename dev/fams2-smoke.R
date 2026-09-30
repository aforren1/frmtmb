# Smoke test: xbeta() and zero_inflated_beta_binomial() reach the post-fit
# methods. Lane build first in the library path.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
try_ <- function(lab, expr) {
  r <- tryCatch({
    v <- withCallingHandlers(expr, warning = function(w) {
      cat("  [warn]", lab, ":", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
    "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", lab, r))
}

set.seed(1)
n <- 600
d <- data.frame(x = rnorm(n), g = gl(20, n / 20))
mu <- plogis(0.3 + 0.6 * d$x)
phi <- 6
kap <- 0.15
z <- rbeta(n, mu * phi, (1 - mu) * phi)
d$y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
cat("xbeta data: zeros", sum(d$y == 0), "ones", sum(d$y == 1), "\n")
fx <- frm(bf(y ~ x), family = xbeta(), data = d)
print(summary(fx))
try_("fitted", f1 <- fitted(fx))
try_("predict", predict(fx))
try_("simulate", simulate(fx, nsim = 2))
try_("resid pearson", residuals(fx, type = "pearson"))
try_("resid osa", residuals(fx, type = "osa"))
try_("conditional_effects", conditional_effects(fx))
try_("emmeans", emmeans::emmeans(fx, ~ x))
try_("default_prior", print(default_prior(bf(y ~ x), family = xbeta(),
                                          data = d)))
try_("kappa ~ x", frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d))
try_("(1|g)", frm(bf(y ~ x + (1 | g)), family = xbeta(), data = d))
try_("mixture", frm(bf(y ~ 1), family = mixture(xbeta, xbeta), data = d))
try_("compat", { cx <- frm_compat("xbeta"); print(table(cx$status)) })

set.seed(2)
d$tr <- 12L
mu <- plogis(-0.4 + 0.5 * d$x)
yb <- rbinom(n, d$tr, rbeta(n, mu * 4, (1 - mu) * 4))
d$yb <- ifelse(runif(n) < 0.25, 0L, yb)
fz <- frm(bf(yb | trials(tr) ~ x), family = zero_inflated_beta_binomial(),
          data = d)
print(summary(fz))
try_("fitted", fitted(fz))
try_("predict", predict(fz))
try_("simulate", simulate(fz, nsim = 2))
try_("resid pearson", residuals(fz, type = "pearson"))
try_("resid osa", residuals(fz, type = "osa"))
try_("conditional_effects", conditional_effects(fz))
try_("emmeans", emmeans::emmeans(fz, ~ x))
try_("zi ~ x", frm(bf(yb | trials(tr) ~ x, zi ~ x),
                   family = zero_inflated_beta_binomial(), data = d))
try_("mixture", frm(bf(yb | trials(tr) ~ 1),
                    family = mixture(zero_inflated_beta_binomial,
                                     beta_binomial), data = d))
try_("default_prior", print(default_prior(bf(yb | trials(tr) ~ x),
                    family = zero_inflated_beta_binomial(), data = d)))
fb <- frm(bf(yb | trials(tr) ~ x), family = beta_binomial(), data = d)
try_("bb resid pearson", residuals(fb, type = "pearson"))
