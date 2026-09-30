# Two shapes where brms 2.23.0's own Stan program fails, measured so the
# likelihood rows can leave them out with a reason. Seed 20260930, the
# data of dev/ordinal-lpcheck.R. Output: dev/ordinal-log-brms-defects.txt
#
# 1. cratio("cloglog") with a modeled disc: the gradient of log_prob.
#    The value is compared as in brms_lp_check(), then the gradient is
#    taken at frmtmb's optimum, at the same point with disc's slope set
#    to 0, and with x's slope moved, to see whether the NaN follows disc.
#    frmtmb's own gradient at its optimum is printed beside it.
# 2. cumulative(threshold = "sum_to_zero") under the logit link with
#    disc at 1: the generated likelihood line and the compile result.
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordinal/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-brms.R", envir = env)
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)

cat("== 1. cratio cloglog, disc ~ 0 + z ==\n")
fit <- frm(bf(y ~ x, disc ~ 0 + z), family = cratio("cloglog"), data = d,
           control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
cat("frmtmb max |gradient| at its optimum:",
    format(max(abs(fit$obj$gr(fit$opt$par))), digits = 3), "\n")
bform <- brms::bf(y ~ x, disc ~ 0 + z)
bfam <- brms::cratio("cloglog")
prior <- env$brms_flat_prior(bform, data = d, family = bfam)
code <- brms::make_stancode(bform, data = d, family = bfam, prior = prior)
sdat <- env$brms_standata(bform, data = d, family = bfam, prior = prior)
sf <- suppressMessages(rstan::sampling(env$brms_stan_model(code),
                                       data = sdat, chains = 0))
pars <- env$stan_pars_from_fit(fit, sdat, code)
at <- function(p, lab) {
  up <- rstan::unconstrain_pars(sf, p)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  cat(sprintf("%-34s log_prob %.10f  gradient: %s\n", lab, lp,
              paste(format(g, digits = 3), collapse = " ")))
}
cat("frmtmb logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")
at(pars, "frmtmb's optimum")
p2 <- pars
p2$b_disc <- p2$b_disc * 0
at(p2, "disc slope 0 (disc = 1)")
p3 <- pars
p3$b <- p3$b * 0.5
at(p3, "x slope halved")
cat("brms's cratio-cloglog q_k line:\n")
cat(grep("log1m_exp|q\\[k\\] =", strsplit(code, "\n")[[1]], value = TRUE),
    sep = "\n")

cat("\n== 2. cumulative logit sum_to_zero, disc at 1 ==\n")
sc <- brms::make_stancode(brms::bf(y ~ x), data = d,
                          family = brms::cumulative(threshold = "sum_to_zero"))
cat(grep("target \\+=", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
r <- tryCatch({
  rstan::stanc(model_code = sc)
  "compiles"
}, error = function(e) paste("stanc ERROR:", conditionMessage(e)))
cat(substr(r, 1, 300), "\n")
sc2 <- brms::make_stancode(brms::bf(y ~ x), data = d,
                           family = brms::cumulative("probit",
                                                     threshold = "sum_to_zero"))
cat("probit link:\n")
cat(grep("target \\+=", strsplit(sc2, "\n")[[1]], value = TRUE), sep = "\n")
