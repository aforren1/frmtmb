## Reviewer, claim 3: name collisions the `b_` argument does not reach.
##   Rscript dev/brmsnames-rev-collide.R base|lane       data seed 12
## For each model: duplicated names in the draws labels and in the
## hypothesis() vocabulary, and which parameter a hypothesis on the
## duplicated name actually reads.
arm <- commandArgs(trailingOnly = TRUE)[1L]
lib <- if (arm == "lane") c("C:/Users/adf44/source/r/brmsnames-lib",
                            "C:/Users/adf44/source/r/rellib-r3") else
  "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, "\n")
set.seed(12)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), x_z = rnorm(n),
                sigma_z = rnorm(n), sigma_Intercept = rnorm(n),
                Intercept = rnorm(n),
                g = factor(rep(1:10, 30)), h2 = factor(rep(1:3, 100)))
d$gh2 <- factor(rep(1:6, 50))
d$y <- 1 + 2 * d$sigma_z - 1.5 * d$sigma_Intercept + 0.7 * d$x_z +
  rnorm(10, 0, 1)[d$g] + rnorm(6, 0, 0.2)[d$gh2] +
  rnorm(30, 0, 2)[interaction(d$g, d$h2)] +
  rnorm(n, 0, exp(0.4 * d$z))
d$y_x <- 0.5 * d$z + rnorm(n)
ms <- list(
  K1_cov_sigma_Intercept = bf(y ~ sigma_Intercept + (1 | g)),
  K2_mu_sigma_z_and_sigma_z = bf(y ~ sigma_z + (1 | g), sigma ~ z),
  K3_group_g_h2_vs_gh2 = bf(y ~ 1 + (1 | g:h2) + (1 | gh2)),
  K4_mv_resp_underscore = mvbf(bf(y ~ x_z), bf(y_x ~ z), rescor = FALSE),
  K5_cov_Intercept = bf(y ~ Intercept + x)
)
for (nm in names(ms)) {
  cat("\n==", nm, "==\n")
  fit <- tryCatch(q(frm(ms[[nm]], family = gaussian(), data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) { cat("frm ERROR:", conditionMessage(fit), "\n"); next }
  cat("fixef:", deparse(lapply(fixef(fit), round, 4)), "\n")
  vo <- frmtmb:::hyp_vals_only(fit)
  env <- frmtmb:::hyp_env_vals(fit, vo$vals, vo$comp)
  cat("variables:", tryCatch(variables(fit), error = function(e) conditionMessage(e)), "\n")
  if (arm == "lane") {
    lab <- brms_par_labels(fit)
    cat("draws labels duplicated:", unique(lab[duplicated(lab)]), "\n")
    cn <- frmtmb:::brms_coef_names(fit)
    cat("brms_coef_names duplicated:", unique(cn[duplicated(cn)]), "\n")
  }
  if (nm == "K3_group_g_h2_vs_gh2") {
    vm <- if (arm == "lane") varcorr_matrices(fit) else VarCorr(fit)
    cat("sqrt block variances:", paste(names(vm), signif(sqrt(sapply(vm, function(m) m[1, 1])), 6)), "\n")
    cat("env sd_gh2__Intercept =", signif(env[["sd_gh2__Intercept"]], 6), "\n")
  }
  hs <- switch(nm,
    K1_cov_sigma_Intercept = "sigma_Intercept = 0",
    K2_mu_sigma_z_and_sigma_z = "sigma_z = 0",
    K3_group_g_h2_vs_gh2 = NULL,
    K4_mv_resp_underscore = "y_x_z = 0",
    K5_cov_Intercept = "Intercept = 0")
  if (!is.null(hs)) {
    h <- tryCatch(q(hypothesis(fit, hs)), error = function(e) e)
    if (inherits(h, "error")) cat("hypothesis ERROR:", conditionMessage(h), "\n") else {
      est <- if (inherits(h, "brmshypothesis")) h$hypothesis$Estimate else h$estimate
      cat("hypothesis('", hs, "') estimate:", signif(est, 6), "\n", sep = "")
    }
  }
}
cat("DONE\n")
