# Reviewer: dev/nanse-rev2-port.R on rellib-r5 (base), generated copy.
# Reviewer, punch round 1: test-portability.R's FN-5 fixture now allows
# a warning on its hypotheses. Is the sd block really flat there, and do
# the hypotheses really move along it?
#   Rscript dev/nanse-rev2-port.R
.libPaths(c(
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(n / 10, 10), y = 0)
d$y <- frm_simulate(bf(y ~ x) + gaussian(), d,
                    newparams = list(b_Intercept = 1, b_x = 2, sigma = 1),
                    nsim = 1, seed = 1001L)[[1]]
fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)), family = gaussian(),
                            data = d))
nm <- ns$outer_par_names(fit)
p <- fit$opt$par
cat("par:", paste(sprintf("%s=%.4g", nm, p), collapse = " "), "\n")
H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
cat("|H| row max:", paste(sprintf("%s=%.2g", nm, apply(abs(H), 1, max)),
                          collapse = " "), "\n")
cat("lost:", paste(sprintf("%s(%s)", names(ns$sdr_of(fit)$se_lost),
                           ns$sdr_of(fit)$se_lost), collapse = " "), "\n")
vc <- VarCorr(fit)
print(vc)
if (requireNamespace("lme4", quietly = TRUE)) {
  m <- lme4::lmer(y ~ x + (1 + x | g), data = d, REML = FALSE)
  cat("lme4 singular:", lme4::isSingular(m), "\n")
  print(lme4::VarCorr(m))
}
for (h in c("sd_g__Intercept - sd_g__x > 0", "sd_g__Intercept > 0")) {
  w <- character()
  r <- withCallingHandlers(hypothesis(fit, h, class = NULL),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  cat(sprintf("%-32s Estimate %.4g Est.Error %s warnings %d\n", h,
              r$hypothesis$Estimate, format(r$hypothesis$Est.Error),
              length(w)))
}
