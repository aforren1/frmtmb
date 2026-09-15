# genrev round 2: does hypothesis() on a frmtmb_draws still give the
# reserved-name shadowing note, now that the arming moved out of core's
# generic and into core's TWO methods?  frmtmb.sample's
# hypothesis.frmtmb_draws() calls hyp_env_vals(), which emits the note
# only when armed.  Same covariate-named-sigma construction as
# tests/testthat/test-naming-collisions.R.  Arm: brms loaded or not.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]; withbrms <- identical(a[2], "brms")
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.sample))
if (withbrms) suppressMessages(loadNamespace("brms"))
options(frmtmb.notices = TRUE)
set.seed(3)
n <- 200
dd <- data.frame(v = rnorm(n), g = factor(rep(1:10, length.out = n)), y = 0)
dd$y <- frmtmb::frm_simulate(frmtmb::bf(y ~ v + (1 | g)) + gaussian(), dd,
                     newparams = list(Intercept = 1, v = 0.7, sigma = 1,
                                      sd_g__Intercept = 0.5),
                     nsim = 1, seed = 1003L)[[1]]
names(dd)[1L] <- "sigma"
fit <- frmtmb::frm(frmtmb::bf(y ~ sigma + (1 | g)) + gaussian(), data = dd)
m_fit <- capture.output(type = "message",
                        h1 <- hypothesis(fit, "sigma = 0"))
dr <- suppressMessages(suppressWarnings(
  frm_sample(fit, chains = 1, iter = 300, warmup = 150, refresh = 0,
             seed = 7)))
m_dr <- tryCatch(capture.output(type = "message",
                                h2 <- hypothesis(dr, "sigma = 0")),
                 error = function(e) paste("ERROR", conditionMessage(e)))
cat(sprintf("LIB %s brms-loaded %s\n", LIB, isNamespaceLoaded("brms")))
cat(sprintf("hypothesis resolves to %s\n",
            environmentName(topenv(environment(get("hypothesis"))))))
cat(sprintf("frmtmb_fit   note lines: %d  %s\n", length(m_fit),
            substr(paste(m_fit, collapse = " "), 1, 70)))
cat(sprintf("frmtmb_draws note lines: %d  %s\n", length(m_dr),
            substr(paste(m_dr, collapse = " "), 1, 70)))
cat("GENREVDONE\n")
