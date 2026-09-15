# Reproduce the regression R CMD check found, and show the fix, in one
# process each rather than by running the whole suite.
#
# `hypothesis()`'s shadowing note fired from the GENERIC. Once frmtmb's
# binding resolves to brms's generic the note is gone. One file per
# process never saw it, because nothing in that file loads brms; the
# single-process check did, because an earlier file had.
#
#   Rscript dev/generics-notice.R <LIB> <brms|nobrms>
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (identical(av[2], "brms")) suppressMessages(loadNamespace("brms"))
suppressMessages(library(frmtmb))

set.seed(3)
n <- 200
dd <- data.frame(v = rnorm(n), g = factor(rep(1:10, length.out = n)),
                 y = 0)
dd$y <- frm_simulate(bf(y ~ v + (1 | g)) + gaussian(), dd,
                     newparams = list(Intercept = 1, v = 0.7, sigma = 1,
                                      sd_g__Intercept = 0.5),
                     nsim = 1, seed = 1003L)[[1]]
names(dd)[1L] <- "sigma"
fit <- frm(bf(y ~ sigma + (1 | g)) + gaussian(), data = dd)

msgs <- testthat::capture_messages(hypothesis(fit, "sigma = 0"))
cat(sprintf("brms loaded            %s\n", isNamespaceLoaded("brms")))
cat(sprintf("hypothesis resolves to %s\n",
            environmentName(topenv(environment(get("hypothesis"))))))
cat(sprintf("shadowing notes        %d\n", length(msgs)))
cat("DONE\n")
