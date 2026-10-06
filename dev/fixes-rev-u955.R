# Reviewer of lane fixes, claim 2: the update of brmsfit-methods:955 on
# frmtmb with fit2's brms priors, which the guard exempts: does the
# default (prior-located) start block it, and does a usable start fit?
#   Rscript dev/fixes-rev-u955.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
d <- brms_fixture_data(2)
fo <- bf(count | weights(AgeSD) ~ a + b, a ~ Age + (1 | ID1 | patient),
         b ~ Age + (1 | ID1 | patient), nl = TRUE)
pr <- c(set_prior("normal(2, 2)", nlpar = "a"),
        set_prior("normal(0, 3)", nlpar = "b"))
show <- function(lab, expr) {
  r <- tryCatch(withCallingHandlers(expr, message = function(m) {
    cat("  [message]", conditionMessage(m)); invokeRestart("muffleMessage")
  }, warning = function(w) {
    cat("  [warning]", substr(conditionMessage(w), 1, 120), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat(lab, "ERROR:", substr(conditionMessage(e), 1, 200), "\n"); NULL
  })
  if (!is.null(r)) {
    cat(lab, "conv", r$opt$convergence, "\n")
    print(round(fixef(r)[, 1:2], 4))
  }
}
show("fit2 priors, default start", frm(fo, data = d, family = Gamma("identity"),
                                       prior = pr))
show("fit2 priors, start intercepts mean/2, slopes 0",
     frm(fo, data = d, family = Gamma("identity"), prior = pr,
         start = list(beta = c(mean(d$count) / 2, 0, mean(d$count) / 2, 0))))
# fit2 itself with its brms priors, then the update: does update() keep
# the priors and what happens
fit2p <- tryCatch(suppressMessages(suppressWarnings(
  frm(bf(count | weights(AgeSD) ~ 1 / (1 + exp(-a)) * exp(b * Trt),
         a ~ Age + (1 | ID1 | patient), b ~ Age + (1 | ID1 | patient),
         nl = TRUE), data = d, family = Gamma("identity"), prior = pr))),
  error = function(e) e)
if (inherits(fit2p, "error")) {
  cat("fit2 with priors ERROR:", conditionMessage(fit2p), "\n")
} else {
  show("update(fit2 with priors, formula. = a + b)",
       update(fit2p, formula. = bf(count ~ a + b, nl = TRUE)))
}
