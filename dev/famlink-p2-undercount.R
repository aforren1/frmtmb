## Punch round 2, MINOR 1: where fit$opt$nonfinite_trials loses trials.
## The reviewer's three designs (dev/famlink-rev2-undercount.R and
## dev/famlink-rev2-silent.R, seed 20260916, n = 200). For each fit: every
## run_optimizer() call with the trials it mapped, whether the call
## returned or raised, and the count the fit reports.
## Usage: Rscript dev/famlink-p2-undercount.R <lane|base>
ARM <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(ARM)) ARM <- "lane"
source("dev/famlink-rev-common.R")
ns <- asNamespace("frmtmb")
calls <- list()
# Count NaN trials at the objective itself, so a call that raises is
# still counted: the exit trace of run_optimizer sees no return value then.
trace("nlminb_trial_fn", where = ns, print = FALSE,
      exit = quote({
        rv <- returnValue()
        calls[[length(calls) + 1L]] <<- rv$count
      }))
mk <- function(n = 200) {
  set.seed(20260916)
  data.frame(x = runif(n, -2, 2), g = factor(rep(1:20, length.out = n)))
}
designs <- list(
  bern_identity_wall = function() {
    d <- mk()
    d$y <- rbinom(200, 1, pmin(0.98, pmax(0.02, 0.5 + 0.45 * d$x)))
    list(y ~ x, d, bernoulli("identity"))
  },
  gamma_inverse_cross = function() {
    d <- mk()
    mu <- 1 / pmax(0.05, 0.6 + 0.35 * d$x)
    d$y <- rgamma(200, 3, 3 / mu)
    list(y ~ x, d, Gamma("inverse"))
  },
  negbin_identity_re = function() {
    d <- mk()
    d$y <- rnbinom(200, mu = pmax(0.1, 2 + 1.5 * d$x +
                                   rnorm(20, 0, 0.4)[d$g]), size = 3)
    list(y ~ x + (1 | g), d, negbinomial("identity"))
  }
)
for (nm in names(designs)) {
  s <- designs[[nm]]()
  calls <- list()
  nanw <- 0L
  fit <- withCallingHandlers(
    frm(s[[1]], s[[2]], family = s[[3]]),
    warning = function(w) {
      if (grepl("NA/NaN", conditionMessage(w))) nanw <<- nanw + 1L
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  per <- vapply(calls, function(f) as.integer(f()), 1L)
  cat(sprintf("%-20s nlminb calls %d; mapped per call %s; sum %d; fit$opt$nonfinite_trials %s; NaN warnings %d\n",
              nm, length(per), paste(per, collapse = " "), sum(per),
              format(fit$opt$nonfinite_trials), nanw))
}
