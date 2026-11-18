# Lane ordstarts, build round: where the time of the swapped start goes
# on cum_sratio (1 | g) n = 1000 seed 13 of dev/ordstarts-run.R. Times
# the fit with starts = 0, the default, the default with the jitter check
# off, and the default with the covariance of the class off, by
# replacing internals in the namespace (dev only).
#   Rscript dev/ordstarts-prof.R [seed]
seed <- as.integer(commandArgs(trailingOnly = TRUE)[1] %||% 13)
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
nobs <- 1000
G <- 50
re <- "int"
run <- readLines("dev/ordstarts-run.R")
eval(parse(text = run[grep("^mkg <- ", run):(grep("^form <- ", run) - 1L)]))
d <- mkg(seed)
q <- function(e) suppressWarnings(suppressMessages(e))
tm <- function(lab, ctl = frmtmb_control()) {
  t0 <- proc.time()[["elapsed"]]
  f <- q(frm(bf(y ~ x + z + (1 | g)), family = mixture(cumulative(), sratio()),
             data = d, control = ctl))
  cat(sprintf("%-28s %.2fs ll %.4f evals %s\n", lab,
              proc.time()[["elapsed"]] - t0, as.numeric(logLik(f)),
              format(f$opt$evals %||% NA)))
}
tm("starts = 0", frmtmb_control(starts = 0))
tm("default")
utils::assignInNamespace("mixture_swap_jitter_runs", 0L, ns = "frmtmb")
tm("default, no jitter check")
cls <- ns$mixture_ord_start_class
utils::assignInNamespace("mixture_ord_thres_se", function(fit, resp,
                                                        force = FALSE) {
  K <- fit$spec$responses[[resp]]$family$mix$ord$K
  out <- rep(NA_real_, K)
  attr(out, "lost") <- rep(FALSE, K)
  out
}, ns = "frmtmb")
tm("default, no jitter, no cov")
