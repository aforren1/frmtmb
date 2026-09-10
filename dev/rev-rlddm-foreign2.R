# REVIEW round 2, attack 2, second design: the case where a foreign
# family that ignores `ndt_floor` fails SILENTLY.
#
#   Rscript dev/rev-rlddm-foreign2.R <lib>
#
# Seed 707. dev/rev-rlddm-foreign.R showed the broken arm dying loudly,
# with "NA/NaN gradient evaluation", and the reason is worth stating:
# ndt_bound_attach() sets the `ndt` starting value to the FRACTION one
# half, so a density that reads the fraction as a time starts at 0.5 s,
# and on a design whose responses are faster than that the very first
# evaluation is a log of a negative number. That is luck, not a guard.
#
# This design removes the luck. Every response is slower than 0.5 s and
# every true shift is inside (0, 1), so the broken family's starting
# value is legal, its parameter space contains the truth, and it
# converges. Same shifted gamma, shape 3, interior maximum.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("lib:", lib, " eam", format(packageVersion("frmtmb.eam")), "\n\n")

shifted_gamma <- function(reads_floor) {
  fam <- custom_family(
    "shiftgam",
    dpars = c("rate", "shape", "ndt"),
    primary_dpars = "rate",
    links = list(rate = "log", shape = "log", ndt = "log"),
    accepts_aterms = c("weights", "ndt_group"),
    lpdf = function(y, dpars, aterms) {
      nd <- if (reads_floor) {
        ndt_apply(dpars, aterms, "shiftgam")
      } else {
        dpars[["ndt"]]
      }
      z <- y - nd
      sh <- dpars[["shape"]]
      rt <- dpars[["rate"]]
      (sh - 1) * log(z) - rt * z + sh * log(rt) - lgamma(sh)
    },
    init_dpars = list(rate = function(y, aterms) 3 / mean(y),
                      shape = function(y, aterms) 3,
                      ndt = function(y, aterms) 0.5 * min(y)),
    family_finalize = function(fam, y, aterms) {
      bd <- ndt_bound(y, aterms, NULL, "shiftgam")
      if (!is.null(ndt_bound_of(fam))) return(fam)
      ndt_bound_attach(fam, bd)
    })
  ndt_bound_attach(fam, ndt_bound_pending(NULL, "shiftgam"))
}

set.seed(707)
nper <- 300L
shift <- c(a = 0.60, b = 0.70, c = 0.80, d = 0.90)
g <- rep(names(shift), each = nper)
y <- shift[g] + stats::rgamma(4L * nper, shape = 3, rate = 20)
d <- data.frame(y = as.numeric(y), g = factor(g))
own <- as.numeric(tapply(d$y, d$g, min))
names(own) <- levels(d$g)
cat("true shifts   :", paste(format(shift, digits = 4), collapse = " "),
    "\n")
cat("group floors  :", paste(format(own, digits = 6), collapse = " "),
    "\n")
cat("fastest response overall:", format(min(d$y), digits = 6),
    " (the broken arm's 0.5 s start is legal here)\n\n")

one <- d[match(levels(d$g), as.character(d$g)), , drop = FALSE]
for (nm in c("ok", "broken")) {
  cat("--", nm, "--\n")
  tryCatch({
    fit <- frm(bf(y | ndt_group(g) ~ 1, shape ~ 1, ndt ~ 0 + g),
               family = shifted_gamma(identical(nm, "ok")), data = d)
    dg <- frmtmb::diagnose(fit, quiet = TRUE)
    fr <- as.numeric(suppressWarnings(stats::predict(
      fit, newdata = one, dpar = "ndt", type = "response")))
    tm <- as.numeric(suppressWarnings(ndt_time(fit, newdata = one)))
    cat("  RETURNED. logLik",
        format(as.numeric(stats::logLik(fit)), digits = 9),
        " conv", fit$opt$convergence,
        " maxgrad", format(dg$max_grad, digits = 4),
        " pdHess", isTRUE(dg$pdHess),
        " bad se", length(dg$bad_se), "\n")
    cat("  predict(dpar='ndt', type='response'):",
        paste(format(fr, digits = 5), collapse = " "), "\n")
    cat("  ndt_time()                         :",
        paste(format(tm, digits = 5), collapse = " "), "\n")
    cat("  the truths                         :",
        paste(format(as.numeric(shift), digits = 5), collapse = " "),
        "\n")
    cat("  worst relative error, ndt_time     :",
        sprintf("%.1f%%", 100 * max(abs(tm - shift) / shift)), "\n")
    cat("  every fitted time below its own floor:",
        all(tm < own), "\n\n")
  }, error = function(e) {
    cat("  REFUSED:", substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                             1, 170), "\n\n")
  })
}
