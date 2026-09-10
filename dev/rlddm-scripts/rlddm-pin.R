# PUNCH ROUND 1, B1. What `bf(ndt = )` on rlddm() actually fitted.
#
#   Rscript dev/rlddm-scripts/rlddm-pin.R <lib>
#
# Seed 4242, six learners by 60 trials, min(rt) = 0.261523214.
#
# THE FIRST VERSION OF THIS SCRIPT WAS WRONG AND ITS NUMBER IS
# RETRACTED. It reconstructed a linear predictor by pushing the pinned
# constant through the family's DECLARED link by hand and reading it
# back through the fitted one, and reported that rlddm() fitted 0.0436 s
# where the user wrote 0.2. Nothing uses that reconstruction. frmtmb
# transforms a constant dpar TWICE and only the second one counts:
#
#   R/parse.R:1266   plain_dpar() calls linkfun(constant) and uses the
#                    result ONLY to decide whether the constant is in
#                    range. That value is discarded and the raw constant
#                    is stored.
#   R/frame.R:2581   par_template[["betad"]][idx] <-
#                      lp[["link"]]$linkfun(lp[["constant"]])
#                    is the transform that reaches the parameter, and by
#                    then family_finalize() has run and lp[["link"]] is
#                    the SETTLED link.
#
# So this reads the number out of the parameter template rather than
# reconstructing one, which is the only place the answer lives.
#
# It also runs the case that IS a defect and that the pending bound
# really does catch: a pinned constant outside the SETTLED link's range.
# There the parse-time check passes on the declared `log` link, and the
# frame's transform then produces NaN with no finiteness check of its
# own.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("lib  :", lib, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n\n")

set.seed(4242)
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 60,
                     seed = 4242)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5),
                       seed = 4242)[[1L]]
cat("min(rt):", format(min(s$rt), digits = 9), "\n\n")

# The pinned constant, read where the fit reads it: out of the parameter
# template, decoded through the link the FITTED family carries.
probe <- function(want, max_ndt = NULL, label) {
  f <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
          bs = 1.6, ndt = want, bias = 0.5)
  fam <- if (is.null(max_ndt)) {
    rlddm(subject = id, trial = trial)
  } else {
    rlddm(subject = id, trial = trial, max_ndt = max_ndt)
  }
  out <- tryCatch({
    o <- frm(f, family = fam, data = s, dry_run = "objective")
    rsp <- frmtmb::single_response(o)
    lk <- rsp[["family"]][["links"]][["ndt"]]
    lps <- o[["frame"]][["linpreds"]]
    hit <- Filter(function(lp) {
      identical(lp[["dpar"]], "ndt") || identical(lp[["name"]], "ndt")
    }, lps)
    eta <- as.numeric(
      o[["frame"]][["par_template"]][["betad"]][hit[[1L]][["idx"]]])
    list(link = lk[["name"]], eta = eta, ndt = lk$linkinv(eta),
         nll = as.numeric(o$obj$fn(o$obj$par)))
  }, error = function(e) conditionMessage(e))
  cat("---", label, "\n")
  if (is.character(out)) {
    cat("  REFUSED:", substr(out, 1, 150), "\n")
    return(invisible(NULL))
  }
  cat("  fitted ndt link       :", out$link, "\n")
  cat("  eta in the template   :", format(out$eta, digits = 9), "\n")
  cat("  the time it decodes to:", format(out$ndt, digits = 9), "\n")
  cat("  relative error        :",
      sprintf("%.2f%%", 100 * abs(out$ndt - want) / want), "\n")
  cat("  negative logLik       :", sprintf("%.9f", out$nll), "\n")
  invisible(NULL)
}

cat("== the case the lane wrongly called a defect ==\n")
probe(0.2, NULL, "bf(ndt = 0.2), no max_ndt")
probe(0.2, 0.26, "bf(ndt = 0.2), max_ndt = 0.26")

cat("\n== the case that IS a defect, and what a pending bound catches ==\n")
probe(0.3, 0.26, "bf(ndt = 0.3), max_ndt = 0.26, OUT OF RANGE")

cat("\n== an ordinary estimated ndt is untouched ==\n")
f_est <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
            bs = 1.6, ndt ~ 1, bias = 0.5)
o <- frm(f_est, family = rlddm(subject = id, trial = trial), data = s,
         dry_run = "objective")
cat("  nll at the start parameter:",
    sprintf("%.9f", as.numeric(o$obj$fn(o$obj$par))), "\n")

cat("\n== what a user SEES for the out-of-range case, on a real fit ==\n")
f_bad <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
            bs = 1.6, ndt = 0.3, bias = 0.5)
# withCallingHandlers, not a warning= handler: the range check itself
# emits `NaNs produced` from log() BEFORE anything refuses, so a
# tryCatch(warning=) returns on the warning and never sees the error.
# A first version of this block did that and reported the two builds as
# identical when they are not.
warns <- character(0)
r <- withCallingHandlers(
  tryCatch({
    fit <- frm(f_bad, family = rlddm(subject = id, trial = trial,
                                     max_ndt = 0.26), data = s)
    paste0("RETURNED A FIT, logLik = ",
           format(as.numeric(stats::logLik(fit))))
  }, error = function(e) paste("REFUSED:", conditionMessage(e))),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
cat("  outcome :", substr(r, 1, 160), "\n")
cat("  warnings:",
    if (length(warns)) paste(unique(warns), collapse = " | ") else "none",
    "\n")

cat("\n== the four eam families already refuse a bare bf(ndt = ) ==\n")
set.seed(1)
dd <- ddm_simulate(200, mu = 1.2, bs = 1.5, ndt = 0.25)
for (nm in c("bare", "max_ndt")) {
  fam <- if (nm == "bare") wiener() else wiener(max_ndt = 0.2)
  m <- tryCatch({
    o <- frm(bf(rt | dec(upper) ~ 1, bs = 1.5, ndt = 0.15, bias = 0.5),
             family = fam, data = dd, dry_run = "objective")
    "ACCEPTED"
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat(sprintf("  wiener(%-8s) bf(ndt = 0.15): %s\n", nm,
              substr(m, 1, 60)))
}

cat("\n== the capability is RESPELLED, not lost ==\n")
# The removed spelling and the supported one must give the same fit, or
# the BREAKING bullet is describing a loss rather than a migration.
f_ok <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
           bs = 1.6, ndt = 0.2, bias = 0.5)
say <- function(lbl, fam) {
  m <- tryCatch(
    sprintf("logLik %.7f", as.numeric(stats::logLik(
      frm(f_ok, family = fam, data = s)))),
    error = function(e) paste("REFUSED:", substr(conditionMessage(e),
                                                 1, 60)))
  cat(sprintf("  %-28s %s\n", lbl, m))
}
say("bf(ndt = 0.2), bare", rlddm(subject = id, trial = trial))
say("bf(ndt = 0.2), max_ndt = 0.26",
    rlddm(subject = id, trial = trial, max_ndt = 0.26))
