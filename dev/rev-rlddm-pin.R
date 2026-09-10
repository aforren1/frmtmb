# REVIEW, item 1.0b, attack 4: the bf(ndt = ) defect, read off the
# OBJECTIVE rather than reconstructed.
#
#   Rscript dev/rev-rlddm-pin.R <lib>
#
# Seed 4242, the same design dev/rlddm-scripts/rlddm-pin.R uses. That
# script reports the non-decision time the density used by pushing the
# constant through the declared link by hand and back through the fitted
# one. This one reads the value the frame actually put in the parameter
# template, because frmtmb transforms a constant dpar TWICE, at
# R/parse.R:1266 and again at R/frame.R:2581, and a reconstruction that
# only knows about the first would report a defect that the second
# repairs.

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
want <- 0.2
cat("min(rt)          :", format(min(s$rt), digits = 9), "\n")

# what the frame PUT in the template for the constant, and what the
# fitted family's own link makes of it
read_pinned <- function(o) {
  fr <- o[["frame"]]
  fam <- frmtmb::single_response(o)[["family"]]
  lk <- fam[["links"]][["ndt"]]
  hit <- NULL
  for (lp in fr[["linpreds"]]) {
    if (identical(lp[["dpar"]], "ndt")) hit <- lp
  }
  if (is.null(hit)) return(list(eta = NA_real_, ndt = NA_real_))
  eta <- fr[["par_template"]][[hit[["par"]]]][hit[["idx"]][1L]]
  list(eta = as.numeric(eta), ndt = as.numeric(lk$linkinv(eta)),
       link = lk$name %||% "(function)")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

probe <- function(lab, fam, form) {
  cat("\n--", lab, "--\n")
  out <- tryCatch({
    o <- frm(form, family = fam, data = s, dry_run = "objective")
    r <- read_pinned(o)
    cat("fitted ndt link  :", r$link, "\n")
    cat("eta in template  :", format(r$eta, digits = 9), "\n")
    cat("ndt the density used:", format(r$ndt, digits = 9), "\n")
    cat("relative error   :",
        sprintf("%.2f%%", 100 * abs(r$ndt - want) / want), "\n")
    cat("nll at par       :",
        sprintf("%.9f", as.numeric(o$obj$fn(o$obj$par))), "\n")
    "ACCEPTED"
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat("outcome          :", substr(out, 1, 260), "\n")
  invisible(out)
}

f_pin <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
            bs = 1.6, ndt = want, bias = 0.5)

probe("bf(ndt = 0.2), no max_ndt",
      rlddm(subject = id, trial = trial), f_pin)

# THE OVER-REFUSAL PROBE. A legitimate pinned non-decision time, with
# the bound stated up front, must still work and must mean what it says.
probe("bf(ndt = 0.2), rlddm(max_ndt = 0.26)",
      rlddm(subject = id, trial = trial, max_ndt = 0.26), f_pin)

# and one that is NOT legitimate: pinned above the stated bound
probe("bf(ndt = 0.3), rlddm(max_ndt = 0.26)",
      rlddm(subject = id, trial = trial, max_ndt = 0.26),
      bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3, bs = 1.6,
         ndt = 0.3, bias = 0.5))

# the ordinary fit, which must be untouched by any of this
cat("\n-- no pin at all, ndt estimated --\n")
o <- tryCatch(frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
                     drift = 3, bs = 1.6, bias = 0.5),
                  family = rlddm(subject = id, trial = trial), data = s,
                  dry_run = "objective"),
              error = function(e) e)
if (inherits(o, "error")) {
  cat("REFUSED:", conditionMessage(o), "\n")
} else {
  cat("nll at start par :",
      sprintf("%.9f", as.numeric(o$obj$fn(o$obj$par))), "\n")
}

# other things a user does with a family object before frm() sees it
cat("\n-- a bare family, before frm() --\n")
fam <- rlddm(subject = id, trial = trial)
cat("class            :", paste(class(fam), collapse = "/"), "\n")
cat("print()          :",
    tryCatch({ capture.output(print(fam)); "ok" },
             error = function(e) paste("REFUSED:", conditionMessage(e))),
    "\n")
cat("links$ndt$name   :",
    tryCatch(fam[["links"]][["ndt"]]$name,
             error = function(e) paste("REFUSED:",
                                       conditionMessage(e))), "\n")
cat("linkinv(0)       :",
    tryCatch(format(fam[["links"]][["ndt"]]$linkinv(0), digits = 9),
             error = function(e) paste("REFUSED:",
                                       substr(conditionMessage(e), 1, 90))),
    "\n")
fam2 <- rlddm(subject = id, trial = trial, max_ndt = 0.26)
cat("max_ndt linkinv(0):",
    tryCatch(format(fam2[["links"]][["ndt"]]$linkinv(0), digits = 9),
             error = function(e) paste("REFUSED:",
                                       substr(conditionMessage(e), 1, 90))),
    "\n")
