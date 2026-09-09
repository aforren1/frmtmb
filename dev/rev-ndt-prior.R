# rev-ndt: what a ported brms prior on `ndt` now means.
#
# SPEC.md section 5 says a distributional parameter's own class is "a
# density on that PARAMETER, carried through its inverse link". The
# inverse link of `ndt` used to give a TIME and now gives a FRACTION, so
# the same prior row is a density on a different quantity. NEWS lists "a
# prior on one" under what stops working; this measures whether it stops
# or merely means something else.
#
# The `class = "Intercept", dpar = "ndt"` spelling is checked beside it,
# because that one is a density on the LINEAR PREDICTOR and the linear
# predictor did not move: log(ndt / (ub - ndt)) and log(f / (1 - f)) are
# the same number.
#
# REV_ARM = "new" | "old". Seed 4242.

arm <- Sys.getenv("REV_ARM", "new")
.libPaths(if (identical(arm, "new")) {
  c(Sys.getenv("REV_LIB",
                "C:/Users/adf44/source/r/rev-ndt-lib2"),
    "C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
} else {
  c("C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
})
suppressMessages({library(frmtmb); library(frmtmb.eam)})

set.seed(4242)
d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
lo <- min(d$rt)
form <- bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)

ndt_of <- function(fit) {
  # the time, whichever build this is
  if (exists("ndt_time", mode = "function")) {
    as.numeric(ndt_time(fit))[1L]
  } else {
    as.numeric(suppressWarnings(
      stats::predict(fit, dpar = "ndt", type = "response")))[1L]
  }
}

# class = "<dpar>" applies only where the dpar has no formula of its
# own, which is brms's rule and frmtmb's, so that spelling gets a model
# that gives `ndt` no predictor.
form_np <- bf(rt | dec(upper) ~ 1, bs ~ 1, bias = 0.5)

run <- function(label, pr, fm = form) {
  r <- tryCatch({
    f <- frm(fm, family = wiener(), data = d, prior = pr)
    sprintf("ndt = %.6f s   logLik %.6f", ndt_of(f),
            as.numeric(logLik(f)))
  }, error = function(e) paste("REFUSED:",
                               substr(conditionMessage(e), 1, 110)))
  cat(sprintf("  %-38s %s\n", label, r))
}

cat("arm", arm, " min(rt)", sprintf("%.6f", lo), "\n")
run("no prior", NULL)
# a brms user's prior, written in SECONDS, which is what brms's own
# wiener ndt is on
run("no prior, ndt with no formula", NULL, form_np)
run("normal(0.30, 0.01) class=ndt",
    prior(normal(0.30, 0.01), class = "ndt"), form_np)
run("normal(0.10, 0.01) class=ndt",
    prior(normal(0.10, 0.01), class = "ndt"), form_np)
run("normal(0.85, 0.01) class=ndt",
    prior(normal(0.85, 0.01), class = "ndt"), form_np)
# and the link-scale spelling, whose meaning did not change
run("normal(0, 0.1) class=Intercept dpar=ndt",
    prior(normal(0, 0.1), class = "Intercept", dpar = "ndt"))
