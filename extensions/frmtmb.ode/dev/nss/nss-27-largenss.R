# NIT 11: the warning went silent at large `n_ss`, in the region
# `?frm_ode` sends a long-half-life user to.
#
# The reviewer found three misses: one compartment, terminal half-life
# 1653 h, `ii` 24, `ss_extrapolate = TRUE`, at `n_ss` 650, 800 and
# 1000, true errors 4.682e-05, 1.035e-05 and 1.385e-06, all silent.
# That `r` is 0.98999, ABOVE `ode_ss_rcap`, so the stand-down declines
# part of the tail; successive extrapolants then agree with each other
# while both sit short of the limit, and a report built only from their
# disagreement cannot see it. The report now adds the part of the tail
# the stand-down declined, which is computable from the same factor.
#
# The trough of a one-compartment intravenous bolus at steady state is
# exact in closed form, so "true error" here is not itself estimated.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-27-largenss.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()
ARM <- Sys.getenv("NSS_ARM", "lane")
SS_TOL <- 1e-6

decay <- function(t, y, p) list(c(-p[1L] * y[1L]))
TAU <- 24
KE <- log(2) / 1653
cat(sprintf("\none compartment iv, half-life %.0f h, ii %g, r = %.6f",
            log(2) / KE, TAU, exp(-KE * TAU)))
cat(sprintf(", rcap = %g\n\n", frmtmb.ode:::ode_ss_rcap))
TR <- 100 * exp(-KE * TAU) / (1 - exp(-KE * TAU))

one <- function(n, ext) {
  args <- list(decay, init = list(0), times = 0, parms = list(KE),
               events = data.frame(time = 0, value = 100, ii = TAU,
                                   ss = TRUE),
               n_ss = n, ss_tol = SS_TOL, atol = 1e-12, rtol = 1e-12)
  if (ARM != "ref") args$ss_extrapolate <- ext
  w <- character(0)
  v <- withCallingHandlers(do.call(frm_ode, args),
                           warning = function(e) {
                             w <<- c(w, conditionMessage(e))
                             invokeRestart("muffleWarning")
                           })
  ww <- w[grepl("run-in", w)]
  list(err = abs(as.numeric(v) - TR) / TR, warned = length(ww) > 0L,
       said = if (!length(ww)) NA_real_ else
         as.numeric(sub("^.*(?:moving by|is still|is about) ([0-9.e+-]+) .*$",
                        "\\1", ww[[1L]])))
}

arms <- if (identical(ARM, "ref")) list(base = NA) else
  list(`ext=TRUE` = TRUE, `ext=FALSE` = FALSE)
for (an in names(arms)) {
  cat(sprintf("== %s ==\n", an))
  cat(sprintf("%7s %12s %8s %12s %10s %s\n", "n_ss", "true error",
              "warned", "said", "said/true", "verdict"))
  miss <- 0L
  fa <- 0L
  for (n in c(20L, 100L, 200L, 400L, 650L, 800L, 1000L, 2000L,
              4000L)) {
    z <- one(n, arms[[an]])
    bad <- z$err > SS_TOL
    v <- if (z$warned && bad) "correct warning"
         else if (!z$warned && !bad) "correct silence"
         else if (z$warned) "FALSE ALARM" else "MISS"
    if (v == "MISS") miss <- miss + 1L
    if (v == "FALSE ALARM") fa <- fa + 1L
    cat(sprintf("%7d %12.3e %8s %12s %10s %s\n", n, z$err, z$warned,
                if (is.na(z$said)) "-" else format(signif(z$said, 3)),
                if (is.na(z$said)) "-" else
                  format(signif(z$said / z$err, 3)), v))
  }
  cat(sprintf("  %d MISS, %d FALSE ALARM\n\n", miss, fa))
}
