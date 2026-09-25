# Item 3.1, follow-up to phase3a-theoph.R: the records spelling of the
# Theoph fit ends at nlminb code 1 (false convergence) where the
# init spelling ends at code 0, with the log likelihoods 6.3e-09 apart.
# Is that the objective's integrator noise? Refit the records spelling
# at tighter tolerances and with restarts, and read the codes.
Sys.setenv(PHASE3A_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.ode)
})
pk_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}
d <- datasets::Theoph
d$Subject <- factor(as.character(d$Subject))
dose_rec <- d[!duplicated(d$Subject), c("Subject", "Time", "Dose")]
rec <- rbind(
  data.frame(Subject = d$Subject, Time = d$Time, evid = 0, amt = 0,
             cmt = 2, conc = d$conc),
  data.frame(Subject = dose_rec$Subject, Time = 0, evid = 1,
             amt = dose_rec$Dose, cmt = 1, conc = NA))
rec <- rec[order(rec$Subject, rec$Time, rec$evid), ]
theo <- frm_ode_records(rec, id = "Subject", time = "Time")
doses <- theo$events
st <- list(beta = c(0.5, log(0.08), log(0.5)))
one <- function(label, tol, ctl) {
  w <- character(0)
  t0 <- proc.time()[["elapsed"]]
  fm <- eval(bquote(bf(
    conc ~ frm_ode(pk_dyn, init = list(0, 0), times = Time,
                   parms = list(exp(lka), exp(lke), exp(lV)),
                   group = Subject, states = c("depot", "central"),
                   output = "central", events = doses, atol = .(tol),
                   rtol = .(tol)),
    lka ~ 1 + (1 | Subject), lke ~ 1 + (1 | Subject), lV ~ 1, nl = TRUE)))
  f <- withCallingHandlers(frm(fm +
      gaussian(), data = theo$data, se = TRUE, start = st, control = ctl),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd))
                              invokeRestart("muffleWarning") })
  cat(sprintf("%-26s code %d  %-24s logLik %.8f  max|grad| %.2e  warnings %d  %.0f s\n",
              label, f$opt$convergence, f$opt$message,
              as.numeric(logLik(f)), diagnose(f, quiet = TRUE)$max_grad,
              length(w), proc.time()[["elapsed"]] - t0))
}
one("default tol, restarts 1", 1e-8, frmtmb_control())
one("default tol, restarts 3", 1e-8, frmtmb_control(restarts = 3))
one("tol 1e-10, restarts 1", 1e-10, frmtmb_control())
