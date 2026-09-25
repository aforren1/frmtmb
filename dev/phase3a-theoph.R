# Item 3.1: the vignette's Theoph fit written through frm_ode_records(),
# against the same model with the dose as the initial condition, which
# is how the vignette wrote it before. Same data, same start, same
# model; the only difference is where the dose enters. frm_lincmt() is
# the third arm, the closed form on the same records.
Sys.setenv(PHASE3A_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.ode)
})
phase3a_where("frmtmb.ode"); phase3a_where("RTMBode")

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
cat("records", nrow(rec), "-> data", nrow(theo$data), "events",
    nrow(doses), "\n")
# the observation frame is Theoph's own rows, in Theoph's order
stopifnot(identical(theo$data$conc, d$conc[order(d$Subject, d$Time)]))
st <- list(beta = c(0.5, log(0.08), log(0.5)))

t0 <- proc.time()[["elapsed"]]
f_init <- frm(bf(
  conc ~ frm_ode(pk_dyn, init = list(Dose, 0), times = Time,
                 parms = list(exp(lka), exp(lke), exp(lV)),
                 group = Subject, states = c("depot", "central"),
                 output = "central"),
  lka ~ 1 + (1 | Subject), lke ~ 1 + (1 | Subject), lV ~ 1, nl = TRUE) +
    gaussian(), data = d, se = TRUE, start = st)
t1 <- proc.time()[["elapsed"]]
f_rec <- frm(bf(
  conc ~ frm_ode(pk_dyn, init = list(0, 0), times = Time,
                 parms = list(exp(lka), exp(lke), exp(lV)),
                 group = Subject, states = c("depot", "central"),
                 output = "central", events = doses),
  lka ~ 1 + (1 | Subject), lke ~ 1 + (1 | Subject), lV ~ 1, nl = TRUE) +
    gaussian(), data = theo$data, se = TRUE, start = st)
t2 <- proc.time()[["elapsed"]]
cat(sprintf("elapsed: init %.1f s, records %.1f s (one run each, not a timing)\n",
            t1 - t0, t2 - t1))
b1 <- unlist(fixef_by_dpar(f_init)); b2 <- unlist(fixef_by_dpar(f_rec))
se1 <- sqrt(diag(vcov(f_init)))[seq_along(b1)]
cat("fixed effects, init:   ", sprintf("%.8f", b1), "\n")
cat("fixed effects, records:", sprintf("%.8f", b2), "\n")
cat(sprintf("max |diff| / se: %.3e\n", max(abs(b1 - b2) / se1)))
cat(sprintf("logLik init %.8f, records %.8f, diff %.3e\n",
            as.numeric(logLik(f_init)), as.numeric(logLik(f_rec)),
            as.numeric(logLik(f_init)) - as.numeric(logLik(f_rec))))
cat("convergence", f_init$opt$convergence, f_rec$opt$convergence, "\n")
# the two objectives at one parameter vector: an identity check on the
# model, free of the optimizer
p <- f_init$obj$env$last.par.best
cat(sprintf("objective at the init fit's optimum: init %.10f, records %.10f\n",
            f_init$obj$fn(f_init$opt$par), f_rec$obj$fn(f_init$opt$par)))
