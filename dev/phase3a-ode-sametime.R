# Item 3.1, punch round 1: events that share a time, and the boundary
# of the reset-inside-infusion refusal.
#
# Part M: the reviewer's constructions M2 to M7 as described in the
# punch, rebuilt here (M5 and M6 are this lane's, replace and multiply
# beside a bolus). One-compartment oral model, ka 1.1, ke 0.23, amounts,
# atol = rtol = 1e-12. For each: frm_ode() and frm_lincmt() on the
# events table, rxode2::rxSolve() on the records in both listed orders,
# and what frm_ode_records() does with the records.
# Part G: 300 random schedules whose times are drawn from a 6-hour grid,
# so that events collide; seeds 20260924 + k. Every schedule is either
# refused by name or agrees with rxSolve.
# Run with PHASE3A_ARM=base to see what the released build's records
# function (absent there) does: part M then reports "absent".
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.ode); library(rxode2)
})
phase3a_where("frmtmb.ode"); phase3a_where("rxode2")
KA <- 1.1; KE <- 0.23; TOL <- 1e-12
mod <- rxode2({
  d/dt(depot) = -ka * depot
  d/dt(central) = ka * depot - ke * central
})
pk_dyn <- function(t, y, p) list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
has_rec <- exists("frm_ode_records", asNamespace("frmtmb.ode"))

rx_at <- function(rec) {
  r <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), rec, atol = TOL,
                             rtol = TOL, maxsteps = 1e6))
  r$central
}
fo_at <- function(ev, tt) {
  tryCatch(as.numeric(frm_ode(pk_dyn, init = list(0, 0), times = tt,
                     parms = list(KA, KE), states = c("depot", "central"),
                     events = ev, atol = TOL, rtol = TOL)[, 2]),
           error = function(e) NA_real_)
}
fl_at <- function(ev, tt) {
  tryCatch(as.numeric(frm_lincmt(parms = list(ka = KA, ke = KE, V = 1),
                                 times = tt, ncmt = 1, depot = TRUE,
                                 output = "central", events = ev)),
           error = function(e) NA_real_)
}
rec_says <- function(rec) {
  if (!has_rec) return("absent")
  r <- tryCatch(frmtmb.ode::frm_ode_records(rec),
                error = function(e) e, warning = function(w) w)
  if (inherits(r, "frmtmb_ode_error")) {
    return(paste("REFUSED:", substr(conditionMessage(r), 19, 95)))
  }
  if (inherits(r, "condition")) return(paste("other:", conditionMessage(r)))
  "accepted"
}
rec <- function(...) {
  d <- data.frame(...)
  for (nm in c("rate", "ii", "addl", "ss")) if (is.null(d[[nm]])) d[[nm]] <- 0
  d$id <- 1
  d
}
show <- function(label, recA, recB, ev, tt) {
  cat(sprintf("\n%s\n", label))
  cat("  frm_ode   ", sprintf("%.2f", fo_at(ev, tt)), "\n")
  cat("  frm_lincmt", sprintf("%.2f", fl_at(ev, tt)), "\n")
  cat("  rxode2 A  ", sprintf("%.2f", rx_at(recA)), "\n")
  if (!is.null(recB)) cat("  rxode2 B  ", sprintf("%.2f", rx_at(recB)), "\n")
  cat("  records A:", rec_says(recA), "\n")
  if (!is.null(recB)) cat("  records B:", rec_says(recB), "\n")
}

cat("== M: constructions, central amount ==\n")
# M2: a bolus into central, then a reset, both at t = 10
m2 <- rec(time = c(10, 10, 11), evid = c(1, 3, 0), amt = c(100, 0, 0),
          cmt = c(2, NA, 2))
show("M2 bolus then evid 3 at t = 10, read at t = 11 (listed order: 0)",
     m2, m2[c(2, 1, 3), ],
     data.frame(time = c(10, 10), state = c(2, NA), value = c(100, 0),
                method = c("add", "reset")), 11)
# M3: a bolus into central, then an ss = 1 dose into the depot, t = 10
m3 <- rec(time = c(10, 10, 11), evid = c(1, 1, 0), amt = c(100, 50, 0),
          cmt = c(2, 1, 2), ii = c(0, 12, 0), ss = c(0, 1, 0))
show("M3 bolus and an ss = 1 dose at t = 10, read at t = 11",
     m3, m3[c(2, 1, 3), ],
     data.frame(time = c(10, 10), state = c(2, 1), value = c(100, 50),
                ii = c(0, 12), ss = c(FALSE, TRUE)), 11)
# M4: infusion into central at 0 every 4 h, three in all, reset at 4
m4 <- rec(time = c(0, 4, 6, 13, 20), evid = c(1, 3, 0, 0, 0),
          amt = c(100, 0, 0, 0, 0), rate = c(25, 0, 0, 0, 0),
          ii = c(4, 0, 0, 0, 0), addl = c(2, 0, 0, 0, 0),
          cmt = c(2, NA, 2, 2, 2))
show("M4 4-hour infusion q4 x3 into central, evid 3 at t = 4 (a repeat starts there)",
     m4, NULL,
     data.frame(time = c(0, 4), state = c(2, NA), value = c(100, 0),
                method = c("add", "reset"), duration = c(4, 0),
                ii = c(4, 0), addl = c(2L, 0L)), c(6, 13, 20))
# M5, M6: a bolus beside a replace, and beside a multiply
m5 <- rec(time = c(10, 10, 11), evid = c(1, 5, 0), amt = c(100, 30, 0),
          cmt = c(2, 2, 2))
show("M5 bolus then replace 30 at t = 10 (listed order: 30 at t = 10)",
     m5, m5[c(2, 1, 3), ],
     data.frame(time = c(10, 10), state = c(2, 2), value = c(100, 30),
                method = c("add", "replace")), 11)
m6 <- rec(time = c(0, 10, 10, 11), evid = c(1, 1, 6, 0),
          amt = c(100, 100, 0.5, 0), cmt = c(2, 2, 2, 2))
show("M6 bolus then multiply 0.5 at t = 10", m6, m6[c(1, 3, 2, 4), ],
     data.frame(time = c(0, 10, 10), state = c(2, 2, 2),
                value = c(100, 100, 0.5),
                method = c("add", "add", "multiply")), 11)
# M7: a bolus listed before an evid 4 at the same time
m7 <- rec(time = c(0, 10, 10, 11), evid = c(1, 1, 4, 0),
          amt = c(100, 100, 20, 0), cmt = c(1, 2, 1, 2))
show("M7 bolus listed before evid 4 at t = 10", m7, NULL,
     data.frame(time = c(0, 10, 10, 10), state = c(1, 2, NA, 1),
                value = c(100, 100, 0, 20),
                method = c("add", "add", "reset", "add")), 11)
# an addl repeat landing on a reset
ma <- rec(time = c(0, 12, 13), evid = c(1, 3, 0), amt = c(100, 0, 0),
          cmt = c(2, NA, 2), ii = c(6, 0, 0), addl = c(3, 0, 0))
show("addl repeat at t = 12 on an evid 3", ma, NULL,
     data.frame(time = c(0, 12), state = c(2, NA), value = c(100, 0),
                method = c("add", "reset"), ii = c(6, 0),
                addl = c(3L, 0L)), 13)
# a sample at an ss record's time, dose into the observed compartment
mo <- rec(time = c(0, 0), evid = c(0, 1), amt = c(0, 100), cmt = c(2, 2),
          ii = c(0, 12), ss = c(0, 1))
show("sample listed before an ss = 1 dose into central at t = 0",
     mo, mo[2:1, ],
     data.frame(time = 0, state = 2, value = 100, ii = 12, ss = TRUE), 0)
# the control: two plain doses at one time pass and agree
mc <- rec(time = c(10, 10, 11), evid = c(1, 1, 0), amt = c(100, 50, 0),
          cmt = c(2, 1, 2))
show("control: two plain doses at t = 10", mc, mc[c(2, 1, 3), ],
     data.frame(time = c(10, 10), state = c(2, 1), value = c(100, 50)), 11)

if (!has_rec) quit(save = "no")

cat("\n== G: 300 random schedules on a 6-hour grid ==\n")
grid_sched <- function(seed) {
  set.seed(seed)
  rows <- list()
  n_id <- sample(1:3, 1)
  for (id in seq_len(n_id)) {
    nd <- sample(1:5, 1)
    tt <- sort(sample(seq(0, 48, by = 6), nd, replace = TRUE))
    for (k in seq_len(nd)) {
      kind <- sample(c("bolus", "bolus", "inf", "ss", "reset", "reset4",
                       "rep", "mul"), 1)
      cmt <- sample(1:2, 1)
      r <- data.frame(id = id, time = tt[k], evid = 1,
                      amt = round(stats::runif(1, 10, 200), 1), cmt = cmt,
                      rate = 0, ii = 0, addl = 0, ss = 0)
      if (kind == "inf") { r$cmt <- 2; r$rate <- round(stats::runif(1, 20, 80)) }
      if (kind == "ss") { r$ii <- 12 }
      if (kind == "ss") r$ss <- 1
      if (kind == "reset") { r$evid <- 3; r$amt <- 0; r$cmt <- NA }
      if (kind == "reset4") r$evid <- 4
      if (kind == "rep") { r$evid <- 5; r$cmt <- 2 }
      if (kind == "mul") { r$evid <- 6; r$cmt <- 2; r$amt <- 0.5 }
      if (kind %in% c("bolus", "inf") && stats::runif(1) < 0.3) {
        r$ii <- 6; r$addl <- sample(1:3, 1)
      }
      rows[[length(rows) + 1L]] <- r
    }
    ot <- sort(unique(round(stats::runif(10, 0.5, 70), 2)))
    rows[[length(rows) + 1L]] <- data.frame(id = id, time = ot, evid = 0,
                                            amt = 0, cmt = 2, rate = 0,
                                            ii = 0, addl = 0, ss = 0)
  }
  d <- do.call(rbind, rows)
  # records at one time: events in drawn order, and an observation at an
  # event time after it, so the ordering refusals are reachable
  d[order(d$id, d$time, d$evid == 0), ]
}
reason <- function(msg) {
  pats <- c(shared = "shares time", infusion = "while the infusion",
            obs_after = "follow.* an event record", ss_obs = "steady-state \\(ss = 1\\) record",
            other = ".")
  names(pats)[vapply(pats, grepl, TRUE, x = msg)][1L]
}
res <- NULL
for (k in 1:300) {
  d <- grid_sched(20260924L + k)
  x <- tryCatch(frm_ode_records(d), error = identity)
  rx <- tryCatch(as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), d,
                                       atol = TOL, rtol = TOL,
                                       maxsteps = 1e6)),
                 error = function(e) NULL)
  rx_min <- if (is.null(rx)) NA else min(c(rx$depot, rx$central))
  if (inherits(x, "error")) {
    res <- rbind(res, data.frame(seed = 20260924L + k, status = "refused",
                                 why = reason(conditionMessage(x)),
                                 classed = inherits(x, "frmtmb_ode_error"),
                                 rel = NA, zero_scale = NA,
                                 rx_min = rx_min))
    next
  }
  fo <- tryCatch(frm_ode(pk_dyn, init = list(0, 0), times = x$data$time,
                parms = list(KA, KE), group = x$data$id,
                states = c("depot", "central"), events = x$events,
                atol = TOL, rtol = TOL), error = identity)
  if (inherits(fo, "error")) {
    # a refusal of frm_ode() itself, by name, on what the split handed it
    res <- rbind(res, data.frame(seed = 20260924L + k,
      status = "frm_ode refused", zero_scale = NA,
      why = if (grepl("more than one row .ss = TRUE", conditionMessage(fo)))
        "two ss rows in one group" else substr(conditionMessage(fo), 1, 40),
      classed = inherits(fo, "frmtmb_ode_error"), rel = NA,
      rx_min = rx_min))
    next
  }
  sc <- max(abs(c(rx$depot, rx$central)))
  res <- rbind(res, data.frame(
    seed = 20260924L + k, status = "accepted", why = "", classed = NA,
    # a schedule whose every amount is zero at the samples is compared
    # absolutely, since its scale is 0
    rel = max(abs(fo[, 2] - rx$central), abs(fo[, 1] - rx$depot)) /
      (if (sc > 0) sc else 1),
    zero_scale = sc == 0,
    rx_min = rx_min))
}
print(table(res$status, res$why))
cat(sprintf("refusals classed frmtmb_ode_error: %d of %d\n",
            sum(res$classed & res$status == "refused", na.rm = TRUE),
            sum(res$status == "refused")))
acc <- res[res$status == "accepted", ]
cat(sprintf("accepted schedules with every sampled amount zero: %d
",
            sum(acc$zero_scale)))
cat(sprintf("accepted %d: worst relative difference %.2e, median %.2e\n",
            nrow(acc), max(acc$rel), stats::median(acc$rel)))
cat(sprintf("accepted with rxode2 below zero anywhere: %d (min %.3g)\n",
            sum(acc$rx_min < -1e-8), min(acc$rx_min)))
cat(sprintf("refused with rxode2 below zero anywhere: %d\n",
            sum(res$status == "refused" & res$rx_min < -1e-8, na.rm = TRUE)))
