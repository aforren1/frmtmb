# False-alarm rate of frm_lincmt()'s refusals, on the schedule space
# the field actually writes.
#
# The vocabulary is NONMEM's data records and rxode2's et(), crossed
# with the model shapes a population pharmacokinetic analysis fits.
# Each row is CALLED, so accept/refuse is the code's answer and not a
# reading of the documentation. A row marked "not linear" is out of
# scope by construction rather than by refusal, and is counted apart.
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.ode)})

tt <- c(0, 0.5, 1, 2, 4, 8, 12, 24, 36, 48)
P1 <- list(ka = 1.1, ke = 0.2, V = 10)
P2 <- list(ka = 1.1, ke = 0.2, k12 = 0.4, k21 = 0.1, V = 10)

rows <- list()
add <- function(record, shape, a) {
  stopifnot(is.list(a), !is.null(a[["times"]]))
  r <- tryCatch({
    do.call(frm_lincmt, a)
    "accepted"
  }, error = function(e) paste0("REFUSED: ", conditionMessage(e)))
  rows[[length(rows) + 1L]] <<- data.frame(
    record = record, shape = shape, result = r,
    stringsAsFactors = FALSE)
}

base1 <- list(parms = P1, times = tt, ncmt = 1, depot = TRUE)
base2 <- list(parms = P2, times = tt, ncmt = 2, depot = TRUE)
baseiv <- list(parms = list(ke = 0.2, k12 = 0.4, k21 = 0.1, V = 10),
               times = tt, ncmt = 2, depot = FALSE)

# EVID = 1, a dose record
add("evid=1 bolus into the depot", "1 cmt oral",
    c(base1, list(events = data.frame(time = 0, state = "depot",
                                      value = 100))))
add("evid=1 bolus into the central compartment", "2 cmt iv",
    c(baseiv, list(events = data.frame(time = 0, state = "central",
                                       value = 100))))
add("evid=1 amt with ii and addl", "1 cmt oral",
    c(base1, list(events = data.frame(time = 0, state = "depot",
                                      value = 100, ii = 12,
                                      addl = 3L))))
add("evid=1 amt with ss=1", "1 cmt oral",
    c(base1, list(events = data.frame(time = 0, state = "depot",
                                      value = 100, ii = 12,
                                      ss = TRUE))))
add("evid=1 ss=1 then ii/addl", "2 cmt oral",
    c(base2, list(events = data.frame(
      time = c(0, 12), state = "depot", value = 100, ii = c(12, 12),
      addl = c(0L, 3L), ss = c(TRUE, FALSE)))))
add("evid=1 rate>0 into the central compartment (iv infusion)",
    "2 cmt iv",
    c(baseiv, list(events = data.frame(time = 0, state = "central",
                                       value = 100, duration = 2))))
add("evid=1 rate>0 with ii/addl (repeated iv infusion)", "2 cmt iv",
    c(baseiv, list(events = data.frame(time = 0, state = "central",
                                       value = 100, duration = 2,
                                       ii = 12, addl = 3L))))
add("evid=1 rate>0 with ss=1 (iv infusion at steady state)",
    "2 cmt iv",
    c(baseiv, list(events = data.frame(time = 0, state = "central",
                                       value = 100, duration = 2,
                                       ii = 12, ss = TRUE))))
add("evid=1 rate>0 into the depot (zero-order absorption)",
    "1 cmt oral",
    c(base1, list(events = data.frame(time = 0, state = "depot",
                                      value = 100, duration = 1))))
add("evid=1 into a peripheral compartment", "2 cmt oral",
    c(base2, list(events = data.frame(time = 0, state = "peripheral1",
                                      value = 100))))
add("evid=3 reset (NONMEM writes zero)", "1 cmt oral",
    c(base1, list(events = data.frame(
      time = c(0, 24), state = c("depot", NA), value = c(100, 0),
      method = c("add", "reset")))))
add("evid=4 reset and dose at one instant", "1 cmt oral",
    c(base1, list(events = data.frame(
      time = c(0, 24, 24), state = c("depot", NA, "depot"),
      value = c(100, 0, 100),
      method = c("add", "reset", "add")))))
add("a dose amount that carries a bioavailability parameter",
    "1 cmt oral",
    c(base1, list(events = data.frame(time = 0, state = "depot",
                                      value = 100),
                  event_scale = 0.8)))
add("a per-subject schedule (events$group)", "1 cmt oral",
    list(parms = P1, times = rep(tt, 2),
         group = rep(c("a", "b"), each = length(tt)), ncmt = 1,
         depot = TRUE,
         events = data.frame(group = c("a", "b"), time = 0,
                             state = "depot", value = c(100, 50))))
add("an amount present at t0 (init)", "1 cmt oral",
    c(base1, list(init = list(depot = 100))))
add("a dosing history starting before the first sample (t0 > 0)",
    "1 cmt oral",
    list(parms = P1, times = tt[tt >= 2], ncmt = 1, depot = TRUE,
         t0 = 2,
         events = data.frame(time = 4, state = "depot", value = 100)))
add("three-compartment disposition", "3 cmt iv",
    list(parms = list(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.05,
                      k31 = 0.01, V = 10),
         times = tt, ncmt = 3, depot = FALSE,
         events = data.frame(time = 0, state = "central",
                             value = 100)))
add("clearance parameterization (CL/V/Q2/V2)", "2 cmt oral",
    list(parms = list(CL = 2, V = 10, Q2 = 4, V2 = 40, ka = 1.1),
         times = tt, ncmt = 2, depot = TRUE,
         events = data.frame(time = 0, state = "depot", value = 100)))
add("method=replace (a state assignment)", "1 cmt oral",
    c(base1, list(events = data.frame(time = 12, state = "depot",
                                      value = 50,
                                      method = "replace"))))
add("method=multiply", "1 cmt oral",
    c(base1, list(events = data.frame(time = 12, state = "depot",
                                      value = 2,
                                      method = "multiply"))))
add("a reset to a non-zero level", "1 cmt oral",
    c(base1, list(events = data.frame(time = 12, state = NA,
                                      value = 5, method = "reset"))))
add("a rate constant that changes with time (tv)", "1 cmt oral",
    c(base1, list(tv = list(0.2), tv_break = rep(1, length(tt)))))
add("an infusion still running when the schedule restarts",
    "2 cmt iv",
    c(baseiv, list(events = data.frame(
      time = c(0, 4), state = c("central", NA), value = c(100, 0),
      duration = c(8, 0), method = c("add", "reset")))))

res <- do.call(rbind, rows)
res$verdict <- ifelse(res$result == "accepted", "accepted", "REFUSED")
cat(sprintf("%-58s %-12s %s\n", "record", "shape", "verdict"))
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-58s %-12s %s\n", res$record[i], res$shape[i],
              res$verdict[i]))
}
cat("\naccepted", sum(res$verdict == "accepted"), "of", nrow(res),
    "\nrefused ", sum(res$verdict == "REFUSED"), "\n")
cat("\nthe refusals, with the reason the code gave:\n")
for (i in which(res$verdict == "REFUSED")) {
  cat("- ", res$record[i], "\n    ",
      substr(res$result[i], 1, 110), "...\n", sep = "")
}
