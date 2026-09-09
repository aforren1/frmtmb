# Targeted steady-state corners the lane's 118 + 9 do not contain, and
# the two places where a difference would be a WRONG READING of the
# grammar rather than arithmetic. Reference: frm_ode() at
# atol = rtol = 1e-12 with a matched n_ss.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
TOL <- 1e-12
NSS <- 20L
P <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1, V = 10)
dyn2 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[[4]] * y[1],
         p[[4]] * y[1] - (p[[1]] + p[[2]]) * y[2] + p[[3]] * y[3],
         p[[2]] * y[2] - p[[3]] * y[3]))
}
cmp <- function(tag, times, ev, t0 = 0, ini = NULL, n_ss = NSS,
                esc = 1) {
  a <- tryCatch(frm_lincmt(parms = P, times = times, ncmt = 2,
                           depot = TRUE, t0 = t0, init = ini,
                           events = ev, event_scale = esc,
                           n_ss = n_ss),
                error = function(e) conditionMessage(e))
  i0 <- list(0, 0, 0)
  if (!is.null(ini)) {
    if (!is.null(ini$depot)) i0[[1]] <- ini$depot
    if (!is.null(ini$central)) i0[[2]] <- ini$central
  }
  b <- tryCatch(frm_ode(dyn2, init = i0, times = times,
                        parms = list(P$ke, P$k12, P$k21, P$ka),
                        states = c("depot", "central", "peripheral1"),
                        output = "central", t0 = t0, events = ev,
                        event_scale = esc, n_ss = n_ss, atol = TOL,
                        rtol = TOL) / P$V,
                error = function(e) conditionMessage(e))
  if (is.character(a) || is.character(b)) {
    cat(sprintf("%-44s lincmt %-8s ode %-8s\n", tag,
                if (is.character(a)) "REFUSED" else "ok",
                if (is.character(b)) "REFUSED" else "ok"))
    if (is.character(a)) cat("      lincmt:", substr(a, 1, 74), "\n")
    if (is.character(b)) cat("      ode   :", substr(b, 1, 74), "\n")
    return(invisible(NULL))
  }
  sc <- max(abs(b))
  cat(sprintf("%-44s err/scale %10.3e   scale %.4g\n", tag,
              max(abs(a - b)) / sc, sc))
}

cat("\n=== steady-state corners ===\n")
cmp("two ss rows in one group",
    c(0, 4, 8, 20, 24, 28, 40),
    data.frame(time = c(0, 24), state = "depot", value = c(100, 60),
               method = "add", duration = 0, ii = c(8, 6), addl = 0L,
               ss = TRUE))
cmp("reset at the same instant as an ss row",
    c(0, 4, 8, 12, 16),
    data.frame(time = c(8, 8, 0), state = c("depot", NA, "depot"),
               value = c(100, 0, 100), method = c("add", "reset",
                                                  "add"),
               duration = 0, ii = c(8, 0, 8), addl = 0L,
               ss = c(TRUE, FALSE, FALSE)))
cmp("ss row at a t0 that is not zero", c(2.5, 6, 10, 18),
    data.frame(time = 2.5, state = "depot", value = 100,
               method = "add", duration = 0, ii = 8, addl = 0L,
               ss = TRUE), t0 = 2.5)
cmp("ss row, observation 30 intervals later",
    c(0, 8, 120, 240),
    data.frame(time = 0, state = "depot", value = 100, method = "add",
               duration = 0, ii = 8, addl = 0L, ss = TRUE))
cmp("ss row and init both present", c(0, 2, 8, 16),
    data.frame(time = 0, state = "depot", value = 100, method = "add",
               duration = 0, ii = 8, addl = 0L, ss = TRUE),
    ini = list(depot = 50, central = 30))
cmp("ss into central with duration == ii", c(0, 2, 4, 8, 12),
    data.frame(time = 0, state = "central", value = 100,
               method = "add", duration = 8, ii = 8, addl = 0L,
               ss = TRUE))
cmp("ss row then a reset then more doses", c(0, 4, 12, 16, 24),
    data.frame(time = c(0, 12, 14), state = c("depot", NA, "depot"),
               value = c(100, 0, 80), method = c("add", "reset",
                                                 "add"),
               duration = 0, ii = c(8, 0, 0), addl = 0L,
               ss = c(TRUE, FALSE, FALSE)))
cmp("ss row with an event_scale", c(0, 3, 8, 11),
    data.frame(time = 0, state = "depot", value = 100, method = "add",
               duration = 0, ii = 8, addl = 0L, ss = TRUE),
    esc = 0.63)
cmp("ss row, n_ss = 1", c(0, 3, 8, 11),
    data.frame(time = 0, state = "depot", value = 100, method = "add",
               duration = 0, ii = 8, addl = 0L, ss = TRUE), n_ss = 1L)
cmp("ss row, n_ss = 200", c(0, 3, 8, 11),
    data.frame(time = 0, state = "depot", value = 100, method = "add",
               duration = 0, ii = 8, addl = 0L, ss = TRUE),
    n_ss = 200L)

cat("\n=== ordinary corners a real dataset produces ===\n")
cmp("dose after the last observation", c(0, 2, 4),
    data.frame(time = c(0, 6), state = "depot", value = 100,
               method = "add", duration = 0, ii = 0, addl = 0L,
               ss = FALSE))
cmp("addl block, 30 doses, observed at the end",
    c(0, 120, 232, 240),
    data.frame(time = 0, state = "depot", value = 100, method = "add",
               duration = 0, ii = 8, addl = 29L, ss = FALSE))
cmp("two infusions overlapping into central", c(0, 1, 2, 3, 6, 12),
    data.frame(time = c(0, 1), state = "central", value = c(100, 60),
               method = "add", duration = c(3, 2), ii = 0, addl = 0L,
               ss = FALSE))
cmp("infusion with ii/addl, duration < ii", c(0, 1, 5, 9, 13),
    data.frame(time = 0, state = "central", value = 100,
               method = "add", duration = 2, ii = 4, addl = 3L,
               ss = FALSE))
cmp("infusion with ii/addl, duration > ii (overlapping)",
    c(0, 1, 5, 9, 13),
    data.frame(time = 0, state = "central", value = 100,
               method = "add", duration = 6, ii = 4, addl = 3L,
               ss = FALSE))
cmp("negative dose value", c(0, 1, 4, 8),
    data.frame(time = c(0, 2), state = c("depot", "central"),
               value = c(100, -20), method = "add", duration = 0,
               ii = 0, addl = 0L, ss = FALSE))
cmp("observation exactly at the end of an infusion",
    c(0, 2, 2.0000001, 4),
    data.frame(time = 0, state = "central", value = 100,
               method = "add", duration = 2, ii = 0, addl = 0L,
               ss = FALSE))
cmp("reset with no dose after it", c(0, 2, 6, 10),
    data.frame(time = c(0, 4), state = c("depot", NA),
               value = c(100, 0), method = c("add", "reset"),
               duration = 0, ii = 0, addl = 0L, ss = FALSE))
