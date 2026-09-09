# Edge cases of the schedule contract, where the closed form and the
# solver could disagree about SEMANTICS rather than about arithmetic.
# A disagreement here is a defect in one of the two, not a tolerance.
source("lincmt-src.R")

dep1 <- function(t, y, p) list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
iv1 <- function(t, y, p) list(c(-p[1] * y[1]))

cmp <- function(label, a, b) {
  s <- max(abs(b), 1e-300)
  cat(sprintf("%-52s  rel/scale %.3e\n", label, max(abs(a - b)) / s))
  invisible(max(abs(a - b)) / s)
}

P <- list(ka = 1.1, ke = 0.2, V = 1)
tol <- 1e-12

# 1. an ss row at a time LATER than t0: everything before it is
#    discarded and the observations before it are not
tt <- c(0, 2, 6, 10, 12, 14, 20, 24)
ev <- data.frame(time = c(0, 12), state = "depot", value = c(50, 100),
                 ii = c(0, 8), ss = c(FALSE, TRUE))
cmp("ss row at t = 12, observations either side",
    frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
               events = ev, n_ss = 20L),
    frm_ode(dep1, init = list(0, 0), times = tt, parms = list(1.1, 0.2),
            states = c("depot", "central"), output = "central",
            events = ev, n_ss = 20L, atol = tol, rtol = tol))

# 2. a reset at a time that is also an observation time: the reading
#    is PRE-reset
ev2 <- data.frame(time = c(0, 12), state = c("depot", NA),
                  value = c(100, 0), method = c("add", "reset"))
cmp("reset at t = 12, with an observation there",
    frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
               events = ev2),
    frm_ode(dep1, init = list(0, 0), times = tt, parms = list(1.1, 0.2),
            states = c("depot", "central"), output = "central",
            events = ev2, atol = tol, rtol = tol))

# 3. a reset at t0 with init present: init is discarded
ev3 <- data.frame(time = 0, state = NA, value = 0, method = "reset")
cmp("reset at t0, init discarded",
    frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
               init = list(depot = 100), events = ev3),
    frm_ode(dep1, init = list(100, 0), times = tt,
            parms = list(1.1, 0.2),
            states = c("depot", "central"), output = "central",
            events = ev3, atol = tol, rtol = tol))

# 4. an observation exactly at every dose of an addl block
ev4 <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                  addl = 3L)
tt4 <- c(0, 8, 16, 24, 4, 12, 20)
cmp("observations exactly on the addl doses",
    frm_lincmt(parms = P, times = tt4, ncmt = 1, depot = TRUE,
               events = ev4),
    frm_ode(dep1, init = list(0, 0), times = tt4,
            parms = list(1.1, 0.2),
            states = c("depot", "central"), output = "central",
            events = ev4, atol = tol, rtol = tol))

# 5. duplicated observation times, and unsorted rows
tt5 <- c(6, 2, 6, 0, 12, 2)
cmp("duplicated and unsorted observation times",
    frm_lincmt(parms = P, times = tt5, ncmt = 1, depot = TRUE,
               events = ev4),
    frm_ode(dep1, init = list(0, 0), times = tt5,
            parms = list(1.1, 0.2),
            states = c("depot", "central"), output = "central",
            events = ev4, atol = tol, rtol = tol))

# 6. two doses at the same instant into the same compartment
ev6 <- data.frame(time = c(4, 4), state = "depot", value = c(60, 40))
cmp("two add rows at one instant",
    frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
               events = ev6),
    frm_ode(dep1, init = list(0, 0), times = tt, parms = list(1.1, 0.2),
            states = c("depot", "central"), output = "central",
            events = ev6, atol = tol, rtol = tol))

# 7. overlapping infusions into the central compartment
ev7 <- data.frame(time = c(0, 2), value = c(100, 50), duration = c(6, 6),
                  state = "central")
cmp("overlapping infusions",
    frm_lincmt(parms = list(ke = 0.2, V = 1), times = tt, ncmt = 1,
               depot = FALSE, events = ev7),
    frm_ode(iv1, init = list(0), times = tt, parms = list(0.2),
            states = "central", events = ev7, atol = tol, rtol = tol))

# 8. an infusion whose duration equals its interval, at steady state
ev8 <- data.frame(time = 0, state = "central", value = 100,
                  duration = 8, ii = 8, ss = TRUE)
cmp("infusion at steady state with duration == ii",
    frm_lincmt(parms = list(ke = 0.2, V = 1), times = tt, ncmt = 1,
               depot = FALSE, events = ev8, n_ss = 20L),
    frm_ode(iv1, init = list(0), times = tt, parms = list(0.2),
            states = "central", events = ev8, n_ss = 20L, atol = tol,
            rtol = tol))

# 9. a group with no events at all beside groups that have them
gd <- data.frame(id = factor(rep(c("a", "b"), each = 5)),
                 time = rep(c(0, 2, 6, 12, 24), 2))
ev9 <- data.frame(group = "a", time = 0, state = "depot", value = 100)
cmp("a group with an empty schedule",
    frm_lincmt(parms = P, times = gd$time, group = gd$id, ncmt = 1,
               depot = TRUE, events = ev9),
    frm_ode(dep1, init = list(0, 0), times = gd$time, group = gd$id,
            parms = list(1.1, 0.2), states = c("depot", "central"),
            output = "central", events = ev9, atol = tol, rtol = tol))

# 10. how far frm_ode()'s n_ss = 20 run-in is from the steady state,
#     which the closed form does not approximate at all
cat("\nfrm_ode()'s run-in against the exact limit, 2 cmt + depot,",
    "ii = 8:\n")
P2 <- list(ka = 1.1, ke = 0.2, k12 = 0.4, k21 = 0.1, V = 1)
dep2 <- function(t, y, p) {
  list(c(-p[4] * y[1],
         p[4] * y[1] - (p[1] + p[2]) * y[2] + p[3] * y[3],
         p[2] * y[2] - p[3] * y[3]))
}
evss <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                   ss = TRUE)
ts <- c(0, 1, 4, 8)
exact <- frm_lincmt(parms = P2, times = ts, ncmt = 2, depot = TRUE,
                    events = evss)
for (n in c(10L, 20L, 40L, 80L)) {
  o <- suppressWarnings(
    frm_ode(dep2, init = list(0, 0, 0), times = ts,
            parms = list(0.2, 0.4, 0.1, 1.1),
            states = c("depot", "central", "p1"), output = "central",
            events = evss, n_ss = n, atol = tol, rtol = tol))
  l <- frm_lincmt(parms = P2, times = ts, ncmt = 2, depot = TRUE,
                  events = evss, n_ss = n)
  cat(sprintf("  n_ss = %3d  ode vs exact %.3e   lincmt(n) vs exact %.3e   ode vs lincmt(n) %.3e\n",
              n, max(abs(o - exact)) / max(abs(exact)),
              max(abs(l - exact)) / max(abs(exact)),
              max(abs(o - l)) / max(abs(exact))))
}
