# Tape build and gradient cost of the learning engine.
#
#   Rscript dev/learn-timing.R
#
# METHOD, and it is the method rather than the numbers that should be
# copied. Warm session, designs and shapes INTERLEAVED (one measurement
# of each per round, not all of one and then all of the next), gc(FALSE)
# before each measurement, median of 7 tape builds, and gradients timed
# over a BATCH of calls with the best batch of 3 reported.
#
# Single-shot system.time() against Windows' roughly 15 ms clock
# granularity is what produced the numbers frmtmb's own dev/rl-findings.md
# had to retract, so none of it is done here.
#
# The build time is frm(dry_run = "objective") with the
# dry_run = "frame" time subtracted, because frame assembly is not what
# this measures.

library(frmtmb.learn)

REPS <- 7L
BATCH <- 50L
BATCHES <- 3L

bench_grad <- function(form, fam, data) {
  o <- frm(form, family = fam, data = data, dry_run = "objective")
  obj <- o$obj
  p <- obj$par
  obj$gr(p)
  best <- Inf
  for (b in seq_len(BATCHES)) {
    gc(FALSE)
    t <- system.time(for (i in seq_len(BATCH)) obj$gr(p))[["elapsed"]]
    best <- min(best, t / BATCH)
  }
  best * 1000
}

build_once <- function(cs) {
  gc(FALSE)
  fr <- system.time(frm(cs$form, family = cs$fam, data = cs$data,
                        dry_run = "frame"))[["elapsed"]]
  gc(FALSE)
  ob <- system.time(frm(cs$form, family = cs$fam, data = cs$data,
                        dry_run = "objective"))[["elapsed"]]
  max(ob - fr, 0)
}

interleaved_builds <- function(cases, reps = REPS) {
  m <- matrix(NA_real_, reps, length(cases))
  for (r in seq_len(reps)) {
    for (j in seq_along(cases)) m[r, j] <- build_once(cases[[j]])
  }
  apply(m, 2, stats::median)
}

delta_case <- function(ns, nt, seed = 1, fam = NULL) {
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.35, tau = 3), seed = seed)[[1L]]$choice
  list(form = bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1),
       fam = if (is.null(fam)) bandit2arm_delta(subject = id,
                                                trial = trial) else fam,
       data = d)
}

## ---- 1. scaling in subjects and trials ----------------------------
grid <- list(c(10, 100), c(40, 25), c(40, 100), c(40, 400), c(160, 100),
             c(40, 800))
cases <- lapply(grid, function(g) delta_case(g[1], g[2]))
scal <- data.frame(
  subjects = vapply(grid, function(g) g[1], numeric(1)),
  trials = vapply(grid, function(g) g[2], numeric(1)),
  rows = vapply(cases, function(cs) nrow(cs$data), numeric(1)),
  build_s = round(interleaved_builds(cases), 3),
  grad_ms = round(vapply(cases, function(cs)
    bench_grad(cs$form, cs$fam, cs$data), numeric(1)), 1))
cat("
== scaling, bandit2arm_delta, one random intercept ==
")
print(scal, row.names = FALSE)

## ---- 2. the elementwise penalty, at four sizes ---------------------
# The engine loops over TRIALS and vectorizes over SUBJECTS. The
# alternative is one iteration per ROW with the value store updated by
# sub-assignment, which is how the model is usually stated. RTMB's
# replacement operator copies the vector it writes into, so the penalty
# grows with the row count. The family below differs from
# bandit2arm_delta() in the loop shape and in nothing else.
ln_slow_loglik <- function(y, dpars, aterms, weights, block, extra) {
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  n <- block[["n"]]
  al <- dpars[["alpha"]]
  tu <- dpars[["tau"]]
  if (length(al) == 1L) al <- al * rep(1, n)
  if (length(tu) == 1L) tu <- tu * rep(1, n)
  p1 <- aterms[["reward1"]]
  p2 <- aterms[["reward2"]]
  q1 <- rep(0, nrow(idx))
  q2 <- q1
  ll <- rep(0, n)
  for (t in seq_len(ncol(idx))) {
    for (s in seq_len(nrow(idx))) {
      if (msk[s, t] == 0) next
      i <- idx[s, t]
      e1 <- tu[i] * q1[s]
      e2 <- tu[i] * q2[s]
      lse <- RTMB::logspace_add(e1, e2)
      ll[i] <- (if (y[i] == 1) e1 else e2) - lse
      if (y[i] == 1) {
        q1[s] <- q1[s] + al[i] * (p1[i] - q1[s])
      } else {
        q2[s] <- q2[s] + al[i] * (p2[i] - q2[s])
      }
    }
  }
  sum(ll)
}

slow_family <- function() {
  f <- bandit2arm_delta(subject = id, trial = trial)
  st <- f$structure
  st$loglik <- frmtmb::frmtmb_ad_overload(ln_slow_loglik)
  f$structure <- st
  f
}

sizes <- list(c(10, 100), c(40, 100), c(50, 200), c(50, 400))
fastc <- lapply(seq_along(sizes), function(j)
  delta_case(sizes[[j]][1], sizes[[j]][2], seed = 10 + j))
slowc <- lapply(fastc, function(cs)
  list(form = cs$form, fam = slow_family(), data = cs$data))
both <- c(fastc, slowc)
bt <- interleaved_builds(both, reps = 3L)
nf <- length(fastc)
vals <- vapply(seq_len(nf), function(j) {
  a <- frm(fastc[[j]]$form, family = fastc[[j]]$fam, data = fastc[[j]]$data,
           dry_run = "objective")
  b <- frm(slowc[[j]]$form, family = slowc[[j]]$fam, data = slowc[[j]]$data,
           dry_run = "objective")
  abs(a$obj$fn(a$obj$par) - b$obj$fn(b$obj$par))
}, numeric(1))
ew <- data.frame(
  rows = vapply(sizes, function(g) g[1] * g[2], numeric(1)),
  vectorized_s = round(bt[seq_len(nf)], 3),
  elementwise_s = round(bt[nf + seq_len(nf)], 3),
  factor = round(bt[nf + seq_len(nf)] / bt[seq_len(nf)], 1),
  value_gap = signif(vals, 3))
cat("
== the elementwise penalty, paid at tape build ==
")
print(ew, row.names = FALSE)

## ---- 3. and NOT paid per gradient ---------------------------------
cs <- delta_case(40, 100, seed = 99)
sc <- list(form = cs$form, fam = slow_family(), data = cs$data)
gf <- bench_grad(cs$form, cs$fam, cs$data)
gs <- bench_grad(sc$form, sc$fam, sc$data)
cat("
== one gradient, 4000 rows, ms ==
")
print(data.frame(shape = c("vectorized", "elementwise"),
                 grad_ms = round(c(gf, gs), 1)), row.names = FALSE)
cat("Once the tape exists it is a node list and the R code that built",
    "it is gone,
so this column is expected to be flat.
")

## ---- 4. every family at one size ----------------------------------
fam_case <- function(nm) {
  switch(nm,
    bandit2arm_delta = delta_case(40, 100, seed = 21),
    bandit2arm_dual = {
      cs <- delta_case(40, 100, seed = 22)
      list(form = bf(choice | reward(pay1, pay2) ~ 1 + (1 | id),
                     Apun ~ 1, tau ~ 1),
           fam = bandit2arm_dual(subject = id, trial = trial),
           data = cs$data)
    },
    prl_fictitious = {
      cs <- delta_case(40, 100, seed = 23)
      list(form = bf(choice | reward(pay1, pay2) ~ 1 + (1 | id),
                     bias ~ 1, tau ~ 1),
           fam = prl_fictitious(subject = id, trial = trial),
           data = cs$data)
    },
    bandit4arm2_kalman_filter = {
      d <- frm_task_design("bandit4arm_restless", n_subject = 40,
                           n_trial = 100, seed = 24)
      fam <- bandit4arm2_kalman_filter(subject = id, trial = trial)
      d$choice <- frm_task_simulate(fam, d, pars = list(
        tau = 0.15, lambda = 0.98, center = 50, mu0 = 50, sigma0 = 10,
        sigmaD = 3), seed = 24)[[1L]]$choice
      list(form = bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1 + (1 | id),
                     lambda ~ 1, center ~ 1, mu0 ~ 1, sigma0 ~ 1,
                     sigmaD ~ 1), fam = fam, data = d)
    },
    igt_pvl_delta = {
      d <- frm_task_design("igt", n_subject = 40, n_trial = 100, seed = 25)
      fam <- igt_pvl_delta(subject = id, trial = trial)
      d$choice <- frm_task_simulate(fam, d, pars = list(
        alpha = 0.3, shape = 0.4, lambda = 1.5, tau = 1),
        seed = 25)[[1L]]$choice
      list(form = bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1 + (1 | id),
                     shape ~ 1, lambda ~ 1, tau ~ 1), fam = fam, data = d)
    },
    ts_par7 = {
      d <- frm_task_design("twostep", n_subject = 40, n_trial = 100,
                           seed = 26)
      fam <- ts_par7(subject = id, trial = trial)
      d <- frm_task_simulate(fam, d, pars = list(
        w = 0.5, alpha1 = 0.4, tau1 = 3, alpha2 = 0.4, tau2 = 3,
        lambda = 0.6, pers = 0.2), seed = 26)[[1L]]
      list(form = bf(choice | stage2(state2, choice2) +
                       payoff(pay1, pay2, pay3, pay4) ~ 1 + (1 | id),
                     alpha1 ~ 1, tau1 ~ 1, alpha2 ~ 1, tau2 ~ 1,
                     lambda ~ 1, pers ~ 1), fam = fam, data = d)
    })
}
fnames <- frm_learn_families()$family
fc <- lapply(fnames, fam_case)
slots <- vapply(fc, function(cs) {
  d1 <- list(mu0 = 0, sigma0 = 1)
  length(cs$fam$learn$spec$init(1L, d1))
}, numeric(1))
per_fam <- data.frame(
  family = fnames,
  state_slots = slots,
  decisions = vapply(fc, function(cs)
    length(cs$fam$learn$spec$n_option), numeric(1)),
  build_s = round(interleaved_builds(fc), 3),
  grad_ms = round(vapply(fc, function(cs)
    bench_grad(cs$form, cs$fam, cs$data), numeric(1)), 1))
cat("
== every family, 40 subjects by 100 trials, 4000 rows ==
")
print(per_fam, row.names = FALSE)
cat("
The cost tracks the number of value stores the rule carries and",
    "the
number of decisions in a trial, which is what the engine",
    "walks.
")
