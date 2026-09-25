# Item 3.1: frm_ode_records() against rxode2 on the same schedule.
#
# A. Expansion: every schedule is built with rxode2::et(), turned into
#    records with as.data.frame(), and split by frm_ode_records(). The
#    dose rows that frm_ode_records() produces, with ii/addl written
#    out, must equal the dose rows of rxode2::etExpand() on the same
#    table: same id, time, compartment, amount and duration. Compared
#    with identical() on the sorted table.
# B. Trajectories: rxode2::rxSolve() on the et() table against
#    frm_ode() and frm_lincmt() on frm_ode_records()' output, same
#    parameters, at the observation times. atol = rtol = 1e-12 on both
#    integrators; reported as the largest absolute difference over the
#    largest absolute state, per schedule.
# C. Real records: nlmixr2data::warfarin and nlmixr2data::theo_md, the
#    latter with rxode2's classic evid 101 rewritten as evid 1 and each
#    time-0 sample moved before its dose (see D), solved at fixed
#    typical parameters by rxSolve on the records and by frm_lincmt()
#    and frm_ode() on the split.
# D. The ordering premise: an observation at a dose time, listed before
#    and after the dose, in rxode2.
#
# Seeds: random schedules use 20260923 + k, k = 1..200.
Sys.setenv(PHASE3A_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.ode); library(rxode2)
})
phase3a_where("frmtmb.ode"); phase3a_where("rxode2"); phase3a_where("RTMBode")

mod <- rxode2({
  d/dt(depot) = -ka * depot
  d/dt(central) = ka * depot - ke * central
})
pk_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
}
KA <- 1.1; KE <- 0.23
TOL <- 1e-12

# dose rows of etExpand(), in frm_ode's vocabulary
rx_doses <- function(e) {
  x <- as.data.frame(rxode2::etExpand(e))
  x <- x[x$evid != 0, , drop = FALSE]
  if (is.null(x$id)) x$id <- 1L
  if (is.null(x$rate)) x$rate <- 0
  if (is.null(x$dur)) x$dur <- 0
  x$rate[is.na(x$rate)] <- 0; x$dur[is.na(x$dur)] <- 0
  dur <- ifelse(x$rate > 0, x$amt / x$rate, x$dur)
  # evid 4 is a reset and then a dose: split it the way frm_ode reads it
  i4 <- which(x$evid == 4)
  if (length(i4)) {
    r <- x[i4, , drop = FALSE]; r$evid <- 3
    x$evid[i4] <- 1
    x <- rbind(x, r)
    dur <- c(dur, rep(0, length(i4)))
  }
  out <- data.frame(id = as.character(x$id), time = x$time,
                    evid = x$evid,
                    cmt = ifelse(x$evid == 3, NA, as.character(x$cmt)),
                    amt = ifelse(x$evid == 3, 0, x$amt),
                    dur = ifelse(x$evid == 3, 0, dur),
                    stringsAsFactors = FALSE)
  out[order(out$id, out$time, out$evid), , drop = FALSE]
}
# frm_ode_records()' events, written out over addl, same vocabulary
fr_doses <- function(ev) {
  if (is.null(ev$group)) ev$group <- 1L
  rows <- lapply(seq_len(nrow(ev)), function(i) {
    k <- 0:ev$addl[i]
    data.frame(id = as.character(ev$group[i]),
               time = ev$time[i] + k * ev$ii[i],
               evid = c(add = 1, reset = 3, replace = 5,
                        multiply = 6)[[ev$method[i]]],
               cmt = as.character(ev$state[i]),
               amt = ev$value[i], dur = ev$duration[i],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$id, out$time, out$evid), , drop = FALSE]
}
same_doses <- function(a, b) {
  rownames(a) <- NULL; rownames(b) <- NULL
  isTRUE(all.equal(a, b, tolerance = 0))
}

solve_all <- function(e, ss_rows = FALSE) {
  recs <- as.data.frame(e)
  x <- frm_ode_records(recs)
  obs <- x$data
  grp <- if (is.null(obs$id)) NULL else obs$id
  rx <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), e, atol = TOL,
                              rtol = TOL, maxsteps = 1e6))
  fo <- frm_ode(pk_dyn, init = list(0, 0), times = obs$time,
                parms = list(KA, KE), group = grp,
                states = c("depot", "central"), events = x$events,
                atol = TOL, rtol = TOL)
  fl <- tryCatch(frm_lincmt(parms = list(ka = KA, ke = KE, V = 1), times = obs$time,
                   group = grp, ncmt = 1, depot = TRUE, output = "central",
                   events = x$events), error = function(err) NA_real_)
  sc <- max(abs(c(rx$depot, rx$central)))
  list(n_obs = nrow(obs), n_ev = nrow(x$events),
       ode_depot = max(abs(fo[, 1] - rx$depot)) / sc,
       ode_central = max(abs(fo[, 2] - rx$central)) / sc,
       lin_central = max(abs(fl - rx$central)) / sc,
       expand = tryCatch(same_doses(fr_doses(x$events), rx_doses(e)),
                         error = function(err) NA))
}

report <- function(label, r) {
  cat(sprintf(paste0("%-28s obs %3d events %3d  expand %-5s  ",
                     "ode depot %.2e  ode central %.2e  lincmt %.2e\n"),
              label, r$n_obs, r$n_ev, format(r$expand), r$ode_depot,
              r$ode_central, r$lin_central))
}

cat("\n== A and B: hand-built schedules ==\n")
obs_t <- c(0.3, 1.7, 5.5, 11.2, 13.4, 25.1, 37.9, 49.3, 61.6, 70.2)
sched <- list(
  "oral q12 x4 (ii/addl)" =
    et(amt = 100, ii = 12, addl = 3, cmt = "depot") |> et(obs_t),
  "IV infusion by rate" =
    et(amt = 80, rate = 40, cmt = "central", ii = 24, addl = 2) |>
    et(obs_t),
  "IV infusion by dur" =
    et(amt = 60, dur = 3.3, cmt = "central", ii = 12, addl = 4) |>
    et(obs_t),
  "oral + IV, two ids" =
    et(amt = 100, cmt = "depot", ii = 12, addl = 2) |>
    et(time = 30, amt = 50, cmt = "central") |> et(obs_t) |> et(id = 1:2),
  "reset (evid 3)" =
    et(amt = 100, cmt = "depot", ii = 12, addl = 4) |>
    et(time = 30, evid = 3) |> et(obs_t),
  "reset and dose (evid 4)" =
    et(amt = 100, cmt = "depot", ii = 12, addl = 1) |>
    et(time = 30, evid = 4, amt = 70, cmt = "depot") |> et(obs_t),
  "replace and multiply (5, 6)" =
    et(amt = 100, cmt = "depot") |>
    et(time = 5, amt = 30, evid = 5, cmt = "central") |>
    et(time = 10, amt = 0.5, evid = 6, cmt = "central") |> et(obs_t),
  "steady state oral q24" =
    et(amt = 100, cmt = "depot", ii = 24, ss = 1) |>
    et(time = 48, amt = 100, cmt = "depot") |> et(obs_t),
  "steady state infusion q12" =
    et(amt = 100, rate = 50, cmt = "central", ii = 12, ss = 1) |>
    et(obs_t)
)
for (nm in names(sched)) report(nm, solve_all(sched[[nm]]))

# Both integrators against the closed form on the first schedule, so
# that an agreement between them to the last bit is read as two LSODA
# runs taking the same steps rather than as one number compared with
# itself.
local({
  e <- sched[[1L]]
  x <- frm_ode_records(as.data.frame(e))
  tt <- x$data$time
  exact <- vapply(tt, function(t) {
    td <- c(0, 12, 24, 36)
    u <- t - td[td <= t]
    sum(100 * KA / (KA - KE) * (exp(-KE * u) - exp(-KA * u)))
  }, 0)
  rx <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), e, atol = TOL,
                              rtol = TOL))
  fo <- frm_ode(pk_dyn, init = list(0, 0), times = tt,
                parms = list(KA, KE), states = c("depot", "central"),
                events = x$events, atol = TOL, rtol = TOL)
  sc <- max(abs(exact))
  cat(sprintf(paste0("closed form, first schedule: rxode2 %.2e, frm_ode %.2e ",
                     "of the scale; rxode2 and frm_ode identical: %s\n"),
              max(abs(rx$central - exact)) / sc,
              max(abs(fo[, 2] - exact)) / sc,
              identical(rx$central, unname(fo[, 2]))))
})

cat("\n== A and B: 200 random schedules ==\n")
rand_sched <- function(seed) {
  set.seed(seed)
  n_id <- sample(1:3, 1)
  e <- NULL
  for (id in seq_len(n_id)) {
    nd <- sample(1:5, 1)
    ti <- sort(round(stats::runif(nd, 0, 48), 3))
    ei <- et(id = id)
    for (k in seq_len(nd)) {
      cmt <- sample(c("depot", "central"), 1)
      amt <- round(stats::runif(1, 10, 200), 2)
      kind <- if (cmt == "central") sample(c("bolus", "rate", "dur"), 1) else
        "bolus"
      rep <- stats::runif(1) < 0.3
      args <- list(time = ti[k], amt = amt, cmt = cmt)
      if (kind == "rate") args$rate <- round(stats::runif(1, 5, 60), 2)
      if (kind == "dur") args$dur <- round(stats::runif(1, 0.5, 4), 2)
      if (rep) { args$ii <- sample(c(6, 8, 12), 1); args$addl <- sample(1:3, 1) }
      ei <- do.call(et, c(list(ei), args))
    }
    if (stats::runif(1) < 0.2) ei <- et(ei, time = round(stats::runif(1, 20, 60), 3),
                                        evid = 3)
    ei <- et(ei, sort(round(stats::runif(12, 0.1, 80), 4)))
    e <- if (is.null(e)) ei else rbind(e, ei)
  }
  e
}
res <- NULL
refused <- NULL
for (k in 1:200) {
  e <- rand_sched(20260923L + k)
  r <- tryCatch(solve_all(e), error = function(err) conditionMessage(err))
  if (is.character(r)) {
    # the one refusal a random schedule can reach: a reset inside a
    # running infusion. What rxode2 does there is recorded beside it.
    rx <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), e, atol = TOL,
                                rtol = TOL, maxsteps = 1e6))
    refused <- rbind(refused, data.frame(
      seed = 20260923L + k,
      reset_in_infusion = grepl("while the infusion from record", r),
      rx_min_amount = min(c(rx$depot, rx$central))))
  } else {
    res <- rbind(res, cbind(seed = 20260923L + k, as.data.frame(r)))
  }
}
cat(sprintf("schedules solved %d of 200, refused %d (reset inside a running infusion: %d)\n",
            nrow(res), nrow(refused), sum(refused$reset_in_infusion)))
if (!is.null(refused)) {
  cat("rxode2's smallest amount on each refused schedule:\n")
  print(refused, row.names = FALSE)
}
cat(sprintf("expansion identical to etExpand() %d of %d (etExpand() errors on %d)\n",
            sum(res$expand, na.rm = TRUE), sum(!is.na(res$expand)),
            sum(is.na(res$expand))))
cat(sprintf("worst relative difference: ode depot %.2e, ode central %.2e, lincmt %.2e\n",
            max(res$ode_depot), max(res$ode_central), max(res$lin_central)))
cat(sprintf("median relative difference: ode central %.2e, lincmt %.2e\n",
            stats::median(res$ode_central), stats::median(res$lin_central)))
cat(sprintf("ode central bitwise equal to rxode2 on %d of %d\n",
            sum(res$ode_central == 0), nrow(res)))
cat(sprintf("observations %d, dose events %d over the solved schedules\n",
            sum(res$n_obs), sum(res$n_ev)))

cat("\n== C: real records ==\n")
real_one <- function(label, recs, id_col) {
  # rxSolve on the records themselves; frm_* on the split
  rx <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), recs,
                              atol = TOL, rtol = TOL, maxsteps = 1e6))
  x <- frm_ode_records(recs)
  obs <- x$data
  g <- obs[[id_col]]
  fo <- frm_ode(pk_dyn, init = list(0, 0), times = obs[[grep("^time$",
                names(obs), ignore.case = TRUE)]], parms = list(KA, KE),
                group = g, events = x$events, atol = TOL, rtol = TOL)
  fl <- frm_lincmt(parms = list(ka = KA, ke = KE, V = 1),
                   times = obs[[grep("^time$", names(obs),
                                     ignore.case = TRUE)]],
                   group = g, ncmt = 1, depot = TRUE, output = "central",
                   events = x$events)
  sc <- max(abs(c(rx$depot, rx$central)))
  cat(sprintf(paste0("%-10s records %4d  obs %4d (rxSolve %4d)  events %3d  ",
                     "ode central %.2e  lincmt %.2e\n"),
              label, nrow(recs), nrow(obs), nrow(rx), nrow(x$events),
              max(abs(fo[, 2] - rx$central)) / sc,
              max(abs(fl - rx$central)) / sc))
}
w <- nlmixr2data::warfarin
# warfarin has no cmt column; rxode2 then doses the first compartment,
# the depot, and frm_ode() refuses a dose with no state on a two-state
# system, so the depot is named
w$cmt <- ifelse(w$evid == 1, 1L, 2L)
# warfarin carries two endpoints per time (dvid cp and pca); rxSolve
# returns one row per observation record, so both are kept
chk <- tryCatch(frm_ode_records(w), error = function(e) conditionMessage(e))
cat("warfarin as shipped:", substr(chk, 1, 110), "...\n")
w2 <- w[order(w$id, w$time, w$evid != 0), ]
real_one("warfarin", w2, "id")
tm <- nlmixr2data::theo_md
tm$EVID[tm$EVID == 101] <- 1
# the time-0 sample follows its dose in the file; in NONMEM order that
# is a post-dose sample, which frm_ode_records() refuses by name
chk <- tryCatch(frm_ode_records(tm), error = function(e) conditionMessage(e))
cat("theo_md as shipped:", substr(chk, 1, 110), "...\n")
o <- order(tm$ID, tm$TIME, tm$EVID != 0)
tm2 <- tm[o, ]
real_one("theo_md", tm2, "ID")

cat("\n== D: the ordering premise in rxode2 ==\n")
d1 <- data.frame(id = 1, time = c(0, 12, 12, 13), evid = c(1, 0, 1, 0),
                 amt = c(100, 0, 100, 0), cmt = c(1, 2, 1, 2))
d2 <- d1[c(1, 3, 2, 4), ]
a <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), d1, atol = TOL,
                           rtol = TOL))
b <- as.data.frame(rxSolve(mod, c(ka = KA, ke = KE), d2, atol = TOL,
                           rtol = TOL))
pre <- frm_ode(pk_dyn, init = list(0, 0), times = 12, parms = list(KA, KE),
               events = data.frame(time = c(0, 12), state = 1, value = 100),
               atol = TOL, rtol = TOL)
cat(sprintf(paste0("depot at t = 12: rxode2 obs listed BEFORE the dose %.6f, ",
                   "AFTER %.6f; frm_ode %.6f\n"),
            a$depot[1], b$depot[1], pre[1, 1]))
