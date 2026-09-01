# Probe D: the edge cases an ODE grammar would have to own, plus cost
# scaling and a stacked-system alternative to the per-subject loop.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
library(frmtmb)
library(RTMB)
library(RTMBode)
source("dev/ode/pk-common.R")

ck <- function(l, e) {
  r <- try(force(e), silent = TRUE)
  cat(if (inherits(r, "try-error")) "[FAIL] " else "[ok]   ", l,
      if (inherits(r, "try-error"))
        paste0(": ", sub("\n$", "", as.character(r))) else "", "\n", sep = "")
  invisible(r)
}
nlform <- bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
             lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
st <- list(beta = c(0, log(0.25), log(8)))
fitit <- function(dd, form = nlform, ...)
  frm(form + gaussian(), data = dd, start = st, ...)

cat("=== 1. edge cases ===\n")
d <- sim_pk()

# 1a. an observation AT t = 0: grid becomes c(0, 0, ...) - duplicated
d0 <- d; d0$time[d0$time == 0.25] <- 0
d0$conc[d0$time == 0] <- 0
ck("observation at t = 0 (duplicated grid point)", {
  f <- fitit(d0); cat("       logLik ", as.numeric(logLik(f)), "\n"); TRUE })

# 1b. duplicate times within a subject (replicate assays)
ddup <- rbind(d, d[d$id == "1" & d$time == 2, ])
ddup$id <- factor(ddup$id)
ck("duplicate time within a subject", {
  f <- fitit(ddup); cat("       logLik ", as.numeric(logLik(f)), "\n"); TRUE })

# 1c. ragged / unsorted design: each subject its own random time set
set.seed(7)
drag <- do.call(rbind, lapply(levels(d$id), function(s) {
  di <- d[d$id == s, ]
  di[sample(nrow(di), sample(4:8, 1)), ]
}))
drag <- drag[sample(nrow(drag)), ]      # scramble row order too
ck("ragged, unsorted, row-shuffled design", {
  f <- fitit(drag)
  cat("       n =", nrow(drag), " logLik ", as.numeric(logLik(f)), "\n"); TRUE })

# 1d. NA rows: na.omit drops them before the body sees the data
dna <- d; dna$conc[c(3, 40, 77)] <- NA
ck("NA responses (na.omit)", {
  f <- fitit(dna)
  cat("       nobs", nobs(f), " predict length", length(predict(f)), "\n")
  TRUE })

# 1e. THE SILENT TRAP: a covariate that varies WITHIN subject. The helper
# reads each subject's parameters off its first row, so a within-subject
# covariate on a dynamics parameter is silently ignored.
d$phase <- factor(ifelse(d$time <= 2, "early", "late"))
ck("within-subject covariate on lke (expected: silently wrong)", {
  f <- frm(bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
              lka ~ 1 + (1 | id), lke ~ 1 + phase + (1 | id), lV ~ 1,
              nl = TRUE) + gaussian(), data = d,
           start = list(beta = c(0, log(0.25), 0, log(8))))
  cat("       phaselate coef:", fixef(f)$lke[[2]],
      "  (a within-subject effect the solve cannot see)\n")
  TRUE })

# 1f. REML
ck("REML = TRUE", {
  f <- fitit(d, REML = TRUE)
  cat("       REML logLik", as.numeric(logLik(f)), "\n"); TRUE })

cat("\n=== 2. cost scaling (per-subject loop) ===\n")
timing <- function(n_id, form = nlform, tag = "loop") {
  dd <- sim_pk(n_id = n_id, seed = 100 + n_id)
  t0 <- proc.time()[["elapsed"]]
  f <- frm(form + gaussian(), data = dd, start = st,
           control = frmtmb_control(verbose = FALSE))
  el <- proc.time()[["elapsed"]] - t0
  data.frame(route = tag, n_id = n_id, n = nrow(dd),
             seconds = round(el, 2), logLik = round(as.numeric(logLik(f)), 3))
}
sc <- do.call(rbind, lapply(c(6, 12, 25, 50), timing))
print(sc, row.names = FALSE)

cat("\n=== 3. ODE overhead vs the closed form ===\n")
# identical model, mu from the analytic solution: isolates what the
# adjoint ODE solve costs relative to plain taped arithmetic.
pk_closed <- function(ka, ke, V, time, dose)
  dose * ka / (V * (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
cform <- bf(conc ~ pk_closed(exp(lka), exp(lke), exp(lV), time, dose),
            lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
sc2 <- do.call(rbind, lapply(c(6, 12, 25, 50), timing, form = cform,
                             tag = "closed form"))
print(sc2, row.names = FALSE)
cat("ODE / closed-form time ratio:",
    format(sc$seconds / sc2$seconds, digits = 3), "\n")

cat("\n=== 4. stacked-system alternative ===\n")
# One ode() call for the whole population: the 2-state system is stacked
# into 2*n_subject states over the union time grid, so the tape carries a
# single adjoint ODE node instead of one per subject.
pk_dyn_stack <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  ns <- length(y) / 2
  A <- y[1:ns]; C <- y[ns + 1:ns]
  ka <- p[1:ns]; ke <- p[ns + 1:ns]; V <- p[2 * ns + 1:ns]
  list(c(-ka * A, ka * A / V - ke * C))
}
pk_ode_stack <- function(ka, ke, V, time, id, dose) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  ids <- as.integer(factor(id))
  ns <- max(ids)
  first <- match(seq_len(ns), ids)
  grid <- sort(unique(c(0, time)))
  sol <- RTMBode::ode(y = c(dose[first], rep(0, ns)), times = grid,
                      func = pk_dyn_stack,
                      parms = c(ka[first], ke[first], V[first]),
                      method = "lsoda", atol = 1e-8, rtol = 1e-8)
  # column 1 is time; state C_s lives in column 1 + ns + s
  sol[cbind(match(time, grid), 1L + ns + ids)]
}
sform <- bf(conc ~ pk_ode_stack(exp(lka), exp(lke), exp(lV), time, id, dose),
            lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
ck("stacked solve fits", {
  f <- fitit(d[, c("id", "time", "dose", "conc")], form = sform)
  cat("       logLik", as.numeric(logLik(f)),
      " (per-subject loop gave -60.462931)\n"); TRUE })
sc3 <- try(do.call(rbind, lapply(c(6, 12, 25, 50), timing, form = sform,
                                 tag = "stacked")), silent = TRUE)
if (!inherits(sc3, "try-error")) print(sc3, row.names = FALSE) else
  cat("stacked scaling FAILED:", as.character(sc3), "\n")
cat("PROBED DONE\n")
