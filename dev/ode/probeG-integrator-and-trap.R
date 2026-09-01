# Probe G: two loose ends on the per-subject-loop route.
# (1) does the integrator choice matter for robustness/cost at 2 states?
# (2) is the within-subject-covariate trap detectable from the fit?
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
library(frmtmb)
library(RTMBode)
source("dev/ode/pk-common.R")
d <- sim_pk()
set.seed(7)
drag <- do.call(rbind, lapply(levels(d$id), function(s) {
  di <- d[d$id == s, ]; di[sample(nrow(di), sample(4:8, 1)), ]
}))
drag <- drag[sample(nrow(drag)), ]

pk_ode_m <- function(ka, ke, V, time, id, dose, method) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  ids <- as.integer(id); out <- numeric(length(ids))
  for (s in unique(ids)) {
    idx <- which(ids == s); idx <- idx[order(time[idx])]
    sol <- RTMBode::ode(y = c(dose[idx[1]], 0), times = c(0, time[idx]),
                        func = pk_dyn,
                        parms = c(ka[idx[1]], ke[idx[1]], V[idx[1]]),
                        method = method, atol = 1e-8, rtol = 1e-8)
    out[idx] <- sol[-1, 3]
  }
  out
}
st <- list(beta = c(0, log(0.25), log(8)))
cat("=== integrator sweep, per-subject loop (2 states) ===\n")
for (m in c("lsoda", "lsode", "adams", "ode45", "rk4")) {
  mm <- m
  pk_m <- function(ka, ke, V, time, id, dose)
    pk_ode_m(ka, ke, V, time, id, dose, mm)
  for (tag in c("regular", "ragged")) {
    dd <- if (tag == "regular") d else drag
    w <- NULL
    t0 <- proc.time()[["elapsed"]]
    f <- withCallingHandlers(
      try(frm(bf(conc ~ pk_m(exp(lka), exp(lke), exp(lV), time, id, dose),
                 lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
                gaussian(), data = dd, start = st), silent = TRUE),
      warning = function(cd) { w <<- c(w, conditionMessage(cd))
                               invokeRestart("muffleWarning") })
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("  %-6s %-8s %6.2fs  logLik %-12s warnings: %d\n", m, tag, el,
                if (inherits(f, "try-error")) "FAILED" else
                  format(as.numeric(logLik(f)), digits = 8),
                length(w)))
    if (length(w)) cat("      first: ", sub("\n.*", "", w[1]), "\n", sep = "")
  }
}

cat("\n=== within-subject covariate trap: is it detectable? ===\n")
d$phase <- factor(ifelse(d$time <= 2, "early", "late"))
f <- frm(bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
            lka ~ 1 + (1 | id), lke ~ 1 + phase + (1 | id), lV ~ 1,
            nl = TRUE) + gaussian(), data = d, se = TRUE,
         start = list(beta = c(0, log(0.25), 0, log(8))))
print(summary(f)$coefficients$lke)
cat("\nlogLik with the ignored term:", as.numeric(logLik(f)),
    " without:", -60.462931, "\n")
cat("(a zero estimate with a huge/NaN SE and an unchanged logLik is the ",
    "only signal)\n", sep = "")
