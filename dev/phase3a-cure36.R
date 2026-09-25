# Item 3.6, punch round 1: cure-fraction designs with a time-varying
# effect. The review measured the widened check firing on 15 of 80
# `gamma1 ~ x` fits at df 3 to 6; this is the lane's own construction of
# the same kind of design, to state the behavior in NEWS and the Rd from
# a run rather than from the review.
#
# Two arms of n/2, a cured fraction that never has the event (40 percent
# in arm 0, 55 percent in arm 1), Weibull 1.3 event times for the rest
# with the arm shifting the scale, administrative censoring at 6 and
# uniform dropout. Fitted as `t | cens(censored) ~ x, gamma1 ~ x` at df
# 3 to 6, n = 500, seeds 20260924 to 20260943 (20 per df, 80 fits).
# For each fit: whether rp_floored() fires, on which rows, whether they
# lie past the last event of their own arm, and how much the fitted
# survival RISES across them (the largest S(t_later) - S(t_earlier)
# over the flagged rows of an arm, from the fitted eta).
Sys.setenv(PHASE3A_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})
phase3a_where("frmtmb.spline")
# What rp_floored(action = "error") does on a fit, measured by calling
# it: "refuse", "warn" or "silent". Punch round 2 split the two.
gate_of <- function(f) {
  w <- FALSE
  e <- tryCatch(withCallingHandlers(rp_floored(f), warning = function(c) {
    w <<- TRUE
    invokeRestart("muffleWarning")
  }), error = function(err) "refuse")
  if (identical(e, "refuse")) "refuse" else if (w) "warn" else "silent"
}

sim_cure <- function(seed, n = 500L) {
  set.seed(seed)
  x <- rep(0:1, each = n / 2)
  cured <- stats::runif(n) < ifelse(x == 1, 0.55, 0.40)
  t <- stats::rweibull(n, 1.3, ifelse(x == 1, 1.6, 1))
  t[cured] <- Inf
  cc <- pmin(6, stats::runif(n, 0, 12))
  data.frame(t = pmin(t, cc), censored = as.integer(t > cc),
             x = factor(x))
}
out <- NULL
for (df in 3:6) {
  for (s in 20260924L + 0:19) {
    d <- sim_cure(s)
    f <- tryCatch(suppressWarnings(frm(
      bf(t | cens(censored) ~ x, gamma1 ~ x),
      family = royston_parmar(df = df), data = d)), error = identity)
    if (inherits(f, "error")) {
      out <- rbind(out, data.frame(df = df, seed = s, fit = "error",
                                   fires = NA, n_rows = NA,
                                   past_last_event = NA, max_rise = NA,
                                   gate = NA, reported_rise = NA))
      next
    }
    r <- rp_floored(f, action = "report")
    rows <- attr(r, "rows")[["nonmonotone_censored"]]
    fx <- frmtmb.spline:::sp_rp_fitted(f, stats::family(f))
    S <- exp(-exp(fx$eta))
    past <- NA; rise <- NA
    if (length(rows)) {
      last_ev <- tapply(d$t[d$censored == 0], d$x[d$censored == 0], max)
      past <- all(d$t[rows] > last_ev[as.character(d$x[rows])])
      rise <- max(vapply(split(rows, d$x[rows]), function(i) {
        i <- i[order(d$t[i])]
        if (length(i) < 2L) return(0)
        max(S[i]) - S[i[1L]]
      }, 0))
    }
    out <- rbind(out, data.frame(
      df = df, seed = s, fit = as.character(f$opt$convergence),
      fires = r[["n_nonmonotone"]] + r[["n_nonmonotone_censored"]] > 0,
      gate = gate_of(f), reported_rise = r[["max_survival_rise"]],
      n_rows = r[["n_nonmonotone_censored"]], past_last_event = past,
      max_rise = rise))
  }
}
cat("fits", nrow(out), " fit errors", sum(out$fit == "error"), "\n")
print(stats::aggregate(cbind(fits = 1, fires = fires) ~ df, data = out,
                       FUN = sum))
fired <- out[isTRUE(TRUE) & out$fires %in% TRUE, ]
cat(sprintf("fired %d of %d; flagged rows all past their arm's last event: %d of %d\n",
            nrow(fired), sum(!is.na(out$fires)),
            sum(fired$past_last_event), nrow(fired)))
if (nrow(fired)) {
  cat(sprintf("rise in fitted S across the flagged rows: %.2e to %.2e\n",
              min(fired$max_rise), max(fired$max_rise)))
}
print(fired, row.names = FALSE)
cat("\nrp_floored() on the 80 fits:\n")
print(table(out$gate, useNA = "ifany"))
cat(sprintf("reported max_survival_rise on the warning fits: %.2e to %.2e\n",
            min(out$reported_rise[out$gate %in% "warn"]),
            max(out$reported_rise[out$gate %in% "warn"])))
