# The diagnostic frm_ode() tells a user to read, against the error it
# is a diagnostic FOR.
#
# ode_run_in() compares the last two cycle-start states and warns when
# they differ by more than ss_tol. After n cycles the run-in holds
# (1 - r^n) of the limit for each mode, r = exp(-lambda ii), so the
# distance to the limit is r^n and the movement the warning reports is
# r^(n-1) (1 - r). Their ratio is (1 - r) / r, which is about
# lambda * ii when lambda * ii is small. The warning therefore
# understates the error by about 1 / (lambda * ii), and it does so
# most where the error is largest.
#
# Script path: dev/rev-lincmt-nsswarn.R.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

dyn1 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[[2]] * y[1], p[[2]] * y[1] - p[[1]] * y[2]))
}
grab <- function(expr) {
  w <- NULL
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  list(v = v, w = w)
}
rel_of <- function(msg) {
  if (is.null(msg)) return(NA_real_)
  m <- regmatches(msg, regexpr("moving by [0-9.e+-]+", msg))
  if (!length(m)) return(NA_real_)
  as.numeric(sub("moving by ", "", m[[1L]]))
}

cat("\nOne compartment with a depot, ka = 1, ii = 12, n_ss = 20.\n")
cat("`warned` is the number frm_ode()'s warning prints; `short` is the",
    "\nmeasured distance to the exact limit, frm_lincmt(n_ss = 20)",
    "against\nfrm_lincmt(n_ss = Inf); `1/(k ii)` is the predicted",
    "understatement.\n\n")
cat(sprintf("%9s %9s %10s %11s %11s %9s %9s\n", "t_half", "k*ii",
            "warned", "short", "ratio", "1/(k ii)", "fires"))
ts <- c(0, 3, 6, 9, 12)
for (kk in c(0.4, 0.15, 0.05, 0.02, 0.01, 0.005)) {
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                   addl = 0L, ss = TRUE)
  g <- grab(frm_ode(dyn1, init = list(0, 0), times = ts,
                    parms = list(kk, 1.0),
                    states = c("depot", "central"),
                    output = "central", events = ev, n_ss = 20L,
                    atol = 1e-12, rtol = 1e-12))
  a <- frm_lincmt(parms = list(ke = kk, ka = 1.0, V = 1), times = ts,
                  ncmt = 1, depot = TRUE, events = ev, n_ss = 20L)
  b <- frm_lincmt(parms = list(ke = kk, ka = 1.0, V = 1), times = ts,
                  ncmt = 1, depot = TRUE, events = ev)
  sh <- max(abs(a - b)) / max(abs(b))
  wr <- rel_of(g$w)
  cat(sprintf("%9.1f %9.4f %10.4g %11.4g %11.1f %9.1f %9s\n",
              log(2) / kk, kk * 12, wr, sh,
              if (is.na(wr)) NA else sh / wr, 1 / (kk * 12),
              !is.null(g$w)))
}

cat("\nWhat ?frm_ode says today, quoted from R/ode.R:\n")
cat("  \"At the default n_ss = 20 that is 1e-21 for a drug eliminated",
    "\n   over its dosing interval and only a percent or two for one",
    "\n   whose half-life is many intervals long.\"\n")
cat("\nA half-life of ten dosing intervals is k*ii = 0.0693, so the",
    "\nshortfall after 20 cycles is exp(-20 * 0.0693) =",
    format(exp(-20 * log(2) / 10), digits = 3), "\n")
cat("The sentence holds at a half-life of about three intervals, not",
    "\n\"many\".\n")
