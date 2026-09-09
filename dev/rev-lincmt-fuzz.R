# Attack 6: an accepted-and-wrong schedule, which is the failure that
# matters. The lane enumerated 23 NONMEM-shaped records by hand; this
# draws schedules at random from a grammar wider than that hand list,
# including combinations the lane's 118 do not contain (two steady-
# state rows in one group, a reset at the same instant as a steady-
# state row, a steady-state row at a t0 that is not zero, a dose after
# the last observation, an observation many intervals after the last
# dose), and compares against frm_ode() at atol = rtol = 1e-12 with a
# MATCHED n_ss.
#
# Seed 20260909. Script path: dev/rev-lincmt-fuzz.R.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

NSS <- 20L
TOL <- 1e-12

mk_dyn <- function(ncmt, depot) function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  ke <- p[[1]]; k12 <- p[[2]]; k21 <- p[[3]]
  k13 <- p[[4]]; k31 <- p[[5]]; ka <- p[[6]]
  ic <- if (depot) 2L else 1L
  d <- numeric(ncmt + as.integer(depot))
  cen <- -(ke + (if (ncmt >= 2L) k12 else 0) +
             (if (ncmt == 3L) k13 else 0)) * y[ic]
  if (depot) {
    d <- c(-ka * y[1], ka * y[1] + cen)
  } else {
    d <- c(cen)
  }
  if (ncmt >= 2L) {
    d[ic] <- d[ic] + k21 * y[ic + 1L]
    d <- c(d, k12 * y[ic] - k21 * y[ic + 1L])
  }
  if (ncmt == 3L) {
    d[ic] <- d[ic] + k31 * y[ic + 2L]
    d <- c(d, k13 * y[ic] - k31 * y[ic + 2L])
  }
  list(d)
}

states_of <- function(ncmt, depot)
  c(if (depot) "depot", "central",
    if (ncmt >= 2L) "peripheral1", if (ncmt == 3L) "peripheral2")

set.seed(20260909)
res <- list()
for (case in 1:140) {
  ncmt <- sample(1:3, 1L)
  depot <- sample(c(TRUE, FALSE), 1L)
  t0 <- sample(c(0, 0, 2.5), 1L)
  p <- list(ke = exp(runif(1, log(0.05), log(1.5))),
            k12 = exp(runif(1, log(0.02), log(2))),
            k21 = exp(runif(1, log(0.02), log(2))),
            k13 = exp(runif(1, log(0.02), log(2))),
            k31 = exp(runif(1, log(0.02), log(2))),
            ka = exp(runif(1, log(0.1), log(3))))
  pl <- p[c(TRUE, ncmt >= 2, ncmt >= 2, ncmt == 3, ncmt == 3, depot)]
  pl$V <- exp(runif(1, log(2), log(40)))
  out <- if (depot) sample(c("conc", "central", "depot"), 1L) else
    sample(c("conc", "central"), 1L)
  times <- sort(t0 + c(0, sort(runif(7, 0, 60))))
  # rows
  rows <- list()
  nadd <- sample(0:3, 1L)
  st_pool <- if (depot) c("depot", "central") else "central"
  for (j in seq_len(nadd)) {
    stt <- sample(st_pool, 1L)
    dur <- if (stt == "central" && runif(1) < 0.35)
      round(runif(1, 0.2, 3), 3) else 0
    hasii <- runif(1) < 0.5
    ii <- if (hasii) round(runif(1, max(dur, 3), 14), 3) else 0
    rows[[length(rows) + 1L]] <- data.frame(
      time = round(t0 + runif(1, 0, 30), 3), state = stt,
      value = round(runif(1, 20, 200), 2),
      method = "add", duration = dur, ii = ii,
      addl = if (hasii) sample(0:4, 1L) else 0L, ss = FALSE)
  }
  nss <- sample(c(0L, 0L, 1L, 2L), 1L)
  for (j in seq_len(nss)) {
    stt <- sample(st_pool, 1L)
    ii <- round(runif(1, 4, 14), 3)
    dur <- if (stt == "central" && runif(1) < 0.3)
      round(runif(1, 0.2, min(3, ii)), 3) else 0
    rows[[length(rows) + 1L]] <- data.frame(
      time = if (j == 1L) t0 else round(t0 + runif(1, 5, 35), 3),
      state = stt, value = round(runif(1, 20, 200), 2),
      method = "add", duration = dur, ii = ii, addl = 0L, ss = TRUE)
  }
  if (runif(1) < 0.25) {
    rows[[length(rows) + 1L]] <- data.frame(
      time = if (runif(1) < 0.3 && nss > 0L) rows[[length(rows)]]$time
             else round(t0 + runif(1, 1, 40), 3),
      state = NA_character_, value = 0, method = "reset",
      duration = 0, ii = 0, addl = 0L, ss = FALSE)
  }
  ev <- if (length(rows)) do.call(rbind, rows) else NULL
  ini <- if (runif(1) < 0.3)
    list(depot = if (depot) round(runif(1, 10, 100), 2) else NULL,
         central = round(runif(1, 10, 100), 2)) else NULL
  if (!is.null(ini)) ini <- ini[!vapply(ini, is.null, TRUE)]
  esc <- if (!is.null(ev) && runif(1) < 0.3) round(runif(1, 0.4, 1.2), 3)
         else 1

  a <- tryCatch(frm_lincmt(parms = pl, times = times, ncmt = ncmt,
                           depot = depot, output = out, t0 = t0,
                           init = ini, events = ev, event_scale = esc,
                           n_ss = NSS),
                error = function(e) structure(conditionMessage(e),
                                              class = "refused"))
  ini_full <- as.list(rep(0, ncmt + as.integer(depot)))
  nmst <- states_of(ncmt, depot)
  if (!is.null(ini)) for (nm in names(ini))
    ini_full[[match(nm, nmst)]] <- ini[[nm]]
  b <- tryCatch(frm_ode(mk_dyn(ncmt, depot), init = ini_full,
                        times = times,
                        parms = list(p$ke, p$k12, p$k21, p$k13, p$k31,
                                     p$ka),
                        states = nmst,
                        output = if (out == "conc") "central" else out,
                        t0 = t0, events = ev, event_scale = esc,
                        n_ss = NSS, atol = TOL, rtol = TOL),
                error = function(e) structure(conditionMessage(e),
                                              class = "refused"))
  if (inherits(a, "refused") || inherits(b, "refused")) {
    res[[length(res) + 1L]] <- list(case = case, kind = "refused",
      lin = inherits(a, "refused"), ode = inherits(b, "refused"),
      msg = if (inherits(a, "refused")) a else b)
    next
  }
  b <- if (out == "conc") b / pl$V else b
  sc <- max(abs(b))
  e <- if (sc > 0) max(abs(a - b)) / sc else max(abs(a - b))
  res[[length(res) + 1L]] <- list(case = case, kind = "ok", err = e,
    ncmt = ncmt, depot = depot, out = out, nss = nss, t0 = t0,
    nadd = nadd, scale = sc)
}

ok <- Filter(function(z) z$kind == "ok", res)
rf <- Filter(function(z) z$kind == "refused", res)
cat("\ncases", length(res), " compared", length(ok),
    " refused by one side", length(rf), "\n")
errs <- vapply(ok, function(z) z$err, 0)
cat("worst difference / trajectory scale:",
    format(max(errs), digits = 4), "\n")
cat("quantiles:", paste(format(quantile(errs, c(0.5, 0.9, 0.99, 1)),
                               digits = 3), collapse = "  "), "\n")
o <- order(errs, decreasing = TRUE)[1:8]
cat("\nworst eight:\n")
for (i in o) {
  z <- ok[[i]]
  cat(sprintf("  case %3d  ncmt %d depot %-5s out %-7s ss %d t0 %.1f",
              z$case, z$ncmt, z$depot, z$out, z$nss, z$t0),
      sprintf(" nadd %d  err %.3e  scale %.3g\n", z$nadd, z$err,
              z$scale))
}
cat("\nrefusals, by side:\n")
lin_only <- Filter(function(z) z$lin && !z$ode, rf)
ode_only <- Filter(function(z) !z$lin && z$ode, rf)
both <- Filter(function(z) z$lin && z$ode, rf)
cat("  frm_lincmt() refused, frm_ode() accepted:", length(lin_only),
    "\n")
cat("  frm_ode() refused, frm_lincmt() accepted:", length(ode_only),
    "\n")
cat("  both refused:", length(both), "\n")
if (length(ode_only)) for (z in ode_only)
  cat("    ACCEPTED BY THE CLOSED FORM ONLY, case", z$case, ":",
      substr(z$msg, 1, 90), "\n")
msgs <- vapply(lin_only, function(z) substr(z$msg, 1, 58), "")
for (m in unique(msgs))
  cat(sprintf("    %3d x %s\n", sum(msgs == m), m))
