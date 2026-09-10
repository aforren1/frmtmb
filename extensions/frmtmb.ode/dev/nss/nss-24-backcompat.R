# Back compatibility of `ss_extrapolate = FALSE`, widened.
#
# Punch round 0 checked five compartment schedules. The reviewer
# widened it to 26 against the PRE-punch build; the code has changed
# since, so it is re-established here against the shipped one. The set
# is deliberately awkward: an oral depot whose differences are exactly
# zero, a state with no steady state, flip-flop absorption, an
# oscillator, infusions, three compartments, `n_ss` of 1, 2, 3, 4 and
# 137, and tolerances from 1e-6 to 1e-12.
#
# Run with NSS_ARM=ref first to write the reference fixture.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-24-backcompat.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()
ARM <- Sys.getenv("NSS_ARM", "lane")
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss"

one <- function(t, y, p) list(c(-p[2L] * y[1L],
                                p[2L] * y[1L] - p[1L] * y[2L]))
one_auc <- function(t, y, p) list(c(-p[2L] * y[1L],
                                    p[2L] * y[1L] - p[1L] * y[2L],
                                    y[2L]))
iv1 <- function(t, y, p) list(c(-p[1L] * y[1L]))
two <- function(t, y, p) list(c(-p[4L] * y[1L],
                                p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] +
                                  p[3L] * y[3L],
                                p[2L] * y[2L] - p[3L] * y[3L]))
three <- function(t, y, p) list(c(-(p[1L] + p[2L] + p[4L]) * y[1L] +
                                    p[3L] * y[2L] + p[5L] * y[3L],
                                  p[2L] * y[1L] - p[3L] * y[2L],
                                  p[4L] * y[1L] - p[5L] * y[3L]))
osc <- function(t, y, p) list(c(y[2L],
                                -p[1L] * p[1L] * y[1L] -
                                  2 * p[2L] * p[1L] * y[2L]))

cases <- list()
add <- function(tag, dyn, ns, pv, ii, out, n, tol, dur = 0,
                addl = 0L) {
  cases[[tag]] <<- list(dyn = dyn, ns = ns, pv = pv, ii = ii,
                        out = out, n = n, tol = tol, dur = dur,
                        addl = addl)
}
P2 <- list(0.15, 0.3, 0.02, 1.0)
P2f <- list(0.0223, 0.0445, 0.0067, 0.0031)
P3 <- list(0.2, 0.4, 0.1, 0.05, 0.01)
add("1cmt oral q12 n20",        one, 2L, list(0.05, 1), 12, 2L, 20L, 1e-8)
add("1cmt oral fast depot q24", one, 2L, list(0.02, 3), 24, 2L, 20L, 1e-8)
add("1cmt oral q12 n1",         one, 2L, list(0.05, 1), 12, 2L, 1L,  1e-8)
add("1cmt oral q12 n2",         one, 2L, list(0.05, 1), 12, 2L, 2L,  1e-8)
add("1cmt oral q12 n3",         one, 2L, list(0.05, 1), 12, 2L, 3L,  1e-8)
add("1cmt oral q12 n4",         one, 2L, list(0.05, 1), 12, 2L, 4L,  1e-8)
add("1cmt oral q12 n137",       one, 2L, list(0.05, 1), 12, 2L, 137L, 1e-8)
add("1cmt oral q12 tol 1e-6",   one, 2L, list(0.05, 1), 12, 2L, 20L, 1e-6)
add("1cmt oral q12 tol 1e-12",  one, 2L, list(0.05, 1), 12, 2L, 20L, 1e-12)
add("1cmt iv bolus q8",         iv1, 1L, list(0.2),      8, 1L, 20L, 1e-8)
add("1cmt iv infusion q24",     iv1, 1L, list(0.02),    24, 1L, 20L, 1e-8, 4)
add("AUC state, read central",  one_auc, 3L, list(0.05, 1), 12, 2L, 20L, 1e-8)
add("AUC state, read the AUC",  one_auc, 3L, list(0.05, 1), 12, 3L, 20L, 1e-8)
add("2cmt oral 107h q24",       two, 3L, P2, 24, 2L, 20L, 1e-8)
add("2cmt oral 107h q12",       two, 3L, P2, 12, 2L, 20L, 1e-8)
add("2cmt oral 107h q24 n5",    two, 3L, P2, 24, 2L, 5L,  1e-8)
add("2cmt oral 107h q24 depot", two, 3L, P2, 24, 1L, 20L, 1e-8)
add("2cmt oral 107h q24 periph", two, 3L, P2, 24, 3L, 20L, 1e-8)
add("2cmt oral 107h q24 addl",  two, 3L, P2, 24, 2L, 20L, 1e-8, 0, 2L)
add("2cmt oral 107h infusion",  two, 3L, P2, 24, 2L, 20L, 1e-8, 3)
add("2cmt FLIP-FLOP q24",       two, 3L, P2f, 24, 2L, 20L, 1e-8)
add("2cmt FLIP-FLOP q24 n40",   two, 3L, P2f, 24, 2L, 40L, 1e-8)
add("3cmt iv q8",               three, 3L, P3, 8, 1L, 20L, 1e-8)
add("3cmt iv q8 infusion",      three, 3L, P3, 8, 1L, 20L, 1e-8, 2)
add("oscillator ii 1",          osc, 2L, list(3.0, 0.05), 1, 1L, 20L, 1e-8)
add("oscillator ii 3",          osc, 2L, list(1.0, 0.10), 3, 1L, 20L, 1e-8)

run <- function(cs) {
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = cs$ii,
                   ss = TRUE)
  if (cs$dur > 0) ev$duration <- cs$dur
  if (cs$addl > 0L) {
    ev <- rbind(ev[, c("time", "state", "value", "ii", "ss")],
                data.frame(time = cs$ii, state = 1L, value = 100,
                           ii = cs$ii, ss = FALSE))
    ev$addl <- c(0L, cs$addl)
  }
  args <- list(cs$dyn, init = rep(list(0), cs$ns),
               times = seq(0, cs$ii, length.out = 9), parms = cs$pv,
               events = ev, output = cs$out, n_ss = cs$n,
               ss_tol = Inf, atol = cs$tol, rtol = cs$tol)
  if (ARM != "ref") args$ss_extrapolate <- FALSE
  suppressWarnings(as.numeric(do.call(frm_ode, args)))
}
got <- lapply(cases, run)
f <- file.path(OUT, paste0("nss-24-", ARM, ".rds"))
saveRDS(got, f)
cat("wrote", basename(f), "with", length(got), "vectors\n")

ref <- file.path(OUT, "nss-24-ref.rds")
if (ARM != "ref" && file.exists(ref)) {
  old <- readRDS(ref)
  ok <- vapply(names(got), function(k) identical(got[[k]], old[[k]]),
               TRUE)
  cat(sprintf("identical() TRUE on %d of %d\n", sum(ok), length(ok)))
  d <- vapply(names(got), function(k)
    max(abs(got[[k]] - old[[k]])), 0)
  cat("largest absolute difference anywhere:", format(max(d)), "\n")
  if (any(!ok)) print(d[!ok])
}
