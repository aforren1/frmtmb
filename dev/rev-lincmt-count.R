# Attack 4: the load-independent count, which is the timing claim that
# does not depend on what else the machine was doing. Claimed: 33
# RTMBode::ode() calls per subject per objective evaluation against 0,
# of which 19 are the ss run-in.
#
# The lane counted on a NUMERIC call. That is only the same number as
# "per objective evaluation" if the tape holds one ADjoint node per
# solve and replays each once, so this script counts BOTH: the numeric
# call, and the tape's own atomic-node census.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

ode_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}
doses <- data.frame(time = c(0, 12), state = "depot", value = 100,
                    ii = c(12, 12), addl = c(0L, 12L),
                    ss = c(TRUE, FALSE))
tt <- 144 + c(0.5, 1, 2, 4, 6, 8, 10, 12)

with_counter <- function(expr) {
  n <- 0L
  ns <- asNamespace("RTMBode")
  f <- get("ode", envir = ns)
  unlockBinding("ode", ns)
  assign("ode", function(...) { n <<- n + 1L; f(...) }, envir = ns)
  on.exit({ assign("ode", f, envir = ns); lockBinding("ode", ns) })
  force(expr)
  n
}

cat("\n=== A. numeric call, one group, eight observations ===\n")
for (n_ss in c(20L, 19L, 10L, 1L)) {
  n <- with_counter(frm_ode(ode_dyn, init = list(0, 0), times = tt,
                            parms = list(1, 0.15, 20),
                            states = c("depot", "central"),
                            output = "central", events = doses,
                            n_ss = n_ss))
  cat(sprintf("  n_ss = %3d : %3d RTMBode::ode() calls\n", n_ss, n))
}

cat("\n=== B. the SAME count on the tape, at build time ===\n")
f <- function(p) sum(frm_ode(ode_dyn, init = list(0, 0), times = tt,
                             parms = list(exp(p[1]), exp(p[2]),
                                          exp(p[3])),
                             states = c("depot", "central"),
                             output = "central", events = doses,
                             n_ss = 20L))
x <- log(c(1, 0.15, 20))
nb <- with_counter(tp <- MakeTape(f, x))
cat("  RTMBode::ode() calls during MakeTape():", nb, "\n")
df <- tp$data.frame()
cat("  tape nodes:", nrow(df), "\n")
op <- table(as.character(df$OpName))
cat("  atomic/adjoint ops on the tape:\n")
print(op[grepl("Adjoint|Atomic|ODE|Tape", names(op),
               ignore.case = TRUE)])

cat("\n=== C. does a REPLAY of the tape solve again? ===\n")
n1 <- with_counter(invisible(tp(x)))
n2 <- with_counter(invisible(tp$jacobian(x)))
cat("  forward replay:", n1, " reverse (jacobian):", n2,
    "  (0 means the solve is inside the atomic, not a fresh call)\n")

cat("\n=== D. the closed form, same group and schedule ===\n")
for (nss in list(20L, Inf)) {
  g <- function(p) sum(frm_lincmt(
    parms = list(ka = exp(p[1]), ke = exp(p[2]), V = exp(p[3])),
    times = tt, ncmt = 1, depot = TRUE, events = doses, n_ss = nss))
  nc <- with_counter(tpl <- MakeTape(g, x))
  cat(sprintf("  n_ss = %-4s : %2d RTMBode::ode() calls, %6d tape",
              format(nss), nc, nrow(tpl$data.frame())), "nodes\n")
}

cat("\n=== E. the ss run-in's share, by construction ===\n")
# the same schedule with the ss row removed, so the difference is the
# run-in and nothing else
d2 <- doses[2L, , drop = FALSE]
n_noss <- with_counter(frm_ode(ode_dyn, init = list(0, 0), times = tt,
                               parms = list(1, 0.15, 20),
                               states = c("depot", "central"),
                               output = "central", events = d2,
                               n_ss = 20L))
n_full <- with_counter(frm_ode(ode_dyn, init = list(0, 0), times = tt,
                               parms = list(1, 0.15, 20),
                               states = c("depot", "central"),
                               output = "central", events = doses,
                               n_ss = 20L))
cat("  with the ss row:", n_full, "  without it:", n_noss,
    "  difference:", n_full - n_noss, "\n")
