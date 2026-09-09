# N1: three compartments with every disposition rate underflowed to
# exactly zero returned NaN in value and gradient. This checks the
# offset that removes it, and checks that the offset changes nothing
# where the rates are ordinary.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src2.R")
tt <- c(0, 1, 4, 12, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                 addl = 2L)
P0 <- list(ke = 0, k12 = 0, k21 = 0, k13 = 0, k31 = 0, ka = 1, V = 1)
v <- frm_lincmt(parms = P0, times = tt, ncmt = 3, depot = TRUE,
                events = ev)
cat("all five disposition rates exactly zero:\n  value ",
    paste(format(v, digits = 8), collapse = " "), "\n")
cat("  finite:", all(is.finite(v)), "\n")
# with no elimination and no transfer, everything absorbed stays in
# the central compartment, so the amount is the dose absorbed so far
want <- vapply(tt, function(t) {
  s <- 0
  for (d in c(0, 12, 24)) if (d < t) s <- s + 100 * (1 - exp(-(t - d)))
  s
}, 0)
cat("  against the analytic answer, max rel",
    format(max(abs(v - want) / pmax(want, 1))), "\n")
g <- RTMB::MakeTape(function(th) {
  sum(frm_lincmt(parms = list(ke = th[1], k12 = th[2], k21 = th[3],
                              k13 = th[4], k31 = th[5], ka = exp(th[6]),
                              V = 1),
                 times = tt, ncmt = 3, depot = TRUE, events = ev))
}, c(0, 0, 0, 0, 0, 0))$jacobian(c(0, 0, 0, 0, 0, 0))
cat("  gradient finite:", all(is.finite(as.numeric(g))), " ",
    paste(format(as.numeric(g), digits = 5), collapse = " "), "\n")

cat("\nordinary rates, the offset must change nothing:\n")
P1 <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.05, k31 = 0.01,
           ka = 1.1, V = 10)
print(frm_lincmt(parms = P1, times = tt, ncmt = 3, depot = TRUE,
                 events = ev), digits = 15)
