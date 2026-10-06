# Punch round 1, B2: can the probit's log-odds form be made to hold past
# |eta| = 38.2 without moving a plain fit? The current form is
# log(pnorm(eta)) - log(pnorm(-eta)) (R/links.R), which is -Inf - 0 and
# then NaN in the densities once pnorm underflows. The candidate reads
# both logs from pnorm(log.p = TRUE), which holds far past that. The
# fix is taken only if it is bitwise the current form wherever the
# current form is finite; this counts the grid points where it is not,
# through RTMB's tape (value and derivative), as a fit reads them.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
eta <- seq(-38, 38, length.out = 2001)
cur <- function(e) log(RTMB::pnorm(e)) - log(RTMB::pnorm(-e))
cand <- function(e) RTMB::pnorm(e, log.p = TRUE) - RTMB::pnorm(-e, log.p = TRUE)
tc <- RTMB::MakeTape(function(e) cur(e), eta)
tn <- RTMB::MakeTape(function(e) cand(e), eta)
vc <- tc(eta)
vn <- tn(eta)
gc <- RTMB::MakeTape(function(e) sum(cur(e)), eta)$jacobian(eta)
gn <- RTMB::MakeTape(function(e) sum(cand(e)), eta)$jacobian(eta)
fin <- is.finite(vc)
cat(sprintf("values: %d of %d finite points differ, max relative %.3g\n",
            sum(vc[fin] != vn[fin]), sum(fin),
            max(abs(vc[fin] - vn[fin]) / pmax(abs(vc[fin]), 1e-300))))
fg <- is.finite(gc)
cat(sprintf("derivatives: %d of %d finite points differ, max relative %.3g\n",
            sum(gc[fg] != gn[fg]), sum(fg),
            max(abs(gc[fg] - gn[fg]) / pmax(abs(gc[fg]), 1e-300))))
far <- c(-60, -40, 40, 60)
cat("past the underflow, eta =", far, "\n  current:  ",
    format(RTMB::MakeTape(cur, far)(far)), "\n  candidate:",
    format(RTMB::MakeTape(cand, far)(far)), "\n")
