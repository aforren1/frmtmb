# A defect found and not fixed: frm_ode() cannot solve a dynamics whose
# derivative for some state does not depend on any state. RTMBode's
# tape reports "Non-consecutive outputs" and the group is recorded as a
# failed solve, so its rows carry the penalty.
#
# Run with NSS_ARM=ref to confirm it is pre-existing.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-16-nonconsec.R
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
cat("arm:", Sys.getenv("NSS_ARM", "lane"), " frmtmb.ode from",
    dirname(find.package("frmtmb.ode")), "\n")

show <- function(tag, dyn, ini) {
  r <- tryCatch(frm_ode(dyn, init = ini, times = c(0, 3),
                        parms = list(0.2), on_error = "error"),
                error = function(e) conditionMessage(e))
  cat(sprintf("  %-34s %s\n", tag,
              if (is.character(r)) paste("ERROR:", r) else
                paste(format(as.numeric(r)), collapse = " ")))
}
show("derivative depends on a state", function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], y[1]))
}, list(100, 0))
show("derivative is 0 * y[2]", function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], 0 * y[2]))
}, list(100, 0))
show("derivative is a constant", function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1]))
}, list(100, 0))
