# What the warning SAYS, base against lane, on the rows the review
# measured the old warning's understatement on. Both arms fire on the
# same rows (nss-12-warn.R); the number they print is what changed.
#
# Run with NSS_ARM=ref for the base commit.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-13-warntext.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
ARM <- Sys.getenv("NSS_ARM", "lane")
cat("arm:", ARM, " frmtmb.ode from", dirname(find.package("frmtmb.ode")),
    "\n\n")

one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
say <- function(ke, ii) {
  tt <- seq(0, ii, length.out = 25)
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = ii,
                   ss = TRUE)
  args <- list(one_oral, init = list(0, 0), times = tt,
               parms = list(ke, 1.0), events = ev, output = 2L,
               n_ss = 20L)
  if (ARM != "ref") args$ss_extrapolate <- FALSE
  msg <- NULL
  v <- withCallingHandlers(do.call(frm_ode, args),
                           warning = function(w) {
                             msg <<- conditionMessage(w)
                             invokeRestart("muffleWarning")
                           })
  ref <- as.numeric(frm_lincmt(parms = list(ke = ke, ka = 1.0, V = 10),
                               times = tt, ncmt = 1, depot = TRUE,
                               events = data.frame(
                                 time = 0, state = "depot",
                                 value = 100, ii = ii, addl = 0L,
                                 ss = TRUE)))
  true <- max(abs(as.numeric(v) / 10 - ref)) / max(abs(ref))
  num <- if (is.null(msg)) NA_real_ else
    as.numeric(sub(".*is about ([0-9.e+-]+) .*|.*moving by ([0-9.e+-]+) .*",
                   "\\1\\2", msg))
  cat(sprintf("t_half %6.1f  ke*ii %5.2f  printed %10.3e  true %10.3e",
              log(2) / ke, ke * ii, num, true))
  cat(sprintf("  understatement %6.2fx\n", true / num))
}
cat(sprintf("%s\n", "one compartment, ka = 1, n_ss = 20, truncated"))
for (ke in c(0.05, 0.02, 0.01, 0.005)) say(ke, 12)
