# The coordinator's wording for merge step 5 (2026-11-14): ?frmtmb_control
# names the two restarts "the sign probe and the pole restart", and says
# that the saddle restart is outside the switch.
f <- "C:/Users/adf44/source/r/frmtmb-wt-rel083/R/fit.R"
x <- readLines(f)
i <- which(x == "#' @param nl_starts Whether a fit with a nonlinear body runs its two")
stopifnot(length(i) == 1, x[i + 10L] == "#'   limit of its body (`frmtmb_nl_limit`) are not affected.")
new <- c(
"#' @param nl_starts Whether a fit with a nonlinear body runs the sign",
"#'   probe and the pole restart. The sign probe, from the default start",
"#'   only, fits from each sign pattern of the body's coefficients (message",
"#'   class `frmtmb_nl_multistart`). The pole restart fits from beyond a",
"#'   pole of the body that holds the fit (message class",
"#'   `frmtmb_nl_pole_start`; from a `start` of your own, a warning of class",
"#'   `frmtmb_nl_pole` that the point is not the maximum). `TRUE` (default)",
"#'   runs both. `FALSE` runs neither, so the fit is the one the optimizer",
"#'   reaches from the start, with no message about a restart. The restart",
"#'   from a saddle point (`frmtmb_nl_restart`) is outside this switch, and",
"#'   the warning that a fit is on a ridge to a limit of its body",
"#'   (`frmtmb_nl_limit`) is not affected.")
x <- c(x[seq_len(i - 1L)], new, x[-seq_len(i + 10L)])
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("done\n")
