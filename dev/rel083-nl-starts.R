# Merge step 5 (user decision 2026-11-13, backlog item 5718):
# frmtmb_control(nl_starts = TRUE). FALSE switches off the two automatic
# restarts of a nonlinear fit, the sign probe (nl_sign_probe()) and the
# restart beyond a pole of the body (nl_pole_escape()). Run once; each
# edit stops if its anchor is not found exactly once.
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel083"
rl <- function(f) readLines(file.path(REL, f), warn = FALSE)
wr <- function(x, f) {
  con <- file(file.path(REL, f), "wb"); on.exit(close(con))
  writeLines(x, con, sep = "\r\n")
}
one <- function(x, pat, fixed = TRUE) {
  i <- grep(pat, x, fixed = fixed)
  if (length(i) != 1L) stop("anchor '", pat, "' found ", length(i), " times")
  i
}

## R/fit.R: frmtmb_control()
f <- "R/fit.R"
x <- rl(f)
stopifnot(!any(grepl("nl_starts =", x, fixed = TRUE)))
i <- one(x, "#' @param quad_max_gb Limit, in GB, on the memory that `frm(quadrature =")
x <- append(x, c(
  "#' @param nl_starts Whether a fit with a nonlinear body runs its two",
  "#'   automatic restarts: from the default start, the fits from each sign",
  "#'   pattern of the body's coefficients (the sign probe, message class",
  "#'   `frmtmb_nl_multistart`), and the fits from beyond a pole of the body",
  "#'   that holds the fit (message class `frmtmb_nl_pole_start`; from a",
  "#'   `start` of your own, a warning of class `frmtmb_nl_pole` that the",
  "#'   point is not the maximum). `TRUE` (default) runs them. `FALSE` runs",
  "#'   neither, so the fit is the one the optimizer reaches from the start,",
  "#'   with no message about a restart. The restart from a saddle point",
  "#'   (`frmtmb_nl_restart`) and the warning that a fit is on a ridge to a",
  "#'   limit of its body (`frmtmb_nl_limit`) are not affected."),
  after = i - 1L)
i <- one(x, "                           quad_max_gb = 4, starts = NULL) {")
x[i] <- "                           quad_max_gb = 4, starts = NULL,"
x <- append(x, "                           nl_starts = TRUE) {", after = i)
i <- one(x, "  check_flag(mo_search, \"mo_search\")")
x <- append(x, "  check_flag(nl_starts, \"nl_starts\")", after = i)
i <- one(x, "       quad_max_gb = quad_max_gb, starts = starts)")
x[i] <- "       quad_max_gb = quad_max_gb, starts = starts,"
x <- append(x, "       nl_starts = isTRUE(nl_starts))", after = i)
# the `start` documentation names the switch
i <- one(x, "#'   (see \"The default starting values\" in [par_template()]). A `start`, or")
stopifnot(identical(x[i + 1L], "#'   `options(frmtmb.nl_sign_probe = FALSE)`, skips that; a located prior"))
x[i] <- "#'   (see \"The default starting values\" in [par_template()]). A `start`,"
x[i + 1L] <- "#'   `options(frmtmb.nl_sign_probe = FALSE)` skips that; a located prior"
x <- append(x, "#'   `frmtmb_control(nl_starts = FALSE)` or", after = i)
wr(x, f)

## R/par-template.R: the sign probe
f <- "R/par-template.R"
x <- rl(f)
i <- one(x, "  if (!isTRUE(getOption(\"frmtmb.nl_sign_probe\", TRUE)) ||")
x <- append(x, "        isFALSE(ctl$nl_starts) ||", after = i)
i <- one(x, "#' `options(frmtmb.nl_sign_probe = FALSE)`.")
x[i] <- "#' `options(frmtmb.nl_sign_probe = FALSE)`."
x <- append(x, "#' `frmtmb_control(nl_starts = FALSE)` or", after = i - 1L)
i <- one(x, "#' `start`, or `options(frmtmb.nl_sign_probe = FALSE)`, skips it; a located")
x[i] <- "#' `options(frmtmb.nl_sign_probe = FALSE)` skips it; a located"
x <- append(x, "#' `start`, `frmtmb_control(nl_starts = FALSE)` or", after = i - 1L)
wr(x, f)

## R/fit-restart.R: the pole restart
f <- "R/fit-restart.R"
x <- rl(f)
i <- one(x, "nl_pole_escape <- function(run, frame, fit) {")
stopifnot(identical(x[i + 3L], "  if (is.null(obj) || is.null(fit$opt) || !is.finite(fit$opt$objective) ||"))
x[i + 3L] <- "  if (isFALSE(ctl$nl_starts) ||"
x <- append(x, "        is.null(obj) || is.null(fit$opt) || !is.finite(fit$opt$objective) ||",
            after = i + 3L)
j <- one(x, "#' Returns `NULL`, or a list as nl_sign_probe() returns it.")
x <- append(x, c(
  "#' Not run under `frmtmb_control(nl_starts = FALSE)`, which asks for the",
  "#' fit from the start alone.",
  "#'"), after = j - 1L)
wr(x, f)
cat("NL_STARTS DONE\n")
