# Lane wt-conditions: prepend the development-version NEWS entries.
# Run once from the worktree root; it refuses to run twice.
core <- c(
"# frmtmb (development version)",
"",
"* **BREAKING: every condition frmtmb raises is classed**, the contract",
"  brms keeps with `brms_error`. An error has the class",
"  `c(\"frmtmb_error\", \"error\", \"condition\")`, a warning",
"  `c(\"frmtmb_warning\", \"warning\", \"condition\")` and a message",
"  `c(\"frmtmb_message\", \"message\", \"condition\")`, so",
"  `tryCatch(frmtmb_error = )` catches any refusal. A condition that an",
"  extension raises also has a subclass named for the package, such as",
"  `frmtmb_eam_error`. See `?frmtmb-conditions`.",
"",
"  What changes for a caller: the class vector no longer contains",
"  `simpleError`, `simpleWarning` or `simpleMessage`, so a handler or",
"  test that catches those classes by name stops matching. A handler for",
"  `error`, `warning` or `message` is not affected. The message text and",
"  the recorded call of every condition are unchanged. A",
"  `frmtmb_fit_error`, raised when the optimizer fails, is now a",
"  `frmtmb_error` and is no longer a `simpleError`.",
"",
"  The class vector does not contain `rlang_error`, as brms's does:",
"  frmtmb builds the conditions with base R and does not import rlang.",
"  Errors from other packages, such as RTMB, TMB and Matrix, keep their",
"  own class, and a failed `stopifnot()` assertion is still a",
"  `simpleError`.",
"",
"* `frm_stop()`, `frm_warning()` and `frm_message()` are exported for",
"  extension authors. They take the arguments of `stop()`, `warning()`",
"  and `message()` and add the classes. Every extension in this",
"  repository now requires them, so the `frmtmb (>= 0.59.0)` floor of",
"  frmtmb.coupling, frmtmb.eam, frmtmb.latent, frmtmb.learn, frmtmb.ode,",
"  frmtmb.sample and frmtmb.spline must move to the release that",
"  includes this change.",
"")

ext <- function(pkg, extra = character(0)) c(
  paste0("# ", pkg, " (development version)"),
  "",
  paste0("* Requires the frmtmb release that exports `frm_stop()`; the"),
  "  `frmtmb (>= 0.59.0)` floor must move to it.",
  paste0("* **BREAKING:** every error, warning and message that ", pkg),
  "  raises is classed. An error has the class",
  paste0("  `c(\"", gsub(".", "_", pkg, fixed = TRUE),
         "_error\", \"frmtmb_error\", \"error\", \"condition\")`,"),
  "  and warnings and messages follow the same pattern, so",
  "  `tryCatch(frmtmb_error = )` catches any refusal. The class vector no",
  "  longer contains `simpleError`, `simpleWarning` or `simpleMessage`.",
  "  The message text is unchanged. See `?frmtmb::frmtmb-conditions`.",
  extra,
  "")

spline_extra <- c(
  "* The `frmtmb_ps_span_warning` that `frm_curve()` and",
  "  `frm_curve_deriv()` raise now also has the classes",
  "  `frmtmb_spline_warning` and `frmtmb_warning`.")
eam_extra <- c(
  "* The `frmtmb_eam_units_warning` now also has the classes",
  "  `frmtmb_eam_warning` and `frmtmb_warning`.")

prepend <- function(path, lines) {
  old <- readLines(path, warn = FALSE, encoding = "UTF-8")
  stopifnot(!grepl("development version", old[1L], fixed = TRUE))
  con <- file(path, open = "wb")
  writeLines(enc2utf8(c(lines, old)), con, sep = "\n", useBytes = TRUE)
  close(con)
}
prepend("NEWS.md", core)
for (p in c("coupling", "eam", "latent", "learn", "ode", "sample",
            "spline")) {
  pkg <- paste0("frmtmb.", p)
  extra <- switch(p, spline = spline_extra, eam = eam_extra, character(0))
  prepend(file.path("extensions", pkg, "NEWS.md"), ext(pkg, extra))
}
cat("done\n")
