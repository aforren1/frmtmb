## Punch round 2, MINOR 2: every call in the test trees and in
## dev/brms-vignettes/ that names a family needing trials() and whose own
## text has no trials() or cbind(). It is a PARSE, not a grep, so a call
## spanning lines is one hit. It over-reports on purpose (a formula held
## in a variable has no trials() in the call text); each hit is read by
## hand, and the decisions are in dev/famlink-findings.src.md.
## Run from the worktree root: Rscript dev/famlink-p2-trials-scan.R
files <- c(list.files("tests/testthat", "[.][Rr]$", full.names = TRUE),
           Sys.glob("extensions/*/tests/testthat/*.R"),
           list.files("dev/brms-vignettes", "[.][Rr]$", full.names = TRUE))
fam_pat <- paste0("(^|[^_a-z])(binomial|multinomial|beta_binomial|",
                  "zero_inflated_binomial)[(]|[\"'](binomial|multinomial|",
                  "beta_binomial|zero_inflated_binomial)[\"']")
top_calls <- function(e, acc) {
  if (!is.call(e)) return(acc)
  fn <- e[[1]]
  nm <- if (is.name(fn)) as.character(fn) else ""
  txt <- paste(deparse(e, width.cutoff = 500L), collapse = " ")
  # the smallest call that fits a model: frm() and its relatives
  if (nm %in% c("frm", "frm_sample", "update", "refit", "get_prior",
                "default_prior", "frm_simulate", "brm", "make_standata",
                "standata", "par_template") &&
      grepl(fam_pat, txt, perl = TRUE)) {
    if (!grepl("trials[(]|cbind[(]", txt)) acc[[length(acc) + 1L]] <- txt
    return(acc)
  }
  args <- as.list(e)[-1]
  for (i in seq_along(args)) {
    # an empty argument, as in x[, 1], is the missing symbol and cannot
    # be bound to a name, so test it before touching it
    if (identical(args[[i]], quote(expr = ))) next
    if (is.call(args[[i]])) acc <- top_calls(args[[i]], acc)
  }
  acc
}
n_hit <- 0L
for (f in files) {
  ex <- tryCatch(parse(f, keep.source = TRUE), error = function(e) NULL)
  if (is.null(ex)) { cat("PARSE FAIL", f, "\n"); next }
  sr <- attr(ex, "srcref")
  for (k in seq_along(ex)) {
    hits <- top_calls(ex[[k]], list())
    for (h in hits) {
      n_hit <- n_hit + 1L
      cat(sprintf("%s (expr at line %d): %s\n", f, sr[[k]][1],
                  substr(h, 1, 240)))
    }
  }
}
cat("files scanned", length(files), " hits", n_hit, "\n")
