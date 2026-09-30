# One-off edit: autocor() joins the brms-owned generics. Kept as the
# record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/"
edit <- function(f, old, new) {
  p <- paste0(wt, f)
  x <- paste(readLines(p), collapse = "\n")
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(f, ": ", n, " matches")
  x <- sub(old, new, x, fixed = TRUE)
  con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
  close(con)
}
edit("R/generic-owners.R", "  expose_functions = \"brms\",\n",
     "  expose_functions = \"brms\",\n  autocor = \"brms\",\n")
edit("tests/testthat/test-generic-collision.R",
     "  \"conditional_effects\", \"conditional_smooths\", \"expose_functions\",\n",
     "  \"conditional_effects\", \"conditional_smooths\", \"expose_functions\",\n  \"autocor\",\n")
