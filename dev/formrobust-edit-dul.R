# One-off edit: frm(drop_unused_levels = ), brms's argument, through
# assemble_frame() and influence(). Kept as the record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/"
edit <- function(f, pairs) {
  p <- paste0(wt, f)
  x <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(x, gregexpr(pr[1], x, fixed = TRUE)))
    if (n != as.integer(pr[3])) stop(f, ": ", n, " matches for ", pr[1])
    x <- gsub(pr[1], pr[2], x, fixed = TRUE)
  }
  con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
  close(con)
}
edit("R/frame.R", list(
  c("                           check_response = TRUE, thres_pin = NULL) {",
    "                           check_response = TRUE, thres_pin = NULL,\n                           drop_unused_levels = TRUE) {", 1),
  c("                             drop.unused.levels = TRUE,",
    "                             drop.unused.levels = drop_unused_levels,", 2)
))
edit("R/fit.R", list(
  c("                data2 = list(), dry_run = NULL, verbose = FALSE) {\n  cl <- match.call()",
    "                data2 = list(), dry_run = NULL, verbose = FALSE,\n                drop_unused_levels = TRUE) {\n  cl <- match.call()", 1),
  c("  check_flag(quadrature, \"quadrature\")\n  check_count(importance, \"importance\", min = 0L)",
    "  check_flag(quadrature, \"quadrature\")\n  check_flag(drop_unused_levels, \"drop_unused_levels\")\n  check_count(importance, \"importance\", min = 0L)", 1),
  c("                          sparse_x = isTRUE(control$sparse_x),\n                          data2 = data2)\n  if (vb) vb_stage(\"frame\"",
    "                          sparse_x = isTRUE(control$sparse_x),\n                          data2 = data2,\n                          drop_unused_levels = drop_unused_levels)\n  if (vb) vb_stage(\"frame\"", 1)
))
edit("R/influence.R", list(
  c("                                data2 = data2, thres_pin = pin)",
    "                                data2 = data2, thres_pin = pin,\n                                drop_unused_levels =\n                                  model$frame[[\"drop_unused_levels\"]] %||%\n                                  TRUE)", 1)
))
