# Lane ceplot punch 2: the mm() divergence under "old_levels" in
# ?fitted.frmtmb_fit and ?predict.frmtmb_fit.
rep1 <- function(f, old, new) {
  s <- paste(readLines(f), collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop(basename(f), ": expected one match, got ", n)
  writeLines(strsplit(sub(old, new, s, fixed = TRUE), "\n", fixed = TRUE)[[1]], f)
}
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/"
rep1(paste0(wt, "predict.R"),
"#'   conditional variance rather than the marginal one. It acts only with
#'   `allow_new_levels = TRUE`, as in brms.",
"#'   conditional variance rather than the marginal one. One choice
#'   serves every response of a multivariate call. Two different unseen
#'   values of one `mm()` term read two seen levels; brms 2.23.0 gives
#'   them one, because it numbers unseen values per member rather than
#'   by value. It acts only with `allow_new_levels = TRUE`, as in brms.")
rep1(paste0(wt, "predict-brms.R"),
"#'   same seen level (within the row's by-level for `gr(by = )`),
#'   from each replicate's draw of the group effects.",
"#'   same seen level (within the row's by-level for `gr(by = )`),
#'   from each replicate's draw of the group effects. One choice serves
#'   every response of a multivariate call. Two different unseen values
#'   of one `mm()` term read two seen levels; brms 2.23.0 gives them
#'   one, because it numbers unseen values per member rather than by
#'   value.")
cat("edited\n")
