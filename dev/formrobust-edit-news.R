# Put the lane's NEWS under a development-version heading at the top of
# NEWS.md and of frmtmb.sample's NEWS.md. Kept as the record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/"
put <- function(f, lines) {
  p <- paste0(wt, f)
  x <- readLines(p)
  con <- file(p, "wb"); writeLines(c(lines, "", x), con, sep = "\r\n")
  close(con)
}
core <- readLines(paste0(wt, "dev/formrobust-news.md"))
put("NEWS.md", core)
put("extensions/frmtmb.sample/NEWS.md", c(
  "# frmtmb.sample (development version)",
  "",
  "Needs the frmtmb that exports `arma_cond_fill_dpars()` and",
  "`response_codes_newdata()` (the next frmtmb release after 0.66.0).",
  "",
  "## Bug fixes",
  "",
  "* **`posterior_predict(newdata = )` of a model with `ar()`, `ma()` or",
  "  `arma()` and `cov = FALSE`** was refused as a structured draw, the",
  "  same refusal as a `cov = TRUE` block, though core's `predict()`",
  "  answers it. It answers now, and with brms's treatment of a missing",
  "  response: a row of newdata whose response is `NA`, or every row when",
  "  the column is absent, is filled with a draw at that draw's",
  "  parameters before the rows after it read its residual.",
  "  `re_formula = ` is no longer refused on such a model either.",
  "* **`predictive_error(newdata = )` of a `bernoulli()` model** codes",
  "  newdata's response as the fit coded its own, so a response held as",
  "  two values other than 0 and 1 (brms's level-order coding, new in",
  "  frmtmb) is compared with the draws on the same scale."))
