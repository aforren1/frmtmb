# One-off edit of frmtmb.sample's posterior_predict(): a cov = FALSE ARMA
# term is drawn row by row, so it is not a structured draw (core's
# predict_simulate() says the same). Kept as the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/extensions/frmtmb.sample/R/methods-draws.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("  if (!is.null(newdata) &&
      sim_is_structured(sim_context(fit, rspec, list(), aterms = av))) {",
"  # brms's cov = FALSE ARMA: brms's posterior_predict() draws each row
  # around its one-step mean, which reads the OBSERVED earlier residuals
  # and which frm_linpred() gives; the rows are then drawn one by one,
  # not by the recursion simulate() runs over its own draws.
  # Core's own predicate, not a read of frame$autocor$cov: two files
  # asking \"is this brms's cov = FALSE form\" must not answer it twice
  arma_cond <- resp %in% arma_cond_resp(fit)
  # so a cov = FALSE term is not a structured draw, on newdata or under
  # re_formula, as core's predict() treats it
  structured <- local({
    ctx0 <- sim_context(fit, rspec, list(), aterms = av)
    if (arma_cond) ctx0[[\"autocor\"]] <- NULL
    sim_is_structured(ctx0)
  })
  if (!is.null(newdata) && structured) {")
rep1("  if (!is.null(re_form) &&
      sim_is_structured(sim_context(fit, rspec, list(), aterms = av))) {",
"  if (!is.null(re_form) && structured) {")
rep1("  # brms's cov = FALSE ARMA: brms's posterior_predict() draws each row
  # around its one-step mean, which reads the OBSERVED earlier residuals
  # and which frm_linpred() gives; the rows are then drawn one by one,
  # not by the recursion simulate() runs over its own draws.
  # Core's own predicate, not a read of frame$autocor$cov: two files
  # asking \"is this brms's cov = FALSE form\" must not answer it twice
  arma_cond <- resp %in% arma_cond_resp(fit)
  # one draw's",
"  # one draw's")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
