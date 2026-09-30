# One-off edit of R/compat.R: the two cov = FALSE rows that said newdata
# must carry the response. Kept as the record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/compat.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(n, " matches")
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("newdata must then carry the response, whose rows are read as their own groups in their own time order, as brms reads new data; without it the call is refused by name.",
     "On newdata the recursion reads newdata's own response, whose rows are read as their own groups in their own time order, as brms reads new data; a response that is NA, or absent, is filled with its expected value.")
rep1("which is brms's posterior_predict(); newdata must carry the response.",
     "which is brms's posterior_predict(). On newdata a response that is NA, or absent, is filled with a draw at its one-step mean before the later rows read its residual, as brms's .predictor_arma() fills it.")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
