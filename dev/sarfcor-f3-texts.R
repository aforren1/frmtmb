# Fix round 1, F3: the compat notes and ?frmtmb-autocor say what
# newdata on the fitted rows does. Rscript <this> <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
same_rows <- paste0(
  "newdata with as many rows as the fit is read as the fitted rows in ",
  "their order, as brms reads it, so a counterfactual (x + 1 everywhere) ",
  "works; other newdata is refused, as brms refuses it without a matrix ",
  "over the new rows (its newdata2).")
edit("R/compat.R", list(
  c("\"The mean is mu, unchanged. fitted() and predict() refuse newdata, as brms does without a matrix over the new rows (its newdata2); frm_linpred() gives mu on newdata.\"",
    paste0("\"The mean is mu, unchanged. ", same_rows,
           " frm_linpred() gives mu on any newdata.\"")),
  c("\"Draws in sample are one draw of the whole response with covariance sigma^2 M. newdata is refused, as in brms without newdata2.\"",
    paste0("\"Draws are one draw of the whole response with covariance ",
           "sigma^2 M. ", same_rows, "\"")),
  c("and the student() draw scales it by one chi-square. newdata is refused: M covers the fitted rows only.\"",
    paste0("and the student() draw scales it by one chi-square. ",
           same_rows, "\""))))
edit("extensions/frmtmb.spatial/R/zzz.R", list(
  c("under \\\"error\\\" it is mu. newdata is refused, as brms refuses it without a matrix over the new rows (its newdata2); frm_linpred() gives mu on newdata.\"",
    paste0("under \\\"error\\\" it is mu. ", same_rows,
           " On such newdata the lag mean is A^-1 mu_new. frm_linpred() ",
           "gives mu on any newdata.\"")),
  c("\"Draws in sample are one draw of the whole response through the model. newdata is refused, as in brms without newdata2.\"",
    paste0("\"Draws are one draw of the whole response through the model. ",
           same_rows, "\"")),
  c("A^-1 (mu + e) or mu + A^-1 e. newdata is refused: M covers the fitted rows only.\"",
    paste0("A^-1 (mu + e) or mu + A^-1 e. ", same_rows, "\""))))
edit("R/autocor.R", list(
  c("#' The term has no parameter, and the mean of the response is `mu`.
#' [fitted()], [predict()] and [simulate()] refuse `newdata`, because
#' `M` covers the fitted rows only; brms refuses the same calls unless a
#' matrix over the new rows comes in its `newdata2`. [frm_linpred()]
#' gives `mu` on new data.",
    "#' The term has no parameter, and the mean of the response is `mu`.
#' [fitted()], [predict()] and [simulate()] take `newdata` with as many
#' rows as the fit, read as the fitted rows in their order, as brms
#' takes it, so a counterfactual such as `x + 1` everywhere works. Other
#' `newdata` is refused, because `M` covers the fitted rows only; brms
#' refuses it too unless a matrix over the new rows comes in its
#' `newdata2`. [frm_linpred()] gives `mu` on any new data.")))
edit("extensions/frmtmb.spatial/R/sar.R", list(
  c("#' frmtmb.sample's default is brms's, flat on rho over its interval. The pointwise
#' `log_lik()` is refused, as for frmtmb's other residual covariance
#' terms;",
    "#' frmtmb.sample's default is brms's, flat on rho over its interval.
#' The pointwise `log_lik()` is refused, as for frmtmb's other residual
#' covariance terms;")))
cat("edited", root, "\n")
