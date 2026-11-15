# One-off edit of the two mi() compat notes (lane arspat2 fix round).
f1 <- "R/compat.R"
f2 <- "R/autocor-seam.R"
new1 <- "mi() on the response makes each missing value a coordinate of the joint draw, brms's Yl. A missing value no other density reads leaves frm()'s fit, which takes the exact margin of the observed rows (the fit over the observed rows with M cut to them), on gaussian() and student(), with random effects or REML beside them too; a value another density reads stays latent under the Laplace approximation. frm_sample() samples them, as brms does. log_lik() is refused, as it is for every mi() model."
new2 <- "mi() on the response makes each missing value a coordinate of the joint draw, brms's Yl. A missing value no other density reads leaves frm()'s fit, which takes the exact margin of the observed rows, on gaussian() and student(), with random effects or REML beside them too; a value another density reads stays latent, which whiten() receives on the tape, under the Laplace approximation. frm_sample() samples them, as brms does. log_lik() is refused, as it is for every mi() model."
fix <- function(f, start, new) {
  x <- readLines(f, warn = FALSE)
  i <- grep(start, x, fixed = TRUE)
  stopifnot(length(i) == 1L)
  lead <- sub("\".*$", "", x[i])
  tail <- sub("^.*\"([)]*)$", "\\1", x[i])
  x[i] <- paste0(lead, "\"", new, "\"", tail)
  writeLines(x, f, useBytes = TRUE)
}
fix(f1, "\"mi() on the response makes each missing value a coordinate of the joint draw, brms's Yl. frm() integrates", new1)
fix(f2, "\"mi() on the response makes each missing value a coordinate of the joint draw, brms's Yl, which whiten()", new2)
