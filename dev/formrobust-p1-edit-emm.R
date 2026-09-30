# Punch round 1, blocker: emmeans() keeps a linear predictor's offset,
# as brms does (emmeans's own .offset. column). Record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/interop.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(n, " matches for ", substr(old, 1, 70))
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("    tt <- stats::delete.response(tg$targets[[1L]]$lp[[\"terms\"]])
    # emmeans adds an offset term of the terms to every grid prediction;
    # brms's emm_basis() predicts with offset = FALSE, and the design
    # route is the linear predictor only
    attr(tt, \"offset\") <- NULL
    return(tt)",
"    # the offset stays in the terms: emmeans then puts it in the grid's
    # .offset. column and adds it at the grid's value, which is how
    # brms's emmeans() includes it (its basis is offset = FALSE)
    return(stats::delete.response(tg$targets[[1L]]$lp[[\"terms\"]]))")
rep1("  # brms's emm_basis() takes the linear predictor with offset = FALSE
  # and the expected response with the offset at the grid's values
  if (!tg$epred) object <- emm_drop_offsets(object)",
"  # a nonlinear body's parameters keep their offsets out, as brms's
  # grid cannot see them; the selected predictor keeps its own
  if (!tg$epred) object <- emm_drop_offsets(object, tg)")
rep1("#' The fit with every linear predictor's `offset()` terms taken out of
#' the prediction. The offset's variables stay in the terms, so a grid
#' still has to carry them, as brms's grid does.
#'
#' @noRd
emm_drop_offsets <- function(object) {
  lps <- object$frame[[\"linpreds\"]]
  for (k in names(lps)) {",
"#' The fit with the `offset()` terms of every linear predictor except
#' the selected ones taken out of the prediction. brms adds a selected
#' predictor's offset back through emmeans's `.offset.` column, which on
#' the grid route is the offset `frm_lp_basis()` already includes; an
#' offset inside a nonlinear parameter's formula never reaches brms's
#' reference grid, so a nonlinear body's emmean leaves it out.
#'
#' @noRd
emm_drop_offsets <- function(object, tg) {
  lps <- object$frame[[\"linpreds\"]]
  keep <- vapply(tg$targets, function(t) linpred_key(t$resp, t$name), \"\")
  for (k in setdiff(names(lps), keep)) {")
rep1("  design <- !epred && !use_re &&
    all(vapply(targets, function(t) emm_is_design(t$lp), NA))",
"  design <- !epred && !use_re &&
    all(vapply(targets, function(t) emm_is_design(t$lp), NA))
  # emmeans has one .offset. column for the whole grid, which would add
  # one response's offset to every response; the grid route adds each
  # response its own
  if (design && length(targets) > 1L &&
        any(vapply(targets, function(t) {
          length(attr(t$lp[[\"terms\"]], \"offset\")) > 0L
        }, NA))) {
    design <- FALSE
  }")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
