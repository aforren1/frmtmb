p <- "R/autocor.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  fam <- resp$family[[\"family\"]]
  if (!fam %in% c(\"gaussian\", \"student\")) {
    frm_stop(ac[[\"fn\"]], \"(): a residual correlation needs a family with real \",",
"  fam <- resp$family[[\"family\"]]
  if (!fam %in% c(\"gaussian\", \"student\") && autocor_is_seam(ac)) {
    frm_stop(ac[[\"fn\"]], \"(): a residual correlation needs a family with \",
             \"real residuals, so gaussian() or student(); '\", fam, \"' has \",
             \"none. brms refuses the same call (\\\"\", toupper(ac[[\"fn\"]]),
             \" terms are not implemented for this family\\\")\", call. = FALSE)
  }
  if (!fam %in% c(\"gaussian\", \"student\")) {
    frm_stop(ac[[\"fn\"]], \"(): a residual correlation needs a family with real \",")
rep1("  if (fam == \"student\") {
    # brms's student_t_time_*_lpdf takes `real nu`: the multivariate t
    # has one shape for the whole group, so a row-varying nu has no
    # multivariate counterpart
    nud <- resp$dpars[[\"nu\"]]
    predicted <- !is.null(nud) && is.null(nud$constant) &&
      (length(nud$re %||% list()) || length(nud$smooth %||% list()) ||
         !identical(deparse1(reformulas::RHSForm(nud$fixed)), \"1\"))
    if (predicted) {",
"  if (autocor_is_seam(ac) && dpar_spec_predicted(resp, \"sigma\") &&
        !isTRUE(ac_block_entry(ac)[[\"sigma_varies\"]])) {
    frm_stop(ac[[\"fn\"]], \"(): this term needs a constant sigma; a \",
             \"predicted 'sigma ~ ...' is refused, as brms refuses it (\\\"\",
             toupper(ac[[\"fn\"]]), \" models are not implemented when \",
             \"predicting 'sigma' or 'nu'\\\")\", call. = FALSE)
  }
  if (fam == \"student\") {
    # brms's student_t_time_*_lpdf takes `real nu`: the multivariate t
    # has one shape for the whole group, so a row-varying nu has no
    # multivariate counterpart
    if (dpar_spec_predicted(resp, \"nu\")) {")
rep1("#' The grouping value of every row, as one label per row.",
"#' Whether a response's spec models `dp` with a formula (a fixed term
#' other than the intercept, a random effect or a smooth) rather than
#' leaving it constant.
#'
#' @noRd
dpar_spec_predicted <- function(resp, dp) {
  d <- resp$dpars[[dp]]
  !is.null(d) && is.null(d$constant) &&
    (length(d$re %||% list()) || length(d$smooth %||% list()) ||
       !identical(deparse1(reformulas::RHSForm(d$fixed)), \"1\"))
}

#' The grouping value of every row, as one label per row.")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
