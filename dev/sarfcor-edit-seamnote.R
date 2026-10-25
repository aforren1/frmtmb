p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/R/autocor-seam.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1('  entry <- list(args = args, parse = parse, build = build, whiten = whiten,',
'  if (!is.null(compat_note) && (!is.character(compat_note) ||
                                  length(compat_note) != 1L ||
                                  is.na(compat_note))) {
    frm_stop(arg("compat_note"), " must be one string", call. = FALSE)
  }
  entry <- list(args = args, parse = parse, build = build, whiten = whiten,')
rep1('  rules <- autocor_seam_compat_rules(feat)', '  rules <- autocor_seam_compat_rules(feat, compat_note)')
rep1('#\' correlation term: the refusals core applies to every such term, with
#\' the notes the `group:autocor` rows give core\'s own terms. What the
#\' term does with the families it takes is its package\'s to say, with
#\' frmtmb_register_compat(), as a registered covariance structure says
#\' what it does (covstruct_compat_rules()).
#\'
#\' @noRd
autocor_seam_compat_rules <- function(feat) {
  force(feat)
  function() {
    b <- compat_rule_builder()
    r <- b$r
    r(feat, "kind:family", "refused",
      "Refused: the term correlates the residuals y - mu, which only gaussian() and student() have.")',
'#\' correlation term: the refusals core applies to every such term, with
#\' the notes the `group:autocor` rows give core\'s own terms, and the
#\' families it takes with the note its package gives (`compat_note`).
#\' One row for the two families, so a package that wrote its own would
#\' tie with it and win no pair.
#\'
#\' @noRd
autocor_seam_compat_rules <- function(feat, note = NULL) {
  force(feat)
  note <- note %||% paste0("The residuals of the response are one joint ",
                           "draw; student() gets the multivariate-t ",
                           "analog, with one nu for the response.")
  function() {
    b <- compat_rule_builder()
    r <- b$r
    r(feat, "kind:family", "refused",
      "Refused: the term correlates the residuals y - mu, which only gaussian() and student() have.")
    r(feat, "group:autocor_families", "works", note)')
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
