# Fix round 1: the block checks (A1), the owner on each block and the
# refusal of a fit whose structure's package is absent (A3), the OLRE
# property (A2), on the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
rep_file <- function(p, old, new) {
  p <- file.path(root, p)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, s, fixed = TRUE)[[1]])
  if (regexpr(old, s, fixed = TRUE) < 0 || n != 1) {
    stop(p, ": count ", n, ": ", substr(old, 1, 70))
  }
  writeLines(sub(old, new, s, fixed = TRUE), p)
}
# A1: npar checked where it is read
rep_file("R/frame.R", '    npar_k <- as.integer(cs_entry(cs_name)$npar(proto))', '    npar_k <- cs_block_npar(cs_name, proto)')
# A3: the block records the package that owns its structure
rep_file("R/frame.R", '      expands = cps[[1]]$expands,
', '      expands = cps[[1]]$expands,
      # the package that owns the structure, so that a fit read where it
      # is absent is refused by its name (require_structure_packages())
      covstruct_package = cs_entry(cs_name)[["package"]],
')
# A1: every block against its entry, once, after theta is laid out
rep_file("R/frame.R", '  par_template <- list(beta = stats::setNames(numeric(length(beta_names)),
                                              beta_names))', '  for (bk in re_blocks) cs_check_block(bk)
  par_template <- list(beta = stats::setNames(numeric(length(beta_names)),
                                              beta_names))')
# A2: the OLRE check reads the property
rep_file("R/fit.R", '    structured <- bk[["covstruct"]] %in% c("gr_cov", "gr_prec", "equalto")', '    structured <- cs_levels_coupled(bk) ||
      identical(bk[["covstruct"]], "equalto")')
rep_file("R/fit.R", '    # A known-structure block is not an OLRE even with one row per
    # level. Its levels are correlated through the fixed relationship
    # matrix (gr(cov = A), gr(prec = Q)) or its covariance is fixed
    # outright (equalto),', '    # A known-structure block is not an OLRE even with one row per
    # level. Its levels are correlated through a fixed matrix
    # (cs_levels_coupled(): gr(cov = A), gr(prec = Q), car(), spde(), a
    # registered field) or its covariance is fixed outright (equalto),')
# A3: the objective and the standard-error machinery check the owners
rep_file("R/objective.R", 'build_objective <- function(frame) {
', 'build_objective <- function(frame) {
  require_structure_packages(frame)
')
rep_file("R/fit.R", '                            "vcov(), confint(), frm_linpred(se.fit =))"))
  cache <- fit$cache
', '                            "vcov(), confint(), frm_linpred(se.fit =))"))
  require_structure_packages(fit)
  cache <- fit$cache
')
cat("fix-a2 applied to", root, "\n")
