# One-off source edit for lane sarfcor: route the residual-block
# helpers of R/autocor.R to the registry for seam terms.
p <- "R/autocor.R"
s <- readLines(p)
txt <- paste(s, collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("autocor_npar <- function(ac) {\n  switch(",
     "autocor_npar <- function(ac) {\n  if (autocor_is_seam(ac)) return(ac[[\"npar\"]])\n  switch(")
rep1("autocor_start <- function(ac) numeric(autocor_npar(ac))",
     "autocor_start <- function(ac) {\n  if (autocor_is_seam(ac)) return(ac_block_entry(ac)[[\"start\"]](ac[[\"aux\"]]))\n  numeric(autocor_npar(ac))\n}")
rep1("autocor_cor <- function(theta, ac) {\n  d <- ac[[\"d\"]]\n",
     "autocor_cor <- function(theta, ac) {\n  d <- ac[[\"d\"]]\n  if (autocor_is_seam(ac)) {\n    # the residual covariance over sigma^2, from the term's map at a\n    # unit sigma; the mean does not enter it\n    B <- autocor_seam_B(ac, theta, numeric(d), 1)\n    return(tcrossprod(B))\n  }\n")
rep1("autocor_natural <- function(theta, ac) {\n",
     "autocor_natural <- function(theta, ac) {\n  if (autocor_is_seam(ac)) {\n    f <- ac_block_entry(ac)[[\"natural\"]]\n    return(if (is.null(f)) numeric(0) else f(theta, ac[[\"aux\"]]))\n  }\n")
rep1("autocor_types <- function(ac) {\n",
     "autocor_types <- function(ac) {\n  # a registered term reports its parameters as they are, on their own\n  # bounded scale, so the interval is the delta method's on that scale\n  if (autocor_is_seam(ac)) return(rep(\"raw\", length(autocor_natural(\n    numeric(ac[[\"npar\"]]), ac))))\n")
rep1("  for (ac in acs) {\n    t0 <- th[ac[[\"theta_idx\"]]]",
     "  for (ac in acs) {\n    # a term with no parameter, fcor(), has no row to report\n    if (!length(ac[[\"theta_idx\"]])) next\n    t0 <- th[ac[[\"theta_idx\"]]]")
rep1("  out <- do.call(rbind, rows)\n  rownames(out) <- NULL\n  out",
     "  if (!length(rows)) return(NULL)\n  out <- do.call(rbind, rows)\n  rownames(out) <- NULL\n  out")
rep1("    is.call(t) && as.character(t[[1L]])[1L] %in% autocor_structs\n",
     "    is.call(t) && as.character(t[[1L]])[1L] %in% autocor_term_names()\n")
rep1("             paste0(autocor_structs, \"()\", collapse = \", \"), \"); '\",",
     "             paste0(autocor_term_names(), \"()\", collapse = \", \"), \"); '\",")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
