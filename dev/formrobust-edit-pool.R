# One-off edit of R/confint.R: the pooling of parameter formulas shared
# by the two update paths, so its message has one template. Kept as the
# record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/confint.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(n, " matches for ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("  pars <- c(old$pforms, old$pfix, new$pforms, new$pfix)
  dup <- duplicated(names(pars), fromLast = TRUE)
  if (any(dup)) {
    frm_message(\"Replacing initial definitions of parameters \",
                paste(unique(names(pars)[dup]), collapse = \", \"))
    pars <- pars[!dup]
  }
  is_form <- vapply(pars, inherits, NA, \"formula\")
  out$pforms <- pars[is_form]
  out$pfix <- pars[!is_form]
  check_dpar_equations(out$pfix,
                       c(names(out$pforms), names(out$nlforms)))
  if (!is.null(new$family)) out$family <- new$family",
"  out <- update_pool_pars(out, old, new)
  if (!is.null(new$family)) out$family <- new$family")
rep1("  pars <- c(old$pforms, old$pfix, new$pforms, new$pfix)
  dup <- duplicated(names(pars), fromLast = TRUE)
  if (any(dup)) {
    frm_message(\"Replacing initial definitions of parameters \",
                paste(unique(names(pars)[dup]), collapse = \", \"))
    pars <- pars[!dup]
  }
  is_form <- vapply(pars, inherits, NA, \"formula\")
  new$pforms <- pars[is_form]
  new$pfix <- pars[!is_form]
  nlf <- c(old$nlforms, new$nlforms)
  new$nlforms <- nlf[!duplicated(names(nlf), fromLast = TRUE)]
  check_dpar_equations(new$pfix,
                       c(names(new$pforms), names(new$nlforms)))
  new
}",
"  nlf <- c(old$nlforms, new$nlforms)
  new$nlforms <- nlf[!duplicated(names(nlf), fromLast = TRUE)]
  update_pool_pars(new, old, new)
}

#' The parameter formulas, constants and equations of a stored model and
#' of an update pooled into `out`, a later one replacing an earlier one
#' of the same name with brms's message, as `update.brmsformula()` pools
#' them.
#'
#' @noRd
update_pool_pars <- function(out, old, new) {
  pars <- c(old$pforms, old$pfix, new$pforms, new$pfix)
  dup <- duplicated(names(pars), fromLast = TRUE)
  if (any(dup)) {
    frm_message(\"Replacing initial definitions of parameters \",
                paste(unique(names(pars)[dup]), collapse = \", \"))
    pars <- pars[!dup]
  }
  is_form <- vapply(pars, inherits, NA, \"formula\")
  out$pforms <- pars[is_form]
  out$pfix <- pars[!is_form]
  check_dpar_equations(out$pfix,
                       c(names(out$pforms), names(out$nlforms)))
  out
}")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
