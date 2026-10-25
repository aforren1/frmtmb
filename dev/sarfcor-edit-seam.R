# One-off source edit for lane sarfcor, part B: frmtmb_register_covstruct
# (full_rank =) and the reader export frmtmb_covstruct_info().
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/R/covstruct-seam.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("#' @param chol `NULL`, or `function(theta, blk)` returning the lower
#'   Cholesky factor of `vcov`, on the tape. With it, and with every
#'   theta entry a standard deviation, `frmtmb.sample::frm_sample()` can
#'   sample the block non-centered.",
"#' @param chol `NULL`, or `function(theta, blk)` returning the lower
#'   Cholesky factor of `vcov`, on the tape. With it, and with every
#'   theta entry a standard deviation, `frmtmb.sample::frm_sample()` can
#'   sample the block non-centered.
#' @param full_rank `TRUE` when the block's covariance is positive
#'   definite at some parameter value, so that an effect a row owns can
#'   reach any value with a finite density; `FALSE` when it cannot. A
#'   weighted model whose every row owns such an effect has no maximum
#'   when the weights sum to more than the rows, and [frm()] refuses it.
#'   `NULL`, the default, reads it off `vcov` at the start: a bar
#'   structure whose levels are not coupled and whose covariance there is
#'   positive definite is full rank, and a field is not.")
rep1("                                      from_natural = NULL, chol = NULL,
                                      replace = FALSE,",
"                                      from_natural = NULL, chol = NULL,
                                      full_rank = NULL, replace = FALSE,")
rep1("  if (!is.null(term)) term <- covstruct_check_term(term, name, arg(\"term\"))",
"  if (!is.null(full_rank) && (!is.logical(full_rank) ||
                                length(full_rank) != 1L ||
                                is.na(full_rank))) {
    frm_stop(arg(\"full_rank\"), \" must be TRUE, FALSE or NULL\",
             call. = FALSE)
  }
  if (!is.null(term)) term <- covstruct_check_term(term, name, arg(\"term\"))")
rep1("                term = term, from_natural = from_natural, chol_L = chol,
                package = package)",
"                term = term, from_natural = from_natural, chol_L = chol,
                full_rank = full_rank, package = package)")
rep1("#' The package a registration came from, or `NULL` for a call from",
"#' What frmtmb knows about a covariance structure
#'
#' Reads one structure's entry in the registry that
#' [frmtmb_register_covstruct()] fills, so that another package, such as
#' frmtmb.sample, asks the registry instead of keeping its own list of
#' structure names. Works for frmtmb's own structures and for a
#' registered one alike.
#'
#' @param name The structure's name, as a block records it
#'   (`fit$frame$re_blocks[[i]]$covstruct`): `\"us\"`, `\"ar1\"`, `\"car\"`.
#' @return A list:
#'   \\describe{
#'     \\item{`name`}{The name.}
#'     \\item{`package`}{The package that registered it, or `NULL` for
#'       one of frmtmb's own.}
#'     \\item{`builtin`}{`TRUE` for one of frmtmb's own structures.}
#'     \\item{`field`}{`TRUE` when the structure reports its own rows
#'       (its `natural` slot), so its parameters have no brms
#'       `sd_<group>__<coef>` or `cor_` names and brms counts no grouping
#'       factor for it: smooths, `gp()`, `car()`.}
#'     \\item{`term`}{`TRUE` when it is written as a predictor special
#'       rather than a bar structure.}
#'     \\item{`levels_coupled`}{`TRUE` when a fixed matrix correlates the
#'       levels.}
#'     \\item{`full_rank`}{`TRUE` or `FALSE` as the structure states it,
#'       or `NA` when it is read off the covariance of each block (see
#'       `full_rank` in [frmtmb_register_covstruct()]).}
#'   }
#'   A name nothing registered is refused.
#' @seealso [frmtmb_register_covstruct()].
#' @examples
#' frmtmb_covstruct_info(\"ar1\")$field
#' frmtmb_covstruct_info(\"smooth\")$field
#' @export
frmtmb_covstruct_info <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    frm_stop(\"frmtmb_covstruct_info() takes one structure name\",
             call. = FALSE)
  }
  e <- cs_entry(name)
  list(name = name, package = e[[\"package\"]],
       builtin = name %in% names(covstruct_core),
       field = !is.null(e[[\"natural\"]]),
       term = !is.null(e[[\"term\"]]),
       levels_coupled = isTRUE(e[[\"levels_coupled\"]]),
       full_rank = if (is.null(e[[\"full_rank\"]])) NA else
         isTRUE(e[[\"full_rank\"]]))
}

#' The package a registration came from, or `NULL` for a call from")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
