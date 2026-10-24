# Fix round 1, A1: the block checks, appended to R/covstruct-seam.R of
# the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
p <- file.path(root, "R/covstruct-seam.R")
x <- readLines(p, warn = FALSE)
if (any(grepl("^cs_check_block <- function", x))) {
  cat("fix-a1 already applied to", root, "\n")
  quit(save = "no")
}
add <- '
#\' The number of theta entries a structure takes for one block, checked:
#\' a whole number, zero or more. A slip there would misalign every
#\' block\'s theta after it.
#\'
#\' @noRd
cs_block_npar <- function(name, blk) {
  np <- cs_entry(name)$npar(blk)
  if (!is.numeric(np) || length(np) != 1L || is.na(np) || np < 0 ||
      np != trunc(np)) {
    frm_stop("The \'", name, "\' structure\'s npar slot gave ",
             paste(format(np), collapse = ", "), " for ",
             blk[["term_label"]] %||% "a block", "; it must be one whole ",
             "number, zero or more", call. = FALSE)
  }
  as.integer(np)
}

#\' Check one block against its structure\'s entry, once, when the frame
#\' is built: the slip a registration can make shows here as a refusal
#\' that names the structure and the slot, and not as a fit that runs
#\' and reports a number the structure never meant.
#\'
#\' Every structure: `start` has one value per theta entry, `sd_idx`
#\' names positions within them, and `vcov` at the start is a finite
#\' `dim x dim` matrix (with the coefficient names as dimnames, for a
#\' structure that is not a field). A structure another package
#\' registers, in addition: `lpdf` at the start is one finite number,
#\' and larger at `b = 0` than far from it, which a log density centered
#\' at zero is and a negative log density is not. The built-in
#\' structures skip the `lpdf` evaluation, which for a dense field costs
#\' a factorization, because their densities are tested directly.
#\'
#\' @noRd
cs_check_block <- function(bk) {
  name <- bk[["covstruct"]]
  e <- cs_entry(name)
  lab <- bk[["term_label"]] %||% name
  bad <- function(slot, what) {
    frm_stop("The \'", name, "\' structure\'s ", slot, " slot ", what,
             " for ", lab, call. = FALSE)
  }
  np <- length(bk[["theta_idx"]])
  st <- e$start(bk)
  if (!is.numeric(st) || length(st) != np || any(!is.finite(st))) {
    bad("start", paste0("gave ", length(st), " value(s) where npar is ",
                        np, ", or a value that is not finite"))
  }
  si <- e$sd_idx(bk)
  if (length(si) && (!is.numeric(si) || anyNA(si) || any(si != trunc(si)) ||
                       any(si < 1 | si > np))) {
    bad("sd_idx", paste0("names position(s) ", paste(si, collapse = ", "),
                         " outside the block\'s ", np, " theta entries"))
  }
  V <- e$vcov(st, bk)
  d <- bk[["dim"]]
  if (!is.matrix(V) || !is.numeric(V) || any(dim(V) != d) ||
      any(!is.finite(V))) {
    bad("vcov", paste0("gave ", if (is.matrix(V)) {
      paste(dim(V), collapse = " x ")
    } else {
      paste("a", class(V)[1L])
    }, " at the start, where the block has dimension ", d, " (a finite ",
    d, " x ", d, " matrix)"))
  }
  if (is.null(e[["natural"]]) && !is.null(bk[["cnms"]]) &&
      !identical(unname(rownames(V)), bk[["cnms"]])) {
    bad("vcov", paste0("gave dimnames ",
                       paste(rownames(V) %||% "NULL", collapse = ", "),
                       " where the coefficients are ",
                       paste(bk[["cnms"]], collapse = ", ")))
  }
  if (name %in% names(covstruct_core)) return(invisible(TRUE))
  nb <- length(bk[["b_idx"]])
  l0 <- e$nll(numeric(nb), st, bk)
  if (!is.numeric(l0) || length(l0) != 1L || !is.finite(l0)) {
    bad("lpdf", "did not give one finite number at b = 0 and the start")
  }
  l1 <- e$nll(rep(1e3, nb), st, bk)
  if (is.numeric(l1) && length(l1) == 1L && !is.na(l1) && l1 >= l0) {
    bad("lpdf", paste0("is no smaller at b = 1000 than at b = 0 (",
                       format(l1, digits = 4), " against ",
                       format(l0, digits = 4), "); it has to be the LOG ",
                       "DENSITY of b, not its negative"))
  }
  invisible(TRUE)
}'
writeLines(c(x, add), p)
cat("fix-a1 applied to", root, "\n")
