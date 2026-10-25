# Fix round 1, part B: F10 (us_t, diag_t are full rank), F11 (an unstated
# full_rank leaves the weights check open) and F14 (the reader's error
# wording and its side effect). Rscript <this>
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
edit <- function(p, pairs) {
  p <- file.path(wt, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
edit("R/covstruct.R", list(
  c("# The built-in structures that are FULL RANK: their covariance is
# positive definite at some parameter value, so an effect a row owns
# can reach any value with a finite density. The weights check
# (R/weights-unbounded.R) reads it; a registered structure declares it
# with frmtmb_register_covstruct(full_rank =), or has it read off its
# covariance at the start (cs_full_rank()). rr() and the structures with
# a fixed or possibly singular covariance are not; the Student-t
# latents are left out, as they were from the weights check before the
# property existed, because no fit measured them.
for (nm_ in names(covstruct_core)) {
  covstruct_core[[nm_]][[\"full_rank\"]] <- nm_ %in% c(
    \"us\", \"diag\", \"homdiag\", \"cs\", \"homcs\", \"ar1\", \"hetar1\", \"ou\", \"toep\",
    \"homtoep\", \"exp\", \"gau\", \"mat\")
}",
"# The built-in structures that are FULL RANK: their covariance is
# positive definite at some parameter value, so an effect a row owns
# can reach any value with a finite density. The weights check
# (R/weights-unbounded.R) reads it; a registered structure states it
# with frmtmb_register_covstruct(full_rank =). rr() and the structures
# with a fixed or possibly singular covariance are not. The Student-t
# latents are: a t density is positive everywhere, and a weighted
# Gamma model with (1 | gr(obs, dist = \"student\")) ran to a logLik of
# 4.7e7 before they were listed here (test-weights-unbounded.R).
for (nm_ in names(covstruct_core)) {
  covstruct_core[[nm_]][[\"full_rank\"]] <- nm_ %in% c(
    \"us\", \"diag\", \"homdiag\", \"cs\", \"homcs\", \"ar1\", \"hetar1\", \"ou\", \"toep\",
    \"homtoep\", \"exp\", \"gau\", \"mat\", \"us_t\", \"diag_t\")
}"),
  c("#' Whether a block's covariance is positive definite at some parameter
#' value, so that an effect a row owns can reach any value with a finite
#' density (the weights check, R/weights-unbounded.R). The entry's
#' `full_rank` when it states one, as every built-in entry does. For a
#' registered structure that does not, it is read off the structure's
#' covariance at the start: a bar structure whose levels are not
#' coupled and whose `vcov` there is positive definite. A field, a
#' structure with coupled levels and a covariance that is singular at
#' the start read FALSE, which leaves the check open, as it was for
#' every registered structure before the property existed.
#'
#' @noRd
cs_full_rank <- function(bk) {
  e <- cs_entry(bk[[\"covstruct\"]])
  fr <- e[[\"full_rank\"]]
  if (!is.null(fr)) return(isTRUE(fr))
  if (!is.null(e[[\"natural\"]]) || isTRUE(e[[\"levels_coupled\"]])) {
    return(FALSE)
  }
  V <- tryCatch(e$vcov(e$start(bk), bk), error = function(err) NULL)
  if (!is.matrix(V) || !is.numeric(V) || any(!is.finite(V)) ||
        nrow(V) != ncol(V) || !nrow(V)) {
    return(FALSE)
  }
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
  isTRUE(min(ev) > sqrt(.Machine$double.eps) * max(abs(ev)))
}",
"#' Whether a block's covariance is positive definite at some parameter
#' value, so that an effect a row owns can reach any value with a finite
#' density (the weights check, R/weights-unbounded.R): the entry's
#' `full_rank`, which every built-in entry states. A registered
#' structure that does not state it reads FALSE, which leaves the check
#' open, as it was for every registered structure before the property
#' existed. It is not inferred: a positive definite covariance says
#' nothing of the latent's support, and a latent with bounded support
#' (a scaled beta) keeps the weighted likelihood bounded, so an
#' inference from the covariance refused a model with a maximum.
#'
#' @noRd
cs_full_rank <- function(bk) {
  isTRUE(cs_entry(bk[[\"covstruct\"]])[[\"full_rank\"]])
}")))
edit("R/covstruct-seam.R", list(
  c("#' @param full_rank `TRUE` when the block's covariance is positive
#'   definite at some parameter value, so that an effect a row owns can
#'   reach any value with a finite density; `FALSE` when it cannot. A
#'   weighted model whose every row owns such an effect has no maximum
#'   when the weights sum to more than the rows, and [frm()] refuses it.
#'   `NULL`, the default, reads it off `vcov` at the start: a bar
#'   structure whose levels are not coupled and whose covariance there is
#'   positive definite is full rank, and a field is not.",
"#' @param full_rank `TRUE` when an effect a row owns can reach any
#'   value with a finite density: the covariance is positive definite at
#'   some parameter value and the density is positive everywhere, as for
#'   a gaussian or Student-t latent. A weighted model whose every row
#'   owns such an effect has no maximum when the weights sum to more
#'   than the rows, and [frm()] refuses it, as it refuses `(1 | obs)`.
#'   `FALSE` when the effects cannot (a singular covariance, a latent
#'   with bounded support). `NULL`, the default, states nothing, and the
#'   check stays open for the structure. A package whose structure is
#'   full rank should state `full_rank = TRUE`, so that a weighted model
#'   with no maximum is refused rather than fitted to a runaway."),
  c("#'     \\item{`full_rank`}{`TRUE` or `FALSE` as the structure states it,
#'       or `NA` when it is read off the covariance of each block (see
#'       `full_rank` in [frmtmb_register_covstruct()]).}
#'   }
#'   A name nothing registered is refused.",
"#'     \\item{`full_rank`}{`TRUE` or `FALSE` as the structure states it,
#'       or `NA` when it states nothing, and the weights check is then
#'       open for it (see `full_rank` in [frmtmb_register_covstruct()]).}
#'   }
#'   A name nothing registered is refused. A name that core lists as an
#'   extension's (`\"car\"`, `\"spde\"`) loads that extension, which
#'   registers it, as a fit that uses the structure would."),
  c("  e <- cs_entry(name)
  list(name = name, package = e[[\"package\"]],",
"  e <- frmtmb_covstruct_registry$reg[[name]]
  if (is.null(e) && !is.null(extension_special_home(name))) {
    e <- tryCatch(cs_entry(name), error = function(err) NULL)
  }
  if (is.null(e)) {
    frm_stop(\"frmtmb_covstruct_info(): no covariance structure '\", name,
             \"' is registered in this session\",
             if (!is.null(extension_special_home(name))) {
               paste0(\"; it is in the extension package \",
                      extension_special_home(name), \", which is not \",
                      \"installed or does not load\")
             }, call. = FALSE)
  }
  list(name = name, package = e[[\"package\"]],")))
cat("ok\n")
