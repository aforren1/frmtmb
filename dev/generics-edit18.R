root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L)
    stop("many in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

# BLOCKER B. Six as_draws* generics and four size generics had no
# frmtmb_fit method at all. On the base build that gave R's own "no
# applicable method"; once the generic is posterior's, a frmtmb_fit is
# a bare list and posterior's own default walks into it and dies with
# "All list elements must be lists themselves", which names neither
# package. The refusal ?as_draws already promises now exists.
blk <- paste0(
'#\' The refusal a `frmtmb_fit` gets from every draws accessor.\n',
'#\'\n',
'#\' Core is maximum likelihood: there is one parameter vector and no\n',
'#\' chains, so there is nothing for these to convert. The method has to\n',
'#\' EXIST rather than be left to fall through, because the exported\n',
'#\' generic is posterior\'s whenever posterior is loaded, and a\n',
'#\' `frmtmb_fit` is a bare list that posterior\'s own default walks into\n',
'#\' and fails inside. Measured before this method existed:\n',
'#\' `as_draws(fit)` gave "All list elements must be lists themselves".\n',
'#\'\n',
'#\' @noRd\n',
'fit_no_draws <- function(fn) {\n',
'  stop(fn, "() needs posterior draws and a frmtmb_fit has none: ",\n',
'       "frm() is maximum likelihood, so it carries one parameter ",\n',
'       "vector rather than chains. Install frmtmb.sample and sample ",\n',
'       "first, with ", fn, "(frmtmb.sample::frm_sample(fit)); for the ",\n',
'       "point estimates and their covariance use fixef(), vcov() or ",\n',
'       "confint() on the fit itself", call. = FALSE)\n',
'}\n',
'\n',
'#\' @rdname as_draws\n',
'#\' @exportS3Method posterior::as_draws\n',
'#\' @export\n',
'as_draws.frmtmb_fit <- function(x, ...) fit_no_draws("as_draws")\n',
'\n',
'#\' @rdname as_draws\n',
'#\' @exportS3Method posterior::as_draws_matrix\n',
'#\' @export\n',
'as_draws_matrix.frmtmb_fit <- function(x, ...) {\n',
'  fit_no_draws("as_draws_matrix")\n',
'}\n',
'\n',
'#\' @rdname as_draws\n',
'#\' @exportS3Method posterior::as_draws_array\n',
'#\' @export\n',
'as_draws_array.frmtmb_fit <- function(x, ...) {\n',
'  fit_no_draws("as_draws_array")\n',
'}\n',
'\n',
'#\' @rdname as_draws\n',
'#\' @exportS3Method posterior::as_draws_df\n',
'#\' @export\n',
'as_draws_df.frmtmb_fit <- function(x, ...) fit_no_draws("as_draws_df")\n',
'\n',
'#\' @rdname as_draws\n',
'#\' @exportS3Method posterior::as_draws_list\n',
'#\' @export\n',
'as_draws_list.frmtmb_fit <- function(x, ...) {\n',
'  fit_no_draws("as_draws_list")\n',
'}\n',
'\n',
'#\' @rdname as_draws\n',
'#\' @exportS3Method posterior::as_draws_rvars\n',
'#\' @export\n',
'as_draws_rvars.frmtmb_fit <- function(x, ...) {\n',
'  fit_no_draws("as_draws_rvars")\n',
'}\n',
'\n')

sub1("R/draws-generics.R",
"#' Size of a draws object\n",
paste0(blk, "#' Size of a draws object\n"))

# the same hazard on the four size generics
sub1("R/draws-generics.R",
"nvariables <- function(x, ...) UseMethod(\"nvariables\")",
paste0("nvariables <- function(x, ...) UseMethod(\"nvariables\")\n",
       "\n",
       "#' @rdname draws-dimensions\n",
       "#' @exportS3Method posterior::ndraws\n",
       "#' @export\n",
       "ndraws.frmtmb_fit <- function(x) fit_no_draws(\"ndraws\")\n",
       "\n",
       "#' @rdname draws-dimensions\n",
       "#' @exportS3Method posterior::nchains\n",
       "#' @export\n",
       "nchains.frmtmb_fit <- function(x) fit_no_draws(\"nchains\")\n",
       "\n",
       "#' @rdname draws-dimensions\n",
       "#' @exportS3Method posterior::niterations\n",
       "#' @export\n",
       "niterations.frmtmb_fit <- function(x) fit_no_draws(\"niterations\")\n",
       "\n",
       "#' @rdname draws-dimensions\n",
       "#' @exportS3Method posterior::nvariables\n",
       "#' @export\n",
       "nvariables.frmtmb_fit <- function(x, ...) {\n",
       "  fit_no_draws(\"nvariables\")\n",
       "}\n"))

# the as_draws page's own example now has a fit to show it on
sub1("R/draws-generics.R",
paste0("#' @examples\n",
       "#' # frm_multiple() pools estimates rather than carrying draws, so it\n",
       "#' # answers with the reason rather than a matrix\n",
       "#' dd <- data.frame(y = rnorm(40), x = rnorm(40))\n",
       "#' fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))\n",
       "#' try(as_draws(fits))\n"),
paste0("#' @examples\n",
       "#' # a maximum-likelihood fit carries no draws, and says so by name\n",
       "#' dd <- data.frame(y = rnorm(40), x = rnorm(40))\n",
       "#' try(as_draws_df(frm(bf(y ~ x) + gaussian(), data = dd)))\n",
       "#'\n",
       "#' # frm_multiple() pools estimates rather than carrying draws, so it\n",
       "#' # answers with the reason rather than a matrix\n",
       "#' fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))\n",
       "#' try(as_draws(fits))\n"))
cat("DONE\n")
