# The public seam for covariance structures.
#
# frmtmb_register_covstruct() is the third registration seam after
# frmtmb_register_aterm() and frmtmb_register_frame_check(), and it is
# built the same way: a contributing package calls it from its own
# .onLoad(), the call checks everything it is given before it commits
# anything, and what it commits is one entry in the registry every
# consumer reads through cs_entry() (R/covstruct.R). The built-in
# structures are entries of the same shape; there is no second path a
# registered structure misses.

#' The names a registered structure may not take: the names the rest of
#' frmtmb, brms and the other extensions own. A structure is a formula
#' special, and the parser aliases every call to a special's name in a
#' formula, so taking one of these would give one spelling two meanings
#' in one formula: a structure named `poly` stopped `poly(x, 2)` from
#' working anywhere.
#'
#' Reserved: core's formula vocabulary and internal structure names, the
#' whole brms vocabulary, core's families, the specials, addition terms
#' and families core lists as other extensions' (only their home package
#' registers those, `covstruct_home_owner()`), and the addition terms
#' and residual correlation terms registered in this session.
#' `covstruct_formula_terms` adds the
#' functions that formulas call as terms.
#'
#' @noRd
covstruct_reserved_names <- function() {
  unique(c(brms_whole_term_specials, "mo", "mi", "me", "gr", "mm",
           "mmc", "propto", "rr", "us_t", "diag_t", "smooth",
           core_aterms,
           unlist(brms_vocabulary, use.names = FALSE),
           names(family_registry), names(extension_family_homes),
           names(extension_aterm_homes), names(extension_special_homes),
           names(frmtmb_aterm_registry$reg),
           names(frmtmb_autocor_registry$reg)))
}

#' The functions that formulas call as terms or around a term's
#' variable, which a structure of the same name would take from every
#' formula: `poly(x, 2)` stops working once `poly` is a structure. It
#' is a fixed list, and not every function of base, stats, splines or
#' mgcv, because a structure conflicts only with a call a formula
#' makes: `toeplitz`, `kernel` or `dist` is free. Being fixed, it also
#' loads no package to decide.
#'
#' @noRd
covstruct_formula_terms <- c(
  # splines, polynomials and mgcv's smooths
  "poly", "polym", "ns", "bs", "s", "t2", "te", "ti",
  # terms with a meaning of their own in a model formula
  "offset", "I", "C", "cbind", "Surv", "strata", "cluster",
  # transformations of a variable
  "log", "log1p", "log2", "log10", "exp", "expm1", "sqrt", "abs",
  "sin", "cos", "tan", "plogis", "qlogis", "pnorm", "qnorm",
  "round", "floor", "ceiling", "pmin", "pmax", "ifelse", "scale",
  # recodings of a variable
  "factor", "ordered", "as.factor", "as.ordered", "as.numeric",
  "as.integer", "as.logical", "interaction", "cut", "relevel"
)

#' The package that owns a reserved name core lists as an extension's,
#' or `NULL` when no extension owns it.
#'
#' @noRd
covstruct_home_owner <- function(name) {
  homes <- c(extension_special_homes, extension_aterm_homes,
             extension_family_homes)
  if (name %in% names(homes)) unname(homes[[name]]) else NULL
}

#' Register a covariance structure from another package
#'
#' Adds a structured covariance, or a precision matrix that depends on
#' parameters, to the structures [frm()] accepts. Call it from the
#' contributing package's `.onLoad()`. A structure is written in a
#' formula in one of two ways:
#'
#' * as a BAR STRUCTURE, `name(terms | g)`, like `ar1(0 + t | g)`: the
#'   block has one coefficient per term at each level of `g`, and the
#'   structure says how the coefficients of one level covary. Give no
#'   `term`.
#' * as a FIELD, a predictor special such as `car(M, gr = g)`: the block
#'   is one intercept per level of the grouping variable `gr`, and the
#'   structure says how the levels covary, usually through a fixed
#'   matrix the special names. Give `term`.
#'
#' Everything else frmtmb does with a block reads the entry this call
#' builds: the objective, the standard errors and their checks,
#' [VarCorr()], [confint_varcorr()], [ranef()], [predict()] with and
#' without new levels, [simulate()], [conditional_effects()],
#' [vcov_cluster()] and `frmtmb.sample::frm_sample()`. The built-in
#' structures are entries of the same shape.
#'
#' @section The block:
#' Each function receives `blk`, the block: a list with `dim` (the
#' number of coefficients per level), `n_levels`, `levels`, `cnms` (the
#' coefficient names), `theta_idx` and `b_idx`, and for a field `aux`,
#' the fixed data `term$build` returned. `b` is level-major: the `dim`
#' coefficients of one level are contiguous, so `matrix(b, nrow = dim)`
#' has one level per column. Read the block with `[[`, never `$`; see
#' [frm_hazard_reads()].
#'
#' `lpdf` runs on the TAPE: `b` and `theta` are RTMB advectors there.
#' Use `RTMB::matrix()` and `RTMB::dmvnorm()` or `RTMB::dgmrf()`, and
#' write the density without loops over levels. The other functions run
#' on numbers.
#'
#' @section A field term:
#' `term` is a list:
#' \describe{
#'   \item{`args`}{The special's argument names in positional order. It
#'     must include `"gr"`, the grouping variable. frmtmb matches a call
#'     against them by position and by exact name, and refuses an
#'     unknown or misspelled argument.}
#'   \item{`parse`}{`function(args, env)`: `args` is the named list of
#'     unevaluated argument expressions (`NULL` for an absent one), and
#'     `env` the formula environment. Check the arguments and return a
#'     list, the spec, which `build` receives. Do not evaluate data
#'     here; the data are not known yet.}
#'   \item{`build`}{`function(spec, gr, lookup)`: `gr` holds the grouping
#'     values of the rows, and `lookup(expr, what)` resolves an argument
#'     expression from `data2`, the data and the formula environment,
#'     as `gr(cov = A)` resolves its matrix. Return
#'     `list(levels, index, aux)`: `levels` the block's level labels,
#'     `index` the level of each row (whole numbers into `levels`), and
#'     `aux` the fixed data the other functions read as `blk[["aux"]]`.
#'     Add `expands = TRUE` when the block's coefficients are a fixed
#'     linear map of its parameters (see `expand`).}
#'   \item{`usage`}{How the special is written, for messages:
#'     `"car(M, gr = g)"`.}
#'   \item{`levels_of`}{Optional. `function(gr)` giving the level labels
#'     of grouping values, when they are not `as.character(gr)`. It is
#'     used on new data too.}
#'   \item{`levels_complete`}{Optional. `TRUE` when `levels` is the whole
#'     index set of the structure (every mesh node), so a value outside
#'     it is not an unseen level.}
#' }
#'
#' @section Rows a field reports:
#' A field's theta is usually not a set of brms `sd_<group>__<coef>`
#' and `cor_` parameters. Give `natural` for such a structure:
#' `function(theta, blk)` returning `list(term, type, est, grad)`, one
#' entry per reported row. `type` is `"sd"` or `"range"` (estimated on
#' the log scale), `"prop"` (logit), `"cor"` (Fisher z) or `"raw"`;
#' `est` is the value on that scale and `grad` its jacobian in `theta`,
#' a matrix with one row per term. [confint_varcorr()] reports the rows
#' with delta-method intervals. A structure with `natural` is a field
#' for brms's names too: [ngrps()], [VarCorr()] and [variables()] do not
#' count it as a grouping factor.
#'
#' @section Testing a structure:
#' Compare the fitted log-likelihood with one you compute independently:
#' for a gaussian response, integrate the block out by hand and evaluate
#' the marginal density with `mvtnorm::dmvnorm()` or a dense
#' `solve()`. A test that registers calls [frmtmb_local_registries()]
#' first, which restores every registry when the test ends.
#'
#' @param name The structure's name, as a formula writes it, without
#'   parentheses. It must be a syntactic name that nothing else in a
#'   formula means: not one of frmtmb's own structures, specials,
#'   families or addition terms, not a brms name, not a name core lists
#'   as another extension's (only that package registers it), and not a
#'   function that formulas call as a term, such as `poly`, `ns`, `log`
#'   or `factor`.
#' @param replace A registered name is replaced when the package that
#'   registered it registers it again, so a reload is not an error. A
#'   name another package, or the session, registered is replaced only
#'   with `replace = TRUE`; without it the call refuses.
#' @param package The package that owns the structure. A fit records it
#'   for each block, so that a fit read where that package is not
#'   installed is refused by the package's name before any method uses
#'   the structure. The default is the package whose namespace calls,
#'   and `NULL` outside a package.
#' @param lpdf `function(b, theta, blk)`: the log density of the block's
#'   `b`, summed over levels. AD-safe. The frame checks only that it is
#'   finite at `b = 0` and the start: a correct density can be flat along
#'   a direction, centered away from zero or bimodal, so its shape tells
#'   nothing, and a density with the wrong sign is for the structure's
#'   own test to catch ("Testing a structure").
#' @param vcov `function(theta, blk)`: the numeric `dim x dim` covariance
#'   within one level, with the coefficient names as dimnames.
#'   [VarCorr()] reports it, and an unseen level is drawn from it unless
#'   `levels_coupled` is `TRUE`.
#' @param npar The number of entries the structure takes in `theta`: a
#'   whole number, or `function(blk)`.
#' @param start Starting values for theta: a numeric vector, or
#'   `function(blk)`. The default is zeros.
#' @param sd_idx The positions within the block's theta that are LOG
#'   standard deviations: an integer vector, or `function(blk)`. They
#'   take the SE checks' boundary reading and class `"sd"` priors. The
#'   default is none.
#' @param natural `NULL`, or `function(theta, blk)` for a field; see
#'   "Rows a field reports".
#' @param draw `NULL`, or `function(theta, blk)` returning one draw of
#'   the block's whole `b`. Give it when the levels are correlated; the
#'   default draws each level from `vcov`.
#' @param levels_coupled `TRUE` when a fixed matrix correlates the
#'   LEVELS of the block. The block then has no law for an unseen level
#'   (prediction gives one no block contribution), it is one unit for
#'   [vcov_cluster()] and the sandwich, and `importance =` refuses it.
#' @param expand `NULL`, or a list for a block whose coefficients are a
#'   fixed linear map `P b` of its parameters: `apply(b, blk)` computes
#'   `P b` (AD-safe), `jacobian(blk)` returns `P` as triplets
#'   `list(i, j, x)` in within-block positions, and `condsd(v, blk)`
#'   gives the conditional standard deviations of `P b` from the
#'   variances `v` of `b` (optional; without it [ranef()] reports `NA`).
#'   It acts on blocks whose `build` returned `expands = TRUE`.
#' @param term `NULL` for a bar structure, or the field term; see "A
#'   field term".
#' @param from_natural `NULL`, or `function(sds, C, blk)` returning the
#'   whole theta segment from natural-scale standard deviations and a
#'   correlation matrix, for [frm_simulate()]'s `newparams` in brms's
#'   names.
#' @param chol `NULL`, or `function(theta, blk)` returning the lower
#'   Cholesky factor of `vcov`, on the tape. With it, and with every
#'   theta entry a standard deviation, `frmtmb.sample::frm_sample()` can
#'   sample the block non-centered.
#' @param full_rank `TRUE` when an effect a row owns can reach any
#'   value with a finite density: the covariance is positive definite at
#'   some parameter value and the density is positive everywhere, as for
#'   a gaussian or Student-t latent. A weighted model whose every row
#'   owns such an effect has no maximum when the weights sum to more
#'   than the rows, and [frm()] refuses it, as it refuses `(1 | obs)`.
#'   `FALSE` when the effects cannot (a singular covariance, a latent
#'   with bounded support). `NULL`, the default, states nothing, and the
#'   check stays open for the structure. A package whose structure is
#'   full rank should state `full_rank = TRUE`, so that a weighted model
#'   with no maximum is refused rather than fitted to a runaway.
#' @return `NULL`, invisibly. Called for the registration.
#' @seealso [frmtmb-extension-api] for the how-to an extension follows,
#'   [frmtmb_register_compat()] for the rows a structure adds to
#'   [frm_compat()] (the structure itself joins the vocabulary as a
#'   feature of kind `"covstruct"` when it is registered),
#'   [frmtmb_register_aterm()] for the seam this one is modeled on, and
#'   [us_sigma()] and [us_theta()] for the parameterization of an
#'   unstructured block, which an `lpdf` can reuse.
#' @examples
#' # an exchangeable block with a FIXED correlation of 0.5 and one
#' # estimated standard deviation: xcs(0 + t | g)
#' xcs_lpdf <- function(b, theta, blk) {
#'   d <- blk[["dim"]]
#'   C <- RTMB::matrix(0.5, d, d) + diag(0.5, d)
#'   Sigma <- exp(2 * theta[1]) * C
#'   dim(b) <- c(d, length(b) %/% d)
#'   sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
#' }
#' xcs_vcov <- function(theta, blk) {
#'   d <- blk[["dim"]]
#'   V <- exp(2 * theta[1]) * (matrix(0.5, d, d) + diag(0.5, d))
#'   dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
#'   V
#' }
#' \dontrun{
#' # in the contributing package's .onLoad()
#' frmtmb_register_covstruct("xcs", lpdf = xcs_lpdf, vcov = xcs_vcov,
#'                           npar = 1L, sd_idx = 1L)
#' }
#' @export
frmtmb_register_covstruct <- function(name, lpdf, vcov, npar, start = NULL,
                                      sd_idx = NULL, natural = NULL,
                                      draw = NULL, levels_coupled = FALSE,
                                      expand = NULL, term = NULL,
                                      from_natural = NULL, chol = NULL,
                                      full_rank = NULL, replace = FALSE,
                                      package = compat_registrant_name(
                                        parent.frame())) {
  arg <- function(nm) paste0("frmtmb_register_covstruct(", nm, " =)")
  package <- seam_check_package(package, arg("package"))
  check_flag(replace, arg("replace"))
  seam_check_name(name, arg("name"))
  if (name %in% names(covstruct_core)) {
    frm_stop("`", name, "` is one of frmtmb's own covariance structures ",
             "and cannot be re-registered", call. = FALSE)
  }
  home <- covstruct_home_owner(name)
  if (!is.null(home) && !identical(package, home)) {
    frm_stop("frmtmb_register_covstruct() cannot register '", name,
             "': it names a term of the extension package ", home,
             ", and only ", home, " registers it", call. = FALSE)
  }
  if (is.null(home) && name %in% covstruct_reserved_names()) {
    frm_stop("frmtmb_register_covstruct() cannot register '", name,
             "', because ", name, "() already means something else in a ",
             "formula. Choose another name.", call. = FALSE)
  }
  if (name %in% covstruct_formula_terms) {
    frm_stop("frmtmb_register_covstruct() cannot register '", name,
             "', because formulas call ", name, "() as a term (as in ",
             "poly(x, 2) or log(x)), and a structure of that name would ",
             "take it from every formula. Choose another name.",
             call. = FALSE)
  }
  old <- frmtmb_covstruct_registry$reg[[name]]
  if (!is.null(old) && !replace && !identical(old[["package"]], package)) {
    frm_stop("frmtmb_register_covstruct(): '", name, "' is already ",
             "registered by ",
             if (is.null(old[["package"]])) "this session" else {
               paste0("the package ", old[["package"]])
             }, ". Only that package replaces it on its own; pass ",
             "replace = TRUE to replace it anyway", call. = FALSE)
  }
  fn_or <- function(x, nm, argn, null_ok = TRUE) {
    seam_fn(x, arg(nm), argn, null_ok)
  }
  lpdf <- fn_or(lpdf, "lpdf", 3L, null_ok = FALSE)
  vcov <- fn_or(vcov, "vcov", 2L, null_ok = FALSE)
  natural <- fn_or(natural, "natural", 2L)
  draw <- fn_or(draw, "draw", 2L)
  from_natural <- fn_or(from_natural, "from_natural", 3L)
  chol <- fn_or(chol, "chol", 2L)
  whole <- function(x, nm, zero_ok = TRUE) {
    if (!is.numeric(x) || anyNA(x) || any(x != trunc(x)) ||
        any(x < if (zero_ok) 0 else 1)) {
      frm_stop(arg(nm), " must be a function of the block or whole ",
               "numbers", call. = FALSE)
    }
    as.integer(x)
  }
  npar_f <- if (is.function(npar)) npar else {
    np <- whole(npar, "npar")
    if (length(np) != 1L) {
      frm_stop(arg("npar"), " must be one whole number or a function of ",
               "the block", call. = FALSE)
    }
    function(blk) np
  }
  start_f <- if (is.null(start)) {
    function(blk) numeric(npar_f(blk))
  } else if (is.function(start)) {
    start
  } else {
    if (!is.numeric(start) || !all(is.finite(start))) {
      frm_stop(arg("start"), " must be finite numbers or a function of the ",
               "block", call. = FALSE)
    }
    st <- as.numeric(start)
    function(blk) st
  }
  sd_f <- if (is.null(sd_idx)) {
    function(blk) integer(0)
  } else if (is.function(sd_idx)) {
    sd_idx
  } else {
    si <- whole(sd_idx, "sd_idx", zero_ok = FALSE)
    function(blk) si
  }
  if (!is.logical(levels_coupled) || length(levels_coupled) != 1L ||
      is.na(levels_coupled)) {
    frm_stop(arg("levels_coupled"), " must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(expand)) {
    if (!is.list(expand) || !is.function(expand[["apply"]]) ||
        !is.function(expand[["jacobian"]]) ||
        !(is.null(expand[["condsd"]]) || is.function(expand[["condsd"]]))) {
      frm_stop(arg("expand"), " must be list(apply = function(b, blk), ",
               "jacobian = function(blk)), with condsd = function(v, blk) ",
               "optional", call. = FALSE)
    }
  }
  if (!is.null(full_rank) && (!is.logical(full_rank) ||
                                length(full_rank) != 1L ||
                                is.na(full_rank))) {
    frm_stop(arg("full_rank"), " must be TRUE, FALSE or NULL",
             call. = FALSE)
  }
  if (!is.null(term)) term <- covstruct_check_term(term, name, arg("term"))
  entry <- list(npar = npar_f, start = start_f, sd_idx = sd_f, nll = lpdf,
                vcov = vcov, natural = natural, draw = draw,
                levels_coupled = levels_coupled, expand = expand,
                term = term, from_natural = from_natural, chol_L = chol,
                full_rank = full_rank, package = package)
  entry <- entry[!vapply(entry, is.null, TRUE)]
  bad <- cs_check_entry(name, entry)
  if (!is.null(bad)) {
    frm_stop("frmtmb_register_covstruct() would register an incomplete ",
             "entry: ", bad, call. = FALSE)
  }
  # The structure joins the compatibility vocabulary as it joins the
  # parser, for the reason an addition term does: a structure frm()
  # accepts and frm_compat() cannot describe is a gap by construction.
  # This is the last step that can refuse, so a refusal leaves both
  # registries as it found them.
  known <- frmtmb_compat_features_tbl()
  add <- compat_new_features(stats::setNames("covstruct", name),
                             "frmtmb_register_covstruct()", known)
  if (length(add)) {
    before <- frmtmb_compat_contrib$features
    frmtmb_compat_contrib$features <- c(before, list(add))
    compat_cache_extend(add, before)
  }
  # Quadrature and importance = refuse every registered structure, and
  # frm_compat() has to say so. The rows are keyed by the structure, so
  # a replacement replaces them rather than leaving a row that wins
  # nothing.
  rules <- covstruct_compat_rules(name, levels_coupled)
  attr(rules, "frmtmb_package") <- package %||% ""
  frmtmb_compat_contrib$rules[[paste0("covstruct:", name)]] <- rules
  frmtmb_covstruct_registry$reg[[name]] <- entry
  invisible(NULL)
}

#' Check a field term's slots and return it in the shape the parser and
#' the frame read.
#'
#' @noRd
covstruct_check_term <- function(term, name, what) {
  bad <- function(msg) frm_stop(what, " ", msg, call. = FALSE)
  if (!is.list(term)) bad("must be a list; see ?frmtmb_register_covstruct")
  a <- term[["args"]]
  if (!is.character(a) || anyNA(a) || !all(nzchar(a)) || anyDuplicated(a) ||
      !"gr" %in% a) {
    bad(paste("needs args, the special's argument names without repeats,",
              "including \"gr\", the grouping variable"))
  }
  for (nm in c("parse", "build")) {
    f <- term[[nm]]
    if (!is.function(f) || length(formals(f)) < c(parse = 2L,
                                                  build = 3L)[[nm]]) {
      bad(paste0("needs ", nm, " = function(",
                 if (nm == "parse") "args, env" else "spec, gr, lookup", ")"))
    }
  }
  if (!is.null(term[["levels_of"]]) && !is.function(term[["levels_of"]])) {
    bad("levels_of must be a function of the grouping values")
  }
  lc <- term[["levels_complete"]] %||% FALSE
  if (!is.logical(lc) || length(lc) != 1L || is.na(lc)) {
    bad("levels_complete must be TRUE or FALSE")
  }
  usage <- term[["usage"]] %||% paste0(name, "(..., gr = g)")
  if (!is.character(usage) || length(usage) != 1L) {
    bad("usage must be one string, how the special is written")
  }
  list(args = a, parse = term[["parse"]], build = term[["build"]],
       levels_of = term[["levels_of"]], levels_complete = lc,
       usage = usage)
}

#' What frmtmb knows about a covariance structure
#'
#' Reads one structure's entry in the registry that
#' [frmtmb_register_covstruct()] fills, so that another package, such as
#' frmtmb.sample, asks the registry instead of keeping its own list of
#' structure names. Works for frmtmb's own structures and for a
#' registered one alike.
#'
#' @param name The structure's name, as a block records it
#'   (`fit$frame$re_blocks[[i]]$covstruct`): `"us"`, `"ar1"`, `"car"`.
#' @return A list:
#'   \describe{
#'     \item{`name`}{The name.}
#'     \item{`package`}{The package that registered it, or `NULL` for
#'       one of frmtmb's own.}
#'     \item{`builtin`}{`TRUE` for one of frmtmb's own structures.}
#'     \item{`field`}{`TRUE` when the structure reports its own rows
#'       (its `natural` slot), so its parameters have no brms
#'       `sd_<group>__<coef>` or `cor_` names and brms counts no grouping
#'       factor for it: smooths, `gp()`, `car()`.}
#'     \item{`term`}{`TRUE` when it is written as a predictor special
#'       rather than a bar structure.}
#'     \item{`levels_coupled`}{`TRUE` when a fixed matrix correlates the
#'       levels.}
#'     \item{`full_rank`}{`TRUE` or `FALSE` as the structure states it,
#'       or `NA` when it states nothing, and the weights check is then
#'       open for it (see `full_rank` in [frmtmb_register_covstruct()]).}
#'   }
#'   A name nothing registered is refused. A name that core lists as an
#'   extension's (`"car"`, `"spde"`) loads that extension, which
#'   registers it, as a fit that uses the structure would. A residual
#'   correlation term (`"fcor"`, frmtmb.spatial's `"sar"`) is no
#'   covariance structure: it is refused with a pointer to
#'   [frmtmb_autocor_pars()].
#' @seealso [frmtmb_register_covstruct()].
#' @examples
#' frmtmb_covstruct_info("ar1")$field
#' frmtmb_covstruct_info("smooth")$field
#' @export
frmtmb_covstruct_info <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    frm_stop("frmtmb_covstruct_info() takes one structure name",
             call. = FALSE)
  }
  # a residual correlation term is registered through another seam, and
  # an extension's (sar()) is in the same table of homes as a structure
  if (name %in% c(ac_entry_names(), extension_autocor_specials)) {
    frm_stop("frmtmb_covstruct_info(): '", name, "' is a residual ",
             "correlation term, registered with frmtmb_register_autocor(), ",
             "not a covariance structure. frmtmb_autocor_pars(fit) lists ",
             "the parameters of a fit's residual correlation terms",
             call. = FALSE)
  }
  e <- frmtmb_covstruct_registry$reg[[name]]
  if (is.null(e) && !is.null(extension_special_home(name))) {
    e <- tryCatch(cs_entry(name), error = keep_missing(NULL))
  }
  if (is.null(e)) {
    frm_stop("frmtmb_covstruct_info(): no covariance structure '", name,
             "' is registered in this session",
             if (!is.null(extension_special_home(name))) {
               paste0("; it is in the extension package ",
                      extension_special_home(name), ", which is not ",
                      "installed or does not load")
             }, call. = FALSE)
  }
  list(name = name, package = e[["package"]],
       builtin = name %in% names(covstruct_core),
       field = !is.null(e[["natural"]]),
       term = !is.null(e[["term"]]),
       levels_coupled = isTRUE(e[["levels_coupled"]]),
       full_rank = if (is.null(e[["full_rank"]])) NA else
         isTRUE(e[["full_rank"]]))
}

#' The checks the registration seams share, one message template each:
#' the owning package, the name as a formula writes it, and a slot that
#' must be a function of so many arguments.
#'
#' @noRd
seam_check_package <- function(package, what) {
  if (identical(package, "frmtmb")) return(NULL)
  if (!is.null(package) && (!is.character(package) ||
                              length(package) != 1L || is.na(package) ||
                              !nzchar(package))) {
    frm_stop(what, " must be one package name, or NULL", call. = FALSE)
  }
  package
}

#' @noRd
seam_check_name <- function(name, what) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name) || !identical(name, make.names(name))) {
    frm_stop(what, " must be one syntactic name, written as it ",
             "appears in a formula and without parentheses", call. = FALSE)
  }
  invisible(name)
}

#' @noRd
seam_fn <- function(x, what, argn, null_ok) {
  if (is.null(x) && null_ok) return(NULL)
  # a function of `...` takes any number of arguments
  if (!is.function(x) || (length(formals(x)) < argn &&
                            !"..." %in% names(formals(x)))) {
    frm_stop(what, " must be a function of ", argn, " argument",
             if (argn > 1L) "s", call. = FALSE)
  }
  x
}

#' The package a registration came from, or `NULL` for a call from
#' outside a package namespace.
#'
#' @noRd
compat_registrant_name <- function(env) {
  top <- topenv(env)
  if (!isNamespace(top)) return(NULL)
  nm <- environmentName(top)
  if (!nzchar(nm)) NULL else nm
}

#' Check one registry entry against the shape every consumer reads.
#' The built-in entries pass it too (test-covstruct-seam.R), which is
#' what says the seam is the path core itself uses.
#'
#' @noRd
cs_check_entry <- function(name, e) {
  need <- c("npar", "start", "sd_idx", "nll", "vcov")
  miss <- need[!vapply(need, function(k) is.function(e[[k]]), TRUE)]
  fns <- c("natural", "draw", "from_natural", "cor_spec", "chol_sd",
           "chol_L", "chol_A")
  badf <- fns[!vapply(fns, function(k) {
    is.null(e[[k]]) || is.function(e[[k]])
  }, TRUE)]
  lc <- e[["levels_coupled"]]
  ok_lc <- is.null(lc) || (is.logical(lc) && length(lc) == 1L && !is.na(lc))
  ok_term <- is.null(e[["term"]]) || all(c("args", "parse", "build", "usage")
                                        %in% names(e[["term"]]))
  ok_ex <- is.null(e[["expand"]]) ||
    all(vapply(c("apply", "jacobian"),
               function(k) is.function(e[["expand"]][[k]]), TRUE))
  if (length(miss) || length(badf) || !ok_lc || !ok_term || !ok_ex) {
    return(paste0(name, ": ", paste(c(
      if (length(miss)) paste("missing", paste(miss, collapse = ", ")),
      if (length(badf)) paste("not functions:", paste(badf, collapse = ", ")),
      if (!ok_lc) "levels_coupled not one logical",
      if (!ok_term) "term without args, parse, build and usage",
      if (!ok_ex) "expand without apply and jacobian"), collapse = "; ")))
  }
  NULL
}

#' The number of theta entries a structure takes for one block, checked:
#' a whole number, zero or more. A slip there would misalign every
#' block's theta after it.
#'
#' @noRd
cs_block_npar <- function(name, blk) {
  np <- cs_entry(name)$npar(blk)
  if (!is.numeric(np) || length(np) != 1L || is.na(np) || np < 0 ||
      np != trunc(np)) {
    frm_stop("The '", name, "' structure's npar slot gave ",
             paste(format(np), collapse = ", "), " for ",
             blk[["term_label"]] %||% "a block", "; it must be one whole ",
             "number, zero or more", call. = FALSE)
  }
  as.integer(np)
}

#' Check one block against its structure's entry, once, when the frame
#' is built: the slip a registration can make shows here as a refusal
#' that names the structure and the slot, and not as a fit that runs
#' and reports a number the structure never meant.
#'
#' Every structure: `start` has one value per theta entry, `sd_idx`
#' names positions within them, and `vcov` at the start is a finite
#' `dim x dim` matrix (with the coefficient names as dimnames, for a
#' structure that is not a field). A structure another package
#' registers, in addition: `lpdf` at `b = 0` and the start is one finite
#' number. Nothing is asked of its shape. A correct density can be flat
#' along a direction (an intrinsic RW1 or ICAR along the constant
#' vector), centered away from zero, or bimodal with a minimum at zero,
#' so no test of its level or curvature at zero separates it from a
#' density with the wrong sign; the author's own test of the structure
#' does. The built-in structures skip the `lpdf` evaluation, which for
#' a dense field costs a factorization, because their densities are
#' tested directly.
#'
#' @noRd
cs_check_block <- function(bk) {
  name <- bk[["covstruct"]]
  e <- cs_entry(name)
  lab <- bk[["term_label"]] %||% name
  bad <- function(slot, what) {
    frm_stop("The '", name, "' structure's ", slot, " slot ", what,
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
                         " outside the block's ", np, " theta entries"))
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
  invisible(TRUE)
}

#' The frm_compat() rows frmtmb writes for a registered structure.
#'
#' @noRd
covstruct_compat_rules <- function(name, levels_coupled) {
  force(name)
  force(levels_coupled)
  function() {
    b <- compat_rule_builder()
    b$r("importance", name, "refused",
        if (isTRUE(levels_coupled)) compat_note_crosslevel_importance
        else compat_note_scope_importance)
    b$r("quadrature", name, "refused", compat_note_wide_quadrature)
    b$rules()
  }
}
