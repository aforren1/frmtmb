# The public seam for residual correlation terms.
#
# A residual correlation term such as brms's fcor() or sar() does not
# add a latent vector to the linear predictor, which is what a
# covariance structure (frmtmb_register_covstruct()) does. It changes
# the density of the response itself: the residuals of all the rows
# are one joint draw, and for brms's lag SAR the mean of y is not even
# mu. So the covariance-structure seam does not fit it, and this seam
# does. A term supplies one affine map, from the response to
# independent standardized innovations ("whiten"), and its inverse
# ("unwhiten"). Core writes the density from the map (a gaussian or a
# multivariate Student-t with the map's log-Jacobian), the mean (the
# inverse at zero innovations) and the draws (the inverse at drawn
# innovations), so a term never writes a density, and the gaussian and
# student() cases cannot drift apart.
#
# fcor() is core's own term of this shape, and the registry's first
# state; frmtmb.spatial registers sar().

#' The names a registered residual correlation term may not take: core's
#' own terms, the names of the covariance-structure registry, and the
#' names core reserves for structures (`covstruct_reserved_names()`).
#' A home that core lists for the name is the one package that may
#' register it.
#'
#' @noRd
autocor_reserved_names <- function() {
  # a term registered here is reserved for structures, not for itself
  setdiff(unique(c(autocor_structs, names(autocor_core), cs_names(),
                   covstruct_reserved_names(), covstruct_formula_terms)),
          setdiff(ac_entry_names(), names(autocor_core)))
}

#' Register a residual correlation term from another package
#'
#' Adds a term that correlates the residuals of a `gaussian()` or
#' `student()` response, such as brms's `sar()`, to the terms [frm()]
#' accepts. Call it from the contributing package's `.onLoad()`. The
#' term is written in the formula of `mu`, `y ~ x + name(...)`, as
#' `ar()` and `fcor()` are.
#'
#' A term is ONE affine map. `whiten` takes the response to innovations
#' `u` that are independent with unit scale; `unwhiten` takes
#' innovations back to a response. frmtmb builds everything else from
#' the two:
#'
#' * the density: for `gaussian()`, \eqn{\sum_i \log \phi(u_i)} plus the
#'   log-Jacobian of the map; for `student()`, the multivariate
#'   Student-t of `u` with scale matrix `I` and the response's `nu`,
#'   plus the same log-Jacobian;
#' * the mean of the response, `unwhiten(0, ...)`: for brms's lag SAR
#'   this is \eqn{(I - \rho W)^{-1} \mu}, which [fitted()] reports;
#' * the draws of [simulate()] and [predict()], `unwhiten()` at drawn
#'   innovations;
#' * the residual covariance divided by `sigma^2` ([autocor_matrix()])
#'   and the marginal variance that pearson residuals divide by.
#'
#' The term holds for all the rows of a response at once, so the
#' likelihood does not factor into rows. frmtmb therefore refuses the
#' term with `weights()`, `cens()`, `trunc()`, `se()`, `mi()`,
#' `rescor = TRUE`, a mixture, a non-linear formula and `quadrature`,
#' as it refuses `ar(cov = TRUE)`. The term's fixed data describe the
#' rows of the data, so [fitted()], [predict()] and [simulate()] take
#' `newdata` with one row per row of the term's matrix, read in their
#' order (a counterfactual on the same units, as brms takes it). When
#' the fit dropped rows for missing values, `newdata` covers them too:
#' the term is built again over every row, its parameters held at
#' their natural values (through `bounds`, for a bounded one), and a
#' row with a missing predictor gives `NA`. Other `newdata` is refused,
#' as brms refuses it without a matrix over the new rows in its
#' `newdata2`.
#'
#' @section The term:
#' `args` names the term's arguments in brms's positional order.
#' `parse(args, env)` receives the named list of unevaluated argument
#' expressions (`NULL` for an absent one) and the formula environment;
#' it checks what it can without data and returns the spec. `build(spec,
#' rows, n_all, lookup)` runs once per fit: `rows` are the rows of the
#' data that the model keeps (missing values and `subset()` removed),
#' `n_all` the number of rows of the data, and `lookup(expr, what)`
#' resolves an argument expression from `data2`, the data and the
#' formula environment, as `gr(cov = A)` resolves its matrix. It returns
#' `aux`, the fixed data every other function receives. A matrix over
#' the observations is subset with `rows`, as brms subsets it.
#'
#' @section The map:
#' `whiten(y, mu, sigma, theta, aux)` runs on the TAPE: `mu`, `sigma`
#' and `theta` are RTMB advectors there, `y` is numeric and `mu` has one
#' value per row. `sigma` is ONE value when `sigma_varies = FALSE`, and
#' one value per row when `sigma_varies = TRUE`. It returns
#' `list(u, logjac)`: the innovations, one per row, and the log of the
#' absolute determinant of the map's Jacobian in `y`, one number. Write
#' it without loops over rows; a sparse `Matrix` times an advector
#' works.
#'
#' `unwhiten(u, mu, sigma, theta, aux)` runs on numbers and returns the
#' response with innovations `u`, with `sigma` as for `whiten`. `u` is a
#' vector, or a matrix with one column per set of innovations, and the
#' result has the same shape: core passes all the columns of the
#' identity at once to read the residual covariance, so a term that
#' factors a matrix factors it once. It must be the inverse of `whiten`
#' at the same parameters.
#'
#' When the frame is built, frmtmb checks the map once on numbers: at
#' the start, `whiten` gives `n` innovations and one finite log-Jacobian,
#' `unwhiten` gives `n` values (and an `n x 2` matrix for two columns),
#' and `whiten(unwhiten(u))` returns `u` to a relative error of 1e-6. A
#' term that fails is refused by its name.
#'
#' @section Bounded parameters and priors:
#' Give `bounds` when each parameter is a bounded quantity mapped as
#' `lo + (hi - lo) / (1 + exp(-theta))`, as `sar()`'s rho is:
#' `function(aux)` returning a matrix with columns `lower` and `upper`,
#' one row per parameter, in the order of `natural`. Then
#' `set_prior(class = "<natural name>")` puts a density on the natural
#' parameter, carried to `theta` with the log-Jacobian of the map, and
#' `frmtmb_autocor_pars()` reports the interval, from which
#' frmtmb.sample takes brms's flat default on it.
#'
#' @section Testing a term:
#' Compare the fitted log-likelihood with a dense multivariate normal
#' density, `mvtnorm::dmvnorm()` at the mean `unwhiten(0, ...)` and the
#' covariance of the map, and check that `whiten(unwhiten(u))` is `u`.
#' A test that registers calls [frmtmb_local_registries()] first.
#'
#' @param name The term's name, as a formula writes it, without
#'   parentheses. It must be a syntactic name that nothing else in a
#'   formula means: not one of frmtmb's own terms or structures, not a
#'   brms name unless core lists it as this package's
#'   (`sar` is frmtmb.spatial's), and not a function formulas call as a
#'   term.
#' @param args The term's argument names, in positional order.
#' @param build `function(spec, rows, n_all, lookup)`; see "The term".
#' @param whiten,unwhiten The map and its inverse; see "The map".
#' @param parse `NULL`, or `function(args, env)`; see "The term". The
#'   default keeps the expressions as they are.
#' @param npar The number of entries the term takes in the parameter
#'   vector `thetaac`: a whole number, or `function(aux)`.
#' @param start Starting values: a numeric vector, or `function(aux)`.
#'   The default is zeros.
#' @param natural `function(theta, aux)` giving the term's parameters
#'   on their natural scale, named as brms names them (`lagsar`). It is
#'   required when `npar` can be positive. [summary()],
#'   [confint_varcorr()], [variables()] and [hypothesis()] read it.
#' @param sigma_varies `TRUE` when the term takes a predicted
#'   `sigma ~ ...`, so that `whiten` and `unwhiten` receive one `sigma`
#'   per row. `FALSE` refuses a predicted `sigma`, and they receive one
#'   `sigma`.
#' @param bounds `NULL`, or `function(aux)`; see "Bounded parameters and
#'   priors".
#' @param usage How the term is written, for messages:
#'   `"sar(M, type = \"lag\")"`.
#' @param compat_note The note of the term's row for `gaussian()` and
#'   `student()` in [frm_compat()], where it works: what the term fits.
#'   `NULL` gives a short default. The rows that refuse the other
#'   families, the addition terms, `rescor`, mixtures, non-linear
#'   formulas, quadrature and a second residual term are written with
#'   the registration.
#' @param replace,package As in [frmtmb_register_covstruct()]. A fit
#'   records the package, so that a fit read where it is not installed
#'   is refused by its name.
#' @return `NULL`, invisibly. Called for the registration.
#' @seealso [frmtmb-autocor] for the terms core holds, `fcor()` among
#'   them; [frmtmb_register_covstruct()] for a structure that adds a
#'   latent vector to the linear predictor instead.
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' # a residual correlation with a FIXED matrix V and no parameter,
#' # myfix(V): y = mu + sigma * L u with L the Cholesky factor of V
#' myfix_build <- function(spec, rows, n_all, lookup) {
#'   V <- lookup(spec$V, "myfix()")[rows, rows, drop = FALSE]
#'   L <- t(chol(V))
#'   list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
#'        logdet = sum(log(diag(L))))
#' }
#' # sigma_varies = FALSE: sigma is one value
#' myfix_whiten <- function(y, mu, sigma, theta, aux) {
#'   u <- as.vector(aux$Linv %*% ((y - mu) / sigma))
#'   list(u = u, logjac = -aux$logdet - length(y) * log(sigma))
#' }
#' myfix_unwhiten <- function(u, mu, sigma, theta, aux) {
#'   mu + sigma * (aux$L %*% u)
#' }
#' # a package calls frmtmb_register_autocor() from its .onLoad(); here
#' # the registration lasts while demo() runs
#' demo <- function() {
#'   frmtmb_local_registries()
#'   frmtmb_register_autocor("myfix", args = "V", build = myfix_build,
#'                           whiten = myfix_whiten,
#'                           unwhiten = myfix_unwhiten,
#'                           sigma_varies = FALSE)
#'   set.seed(1)
#'   V <- 0.5^abs(outer(1:20, 1:20, "-"))
#'   d <- data.frame(x = rnorm(20))
#'   d$y <- 1 + d$x + as.vector(t(chol(V)) %*% rnorm(20))
#'   fit <- frm(y ~ x + myfix(V), data = d, data2 = list(V = V))
#'   # the same model as core's own fcor(V)
#'   c(myfix = logLik(fit),
#'     fcor = logLik(frm(y ~ x + fcor(V), data = d, data2 = list(V = V))))
#' }
#' demo()
#' @export
frmtmb_register_autocor <- function(name, args, build, whiten, unwhiten,
                                    parse = NULL, npar = 0L, start = NULL,
                                    natural = NULL, sigma_varies = TRUE,
                                    bounds = NULL, usage = NULL,
                                    compat_note = NULL,
                                    replace = FALSE,
                                    package = compat_registrant_name(
                                      parent.frame())) {
  arg <- function(nm) paste0("frmtmb_register_autocor(", nm, " =)")
  package <- seam_check_package(package, arg("package"))
  check_flag(replace, arg("replace"))
  check_flag(sigma_varies, arg("sigma_varies"))
  seam_check_name(name, arg("name"))
  if (name %in% c(autocor_structs, names(autocor_core))) {
    frm_stop("`", name, "` is one of frmtmb's own residual correlation ",
             "terms and cannot be re-registered", call. = FALSE)
  }
  home <- covstruct_home_owner(name)
  if (!is.null(home) && !identical(package, home)) {
    frm_stop("frmtmb_register_autocor() cannot register '", name,
             "': it names a term of the extension package ", home,
             ", and only ", home, " registers it", call. = FALSE)
  }
  if (is.null(home) && name %in% autocor_reserved_names()) {
    frm_stop("frmtmb_register_autocor() cannot register '", name,
             "', because ", name, "() already means something else in a ",
             "formula. Choose another name.", call. = FALSE)
  }
  old <- frmtmb_autocor_registry$reg[[name]]
  if (!is.null(old) && !replace && !identical(old[["package"]], package)) {
    frm_stop("frmtmb_register_autocor(): '", name, "' is already ",
             "registered by ",
             if (is.null(old[["package"]])) "this session" else {
               paste0("the package ", old[["package"]])
             }, ". Only that package replaces it on its own; pass ",
             "replace = TRUE to replace it anyway", call. = FALSE)
  }
  if (!is.character(args) || anyNA(args) || !all(nzchar(args)) ||
      anyDuplicated(args)) {
    frm_stop(arg("args"), " must be the term's argument names, without ",
             "repeats", call. = FALSE)
  }
  fn_or <- function(x, nm, argn, null_ok = FALSE) {
    seam_fn(x, arg(nm), argn, null_ok)
  }
  build <- fn_or(build, "build", 4L)
  whiten <- fn_or(whiten, "whiten", 5L)
  unwhiten <- fn_or(unwhiten, "unwhiten", 5L)
  parse <- fn_or(parse, "parse", 2L, null_ok = TRUE) %||%
    function(args, env) args
  natural <- fn_or(natural, "natural", 2L, null_ok = TRUE)
  bounds <- fn_or(bounds, "bounds", 1L, null_ok = TRUE)
  npar_f <- if (is.function(npar)) npar else {
    if (!is.numeric(npar) || length(npar) != 1L || is.na(npar) ||
        npar < 0 || npar != trunc(npar)) {
      frm_stop(arg("npar"), " must be one whole number or a function of ",
               "aux", call. = FALSE)
    }
    np <- as.integer(npar)
    function(aux) np
  }
  if (is.null(natural) && !(is.numeric(npar) && npar == 0)) {
    frm_stop(arg("natural"), " is needed for a term with parameters: ",
             "summary(), confint_varcorr() and variables() report them ",
             "through it", call. = FALSE)
  }
  start_f <- if (is.null(start)) {
    function(aux) numeric(npar_f(aux))
  } else if (is.function(start)) {
    start
  } else {
    if (!is.numeric(start) || !all(is.finite(start))) {
      frm_stop(arg("start"), " must be finite numbers or a function of ",
               "aux", call. = FALSE)
    }
    st <- as.numeric(start)
    function(aux) st
  }
  usage <- usage %||% paste0(name, "(", paste(args, collapse = ", "), ")")
  if (!is.character(usage) || length(usage) != 1L) {
    frm_stop(arg("usage"), " must be one string, how the term is written",
             call. = FALSE)
  }
  if (!is.null(compat_note) && (!is.character(compat_note) ||
                                  length(compat_note) != 1L ||
                                  is.na(compat_note))) {
    frm_stop(arg("compat_note"), " must be one string", call. = FALSE)
  }
  entry <- list(args = args, parse = parse, build = build, whiten = whiten,
                unwhiten = unwhiten, npar = npar_f, start = start_f,
                natural = natural, sigma_varies = sigma_varies,
                bounds = bounds, usage = usage, package = package)
  entry <- entry[!vapply(entry, is.null, TRUE)]
  # the term joins the compatibility vocabulary as it joins the parser,
  # so frm_compat() can describe every term frm() accepts; the last step
  # that can refuse, so a refusal leaves both registries as they were
  feat <- paste0(name, "()")
  known <- frmtmb_compat_features_tbl()
  add <- compat_new_features(stats::setNames("autocor", feat),
                             "frmtmb_register_autocor()", known)
  if (length(add)) {
    before <- frmtmb_compat_contrib$features
    frmtmb_compat_contrib$features <- c(before, list(add))
    compat_cache_extend(add, before)
  }
  rules <- autocor_seam_compat_rules(feat, compat_note)
  attr(rules, "frmtmb_package") <- package %||% ""
  frmtmb_compat_contrib$rules[[paste0("autocor:", name)]] <- rules
  frmtmb_autocor_registry$reg[[name]] <- entry
  invisible(NULL)
}

#' The frm_compat() rows frmtmb writes for a registered residual
#' correlation term: the refusals core applies to every such term, with
#' the notes the `group:autocor` rows give core's own terms, and the
#' families it takes with the note its package gives (`compat_note`).
#' One row for the two families, so a package that wrote its own would
#' tie with it and win no pair.
#'
#' @noRd
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
    r(feat, "group:autocor_families", "works", note)
    r(feat, "kind:aterm", "refused",
      "Refused: the density is joint over the rows, so there is no per-row contribution for a weight, a censoring indicator, a truncation bound or a known standard error to act on.")
    r(feat, "rescor", "refused",
      "Refused: rescor = TRUE is another residual covariance, across responses.")
    r(feat, "mixture", "refused",
      "Refused: a mixture likelihood has no single residual to correlate.")
    r(feat, "nl", "refused",
      "Refused: the term is written in the formula of mu.")
    r(feat, "quadrature", "refused",
      "Refused: quadrature integrates a random effect against per-row densities, and this density is joint over the rows.")
    r(feat, "kind:autocor", "refused",
      "Refused: a response has one residual covariance, so it carries one such term.")
    r(feat, "residuals_osa", "refused",
      "Refused: the tape holds one joint density, not a density per row.")
    b$rules()
  }
}

#' A registered residual correlation term's entry, or `NULL` when the
#' name is not one.
#'
#' @noRd
ac_entry <- function(name) frmtmb_autocor_registry$reg[[name]]

#' The names of every registered residual correlation term, core's
#' `fcor` first.
#'
#' @noRd
ac_entry_names <- function() names(frmtmb_autocor_registry$reg)

#' Every residual correlation term this session resolves: the ar()
#' family of R/autocor.R and the registry's terms.
#'
#' @noRd
autocor_term_names <- function() c(autocor_structs, ac_entry_names())

#' Whether a residual block is a registry term rather than one of the
#' ar() family that R/autocor.R implements directly.
#'
#' @noRd
autocor_is_seam <- function(ac) !is.null(ac) && isTRUE(ac[["seam"]])

#' The parts of a fit's residual correlation terms that this session
#' cannot supply: one `list(kind, pkg, what)` per registry term with no
#' entry, `kind` "session" (registered in an R session, not by a
#' package), "missing" (its package is not installed or does not load)
#' or "stale" (the package loads and does not register it). Loading the
#' package is tried first, which registers the term. Kept apart from the
#' covariance structures' parts so that one refusal lists both
#' (`seam_missing_parts()`, R/extension-home.R).
#'
#' @noRd
autocor_missing_parts <- function(x) {
  frame <- if (!is.null(x[["frame"]])) x[["frame"]] else x
  parts <- list()
  for (ac in frame[["autocor"]] %||% list()) {
    if (!autocor_is_seam(ac) || !is.null(ac_entry(ac[["struct"]]))) next
    nm <- ac[["struct"]]
    pkg <- ac[["package"]] %||% extension_special_home(nm)
    usable <- !is.null(pkg) && (isNamespaceLoaded(pkg) ||
                                  requireNamespace(pkg, quietly = TRUE))
    if (usable && !is.null(ac_entry(nm))) next
    parts[[length(parts) + 1L]] <- list(
      kind = if (is.null(pkg)) "session" else if (!usable) "missing" else
        "stale",
      pkg = pkg, what = paste0("The ", nm, "() term of this model"),
      register = "frmtmb_register_autocor()")
  }
  parts
}

#' A block's entry, refused by its package's name when that package is
#' not loaded: a saved fit with a registered term, read without the
#' package. The refusal is the one every missing part of a model gets
#' (`stop_extension_parts()`, R/extension-home.R).
#'
#' @noRd
ac_block_entry <- function(ac) {
  e <- ac_entry(ac[["struct"]])
  if (!is.null(e)) return(e)
  parts <- autocor_missing_parts(list(autocor = list(ac)))
  if (length(parts)) stop_extension_parts(parts)
  ac_entry(ac[["struct"]])
}

#' Parse a registered term's call into the spec entry `pull_autocor()`
#' carries to the frame.
#'
#' @noRd
parse_autocor_seam_call <- function(tm, env) {
  fn <- as.character(tm[[1L]])[1L]
  e <- ac_entry(fn)
  a <- match_special_args(tm, e[["args"]], fn)
  spec <- e[["parse"]](a, env)
  list(fn = fn, struct = fn, seam = TRUE, spec = spec, cov = TRUE,
       time_expr = NULL, gr_expr = NULL, gr_vars = character(0),
       p = 0L, q = 0L, label = deparse1(tm))
}

#' The rows of the data a registered term's block covers: every row,
#' less the ones `na.action` dropped and the ones a univariate
#' `subset()` leaves out, as brms subsets a matrix over the
#' observations. A multivariate `subset()` gives each response its own
#' rows, which a term over one response's rows could read, but nothing
#' tests it, so it is refused.
#'
#' @noRd
autocor_seam_rows <- function(ac, n_all, mf, uni_rows, sub_rows, n) {
  if (length(sub_rows)) {
    frm_stop(ac[["label"]], " cannot be combined with subset() in a ",
             "multivariate model", call. = FALSE)
  }
  rows <- seq_len(n_all)
  na <- attr(mf, "na.action")
  if (length(na)) rows <- rows[-as.integer(na)]
  if (!is.null(uni_rows)) rows <- rows[uni_rows]
  if (length(rows) != n) {
    frm_stop("internal: ", ac[["label"]], " covers ", length(rows),
             " rows where the model has ", n, call. = FALSE)
  }
  rows
}

#' The sigma a registered term's map receives: one value under
#' `sigma_varies = FALSE`, as its documentation says, and one per row
#' otherwise. Advector-safe.
#'
#' @noRd
autocor_seam_sigma <- function(e, sigma, n) {
  if (isTRUE(e[["sigma_varies"]])) {
    if (length(sigma) == 1L) sigma + numeric(n) else sigma
  } else {
    sigma[1L]
  }
}

#' Check a registered term's map once, on numbers, when its block is
#' built: the lengths `whiten` and `unwhiten` give, a finite log-Jacobian
#' and a round trip `whiten(unwhiten(u)) = u` at the start, `mu = 0` and
#' `sigma = 1`. A slip there would otherwise fit a different model
#' without a word (one innovation short), or report a covariance the
#' term never meant (an `unwhiten` that is not the inverse). The
#' tolerance is relative to the size of `u`.
#'
#' @noRd
autocor_seam_check_map <- function(e, ac, aux, st, n) {
  bad <- function(what) {
    frm_stop("The '", ac[["struct"]], "' term's map ", what, " for ",
             ac[["label"]], call. = FALSE)
  }
  u <- sin(seq_len(n))
  sg <- autocor_seam_sigma(e, 1, n)
  y <- tryCatch(e[["unwhiten"]](u, numeric(n), sg, st, aux),
                error = function(err) bad(paste0("could not be evaluated: ",
                                                 conditionMessage(err))))
  if (!is.numeric(y) || length(y) != n || any(!is.finite(as.numeric(y)))) {
    bad(paste0("is wrong: unwhiten() gave ", length(y), " value(s) for ",
               n, " innovations, or a value that is not finite"))
  }
  U <- cbind(u, rev(u))
  Y <- tryCatch(e[["unwhiten"]](U, numeric(n), sg, st, aux),
                error = function(err) NULL)
  # each column must be the vector call on that column alone
  y2 <- tryCatch(as.numeric(e[["unwhiten"]](rev(u), numeric(n), sg, st, aux)),
                 error = function(err) NULL)
  Ym <- if (is.null(Y)) NULL else as.matrix(Y)
  if (is.null(Ym) || is.null(y2) || !identical(dim(Ym), c(n, 2L)) ||
        any(abs(Ym - cbind(as.numeric(y), y2)) >
              1e-8 * max(1, abs(as.numeric(y)), abs(y2)))) {
    bad(paste0("is wrong: unwhiten() of a matrix of innovations must ",
               "give the matrix of responses, one column each"))
  }
  w <- tryCatch(e[["whiten"]](as.numeric(y), numeric(n), sg, st, aux),
                error = function(err) bad(paste0("could not be evaluated: ",
                                                 conditionMessage(err))))
  uu <- as.numeric(w[["u"]])
  lj <- w[["logjac"]]
  if (length(uu) != n) {
    bad(paste0("is wrong: whiten() gave ", length(uu), " innovation(s) ",
               "for ", n, " rows"))
  }
  if (!is.numeric(lj) || length(lj) != 1L || !is.finite(lj)) {
    bad("is wrong: whiten() must give one finite log-Jacobian")
  }
  if (max(abs(uu - u)) > 1e-6 * max(abs(u))) {
    bad(paste0("is wrong: whiten(unwhiten(u)) is not u (largest ",
               "difference ", format(max(abs(uu - u)), digits = 3),
               "), so unwhiten() is not the inverse of whiten()"))
  }
  b <- e[["bounds"]]
  if (!is.null(b)) {
    B <- b(aux)
    nat <- e[["natural"]](st, aux)
    ok <- is.matrix(B) && is.numeric(B) && nrow(B) == length(st) &&
      ncol(B) == 2L && all(is.finite(B)) && all(B[, 1L] < B[, 2L])
    if (!ok) {
      bad("is wrong: bounds() must give one row (lower, upper) per parameter")
    }
    want <- B[, 1L] + (B[, 2L] - B[, 1L]) / (1 + exp(-st))
    if (any(abs(as.numeric(nat) - want) >
              1e-8 * pmax(1, abs(B[, 2L] - B[, 1L])))) {
      bad(paste0("is wrong: natural() is not the bounded map of bounds() ",
                 "at the start"))
    }
  }
  invisible(TRUE)
}

#' Build a registered term's block once, when the frame is assembled:
#' the term's fixed data from `build`, the checks of its slots, and the
#' bookkeeping the other residual blocks carry (one group holding every
#' row, in row order).
#'
#' @noRd
autocor_seam_block <- function(ac, resp, rows, n_all, lookup) {
  e <- ac_entry(ac[["struct"]])
  n <- length(rows)
  aux <- e[["build"]](ac[["spec"]], rows, n_all, lookup)
  np <- e[["npar"]](aux)
  if (!is.numeric(np) || length(np) != 1L || is.na(np) || np < 0 ||
      np != trunc(np)) {
    frm_stop("The '", ac[["struct"]], "' term's npar slot gave ",
             paste(format(np), collapse = ", "), " for ", ac[["label"]],
             "; it must be one whole number, zero or more", call. = FALSE)
  }
  st <- e[["start"]](aux)
  if (!is.numeric(st) || length(st) != np || any(!is.finite(st))) {
    frm_stop("The '", ac[["struct"]], "' term's start slot gave ",
             length(st), " value(s) where npar is ", np, ", or a value ",
             "that is not finite, for ", ac[["label"]], call. = FALSE)
  }
  autocor_seam_check_map(e, ac, aux, st, n)
  ac[["aux"]] <- aux
  ac[["d"]] <- n
  # the rows of the data the term's matrix covers, the dropped ones
  # among them: brms reads newdata of this many rows against it
  ac[["n_all"]] <- n_all
  ac[["time_levels"]] <- as.character(seq_len(n))
  ac[["group_levels"]] <- "1"
  ac[["patterns"]] <- list()
  ac[["npar"]] <- as.integer(np)
  ac[["resp"]] <- resp
  ac[["n_groups"]] <- 1L
  ac[["package"]] <- e[["package"]]
  ac
}

#' The log density of a response under a registered term, on the tape:
#' the gaussian or multivariate Student-t of the term's innovations,
#' plus the log-Jacobian of its map.
#'
#' @noRd
autocor_seam_loglik <- function(e, y, mu, sigma, theta, ac, nu = NULL) {
  n <- length(y)
  w <- e[["whiten"]](y, mu, autocor_seam_sigma(e, sigma, n), theta,
                     ac[["aux"]])
  u <- w[["u"]]
  q <- sum(u * u)
  if (is.null(nu)) {
    return(-0.5 * n * log(2 * pi) - 0.5 * q + w[["logjac"]])
  }
  # brms's multi_student_t(nu, mu, Sigma) at innovations of unit scale;
  # lgamma_shift_diff() keeps the head accurate for a large nu, as the
  # residual correlation of the ar() family does (mvt_std_loglik())
  lgamma_shift_diff(nu / 2, n / 2) - 0.5 * n * (log(nu) + log(pi)) -
    0.5 * (nu + n) * log1p(q / nu) + w[["logjac"]]
}

#' The natural names of a registered block's parameters, read off its
#' `natural` slot at the start.
#'
#' @noRd
autocor_seam_par_names <- function(ac) {
  e <- ac_block_entry(ac)
  if (!ac[["npar"]] || is.null(e[["natural"]])) return(character(0))
  names(e[["natural"]](e[["start"]](ac[["aux"]]), ac[["aux"]])) %||%
    character(0)
}

#' A registered block's BOUNDED parameters, the ones its `bounds` slot
#' maps by `lo + (hi - lo) plogis(theta)`: one `list(name, idx, lower,
#' upper)` each, `idx` within the block's theta. Empty for a term with
#' no `bounds`.
#'
#' @noRd
autocor_seam_bounded <- function(ac) {
  if (!autocor_is_seam(ac)) return(list())
  e <- ac_block_entry(ac)
  if (is.null(e[["bounds"]])) return(list())
  B <- e[["bounds"]](ac[["aux"]])
  nm <- autocor_seam_par_names(ac)
  lapply(seq_along(nm), function(k) {
    list(name = nm[[k]], idx = k, lower = B[k, 1L], upper = B[k, 2L])
  })
}

#' The internal names (`thetaac_<k>`) of a frame's bounded registry
#' parameters, named by their class (`lagsar`): such a parameter has
#' one prior slot, its class, so the internal spelling is not listed
#' in the prior table and a prior on it is refused.
#'
#' @noRd
autocor_seam_bounded_coefs <- function(frame) {
  out <- character(0)
  for (ac in frame[["autocor"]] %||% list()) {
    for (pr in autocor_seam_bounded(ac)) {
      out[paste0("thetaac_", ac[["theta_idx"]][pr$idx])] <- pr$name
    }
  }
  out
}

#' The class names `set_prior()` reaches on a frame's registry terms:
#' `bounded = TRUE` the bounded parameters, `FALSE` the others.
#'
#' @noRd
autocor_seam_classes <- function(frame, bounded = TRUE) {
  out <- character(0)
  for (ac in frame[["autocor"]] %||% list()) {
    if (!autocor_is_seam(ac)) next
    b <- vapply(autocor_seam_bounded(ac), `[[`, "", "name")
    out <- c(out, if (bounded) b else setdiff(autocor_seam_par_names(ac), b))
  }
  unique(out)
}

#' The parameters of a fit's residual correlation terms
#'
#' Lists the parameters of the residual correlation terms a fit takes
#' from the registry of [frmtmb_register_autocor()] (frmtmb.spatial's
#' `sar()`): their names as brms names them, the internal parameter each
#' is, and the interval a bounded one lives on. frmtmb.sample reads it
#' for the names of its draws and for brms's flat default prior on a
#' bounded parameter, so that it keeps no list of names of its own.
#'
#' @param fit A `frmtmb_fit`.
#' @return A data frame with one row per parameter: `resp`, `term` (the
#'   term as written), `name` (`lagsar`), `coef` (the internal
#'   parameter, `thetaac_1`), `lower` and `upper` (`NA` when the term
#'   declares no bounds). No rows when the fit has no such term.
#' @seealso [frmtmb_register_autocor()].
#' @examples
#' set.seed(1)
#' V <- 0.5^abs(outer(1:20, 1:20, "-"))
#' d <- data.frame(x = rnorm(20), y = rnorm(20))
#' fit <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
#' frmtmb_autocor_pars(fit)   # fcor() has no parameter
#' @export
frmtmb_autocor_pars <- function(fit) {
  out <- data.frame(resp = character(0), term = character(0),
                    name = character(0), coef = character(0),
                    lower = numeric(0), upper = numeric(0),
                    stringsAsFactors = FALSE)
  if (length(fit$frame[["autocor"]])) require_fit_packages(fit)
  for (r in names(fit$frame[["autocor"]] %||% list())) {
    ac <- fit$frame[["autocor"]][[r]]
    if (!autocor_is_seam(ac) || !ac[["npar"]]) next
    nm <- autocor_seam_par_names(ac)
    b <- autocor_seam_bounded(ac)
    lo <- rep(NA_real_, length(nm))
    hi <- lo
    for (pr in b) {
      lo[pr$idx] <- pr$lower
      hi[pr$idx] <- pr$upper
    }
    out <- rbind(out, data.frame(
      resp = r, term = ac[["label"]], name = nm,
      coef = paste0("thetaac_", ac[["theta_idx"]]), lower = lo,
      upper = hi, stringsAsFactors = FALSE))
  }
  out
}

#' The parameters of one registered block at a fit's estimates.
#'
#' @noRd
autocor_seam_theta <- function(fit, ac) {
  fit$estimates[["thetaac"]][ac[["theta_idx"]]]
}

#' The mean of the response under a registered term, `unwhiten()` at
#' zero innovations. For brms's lag SAR this is `(I - rho W)^-1 mu`,
#' brms's `posterior_epred()`; for fcor() and the error SAR it is `mu`.
#'
#' @noRd
autocor_seam_mean <- function(ac, theta, mu, sigma) {
  theta <- autocor_seam_theta_on(ac, theta)
  e <- ac_block_entry(ac)
  n <- ac[["d"]]
  as.numeric(e[["unwhiten"]](numeric(n), rep_len(as.numeric(mu), n),
                             autocor_seam_sigma(e, as.numeric(sigma), n),
                             theta, ac[["aux"]]))
}

#' The linear map `B` from innovations to response, `y = m + B u`: one
#' call of `unwhiten()` on the identity, so a term that factors a matrix
#' factors it once. `B B'` is the covariance of the response for
#' unit-scale innovations; `sigma = 1` gives it over `sigma^2`.
#'
#' @noRd
autocor_seam_B <- function(ac, theta, mu, sigma) {
  theta <- autocor_seam_theta_on(ac, theta)
  e <- ac_block_entry(ac)
  n <- ac[["d"]]
  mu <- rep_len(as.numeric(mu), n)
  sigma <- autocor_seam_sigma(e, as.numeric(sigma), n)
  m <- as.numeric(e[["unwhiten"]](numeric(n), mu, sigma, theta,
                                  ac[["aux"]]))
  as.matrix(e[["unwhiten"]](diag(n), mu, sigma, theta, ac[["aux"]])) - m
}

#' The marginal variance of each row over `sigma^2`, `diag(B B')`,
#' without forming `B B'`.
#'
#' @noRd
autocor_seam_marginal <- function(ac, theta) {
  B <- autocor_seam_B(ac, theta, numeric(ac[["d"]]), 1)
  rowSums(B * B)
}

#' One draw of the response under a registered term: gaussian
#' innovations, or one chi-square scaling of them for the multivariate
#' Student-t, through `unwhiten()`.
#'
#' @noRd
autocor_seam_draw <- function(ac, theta, mu, sigma, nu = NULL) {
  theta <- autocor_seam_theta_on(ac, theta)
  e <- ac_block_entry(ac)
  n <- ac[["d"]]
  u <- stats::rnorm(n)
  if (!is.null(nu)) u <- u * sqrt(nu / stats::rchisq(1L, nu))
  as.numeric(e[["unwhiten"]](u, rep_len(as.numeric(mu), n),
                             autocor_seam_sigma(e, as.numeric(sigma), n),
                             theta, ac[["aux"]]))
}

#' The refusal of `newdata` whose rows are not the rows of the term's
#' matrix. brms 2.23.0 reads newdata against the matrix of `data2`, so
#' it takes exactly `nrow(M)` rows, the rows the fit dropped for missing
#' values among them, and refuses any other number unless `newdata2`
#' holds a matrix over the new rows (`dev/sarfcor-fr2-rows.R`).
#'
#' @noRd
stop_autocor_seam_newdata <- function(ac, m, what) {
  n_all <- ac[["n_all"]] %||% ac[["d"]]
  drop <- n_all - ac[["d"]]
  frm_stop(what, " with ", ac[["label"]], " needs newdata with one ",
           "row per row of the term's matrix: ", n_all, " rows, read as ",
           "the rows of the data in their order (a counterfactual on the ",
           "same units), where newdata has ", m, ".",
           if (drop > 0L) {
             paste0(" The fit dropped ", drop, " of those rows for ",
                    "missing values; newdata covers them too, as brms ",
                    "reads it, and a row with a missing predictor gives NA.")
           },
           " brms refuses any other number of rows unless newdata2 holds ",
           "the term's matrix over the new rows, which frmtmb does not ",
           "take. frm_linpred() gives the linear predictor mu on any ",
           "newdata, without the term", call. = FALSE)
}

#' Whether a block covers `newdata`'s rows: a block's own rows, for a
#' draw on the rows it was built for.
#'
#' @noRd
check_autocor_seam_newdata <- function(ac, newdata, what) {
  if (is.null(newdata) || NROW(newdata) == ac[["d"]]) {
    return(invisible(TRUE))
  }
  stop_autocor_seam_newdata(ac, NROW(newdata), what)
}

#' A registered block for `newdata`, read as brms reads it: one row per
#' row of the term's matrix, in order. When the fit dropped no row the
#' fitted block serves. When it dropped rows for missing values, the
#' block is built again over every row of the matrix, from the fit's
#' `data2`, and its parameters are carried over at their natural
#' values: brms's `lagsar` is rho itself, while frmtmb's internal
#' parameter maps onto an interval that depends on the rows
#' (`autocor_seam_theta_map()`). Any other number of rows is refused.
#'
#' @noRd
autocor_seam_for_newdata <- function(fit, ac, rspec, newdata, what) {
  if (is.null(newdata)) return(ac)
  n_all <- ac[["n_all"]] %||% ac[["d"]]
  m <- NROW(newdata)
  if (m != n_all) stop_autocor_seam_newdata(ac, m, what)
  if (n_all == ac[["d"]]) return(ac)
  e <- ac_block_entry(ac)
  env <- rspec$formula_env
  data2 <- fit$data2 %||% list()
  data <- fit$data
  aux <- e[["build"]](ac[["spec"]], seq_len(n_all), n_all,
                      function(expr, what) {
                        lookup_structural(expr, data2, data, env, what)
                      })
  acn <- ac
  acn[["aux"]] <- aux
  acn[["d"]] <- n_all
  acn[["time_levels"]] <- as.character(seq_len(n_all))
  acn[["theta_map"]] <- autocor_seam_theta_map(e, ac, ac[["aux"]], aux)
  acn
}

#' The map from a block's fitted internal parameters to those of the
#' same term built over other rows, holding the natural values fixed.
#' A term with `bounds` is mapped through them; a value outside the new
#' interval is refused (for sar(), `I - rho M` is then singular or has
#' a negative determinant over the new rows). A term with parameters
#' and no `bounds` is carried over only when its `natural` slot gives
#' the same values on both sets of rows.
#'
#' @noRd
autocor_seam_theta_map <- function(e, ac, aux0, aux1) {
  if (!ac[["npar"]]) return(NULL)
  b <- e[["bounds"]]
  if (!is.null(b)) {
    B0 <- b(aux0)
    B1 <- b(aux1)
    return(function(theta) {
      nat <- B0[, 1L] + (B0[, 2L] - B0[, 1L]) / (1 + exp(-theta))
      if (any(nat <= B1[, 1L] | nat >= B1[, 2L])) {
        frm_stop("A prediction with ", ac[["label"]], " on every row ",
                 "of its matrix: the fitted value ",
                 paste(format(nat, digits = 4), collapse = ", "),
                 " is outside the interval the term has over those ",
                 "rows, (", paste(format(B1[, 1L], digits = 4),
                                    format(B1[, 2L], digits = 4),
                                    sep = ", ", collapse = "; "),
                 ")", call. = FALSE)
      }
      stats::qlogis((nat - B1[, 1L]) / (B1[, 2L] - B1[, 1L]))
    })
  }
  function(theta) {
    n0 <- e[["natural"]](theta, aux0)
    if (!isTRUE(all.equal(n0, e[["natural"]](theta, aux1)))) {
      frm_stop("A prediction with ", ac[["label"]], " on every row of ",
               "its matrix needs the term's parameters on those rows, ",
               "and the term declares no bounds to carry them over by",
               call. = FALSE)
    }
    theta
  }
}

#' A block's internal parameters as its own rows read them: the fitted
#' ones, or their image under the block's `theta_map` for a block built
#' over other rows.
#'
#' @noRd
autocor_seam_theta_on <- function(ac, theta) {
  f <- ac[["theta_map"]]
  if (is.null(f)) theta else f(theta)
}

# ------------------------------------------------------------- fcor()
#
# brms 2.23.0's fcor(M): the residual covariance is sigma^2 M for a
# constant sigma, and D L L' D with D = diag(sigma) and L the Cholesky
# factor of M for a predicted sigma (its normal_fcor_het_lpdf). M is a
# known matrix over the observations: a correlation, or a covariance
# whose scale sigma then multiplies. A vector is the diagonal.

#' @noRd
fcor_build <- function(spec, rows, n_all, lookup) {
  M <- lookup(spec[["M"]], "fcor()")
  if (inherits(M, "Matrix")) M <- as.matrix(M)
  if (is.null(dim(M)) && is.numeric(M)) M <- diag(as.numeric(M), length(M))
  if (!is.matrix(M) || !is.numeric(M) || nrow(M) != ncol(M)) {
    frm_stop("'M' for fcor() must be a square numeric matrix, or a ",
             "vector of variances, as in brms", call. = FALSE)
  }
  if (anyNA(M) || any(!is.finite(M))) {
    frm_stop("'M' for fcor() has missing or infinite values", call. = FALSE)
  }
  if (nrow(M) != n_all) {
    frm_stop("Dimensions of 'M' for fcor() must be equal to the number ",
             "of observations: M is ", nrow(M), " x ", ncol(M), " and the ",
             "data have ", n_all, " rows. As in brms, M covers every row ",
             "of the data, and the rows the model drops are dropped from ",
             "M", call. = FALSE)
  }
  M <- M[rows, rows, drop = FALSE]
  dimnames(M) <- NULL
  if (!isSymmetric(M, tol = 100 * .Machine$double.eps)) {
    frm_stop("'M' for fcor() must be symmetric", call. = FALSE)
  }
  # brms 2.23.0 tests min(eigen(M)$values <= 0), the smallest of a
  # logical vector, so it refuses only a matrix whose eigenvalues are ALL
  # non-positive, and an indefinite M reaches its Stan program
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= length(ev) * .Machine$double.eps * max(abs(ev))) {
    frm_stop("'M' for fcor() must be positive definite; its smallest ",
             "eigenvalue is ", format(min(ev), digits = 4),
             call. = FALSE)
  }
  L <- t(chol(M))
  list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
       logdet = sum(log(diag(L))))
}

#' @noRd
fcor_whiten <- function(y, mu, sigma, theta, aux) {
  n <- length(y)
  u <- as.vector(aux[["Linv"]] %*% ((y - mu) / sigma))
  lsig <- if (length(sigma) == 1L) n * log(sigma) else sum(log(sigma))
  list(u = u, logjac = -aux[["logdet"]] - lsig)
}

#' @noRd
fcor_unwhiten <- function(u, mu, sigma, theta, aux) {
  mu + sigma * (aux[["L"]] %*% u)
}

# The built-in terms: the registry's first state, entries of the shape
# frmtmb_register_autocor() builds.
autocor_core <- list(
  fcor = list(
    args = "M",
    parse = function(args, env) {
      if (is.null(args[["M"]])) {
        frm_stop("Argument 'M' is missing in fcor().", call. = FALSE)
      }
      args
    },
    build = fcor_build, whiten = fcor_whiten, unwhiten = fcor_unwhiten,
    npar = function(aux) 0L, start = function(aux) numeric(0),
    sigma_varies = TRUE, usage = "fcor(M)")
)

frmtmb_autocor_registry <- new.env(parent = emptyenv())
frmtmb_autocor_registry$reg <- autocor_core
