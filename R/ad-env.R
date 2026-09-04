# Tape-safe scope for user-written functions and formula bodies.

#' Make a user-written function tape-safe
#'
#' Returns `f` with RTMB's tape-safe replacements for `c()`, `[<-` and
#' `diag<-` bound at the top of its body. Call it on any function a user
#' supplies that frmtmb will later put on an AD tape, so that the user
#' does not have to carry the binding boilerplate. frmtmb applies it to
#' every function slot of [frmtmb_family()] and [frmtmb_structure()]; a
#' package that takes tape-side functions of its own applies it to
#' those. It is safe to call on any object: a non-function and a
#' primitive come back unchanged.
#'
#' RTMB's replacements are LEXICAL. frmtmb's own objective code binds
#' them where it needs them (`"c" <- RTMB::ADoverload("c")`), but a
#' function the USER wrote inherits nothing from those bindings, and
#' base `c()` silently strips the advector class when a plain numeric
#' comes first. The result is a wrong gradient with no error.
#'
#' A body rewrite, not an environment injection, on purpose: the
#' byte-code compiler may inline base builtins unless the shadowing is
#' visible in the body it compiles, so an environment binding is not
#' reliably honored while these assignments always are. `ADoverload()`
#' resolves at CALL time - base versions off the tape, AD versions on
#' it - so the wrapped function stays correct in both worlds, and the
#' cost is three lookups per call of a function that is taped once.
#'
#' A function that already names `ADoverload` manages its own scope and
#' is returned untouched, so a hand-placed binding is respected and a
#' double wrap never stacks. Helpers the wrapped function CALLS still
#' use their own environments - lexical scope cannot be injected
#' transitively - so a helper that builds advectors needs the explicit
#' spelling of its own.
#'
#' @param f A function, or any object. Only a non-primitive function is
#'   changed.
#' @return `f`, with the three bindings prepended to its body when it is
#'   a function that needs them, and unchanged otherwise.
#' @seealso [frmtmb_family()] and [frmtmb_structure()], whose function
#'   slots frmtmb wraps with this
#' @examples
#' dyn <- function(t, y, p) list(c(-p[1] * y[1], p[1] * y[1]))
#' body(frmtmb_ad_overload(dyn))[[2]]
#' @export
frmtmb_ad_overload <- function(f) {
  if (!is.function(f) || is.primitive(f)) return(f)
  if ("ADoverload" %in% all.names(body(f))) return(f)
  body(f) <- bquote({
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")
    .(body(f))
  })
  f
}

#' The RTMB exports that a nonlinear body sees before base and stats.
#'
#' A nonlinear formula body is a DSL for the tape, not ordinary R, so a
#' bare `pnorm()` in one has to mean RTMB's `pnorm()`. Without that the
#' name resolves lexically to `stats::pnorm()`, which either fails on an
#' advector or returns a number that the tape then differentiates as a
#' constant, and the only cure was to write `RTMB::pnorm()` by hand.
#'
#' The set is mechanical: the RTMB exports whose names collide with
#' `base` or `stats`, minus the ones that are not numerically
#' transparent. Transparency is what makes the shadowing safe on the
#' numeric paths - `simulate()`, `predict()` and the plain-numeric
#' re-evaluation of the objective all run the same body off the tape -
#' because a shadowed name that returned a different number off the tape
#' would silently split the two paths. `test-nl-rtmb-scope.R` re-derives
#' the collision set from the installed RTMB and probes every member
#' against its `base`/`stats` original, so a new RTMB export is caught
#' rather than assumed.
#'
#' `qchisq` is the one collision left out. RTMB reaches it through its
#' own gamma quantile and lands about 1e-14 relative away from
#' `stats::qchisq()`, which makes it a reimplementation rather than the
#' same function, so the taped and numeric paths would not agree bit for
#' bit. `stats::qchisq()` keeps working in a body, and `RTMB::qchisq()`
#' is there for anyone who wants the AD-capable one.
#'
#' @noRd
nl_rtmb_shadow <- c(
  "apply", "besselI", "besselJ", "besselK", "besselY", "colSums",
  "cov2cor", "dbeta", "dbinom", "dcauchy", "dchisq", "dexp", "df",
  "dgamma", "diag", "dlnorm", "dlogis", "dmultinom", "dnbinom",
  "dnorm", "dpois", "dt", "dweibull", "eigen", "fft", "findInterval",
  "ifelse", "integrate", "lbeta", "matrix", "order", "pbeta", "pbinom",
  "pchisq", "pexp", "pgamma", "plogis", "pnbinom", "pnorm", "ppois",
  "pweibull", "qbeta", "qexp", "qgamma", "qlogis", "qnorm", "qweibull",
  "rowSums", "sapply", "solve", "sort", "splinefun", "svd", "uniroot",
  "Vectorize"
)

# resolved once per session: getExportedValue() is not free and the
# enclosure is rebuilt on every numeric re-evaluation of a body
nl_shadow_cache <- new.env(parent = emptyenv())

nl_shadow_fun <- function(nm) {
  f <- nl_shadow_cache[[nm]]
  if (is.null(f)) {
    # a name a later RTMB drops is skipped, not fatal: the body then
    # means what it meant before the shadowing existed
    f <- tryCatch(getExportedValue("RTMB", nm), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    nl_shadow_cache[[nm]] <- f
  }
  f
}

#' Which shadowed names does this body CALL?
#'
#' Only a name the body calls is shadowed, never a name it reads. That
#' difference is what keeps the shadowing clear of the other lexical
#' rule of a nonlinear body: a name used as a value is a request for a
#' data column, or an object of the formula environment that
#' `drop_nl_lexical_datavars()` left to resolve there, and plenty of
#' those objects are called `df`, `matrix` or `order`. Binding every
#' shadowed name unconditionally would turn a user's `df` data frame
#' into `RTMB::df` the function and break a body that works today.
#'
#' `all.vars()` reports the symbols used as values and `all.names()`
#' every symbol, so the difference is the names that appear only in call
#' position. A name used both ways in one body - `df(x, 4, 9) + df$k` -
#' falls out of that difference and is left alone, which is the
#' conservative half of an ambiguity no formula should contain.
#'
#' @noRd
nl_shadow_for <- function(expr) {
  if (is.null(expr)) return(character(0))
  intersect(nl_rtmb_shadow, setdiff(all.names(expr), all.vars(expr)))
}

#' An evaluation enclosure carrying the AD overloads and RTMB's math.
#'
#' For `eval(expr, values, enclos)` sites - the nonlinear formula
#' bodies - where the expression is interpreted, not byte-compiled, so
#' active bindings are honored: each access to `c`, `[<-` and `diag<-`
#' resolves through `RTMB::ADoverload()` at call time.
#'
#' The enclosure also carries RTMB's replacements for the base and stats
#' functions the body calls - see `nl_rtmb_shadow` - and it sits BELOW
#' the formula environment in the chain, so those win over the search
#' path and over a same-named function the user defined. Every other
#' name still reaches the formula environment unchanged, so a helper the
#' user wrote is found exactly as before.
#'
#' @param parent The formula environment of the body.
#' @param expr The body itself, which decides the shadowed names it
#'   calls. `NULL` shadows nothing.
#' @noRd
ad_overload_env <- function(parent, expr = NULL) {
  e <- new.env(parent = parent %||% baseenv())
  for (nm in c("c", "[<-", "diag<-")) {
    local({
      n <- nm
      makeActiveBinding(n, function() RTMB::ADoverload(n), e)
    })
  }
  for (nm in nl_shadow_for(expr)) {
    f <- nl_shadow_fun(nm)
    if (!is.null(f)) assign(nm, f, envir = e)
  }
  e
}
