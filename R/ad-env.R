# Tape-safe scope for user-written functions and formula bodies.

#' Prepend RTMB's AD overloads to a user function's body.
#'
#' RTMB's tape-safe replacements for `c()`, `[<-` and `diag<-` are
#' LEXICAL: frmtmb's own objective code binds them where it needs them
#' (`"c" <- RTMB::ADoverload("c")`), but a function the USER wrote - an
#' ODE `dynamics`, a `custom_family()` lpdf - inherits nothing from our
#' bindings, and base `c()` silently strips the advector class when a
#' plain numeric comes first. Rather than asking users to carry that
#' boilerplate, the same three bindings are spliced onto the front of
#' the function's body, which is exactly the spelling RTMB documents.
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
#' is returned untouched, so frmtmb's internal families keep their
#' hand-placed bindings and a double wrap never stacks. Helpers the
#' wrapped function CALLS still use their own environments - lexical
#' scope cannot be injected transitively - so a helper that builds
#' advectors keeps the explicit spelling.
#'
#' @noRd
ad_overload_fn <- function(f) {
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

#' An evaluation enclosure carrying the AD overloads.
#'
#' For `eval(expr, values, enclos)` sites - the nonlinear formula
#' bodies - where the expression is interpreted, not byte-compiled, so
#' active bindings are honored: each access resolves through
#' `RTMB::ADoverload()` at call time, and every other name still
#' reaches the formula environment unchanged.
#'
#' @noRd
ad_overload_env <- function(parent) {
  e <- new.env(parent = parent %||% baseenv())
  for (nm in c("c", "[<-", "diag<-")) {
    local({
      n <- nm
      makeActiveBinding(n, function() RTMB::ADoverload(n), e)
    })
  }
  e
}
