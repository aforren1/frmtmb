# Tape-safe scope for user-written functions and formula bodies.

#' Make a user-written function tape-safe
#'
#' Returns `f` with RTMB's tape-safe replacements for `c()`, `[<-` and
#' `diag<-` bound at the top of its body. Call it on any function a user
#' supplies that frmtmb will later put on an AD tape, so that the user
#' does not have to carry the binding boilerplate. frmtmb applies it to
#' the tape-side slots and to those only: [frmtmb_family()]'s `lpdf`
#' and `lcdf`, and [frmtmb_structure()]'s `loglik`. `valid_y`, `sim`,
#' `post$mean_fn` and the rest run on doubles, off the tape, where the
#' bindings would buy nothing. A package that takes tape-side functions
#' of its own applies this to those. It is safe to call on any object:
#' a non-function and a primitive come back unchanged.
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
#' @seealso [frmtmb_family()] and [frmtmb_structure()], whose tape-side
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
#' `atan2` arrived with RTMB 2.0 and is shadowed although this machine's
#' RTMB 1.9 does not export it: `nl_shadow_fun()` skips a name the
#' installed RTMB lacks, so listing it early is free, and CI's R-devel
#' already has 2.0. It is transparent by construction rather than by
#' probe: RTMB defines it as an S4 generic with methods for `advector`
#' signatures only, so numeric arguments reach `base::atan2()` itself.
#'
#' @noRd
nl_rtmb_shadow <- c(
  "apply", "atan2", "besselI", "besselJ", "besselK", "besselY", "colSums",
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

#' Stan's functions a nonlinear body can call.
#'
#' brms writes a body into its Stan program, so a body there can call any
#' function of the Stan language that takes one row's values, and brms
#' mirrors a list of them in R for its post-processing
#' (`brms:::names_stan_functions()`). These are the members of that list
#' that Stan 2.32 parses with real arguments (dev/nlfix-stanc.R), which
#' is what a row-wise body passes. `softmax()` and `log_softmax()` take a
#' vector only, so no brms body can call them on a row. Each is written
#' for the tape and for doubles alike, so the fitted and the post-fit
#' paths read the same number, and through RTMB's log-space sums where
#' the plain formula loses its derivative at the ends (`logspace_add()`
#' has a finite gradient at x = 800, where `log1p(exp(x))` has none).
#' `test-nl-stan-funs.R` holds the values against Stan's.
#'
#' `step()` is Stan's: 0 below zero and 1 otherwise, with derivative 0,
#' written through `sign()` because a comparison cannot go on the tape.
#'
#' @noRd
nl_stan_funs <- list(
  logit = function(x) log(x) - log1p(-x),
  inv_logit = function(x) RTMB::plogis(x),
  # through the log-space difference: expm1's derivative rounds to 0 long
  # before the true one, 5e-63 at x = 5
  inv_cloglog = function(x) {
    u <- exp(x)
    exp(RTMB::logspace_sub(0 * u, -u))
  },
  Phi = function(x) RTMB::pnorm(x),
  square = function(x) x * x,
  cbrt = function(x) sign(x) * abs(x)^(1 / 3),
  exp2 = function(x) 2^x,
  pow = function(x, y) nl_pow(x, y),
  inv = function(x) 1 / x,
  inv_sqrt = function(x) 1 / sqrt(x),
  # x^(-2) and not 1 / (x * x), whose derivative is 0 / 0 at x = 0
  inv_square = function(x) x^(-2),
  # scaled by the larger side, so that neither square overflows or
  # underflows; sign(m) keeps both sides at 0 from dividing by 0
  hypot = function(x, y) {
    a <- abs(x)
    b <- abs(y)
    m <- 0.5 * (a + b + abs(a - b))
    r <- (a + b - m) / (m + (1 - sign(m)))
    m * sqrt(1 + r * r)
  },
  log1m = function(x) log1p(-x),
  step = function(x) {
    s <- sign(x)
    1 + s / 2 - s * s / 2
  },
  log1p_exp = function(x) RTMB::logspace_add(0 * x, x),
  log1m_exp = function(x) RTMB::logspace_sub(0 * x, x),
  log_diff_exp = function(x, y) {
    RTMB::logspace_sub(x + 0 * y, y + 0 * x)
  },
  log_sum_exp = function(x, y) {
    RTMB::logspace_add(x + 0 * y, y + 0 * x)
  },
  log_inv_logit = function(x) -RTMB::logspace_add(0 * x, -x),
  log1m_inv_logit = function(x) -RTMB::logspace_add(0 * x, x)
)

#' `x^y` with Stan's value and gradient where a data `x` is 0.
#'
#' On the tape, `x^y` differentiates to `x^y * log(x)` in `y`, which is
#' `0 * -Inf` at `x = 0`, so `a * x^b` with a zero row stops with a NaN
#' gradient where Stan's is finite (`stan::math::pow` skips the log term
#' there). The value at 0 is a function of `y` alone: 1 at `y = 0`, 0 above
#' it, Inf below. With `s = sign(y)` and `t = sign(-y)`, which the tape takes
#' as 1 at 0, `(1 - s) / (1 + s)` gives the two sides and
#' `(1 + s) * (1 + t) / 4` the point, and their derivative in `y` is 0. Only a
#' `y` on the tape with a data `x` that holds a 0 takes this branch: everything
#' else is `x^y`, the same operation.
#'
#' @noRd
nl_pow <- function(x, y) {
  if (!inherits(y, "advector") || inherits(x, "advector") || !any(x == 0)) {
    return(x^y)
  }
  "[<-" <- RTMB::ADoverload("[<-")
  n <- max(length(x), length(y))
  x <- rep_len(x, n)
  i0 <- which(x == 0)
  x[i0] <- 1
  s <- sign(y)
  f <- (1 - s) / (1 + s) + (1 + s) * (1 + sign(-y)) / 4
  out <- x^y
  out[i0] <- if (length(f) == 1L) f else f[i0]
  out
}

#' The Stan functions to bind for this body: the names it calls, less one
#' its own environment already resolves. A helper the user wrote under
#' one of these names is found as before. `stats::step()` does not count:
#' it reaches a body only through the search path, and it fits a stepwise
#' regression, which no body means.
#'
#' @noRd
nl_stan_for <- function(expr, parent) {
  if (is.null(expr)) return(character(0))
  # a name in call position, even if the same body also reads a column of
  # that name: eval() looks for a function there and skips the column
  called <- called_functions(expr, names(nl_stan_funs))
  Filter(function(nm) {
    f <- get0(nm, envir = parent %||% baseenv(), mode = "function")
    is.null(f) ||
      (identical(nm, "step") && identical(f, stats::step))
  }, called)
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
  for (nm in nl_stan_for(expr, parent)) {
    assign(nm, nl_stan_funs[[nm]], envir = e)
  }
  # the same zero row for the operator, unless the user's environment
  # redefined it
  if ("^" %in% all.names(expr %||% quote(x)) &&
        identical(get0("^", envir = parent %||% baseenv(), mode = "function"),
                  base::`^`)) {
    assign("^", nl_pow, envir = e)
  }
  e
}
