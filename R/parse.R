# Parse the left-hand side of a formula into the response expression and
# addition terms (brms aterms): `y | weights(w) + trials(n) ~ ...`.
parse_response <- function(formula) {
  lhs <- formula[[2]]
  aterms <- list()
  resp <- lhs
  if (is.call(lhs) && identical(lhs[[1]], as.name("|"))) {
    resp <- lhs[[2]]
    for (tm in split_plus(lhs[[3]])) {
      if (!is.call(tm)) {
        stop("Malformed addition term: ", deparse1(tm), call. = FALSE)
      }
      nm <- as.character(tm[[1]])
      supported <- c("weights", "trials", "cens", "trunc")
      if (!nm %in% supported) {
        stop("Addition term `", nm, "()` is not supported yet ",
             "(currently supported: ",
             paste0(supported, "()", collapse = ", "), ")", call. = FALSE)
      }
      if (nm %in% names(aterms) ||
          (nm == "trunc" && any(c("trunc_lb", "trunc_ub") %in%
                                  names(aterms)))) {
        stop("Duplicated addition term `", nm, "()`", call. = FALSE)
      }
      if (nm == "trunc") {
        args <- as.list(tm)[-1]
        if (!length(args) || is.null(names(args)) ||
            !all(names(args) %in% c("lb", "ub"))) {
          stop("trunc() takes named bounds: trunc(lb = ...), ",
               "trunc(ub = ...), or both", call. = FALSE)
        }
        if (!is.null(args$lb)) aterms$trunc_lb <- args$lb
        if (!is.null(args$ub)) aterms$trunc_ub <- args$ub
      } else {
        if (length(tm) != 2) {
          stop("`", nm, "()` takes exactly one argument", call. = FALSE)
        }
        aterms[[nm]] <- tm[[2]]
      }
    }
  }
  list(resp = resp, aterms = aterms)
}

# Split one linear-predictor RHS (a one-sided formula) into a parametric
# fixed formula, random-effect terms (reformulas), and mgcv smooth
# specifications.
parse_linpred <- function(rhs_form, env) {
  environment(rhs_form) <- env
  rhs_form <- reformulas::expandDoubleVerts(rhs_form)

  # Pull smooth terms out at the top level before splitForm sees them
  # (splitForm strips s()/t2() regardless of its specials argument).
  # Evaluating the calls with mgcv's own constructors parses k=, by=,
  # bs=, and multi-variable smooths for free.
  terms_list <- split_plus(reformulas::RHSForm(rhs_form))
  is_smooth_call <- function(tm) {
    is.call(tm) && as.character(tm[[1]])[1] %in% c("s", "t2", "te", "ti")
  }
  smooth <- list()
  rest <- list()
  for (tm in terms_list) {
    if (is_smooth_call(tm)) {
      fn <- as.character(tm[[1]])[1]
      if (fn %in% c("te", "ti")) {
        stop("te() and ti() smooths are not supported (no random-effect ",
             "representation); use t2() instead", call. = FALSE)
      }
      smooth[[length(smooth) + 1L]] <-
        eval(tm, list(s = mgcv::s, t2 = mgcv::t2), enclos = env)
    } else {
      rest[[length(rest) + 1L]] <- tm
    }
  }
  if (length(rest)) {
    rhs_expr <- Reduce(function(a, b) call("+", a, b), rest)
  } else {
    rhs_expr <- 1
  }
  bare_form <- stats::as.formula(call("~", rhs_expr), env = env)
  sf <- reformulas::splitForm(bare_form, defaultTerm = "us")

  supported_cs <- setdiff(names(covstruct_registry), "smooth")
  bad <- setdiff(sf$reTrmClasses, supported_cs)
  if (length(bad)) {
    stop("Covariance structure(s) not supported yet: ",
         paste(unique(bad), collapse = ", "),
         " (currently supported: ",
         paste(supported_cs, collapse = ", "), ")", call. = FALSE)
  }
  re <- Map(function(bar, cls) {
    # brms |ID| syntax: (x | p | g) parses as ((x | p) | g); the middle
    # element keys random-effect correlation across formulas
    id <- NULL
    if (is.call(bar[[2]]) && identical(bar[[2]][[1]], as.name("|"))) {
      if (cls != "us") {
        stop("|ID| correlation is only supported for default (us) ",
             "random-effect terms", call. = FALSE)
      }
      id <- paste0(deparse1(bar[[2]][[3]]), "|", deparse1(bar[[3]]))
      bar <- call("|", bar[[2]][[2]], bar[[3]])
    }
    list(bar = bar, group = bar[[3]], covstruct = cls, id = id)
  }, sf$reTrmFormulas, sf$reTrmClasses)
  names(re) <- vapply(re, function(z) deparse1(z$bar), "")

  fixed <- sf$fixedFormula
  environment(fixed) <- env
  list(fixed = fixed, re = re, smooth = smooth, rhs = rhs_form)
}

# Default (intercept-only) or constant dpar spec.
plain_dpar <- function(dp, fam, constant = NULL) {
  if (!is.null(constant)) {
    lv <- fam$links[[dp]]$linkfun(constant)
    if (!is.finite(lv)) {
      stop("Constant ", dp, " = ", constant, " is not in the range of ",
           "the ", fam$links[[dp]]$name, " link", call. = FALSE)
    }
  }
  list(name = dp, link = fam$links[[dp]], fixed = ~1, re = list(),
       rhs = ~1, smooth = list(), constant = constant)
}

# One bf() -> one response entry of the spec.
parse_one_response <- function(bform) {
  fam <- bform$family
  if (is.null(fam)) {
    stop("No family specified. Attach one with `bf(...) + gaussian()` or ",
         "the `family` argument of frm()", call. = FALSE)
  }
  f <- bform$formula
  if (length(f) != 3L) {
    stop("The model formula needs a response (left-hand side)", call. = FALSE)
  }
  env <- environment(f) %||% globalenv()
  ri <- parse_response(f)

  if (isTRUE(bform$nl)) {
    if (!identical(fam$primary_dpars %||% "mu", "mu")) {
      stop("nl = TRUE requires a family with a single 'mu' location ",
           "parameter", call. = FALSE)
    }
    body <- reformulas::RHSForm(f)
    nlpars <- setdiff(names(bform$pforms), fam$dpars)
    if (!length(nlpars)) {
      stop("nl = TRUE needs at least one nonlinear-parameter formula ",
           "whose name appears in the model formula", call. = FALSE)
    }
    body_vars <- all.vars(body)
    miss <- setdiff(nlpars, body_vars)
    if (length(miss)) {
      stop("Nonlinear parameter(s) not used in the model formula: ",
           paste(miss, collapse = ", "), call. = FALSE)
    }
    datavars <- setdiff(body_vars, nlpars)

    dpars <- list()
    for (np in nlpars) {
      pf <- bform$pforms[[np]]
      lp <- parse_linpred(reformulas::RHSForm(pf, as.form = TRUE),
                          environment(pf) %||% env)
      dpars[[np]] <- c(list(name = np, link = get_link("identity"),
                            constant = NULL), lp)
    }
    dpars$mu <- list(name = "mu", link = fam$links$mu, nl_body = body,
                     datavars = datavars, nl_env = env,
                     fixed = NULL, re = list(), rhs = NULL,
                     smooth = list(), constant = NULL)
    for (dp in setdiff(fam$dpars, "mu")) {
      if (dp %in% names(bform$pforms)) {
        pf <- bform$pforms[[dp]]
        lp <- parse_linpred(reformulas::RHSForm(pf, as.form = TRUE),
                            environment(pf) %||% env)
        dpars[[dp]] <- c(list(name = dp, link = fam$links[[dp]],
                              constant = NULL), lp)
      } else {
        dpars[[dp]] <- plain_dpar(dp, fam, bform$pfix[[dp]])
      }
    }
    return(list(
      resp_name = deparse1(ri$resp),
      resp_expr = ri$resp,
      family = fam,
      aterms = ri$aterms,
      dpars = dpars,
      primary_dpars = nlpars,   # REML integrates the nlpar coefficients
      nlpars = nlpars,
      formula_env = env
    ))
  }

  primaries <- fam$primary_dpars %||% "mu"
  extra <- c(names(bform$pforms), names(bform$pfix))
  allowed <- setdiff(fam$dpars, primaries[1])
  unknown <- setdiff(extra, allowed)
  if (length(unknown)) {
    stop("dpar(s) not available for family '", fam$family, "': ",
         paste(unknown, collapse = ", "),
         " (available: ", paste(allowed, collapse = ", "), ")",
         call. = FALSE)
  }

  main_lp <- parse_linpred(reformulas::RHSForm(f, as.form = TRUE), env)

  dpars <- list()
  for (dp in fam$dpars) {
    if (dp %in% names(bform$pforms)) {
      pf <- bform$pforms[[dp]]
      lp <- parse_linpred(reformulas::RHSForm(pf, as.form = TRUE),
                          environment(pf) %||% env)
      dpars[[dp]] <- c(list(name = dp, link = fam$links[[dp]],
                            constant = NULL), lp)
    } else if (dp %in% primaries) {
      dpars[[dp]] <- c(list(name = dp, link = fam$links[[dp]],
                            constant = NULL), main_lp)
    } else {
      dpars[[dp]] <- plain_dpar(dp, fam, bform$pfix[[dp]])
    }
  }

  list(
    resp_name = deparse1(ri$resp),
    resp_expr = ri$resp,
    family = fam,
    aterms = ri$aterms,
    dpars = dpars,
    primary_dpars = primaries,
    nlpars = character(0),
    formula_env = env
  )
}

# frmtmb_formula / frmtmb_mvformula -> frmtmb_spec.
parse_spec <- function(bform) {
  if (inherits(bform, "frmtmb_mvformula")) {
    resps <- lapply(bform$forms, parse_one_response)
    names(resps) <- vapply(resps, `[[`, "", "resp_name")
    if (anyDuplicated(names(resps))) {
      stop("Duplicated response in mvbf(): ",
           names(resps)[duplicated(names(resps))][1], call. = FALSE)
    }
    rescor <- isTRUE(bform$rescor)
    if (rescor) {
      fams <- vapply(resps, function(r) r$family$family, "")
      if (!all(fams == "gaussian")) {
        stop("rescor = TRUE requires all responses to be gaussian ",
             "(got: ", paste(unique(fams), collapse = ", "), ")",
             call. = FALSE)
      }
    }
    return(structure(
      list(responses = resps, rescor = rescor, re_ids = list()),
      class = "frmtmb_spec"
    ))
  }
  stopifnot(inherits(bform, "frmtmb_formula"))
  resp <- parse_one_response(bform)
  structure(
    list(responses = stats::setNames(list(resp), resp$resp_name),
         rescor = FALSE, re_ids = list()),
    class = "frmtmb_spec"
  )
}

#' @export
print.frmtmb_spec <- function(x, ...) {
  cat("<frmtmb spec>\n")
  for (r in x$responses) {
    cat("Response: ", r$resp_name, "  [", r$family$family, "]\n", sep = "")
    if (length(r$aterms)) {
      cat("  aterms: ",
          paste0(names(r$aterms), "(", vapply(r$aterms, deparse1, ""), ")",
                 collapse = ", "), "\n", sep = "")
    }
    for (dp in r$dpars) {
      if (!is.null(dp$constant)) {
        cat("  ", dp$name, " = ", dp$constant, " (fixed)\n", sep = "")
        next
      }
      cat("  ", dp$name, " (", dp$link$name, "): ",
          deparse1(dp$fixed), sep = "")
      if (length(dp$smooth)) {
        cat(" + smooths: ",
            paste(vapply(dp$smooth, function(s) s$label %||%
                           paste0("s(", paste(s$term, collapse = ","), ")"),
                         ""), collapse = ", "), sep = "")
      }
      if (length(dp$re)) {
        cat(" + RE: ",
            paste0(vapply(dp$re, function(z) deparse1(z$bar), ""),
                   " [", vapply(dp$re, `[[`, "", "covstruct"),
                   ifelse(vapply(dp$re, function(z) is.null(z$id), TRUE),
                          "", ", ID"), "]",
                   collapse = ", "), sep = "")
      }
      cat("\n")
    }
  }
  if (x$rescor) cat("Residual correlation: yes\n")
  invisible(x)
}
