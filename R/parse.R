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
      supported <- c("weights", "trials", "cens", "trunc", "se",
                     "vint", "vreal", "mi")
      if (!nm %in% supported) {
        stop("Addition term `", nm, "()` is not supported yet ",
             "(currently supported: ",
             paste0(supported, "()", collapse = ", "), ")", call. = FALSE)
      }
      if (nm %in% names(aterms) ||
          (nm == "trunc" && any(c("trunc_lb", "trunc_ub") %in%
                                  names(aterms))) ||
          (nm %in% c("vint", "vreal") &&
             paste0(nm, "1") %in% names(aterms))) {
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
      } else if (nm == "cens") {
        args <- as.list(tm)[-1]
        if (length(args) < 1 || length(args) > 2) {
          stop("cens() takes the censoring code and optionally interval ",
               "upper bounds: cens(c) or cens(c, y2)", call. = FALSE)
        }
        aterms$cens <- args[[1]]
        if (length(args) == 2) aterms$cens_y2 <- args[[2]]
      } else if (nm == "mi") {
        # x | mi() ~ ...: the response may contain NAs; missing entries
        # become latent parameters (one-step imputation). With known
        # measurement SDs - x | mi(sdx) - every value is latent and the
        # observed ones get a measurement model (brms me()).
        if (length(tm) > 2L) {
          stop("mi() on the response side takes at most one argument ",
               "(known measurement SDs)", call. = FALSE)
        }
        aterms$mi <- TRUE
        if (length(tm) == 2L) aterms$mi_sd <- tm[[2]]
      } else if (nm %in% c("vint", "vreal")) {
        # custom-family data vectors (brms vint()/vreal()): each
        # argument becomes aterms$vint1, vint2, ... for the lpdf
        args <- as.list(tm)[-1]
        if (!length(args)) {
          stop(nm, "() needs at least one variable", call. = FALSE)
        }
        for (i in seq_along(args)) {
          aterms[[paste0(nm, i)]] <- args[[i]]
        }
      } else if (nm == "se") {
        args <- as.list(tm)[-1]
        nms <- names(args) %||% rep("", length(args))
        if (sum(nms == "") != 1L || !all(nms %in% c("", "sigma"))) {
          stop("se() takes the known SDs and optionally sigma = TRUE: ",
               "se(x) or se(x, sigma = TRUE)", call. = FALSE)
        }
        aterms$se <- args[[which(nms == "")]]
        if ("sigma" %in% nms) {
          aterms$se_sigma <- isTRUE(eval(args[["sigma"]], baseenv()))
        }
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
# specifications. `shared` is the response-level environment that keeps
# protected-function aliases visible to the combined model frame.
parse_linpred <- function(rhs_form, env, shared = NULL) {
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
  supported_cs <- setdiff(names(covstruct_registry), "smooth")
  cs_specials <- unique(c(supported_cs, "rr", "propto"))
  smooth <- list()
  mo <- list()
  miterms <- list()
  csterms <- list()
  gpterms <- list()
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
    } else if (is.call(tm) && identical(tm[[1]], as.name("mo"))) {
      if (length(tm) != 2L) {
        stop("mo() takes exactly one variable", call. = FALSE)
      }
      mo[[length(mo) + 1L]] <- list(expr = tm[[2]], mult = NULL)
    } else if (is.call(tm) && identical(tm[[1]], as.name("mi"))) {
      if (length(tm) != 2L || !is.name(tm[[2]])) {
        stop("mi() in a predictor takes one variable name: mi(x)",
             call. = FALSE)
      }
      miterms[[length(miterms) + 1L]] <- list(expr = tm[[2]],
                                              mult = NULL)
    } else if (is.call(tm) &&
               (identical(tm[[1]], as.name(":")) ||
                  identical(tm[[1]], as.name("*"))) &&
               any(c("mo", "mi") %in% all.names(tm))) {
      # mo(x):z / mo(x)*z (and mi versions): the special call on one
      # side, a plain multiplier term on the other. `*` also emits the
      # main effects; mo() interactions share their variable's simplex.
      op_star <- identical(tm[[1]], as.name("*"))
      sides <- list(tm[[2]], tm[[3]])
      is_sp <- vapply(sides, function(s) {
        is.call(s) && as.character(s[[1]])[1] %in% c("mo", "mi")
      }, TRUE)
      if (sum(is_sp) != 1L) {
        stop("mo()/mi() interactions need the special on exactly one ",
             "side of ':' or '*': ", deparse1(tm), call. = FALSE)
      }
      sp <- sides[[which(is_sp)]]
      other <- sides[[which(!is_sp)]]
      if (any(c("mo", "mi") %in% all.names(other))) {
        stop("mo()/mi() cannot interact with another mo()/mi() term: ",
             deparse1(tm), call. = FALSE)
      }
      spn <- as.character(sp[[1]])[1]
      entry <- list(expr = sp[[2]], mult = other)
      if (spn == "mo") {
        mo[[length(mo) + 1L]] <- entry
        if (op_star) mo[[length(mo) + 1L]] <- list(expr = sp[[2]],
                                                   mult = NULL)
      } else {
        if (!is.name(sp[[2]])) {
          stop("mi() in a predictor takes one variable name: mi(x)",
               call. = FALSE)
        }
        miterms[[length(miterms) + 1L]] <- entry
        if (op_star) miterms[[length(miterms) + 1L]] <-
          list(expr = sp[[2]], mult = NULL)
      }
      if (op_star) rest[[length(rest) + 1L]] <- other
    } else if (is.call(tm) && identical(tm[[1]], as.name("gp"))) {
      aa <- as.list(tm)[-1]
      nms <- names(aa) %||% rep("", length(aa))
      vars <- aa[nms == ""]
      if (length(vars) < 1L || !all(nms %in% c("", "k", "c", "iso"))) {
        stop("gp() takes 1-3 variables plus optional k = (basis size ",
             "per dimension), c = (boundary factor), and iso = ",
             "(shared lengthscale): gp(x), gp(x, k = 30), gp(x1, x2)",
             call. = FALSE)
      }
      if (length(vars) > 3L) {
        stop("gp() supports at most 3 dimensions (got ", length(vars),
             ")", call. = FALSE)
      }
      gpterms[[length(gpterms) + 1L]] <- list(
        exprs = vars,
        k = if (!is.null(aa$k)) as.integer(eval(aa$k, baseenv())),
        c = if (!is.null(aa$c)) as.numeric(eval(aa$c, baseenv()))
            else 1.25,
        iso = if (!is.null(aa$iso)) isTRUE(eval(aa$iso, baseenv()))
              else FALSE
      )
    } else if (is.call(tm) && identical(tm[[1]], as.name("cs")) &&
               !("|" %in% all.names(tm))) {
      # barless cs(x): category-specific ordinal effect (the bar form
      # cs(x | g) stays a compound-symmetry covariance structure)
      if (length(tm) != 2L) {
        stop("cs() takes one variable: cs(x)", call. = FALSE)
      }
      csterms[[length(csterms) + 1L]] <- tm[[2]]
    } else {
      for (sp_nm in c("mo", "mi")) {
        if (sp_nm %in% all.names(tm)) {
          stop(sp_nm, "() is only supported as a standalone additive ",
               "term or a two-way ':'/'*' interaction: ",
               deparse1(tm), call. = FALSE)
        }
      }
      rest[[length(rest) + 1L]] <- tm
    }
  }
  # splitForm silently strips ANY term whose expression mentions a
  # special name - even inside I() or nested calls - so exp(x) in a
  # fixed formula would vanish. Protect by rewriting special-named
  # CALLS (never bar terms, which were routed to splitForm above) to
  # .frm_<name> aliases bound in a child environment.
  prot <- new.env(parent = env)
  protected <- FALSE
  sub_specials <- function(e) {
    if (is.call(e)) {
      hd <- e[[1]]
      if (is.name(hd) && as.character(hd) %in% cs_specials &&
          !("|" %in% all.names(e))) {
        nm <- as.character(hd)
        pn <- paste0(".frm_", nm)
        fn <- tryCatch(eval(hd, env), error = function(err) NULL)
        if (is.function(fn)) {
          assign(pn, fn, envir = prot)
          if (!is.null(shared)) assign(pn, fn, envir = shared)
          e[[1]] <- as.name(pn)
          protected <<- TRUE
        }
      }
      if (length(e) > 1L) {
        for (i in seq_along(e)[-1]) {
          ei <- e[[i]]
          if (!(is.symbol(ei) && !nzchar(as.character(ei)))) {
            e[[i]] <- sub_specials(ei)
          }
        }
      }
    }
    e
  }
  rest <- lapply(rest, sub_specials)
  env_lp <- if (protected) prot else env

  if (length(rest)) {
    rhs_expr <- Reduce(function(a, b) call("+", a, b), rest)
  } else {
    rhs_expr <- 1
  }
  bare_form <- stats::as.formula(call("~", rhs_expr), env = env_lp)
  # newer structure names (hetar1, homcs, ...) are not in reformulas's
  # default specials list; pass the full set explicitly (rr/propto so
  # they reach the informative not-supported error below)
  sf <- reformulas::splitForm(bare_form, defaultTerm = "us",
                              specials = cs_specials)
  bad <- setdiff(sf$reTrmClasses, supported_cs)
  if (length(bad)) {
    stop("Covariance structure(s) not supported yet: ",
         paste(unique(bad), collapse = ", "),
         " (currently supported: ",
         paste(supported_cs, collapse = ", "), ")", call. = FALSE)
  }
  re <- Map(function(bar, cls, addargs) {
    # brms |ID| syntax: (x | p | g) parses as ((x | p) | g); the middle
    # element keys random-effect correlation across formulas
    id <- NULL
    cov_expr <- NULL
    if (cls %in% c("gp", "hsgp")) {
      stop("gp() is not a bar term; write gp(x) or gp(x, k = 30)",
           call. = FALSE)
    }
    rank <- NULL
    if (cls == "rr") {
      aa <- as.list(addargs)[-1]
      rank <- as.integer(eval(aa$d %||% 2, baseenv()))
      if (is.na(rank) || rank < 1L) {
        stop("rr() needs a positive integer rank: rr(x | g, d = 2)",
             call. = FALSE)
      }
    }
    if (cls == "equalto") {
      aa <- as.list(addargs)[-1]
      aa <- aa[!nzchar(names(aa) %||% rep("", length(aa)))]
      if (length(aa) != 1L) {
        stop("equalto() needs the fixed covariance matrix: ",
             "equalto(x + 0 | g, V)", call. = FALSE)
      }
      cov_expr <- aa[[1L]]
    }
    if (is.call(bar[[2]]) && identical(bar[[2]][[1]], as.name("|"))) {
      if (cls != "us") {
        stop("|ID| correlation is only supported for default (us) ",
             "random-effect terms", call. = FALSE)
      }
      id <- paste0(deparse1(bar[[2]][[3]]), "|", deparse1(bar[[3]]))
      bar <- call("|", bar[[2]][[2]], bar[[3]])
    }
    # brms (x | gr(g, cov = A)): known covariance over the levels;
    # gr(g, prec = Q) takes a (sparse) precision matrix instead
    if (is.call(bar[[3]]) && identical(bar[[3]][[1]], as.name("gr"))) {
      ga <- as.list(bar[[3]])[-1]
      nms <- names(ga) %||% rep("", length(ga))
      gvar <- ga[nms == ""]
      has_cov <- !is.null(ga$cov)
      has_prec <- !is.null(ga$prec)
      if (length(gvar) != 1 || (has_cov + has_prec) != 1 ||
          !all(nms %in% c("", "cov", "prec"))) {
        stop("gr() supports (x | gr(g, cov = A)) or ",
             "(1 | gr(g, prec = Q))", call. = FALSE)
      }
      cov_expr <- ga$cov %||% ga$prec
      bar <- call("|", bar[[2]], gvar[[1]])
      cls <- if (has_cov) "gr_cov" else "gr_prec"
    }
    list(bar = bar, group = bar[[3]], covstruct = cls, id = id,
         cov_expr = cov_expr, rank = rank)
  }, sf$reTrmFormulas, sf$reTrmClasses, sf$reTrmAddArgs)
  names(re) <- vapply(re, function(z) deparse1(z$bar), "")

  fixed <- sf$fixedFormula
  environment(fixed) <- env_lp
  list(fixed = fixed, re = re, smooth = smooth, mo = mo,
       miterms = miterms, csterms = csterms, gpterms = gpterms,
       rhs = rhs_form)
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
  # response-level alias environment: parse_linpred registers protected
  # special-named functions here so the combined model frame (which
  # evaluates with formula_env) can resolve them
  shared_env <- new.env(parent = env)
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
                          environment(pf) %||% env, shared_env)
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
                            environment(pf) %||% env, shared_env)
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
      formula_env = shared_env
    ))
  }

  primaries <- fam$primary_dpars %||% "mu"
  pfix <- bform$pfix
  if (!is.null(ri$aterms$se) && !isTRUE(ri$aterms$se_sigma) &&
      "sigma" %in% fam$dpars &&
      !"sigma" %in% c(names(bform$pforms), names(pfix))) {
    # se() without sigma = TRUE: the residual SD is the known se alone,
    # so the sigma dpar is mapped out (its value is unused by the lpdf)
    pfix$sigma <- 1
  }
  extra <- c(names(bform$pforms), names(pfix))
  allowed <- setdiff(fam$dpars, primaries[1])
  unknown <- setdiff(extra, allowed)
  if (length(unknown)) {
    stop("dpar(s) not available for family '", fam$family, "': ",
         paste(unknown, collapse = ", "),
         " (available: ", paste(allowed, collapse = ", "), ")",
         call. = FALSE)
  }

  main_lp <- parse_linpred(reformulas::RHSForm(f, as.form = TRUE), env,
                           shared_env)

  dpars <- list()
  for (dp in fam$dpars) {
    if (dp %in% names(bform$pforms)) {
      pf <- bform$pforms[[dp]]
      lp <- parse_linpred(reformulas::RHSForm(pf, as.form = TRUE),
                          environment(pf) %||% env, shared_env)
      dpars[[dp]] <- c(list(name = dp, link = fam$links[[dp]],
                            constant = NULL), lp)
    } else if (dp %in% primaries) {
      dpars[[dp]] <- c(list(name = dp, link = fam$links[[dp]],
                            constant = NULL), main_lp)
    } else {
      dpars[[dp]] <- plain_dpar(dp, fam, pfix[[dp]])
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
    formula_env = shared_env
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
