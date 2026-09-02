#' Parse the left-hand side of a formula into the response expression and
#' addition terms (brms aterms): `y | weights(w) + trials(n) ~ ...`.
#'
#' @noRd
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
          aterms$se_sigma <- eval_spec_arg(args[["sigma"]], "sigma",
                                           environment(formula),
                                           fn = "se")
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

#' glm(), lme4 and glmmTMB all spell a two-column binomial response
#' `cbind(successes, failures)`, so ported code lands here before it
#' reaches anything else. Rewrite it into the internal
#' `successes | trials(successes + failures)` form: every later stage
#' (frame assembly, valid_y, fitted, predict, simulate) then sees an
#' ordinary count response and needs no matrix branch of its own.
#' Other matrix-response families (multinomial) keep their own spelling.
#' `[glmmTMB#1319, #1325]`
#'
#' @noRd
rewrite_cbind_response <- function(ri, fam) {
  resp <- ri$resp
  if (!is.call(resp) || !identical(resp[[1L]], as.name("cbind")) ||
      !fam$family %in% c("binomial", "beta_binomial")) {
    return(ri)
  }
  args <- as.list(resp)[-1L]
  nms <- names(args) %||% rep("", length(args))
  if (length(args) != 2L || any(nzchar(nms))) {
    stop("A cbind() ", fam$family, " response takes exactly two unnamed ",
         "columns: cbind(successes, failures)", call. = FALSE)
  }
  if (!is.null(ri$aterms$trials)) {
    stop("cbind(successes, failures) already carries the number of ",
         "trials; drop the trials() addition term (write either ",
         "cbind(s, f) ~ ... or s | trials(n) ~ ...)", call. = FALSE)
  }
  ri$resp <- args[[1L]]
  ri$aterms$trials <- call("+", args[[1L]], args[[2L]])
  ri$cbind_resp <- TRUE
  ri
}

#' Evaluate one special-term tuning argument (gp's k/c/iso, rr's d,
#' se's sigma flag). These are user code, so they resolve in the
#' environment the formula was written in (brms does the same); only the
#' evaluated value reaches the spec, so the frame and predict() never
#' re-evaluate them. gp's c is per-dimension (brms convention), so it
#' alone may be a vector; its length-vs-D check lives in the frame,
#' which knows D.
#'
#' @noRd
eval_spec_arg <- function(expr, nm, env, fn = "gp") {
  val <- tryCatch(eval(expr, env), error = function(e) {
    stop(fn, "(): cannot evaluate ", nm, " = ", deparse1(expr), ": ",
         conditionMessage(e), call. = FALSE)
  })
  if (length(val) != 1L && nm != "c") {
    stop(fn, "(): ", nm, " = ", deparse1(expr),
         " must be a single value (got length ", length(val), ")",
         call. = FALSE)
  }
  if (nm %in% c("iso", "sigma")) {
    if (!is.logical(val) || is.na(val)) {
      stop(fn, "(): ", nm, " = ", deparse1(expr),
           " must be TRUE or FALSE", call. = FALSE)
    }
    return(isTRUE(val))
  }
  if (!is.numeric(val) || !length(val) || !all(is.finite(val))) {
    stop(fn, "(): ", nm, " = ", deparse1(expr),
         " must be finite and numeric", call. = FALSE)
  }
  if (nm %in% c("k", "d")) {
    if (val < 1 || val != trunc(val)) {
      stop(fn, "(): ", nm, " = ", deparse1(expr),
           " must be a positive whole number", call. = FALSE)
    }
    return(as.integer(val))
  }
  if (any(val <= 0)) {
    stop(fn, "(): ", nm, " = ", deparse1(expr), " must be positive",
         call. = FALSE)
  }
  as.numeric(val)
}

#' Positional-or-named argument matching for the predictor specials that
#' take a fixed argument list (brms's car() and our spde()). R's own
#' partial matching is deliberately not reproduced: a misspelled name
#' should be an error, not a silent match.
#'
#' @noRd
match_special_args <- function(call_expr, argn, fn) {
  aa <- as.list(call_expr)[-1]
  nms <- names(aa) %||% rep("", length(aa))
  bad <- setdiff(nms[nzchar(nms)], argn)
  if (length(bad)) {
    stop(fn, "(): unknown argument(s) ", paste(bad, collapse = ", "),
         " (takes ", paste(argn, collapse = ", "), ")", call. = FALSE)
  }
  out <- stats::setNames(vector("list", length(argn)), argn)
  for (i in which(nzchar(nms))) out[[nms[i]]] <- aa[[i]]
  free <- argn[vapply(out, is.null, TRUE)]
  pos <- which(!nzchar(nms))
  if (length(pos) > length(free)) {
    stop(fn, "(): too many arguments (takes ",
         paste(argn, collapse = ", "), ")", call. = FALSE)
  }
  for (i in seq_along(pos)) out[[free[i]]] <- aa[[pos[i]]]
  out
}

#' brms's car(M, gr, type): M is an adjacency matrix (looked up in the
#' data or the formula environment, as gr(cov = A) is), gr the grouping
#' variable naming the locations, type one of brms's four.
#'
#' @noRd
parse_car_call <- function(tm, env) {
  a <- match_special_args(tm, c("M", "gr", "type", "con_sd"), "car")
  if (is.null(a$M) || is.null(a$gr)) {
    stop("car() needs an adjacency matrix and a grouping variable: ",
         "car(M, gr = g, type = \"icar\")", call. = FALSE)
  }
  # brms's gr = NA default (one location per observation) is deprecated
  # there and never supported here: the field's levels come from the
  # factor, so the call has to name one
  if (identical(a$gr, NA) || identical(a$gr, quote(NA))) {
    stop("car(): gr must name a grouping variable. brms's gr = NA ",
         "(one location per observation) is deprecated; build the ",
         "location factor and pass it", call. = FALSE)
  }
  type <- if (is.null(a$type)) "escar" else {
    tp <- eval(a$type, env)
    if (!is.character(tp) || length(tp) != 1L || !(tp %in% car_types)) {
      stop("car(): type must be one of ",
           paste0("\"", car_types, "\"", collapse = ", "), call. = FALSE)
    }
    tp
  }
  con_sd <- if (is.null(a$con_sd)) car_con_sd_default else {
    eval_spec_arg(a$con_sd, "con_sd", env, fn = "car")
  }
  list(M_expr = a$M, gr_expr = a$gr, type = type, con_sd = con_sd,
       label = deparse1(tm))
}

#' spde(fem, gr): fem is a list of the mesh's finite-element matrices,
#' gr the factor mapping observations to mesh nodes.
#'
#' @noRd
parse_spde_call <- function(tm, env) {
  a <- match_special_args(tm, c("fem", "gr"), "spde")
  if (is.null(a$fem) || is.null(a$gr)) {
    stop("spde() needs the mesh matrices and a grouping variable: ",
         "spde(fm_fem(mesh), gr = node)", call. = FALSE)
  }
  list(fem_expr = a$fem, gr_expr = a$gr, label = deparse1(tm))
}

#' Drop redundant parentheses so a term can be inspected and re-split.
#'
#' @noRd
strip_parens <- function(e) {
  while (is.call(e) && identical(e[[1]], as.name("("))) e <- e[[2]]
  e
}

#' `(x || g)` promises uncorrelated effects, but lme4's expansion is
#' purely syntactic: a FACTOR contributes a single term carrying all of
#' its contrast columns, which then lands in one default (`us`) block and
#' comes back fully correlated - the opposite of what was asked for
#' `[lme4#818]`. Expand each `||` term the way lme4 does, so the block
#' structure users already rely on is preserved, then tag every piece
#' `diag()`, the independent-variance structure for any column count. A
#' one-column piece has no correlation to lose, so `diag` and `us`
#' coincide there and numeric double bars are bit-identical to before.
#'
#' @noRd
expand_double_verts <- function(form) {
  raw <- split_plus(reformulas::RHSForm(form))
  out <- list()
  for (tm in raw) {
    inner <- strip_parens(tm)
    if (!(is.call(inner) && identical(inner[[1]], as.name("||")))) {
      out[[length(out) + 1L]] <- tm
      next
    }
    ex <- reformulas::expandDoubleVerts(
      stats::as.formula(call("~", inner), env = environment(form)))
    for (p in split_plus(strip_parens(reformulas::RHSForm(ex)))) {
      out[[length(out) + 1L]] <- call("diag", strip_parens(p))
    }
  }
  reformulas::RHSForm(form) <- Reduce(function(a, b) call("+", a, b), out)
  form
}

#' Split one linear-predictor RHS (a one-sided formula) into a parametric
#' fixed formula, random-effect terms (reformulas), and mgcv smooth
#' specifications. `shared` is the response-level environment that keeps
#' protected-function aliases visible to the combined model frame.
#'
#' This is the first stage of the formula-to-design-matrix pipeline and
#' the only place that reads predictor syntax. It takes the RHS of one
#' dpar formula plus the environment it was written in, and returns a
#' list with one slot per term family: `fixed` (the parametric formula),
#' `re` (bar terms with their covariance structure, `|ID|` key, and
#' known covariance or rank), `smooth`, `mo`, `miterms`, `csterms`,
#' `gpterms`, `carterms`, `spdeterms`, and `rhs` (the expanded formula).
#' Every entry is an unevaluated expression plus its tuning values; the
#' data is never touched here. `assemble_frame()` turns these slots into
#' design matrices and random-effect blocks.
#'
#' @noRd
parse_linpred <- function(rhs_form, env, shared = NULL) {
  environment(rhs_form) <- env
  rhs_form <- expand_double_verts(rhs_form)

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
  carterms <- list()
  spdeterms <- list()
  acterms <- list()
  rest <- list()
  for (tm in terms_list) {
    # `x * (1 | g)` and `x:(1 | g)` are almost always a typo for `+`.
    # splitForm hoists the bar out and silently fits the `+` model, so
    # the user never learns the interaction was ignored. [lme4#196]
    if (is.call(tm) &&
        as.character(tm[[1]])[1] %in% c(":", "*", "/") &&
        any(c("|", "||") %in% all.names(tm))) {
      stop("A random-effect term cannot be crossed with '",
           as.character(tm[[1]])[1], "': ", deparse1(tm),
           ". Did you mean '+'?", call. = FALSE)
    }
    # same trap for the autocor terms: without this, ar(...) inside an
    # interaction reaches the model matrix and dies inside stats::ar
    # with a message that points nowhere near the mistake
    if (is.call(tm) &&
        as.character(tm[[1]])[1] %in% c(":", "*", "/") &&
        any(autocor_structs %in% all.names(tm))) {
      stop("An autocorrelation term cannot be crossed with '",
           as.character(tm[[1]])[1], "': ", deparse1(tm),
           ". Write it as a separate term: y ~ ... + ",
           intersect(autocor_structs, all.names(tm))[1L], "(...)",
           call. = FALSE)
    }
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
      # the interaction branch below repeats this check with its own
      # message, because there the offending mi() sits inside a `:`
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
          stop("mi() in an interaction takes one variable name: ",
               "mi(x):z", call. = FALSE)
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
        k = if (!is.null(aa$k)) eval_spec_arg(aa$k, "k", env),
        c = if (!is.null(aa$c)) eval_spec_arg(aa$c, "c", env) else 1.25,
        iso = if (!is.null(aa$iso)) eval_spec_arg(aa$iso, "iso", env)
              else FALSE
      )
    } else if (is.call(tm) && identical(tm[[1]], as.name("car"))) {
      carterms[[length(carterms) + 1L]] <- parse_car_call(tm, env)
    } else if (is.call(tm) && identical(tm[[1]], as.name("spde"))) {
      spdeterms[[length(spdeterms) + 1L]] <- parse_spde_call(tm, env)
    } else if (is.call(tm) &&
               as.character(tm[[1]])[1] %in% autocor_structs &&
               !("|" %in% all.names(tm))) {
      # brms R-side autocorrelation: ar(), ma(), arma(), cosy(), unstr()
      acterms[[length(acterms) + 1L]] <- parse_autocor_call(tm, env)
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
  # default specials list; pass the full set explicitly (propto so it
  # reaches the informative not-supported error below)
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
    id_label <- NULL
    id_group <- NULL
    cov_expr <- NULL
    if (cls %in% c("gp", "hsgp")) {
      stop("gp() is not a bar term; write gp(x) or gp(x, k = 30)",
           call. = FALSE)
    }
    if (cls %in% c("car", "spde")) {
      stop(cls, "() is not a bar term; write ",
           if (cls == "car") "car(M, gr = g)" else "spde(fem, gr = g)",
           " as a predictor term, not on the left of ( | )", call. = FALSE)
    }
    rank <- NULL
    if (cls == "rr") {
      aa <- as.list(addargs)[-1]
      rank <- if (is.null(aa$d)) 2L else {
        eval_spec_arg(aa$d, "d", env, fn = "rr")
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
      id_label <- deparse1(bar[[2]][[3]])
      # the merge key carries the grouping expression as well, so two
      # terms only ever merge when they name the same factor (and, for
      # gr(), the same relationship matrix). check_id_covstructs()
      # refuses one label spread over several grouping expressions
      # rather than letting them drift into separate blocks.
      id_group <- bar[[3]]
      id <- paste0(id_label, "|", deparse1(id_group))
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
    if (is.call(bar[[3]]) &&
        as.character(bar[[3]][[1]])[1] %in% c("car", "spde")) {
      nm <- as.character(bar[[3]][[1]])[1]
      stop(nm, "() is not a bar term; write ",
           if (nm == "car") "car(M, gr = g)" else "spde(fem, gr = g)",
           " as a predictor term, not as the grouping factor of ( | )",
           call. = FALSE)
    }
    list(bar = bar, group = bar[[3]], covstruct = cls, id = id,
         id_label = id_label, id_group = id_group,
         cov_expr = cov_expr, rank = rank)
  }, sf$reTrmFormulas, sf$reTrmClasses, sf$reTrmAddArgs)
  names(re) <- vapply(re, function(z) deparse1(z$bar), "")

  fixed <- sf$fixedFormula
  environment(fixed) <- env_lp
  list(fixed = fixed, re = re, smooth = smooth, mo = mo,
       miterms = miterms, csterms = csterms, gpterms = gpterms,
       carterms = carterms, spdeterms = spdeterms, acterms = acterms,
       rhs = rhs_form)
}

#' Lift the response's residual-correlation term out of its linear
#' predictors.
#'
#' An R-side term changes the shape of the LIKELIHOOD, not of a linear
#' predictor, so it belongs to the response and not to a dpar. brms
#' takes the same view and refuses one written anywhere but `mu`
#' ("Explicit covariance terms can only be specified on 'mu'"), which
#' is also what rules it out of mixture models, where every location
#' parameter is `mu1`, `mu2`, ... rather than `mu`.
#'
#' @noRd
pull_autocor <- function(dpars, resp_name) {
  found <- list()
  for (nm in names(dpars)) {
    ats <- dpars[[nm]]$acterms %||% list()
    if (!length(ats)) next
    if (!identical(nm, "mu")) {
      stop("Residual correlation terms can only be written on 'mu'; ",
           ats[[1L]]$label, " appears in the formula for '", nm,
           "'. The term changes the residual density of the response, ",
           "not a linear predictor",
           if (grepl("^mu[0-9]", nm)) {
             paste0(". '", nm, "' is a mixture component, and a mixture ",
                    "likelihood has no single residual to correlate; ",
                    "an ar1()/toep() random effect over the time factor ",
                    "is the available alternative there")
           } else "", call. = FALSE)
    }
    found <- c(found, ats)
    dpars[[nm]]$acterms <- list()
  }
  if (length(found) > 1L) {
    stop("Response '", resp_name, "' carries ", length(found),
         " residual correlation terms (",
         paste(vapply(found, `[[`, "", "label"), collapse = ", "),
         "); a response has one residual covariance, so keep one. ",
         "Structures can be nested through random effects instead, ",
         "e.g. ar(week, subj, cov = TRUE) + (1 | site)", call. = FALSE)
  }
  list(dpars = dpars, autocor = if (length(found)) found[[1L]])
}

#' Default (intercept-only) or constant dpar spec.
#'
#' @noRd
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

#' One bf() -> one response entry of the spec.
#'
#' Takes one `frmtmb_formula` (a formula plus its family, dpar formulas
#' and fixed dpar values) and returns the response entry the frame
#' builds from: `resp_name`, `resp_expr`, `family`, the addition terms
#' `aterms`, one parsed linear predictor per `dpars` entry, the primary
#' (location) dpar names, the nonlinear parameter names, and the shared
#' `formula_env`. It runs `parse_linpred()` once per dpar, gives every
#' dpar without a formula an intercept or a constant, and handles the
#' nonlinear (`nl = TRUE`) case where mu is a body expression instead of
#' a design. It checks the formula, never the data.
#'
#' @noRd
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
  ri <- rewrite_cbind_response(parse_response(f), fam)

  if (isTRUE(bform$nl)) {
    if (!identical(fam$primary_dpars %||% "mu", "mu")) {
      stop("nl = TRUE requires a family with a single 'mu' location ",
           "parameter", call. = FALSE)
    }
    body <- reformulas::RHSForm(f)
    # A nonlinear mu is arbitrary R code, so an ar() written there is
    # EVALUATED, not parsed, and fails deep inside the objective with a
    # message about the body. Say what is wrong instead.
    ac_in_body <- intersect(autocor_structs, all.names(body))
    if (length(ac_in_body)) {
      stop("Residual correlation terms are not supported in a ",
           "nonlinear (nl = TRUE) formula; '", ac_in_body[1L],
           "()' appears in the model body, where it would be evaluated ",
           "as ordinary R code rather than read as a term. brms reaches ",
           "the same model through acformula(), which has no analog ",
           "here", call. = FALSE)
    }
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
    pa <- pull_autocor(dpars, deparse1(ri$resp))
    return(list(
      resp_name = deparse1(ri$resp),
      resp_expr = ri$resp,
      family = fam,
      aterms = ri$aterms,
      dpars = pa$dpars,
      autocor = pa$autocor,
      primary_dpars = nlpars,   # REML integrates the nlpar coefficients
      nlpars = nlpars,
      cbind_resp = isTRUE(ri$cbind_resp),
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

  # A family may ship a DEFAULT formula for some of its own dpars, which
  # stands in wherever the user wrote neither a formula nor a fixed
  # value: hmm(trans = ~x) gives every transition cell the same
  # predictor without making the user spell out K * (K - 1) of them. It
  # is a default, so an explicit `tr12 ~ z` or `tr12 = 0` still wins.
  for (dp in names(fam[["default_forms"]] %||% list())) {
    if (dp %in% fam$dpars && !dp %in% extra) {
      bform$pforms[[dp]] <- fam[["default_forms"]][[dp]]
    }
  }

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

  # A mixture family has no 'mu': its location dpars are mu1, mu2, ...,
  # and main_lp is copied into each. pull_autocor() therefore refuses a
  # residual correlation term on a mixture, which is what brms does and
  # for the same reason - a mixture likelihood has no single residual.
  pa <- pull_autocor(dpars, deparse1(ri$resp))

  list(
    resp_name = deparse1(ri$resp),
    resp_expr = ri$resp,
    family = fam,
    aterms = ri$aterms,
    dpars = pa$dpars,
    autocor = pa$autocor,
    primary_dpars = primaries,
    nlpars = character(0),
    cbind_resp = isTRUE(ri$cbind_resp),
    formula_env = shared_env
  )
}

#' frmtmb_formula / frmtmb_mvformula -> frmtmb_spec.
#'
#' The entry point of the parsing stage. It calls `parse_one_response()`
#' for each response of a univariate or multivariate formula and returns
#' a `frmtmb_spec`: the named `responses` list plus the `rescor` flag.
#' The spec is pure formula information, so `frm()` can inspect it
#' before any data arrives; `assemble_frame(spec, data)` is the next
#' step and produces the design matrices.
#'
#' Spec-level consistency of the |ID| keys.
#'
#' Two rules, both properties of the whole spec rather than of one
#' linear predictor, which is why they live here and not in
#' `parse_linpred()`.
#'
#' 1. One `|ID|` label, one grouping specification. The merge key is the
#'    label PLUS the deparsed grouping expression, so terms that write
#'    the same label over different grouping expressions land in
#'    different keys and quietly fail to correlate at all - the user
#'    asked for a link and got none. Refusing is the only honest answer:
#'    `(1 | q | g1)` with `(1 | q | g2)`, or `(1 | q | gr(g, cov = A))`
#'    with `(1 | q | gr(g, cov = B))`, name one link over two different
#'    structures.
#'
#' 2. One `|ID|` key, one covariance structure. The guard inside
#'    `parse_linpred()` runs before the `gr()` rewrite turns `cls` from
#'    `"us"` into `"gr_cov"`/`"gr_prec"`, so those two are the only
#'    structures that reach a shared key. Since v0.32 a shared key whose
#'    terms are ALL `gr(cov = )` (or all `gr(prec = )`) over the same
#'    grouping factor and the same matrix is supported: the merged block
#'    is built as one `gr_cov`/`gr_prec` block of the total merged
#'    dimension, so its covariance is `A (x) Sigma` with `Sigma` the
#'    unstructured covariance across the merged coefficients - the same
#'    joint density as the long-format spelling. Rule 1 already forces
#'    the grouping expressions to agree; this rule catches the residue
#'    (a key whose terms disagree on the structure itself), and frame
#'    assembly re-checks that the two `cov =` expressions RESOLVE to the
#'    same matrix, which formula environments can make them not do.
#'
#' A key used by a single term is a no-op either way: it takes the
#' length-1 branch of the phase-2 loop and keeps its own structure.
#'
#' @noRd
check_id_covstructs <- function(spec) {
  ids <- character(0)
  id_labels <- character(0)
  groups <- character(0)
  cls <- character(0)
  labs <- character(0)
  for (resp in spec$responses) {
    for (dp in resp$dpars) {
      for (z in dp$re %||% list()) {
        if (is.null(z$id)) next
        ids <- c(ids, z$id)
        id_labels <- c(id_labels, z$id_label)
        groups <- c(groups, deparse1(z$id_group))
        cls <- c(cls, z$covstruct)
        labs <- c(labs, paste0(resp$resp_name, " ", dp$name, ": (",
                               deparse1(z$bar), ") [", z$covstruct, "]"))
      }
    }
  }
  if (!length(ids)) return(invisible(NULL))
  for (lb in unique(id_labels)) {
    at <- which(id_labels == lb)
    gs <- unique(groups[at])
    if (length(gs) > 1L) {
      stop("The |", lb, "| key is used over more than one grouping ",
           "specification (", paste(gs, collapse = ", "), "): ",
           paste(labs[at], collapse = "; "),
           ". Terms keyed by the same |ID| merge into one covariance ",
           "block, which needs a single grouping factor and, for ",
           "gr(cov = ) / gr(prec = ), a single relationship matrix. ",
           "Use one spelling for all of them, or give the terms ",
           "different |ID| labels. Both spellings of a multi-trait ",
           "model over ONE relationship matrix are supported and give ",
           "the same fit:\n",
           "  mvbf(bf(y1 ~ (1 | q | gr(id, cov = A))), ",
           "bf(y2 ~ (1 | q | gr(id, cov = A))))\n",
           "  bf(value ~ 0 + trait + (0 + trait | gr(id, cov = A)), ",
           "sigma ~ 0 + trait)", call. = FALSE)
    }
  }
  shared <- ids %in% ids[duplicated(ids)]
  for (k in unique(ids[shared])) {
    at <- which(ids == k)
    cs <- unique(cls[at])
    if (length(cs) == 1L) next
    stop("Terms sharing the |ID| key '", k, "' mix covariance ",
         "structures (", paste(cs, collapse = ", "), "): ",
         paste(labs[at], collapse = "; "),
         ". A merged block has one structure. Either give them all the ",
         "same structure, or split the |ID| label. A key whose terms ",
         "are all gr(cov = ) (or all gr(prec = )) over the same factor ",
         "and the same matrix is supported and fits the same model as ",
         "the long-format spelling, ",
         "bf(value ~ 0 + trait + (0 + trait | gr(id, cov = A)), ",
         "sigma ~ 0 + trait).", call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
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
    out <- structure(
      list(responses = resps, rescor = rescor),
      class = "frmtmb_spec"
    )
    check_id_covstructs(out)
    return(out)
  }
  stopifnot(inherits(bform, "frmtmb_formula"))
  resp <- parse_one_response(bform)
  out <- structure(
    list(responses = stats::setNames(list(resp), resp$resp_name),
         rescor = FALSE),
    class = "frmtmb_spec"
  )
  check_id_covstructs(out)
  out
}

#' @export
print.frmtmb_spec <- function(x, ...) {
  cat("<frmtmb spec>\n")
  for (r in x$responses) {
    cat("Response: ", r$resp_name, "  [", r$family$family, "]\n", sep = "")
    if (!is.null(r$autocor)) {
      cat("  residual correlation: ", r$autocor$label, "\n", sep = "")
    }
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
