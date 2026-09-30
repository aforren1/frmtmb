brmsformula <- 
function (formula, ..., flist = NULL, family = NULL, autocor = NULL, 
    nl = NULL, loop = NULL, center = NULL, cmc = NULL, sparse = NULL, 
    decomp = NULL, unused = NULL) 
{
    if (is.brmsformula(formula)) {
        out <- formula
    }
    else {
        out <- list(formula = as_formula(formula))
        class(out) <- "brmsformula"
    }
    dots <- c(out$pforms, out$pfix, list(...), flist)
    dots <- lapply(dots, function(x) if (is.list(x)) 
        x
    else list(x))
    dots <- unlist(dots, recursive = FALSE)
    forms <- list()
    for (i in seq_along(dots)) {
        c(forms) <- validate_par_formula(dots[[i]], par = names(dots)[i])
    }
    is_dupl_pars <- duplicated(names(forms), fromLast = TRUE)
    if (any(is_dupl_pars)) {
        dupl_pars <- collapse_comma(names(forms)[is_dupl_pars])
        message("Replacing initial definitions of parameters ", 
            dupl_pars)
        forms[is_dupl_pars] <- NULL
    }
    not_form <- ulapply(forms, function(x) !is.formula(x))
    fix <- forms[not_form]
    forms[names(fix)] <- NULL
    out$pforms <- forms
    fix_theta <- fix[dpar_class(names(fix)) %in% "theta"]
    if (length(fix_theta)) {
        sum_theta <- sum(unlist(fix_theta))
        fix_theta <- lapply(fix_theta, "/", sum_theta)
        fix[names(fix_theta)] <- fix_theta
    }
    out$pfix <- fix
    for (dp in names(out$pfix)) {
        if (is.character(out$pfix[[dp]])) {
            if (identical(dp, out$pfix[[dp]])) {
                stop2("Equating '", dp, "' with itself is not meaningful.")
            }
            ap_class <- dpar_class(dp)
            if (ap_class == "mu") {
                stop2("Equating parameters of class 'mu' is not allowed.")
            }
            if (!identical(ap_class, dpar_class(out$pfix[[dp]]))) {
                stop2("Can only equate parameters of the same class.")
            }
            if (out$pfix[[dp]] %in% names(out$pfix)) {
                stop2("Cannot use fixed parameters on ", "the right-hand side of an equation.")
            }
            if (out$pfix[[dp]] %in% names(out$pforms)) {
                stop2("Cannot use predicted parameters on ", 
                  "the right-hand side of an equation.")
            }
        }
    }
    if (!is.null(nl)) {
        attr(out$formula, "nl") <- as_one_logical(nl)
    }
    else if (!is.null(out[["nl"]])) {
        attr(out$formula, "nl") <- out[["nl"]]
        out[["nl"]] <- NULL
    }
    if (is.null(attr(out$formula, "nl"))) {
        attr(out$formula, "nl") <- FALSE
    }
    if (!is.null(loop)) {
        attr(out$formula, "loop") <- as_one_logical(loop)
    }
    if (is.null(attr(out$formula, "loop"))) {
        attr(out$formula, "loop") <- TRUE
    }
    if (!is.null(center)) {
        attr(out$formula, "center") <- as_one_logical(center)
    }
    if (!is.null(cmc)) {
        attr(out$formula, "cmc") <- as_one_logical(cmc)
    }
    if (!is.null(sparse)) {
        attr(out$formula, "sparse") <- as_one_logical(sparse)
    }
    if (!is.null(decomp)) {
        attr(out$formula, "decomp") <- match.arg(decomp, decomp_opts())
    }
    if (!is.null(unused)) {
        attr(out$formula, "unused") <- as.formula(unused)
    }
    if (!is.null(autocor)) {
        attr(out$formula, "autocor") <- validate_autocor(autocor)
    }
    else if (!is.null(out$autocor)) {
        attr(out$formula, "autocor") <- validate_autocor(out$autocor)
        out$autocor <- NULL
    }
    if (!is.null(family)) {
        out$family <- validate_family(family)
    }
    if (!is.null(lhs(formula))) {
        out$resp <- terms_resp(formula)
    }
    defs <- list(pforms = list(), pfix = list(), family = NULL, 
        resp = NULL)
    defs <- defs[setdiff(names(defs), names(rmNULL(out, FALSE)))]
    out[names(defs)] <- defs
    class(out) <- c("brmsformula", "bform")
    split_bf(out)
}
