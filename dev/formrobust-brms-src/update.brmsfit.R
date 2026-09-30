update.brmsfit <- 
function (object, formula., newdata = NULL, recompile = NULL, 
    ...) 
{
    dots <- list(...)
    testmode <- isTRUE(dots[["testmode"]])
    dots$testmode <- NULL
    if ("silent" %in% names(dots)) {
        dots$silent <- validate_silent(dots$silent)
    }
    else {
        dots$silent <- object$stan_args$silent %||% 1
    }
    silent <- dots$silent
    object <- restructure(object)
    if (isTRUE(object$version$brms < "2.0.0")) {
        warning2("Updating models fitted with older versions of brms may fail.")
    }
    object$file <- NULL
    if ("data" %in% names(dots)) {
        stop2("Please use argument 'newdata' to update the data.")
    }
    if (!is.null(newdata)) {
        dots$data <- newdata
        data_name <- substitute_name(newdata)
    }
    else {
        dots$data <- object$data
        data_name <- get_data_name(object$data)
    }
    if (missing(formula.) || is.null(formula.)) {
        dots$formula <- object$formula
        if (!is.null(dots[["family"]])) {
            dots$formula <- bf(dots$formula, family = dots$family)
        }
        if (!is.null(dots[["autocor"]])) {
            dots$formula <- bf(dots$formula, autocor = dots$autocor)
        }
    }
    else {
        if (is.mvbrmsformula(formula.) || is.mvbrmsformula(object$formula)) {
            stop2("Updating formulas of multivariate models is not yet possible.")
        }
        if (is.brmsformula(formula.)) {
            nl <- get_nl(formula.)
        }
        else {
            formula. <- as.formula(formula.)
            nl <- get_nl(formula(object))
        }
        family <- get_arg("family", formula., dots, object)
        autocor <- get_arg("autocor", formula., dots, object)
        dots$formula <- bf(formula., family = family, autocor = autocor, 
            nl = nl)
        if (is_nonlinear(object)) {
            if (length(setdiff(all.vars(dots$formula$formula), 
                ".")) == 0) {
                dots$formula <- update(object$formula, dots$formula, 
                  mode = "keep")
            }
            else {
                dots$formula <- update(object$formula, dots$formula, 
                  mode = "replace")
                if (silent < 2) {
                  message("Argument 'formula.' will completely replace the ", 
                    "original formula in non-linear models.")
                }
            }
        }
        else {
            mvars <- all.vars(dots$formula$formula)
            mvars <- setdiff(mvars, c(names(object$data), "."))
            if (length(mvars) && is.null(newdata)) {
                stop2("New variables found: ", collapse_comma(mvars), 
                  "\nPlease supply your data again via argument 'newdata'.")
            }
            dots$formula <- update(formula(object), dots$formula)
        }
    }
    dots$formula <- validate_formula(dots$formula, data = dots$data)
    if (is.null(dots$prior)) {
        dots$prior <- object$prior
    }
    else {
        if (!is.brmsprior(dots$prior)) {
            stop2("Argument 'prior' needs to be a 'brmsprior' object.")
        }
        old_user_prior <- subset2(object$prior, source = "user")
        dots$prior <- rbind(dots$prior, old_user_prior)
        dupl_priors <- duplicated(dots$prior[, rcols_prior()])
        dots$prior <- dots$prior[!dupl_priors, ]
    }
    attr(dots$prior, "allow_invalid_prior") <- TRUE
    if (!"sample_prior" %in% names(dots)) {
        dots$sample_prior <- attr(object$prior, "sample_prior")
        if (is.null(dots$sample_prior)) {
            has_prior_pars <- any(grepl("^prior_", variables(object)))
            dots$sample_prior <- if (has_prior_pars) 
                "yes"
            else "no"
        }
    }
    if (!"data2" %in% names(dots)) {
        dots$data2 <- object$data2
    }
    if (!"stanvars" %in% names(dots)) {
        dots$stanvars <- object$stanvars
    }
    if (!"algorithm" %in% names(dots)) {
        dots$algorithm <- object$algorithm
    }
    if (!"backend" %in% names(dots)) {
        dots$backend <- object$backend
    }
    if (!"threads" %in% names(dots)) {
        dots$threads <- object$threads
    }
    if (!"save_pars" %in% names(dots)) {
        dots$save_pars <- object$save_pars
    }
    if (!"knots" %in% names(dots)) {
        dots$knots <- get_knots(object$data)
    }
    if (!"drop_unused_levels" %in% names(dots)) {
        dots$drop_unused_levels <- get_drop_unused_levels(object$data)
    }
    if (!"normalize" %in% names(dots)) {
        dots$normalize <- is_normalized(object$model)
    }
    dots$algorithm <- match.arg(dots$algorithm, algorithm_choices())
    dots$backend <- match.arg(dots$backend, backend_choices())
    same_algorithm <- is_equal(dots$algorithm, object$algorithm)
    same_backend <- is_equal(dots$backend, object$backend)
    if (same_algorithm) {
        if (is.null(dots$iter)) {
            dots$warmup <- first_not_null(dots$warmup, object$fit@sim$warmup)
        }
        dots$iter <- first_not_null(dots$iter, object$fit@sim$iter)
        dots$chains <- first_not_null(dots$chains, object$fit@sim$chains)
        dots$thin <- first_not_null(dots$thin, object$fit@sim$thin)
        if (same_backend) {
            control <- attr(object$fit@sim$samples[[1]], "args")$control
            control <- control[setdiff(names(control), names(dots$control))]
            dots$control[names(control)] <- control
            names_old_stan_args <- setdiff(names(object$stan_args), 
                names(dots))
            dots[names_old_stan_args] <- object$stan_args[names_old_stan_args]
        }
    }
    if (is.null(recompile)) {
        new_stancode <- suppressMessages(do_call(make_stancode, 
            dots))
        new_stancode <- sub("^[^\n]+\n", "", new_stancode)
        old_stancode <- stancode(object, version = FALSE)
        recompile <- needs_recompilation(object) || !same_backend || 
            !is_equal(new_stancode, old_stancode)
        if (recompile && silent < 2) {
            message("The desired updates require recompiling the model")
        }
    }
    recompile <- as_one_logical(recompile)
    if (recompile) {
        dots$fit <- NA
        if (!testmode) {
            object <- do_call(brm, dots)
        }
    }
    else {
        if (!is.null(dots$formula)) {
            object$formula <- dots$formula
            dots$formula <- NULL
        }
        bterms <- brmsterms(object$formula)
        object$data2 <- validate_data2(dots$data2, bterms = bterms)
        object$data <- validate_data(dots$data, bterms = bterms, 
            data2 = object$data2, knots = dots$knots, drop_unused_levels = dots$drop_unused_levels)
        bframe <- brmsframe(bterms, data = object$data)
        object$prior <- .validate_prior(dots$prior, bframe = bframe, 
            sample_prior = dots$sample_prior)
        object$family <- get_element(object$formula, "family")
        object$autocor <- get_element(object$formula, "autocor")
        object$ranef <- frame_re(bterms, data = object$data)
        object$stanvars <- validate_stanvars(dots$stanvars)
        object$threads <- validate_threads(dots$threads)
        if ("sample_prior" %in% names(dots)) {
            dots$sample_prior <- validate_sample_prior(dots$sample_prior)
            attr(object$prior, "sample_prior") <- dots$sample_prior
        }
        object$save_pars <- validate_save_pars(save_pars = dots$save_pars, 
            save_ranef = dots$save_ranef, save_mevars = dots$save_mevars, 
            save_all_pars = dots$save_all_pars)
        object$basis <- frame_basis(bframe, data = object$data)
        algorithm <- match.arg(dots$algorithm, algorithm_choices())
        dots$algorithm <- object$algorithm <- algorithm
        dots$backend <- object$backend
        if (!testmode) {
            dots$fit <- object
            object <- do_call(brm, dots)
        }
    }
    attr(object$data, "data_name") <- data_name
    object
}
