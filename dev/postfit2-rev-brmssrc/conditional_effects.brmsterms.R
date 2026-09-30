function (x, fit, cond_data, int_conditions, method, surface, 
    spaghetti, categorical, ordinal, probs, robust, dpar = NULL, 
    nlpar = NULL, resp = NULL, ...) 
{
    stopifnot(is.brmsfit(fit))
    effects <- attr(cond_data, "effects")
    types <- attr(cond_data, "types")
    catscale <- NULL
    pred_args <- list(fit, newdata = cond_data, allow_new_levels = TRUE, 
        dpar = dpar, nlpar = nlpar, resp = if (nzchar(x$resp)) x$resp, 
        incl_autocor = FALSE, ...)
    if (method != "posterior_predict") {
        pred_args$transform <- NULL
    }
    out <- do_call(method, pred_args)
    rownames(cond_data) <- NULL
    if (categorical || ordinal) {
        if (method != "posterior_epred") {
            stop2("Can only use 'categorical' with method = 'posterior_epred'.")
        }
        if (!is_polytomous(x)) {
            stop2("Argument 'categorical' may only be used ", 
                "for categorical or ordinal models.")
        }
        if (categorical && ordinal) {
            stop2("Please use argument 'categorical' instead of 'ordinal'.")
        }
        catscale <- str_if(is_multinomial(x), "Count", "Probability")
        cats <- dimnames(out)[[3]]
        if (is.null(cats)) 
            cats <- seq_dim(out, 3)
        cond_data <- repl(cond_data, length(cats))
        cond_data <- do_call(rbind, cond_data)
        cond_data$cats__ <- factor(rep(cats, each = ncol(out)), 
            levels = cats)
        effects[2] <- "cats__"
        types[2] <- "factor"
    }
    else {
        if (conv_cats_dpars(x$family) && is.null(dpar)) {
            stop2("Please set 'categorical' to TRUE.")
        }
        if (is_ordinal(x$family) && is.null(dpar) && method != 
            "posterior_linpred") {
            warning2("Predictions are treated as continuous variables in ", 
                "'conditional_effects' by default which is likely invalid ", 
                "for ordinal families. Please set 'categorical' to TRUE.")
            if (method == "posterior_epred") {
                out <- ordinal_probs_continuous(out)
            }
        }
    }
    cond_data <- add_effects__(cond_data, effects)
    first_numeric <- types[1] %in% "numeric"
    second_numeric <- types[2] %in% "numeric"
    both_numeric <- first_numeric && second_numeric
    if (second_numeric && !surface) {
        mde2 <- round(cond_data[[effects[2]]], 2)
        levels2 <- sort(unique(mde2), TRUE)
        cond_data$effect2__ <- factor(mde2, levels = levels2)
        labels2 <- names(int_conditions[[effects[2]]])
        if (length(labels2) == length(levels2)) {
            levels(cond_data$effect2__) <- labels2
        }
    }
    spag <- NULL
    if (first_numeric && spaghetti) {
        if (surface) {
            stop2("Cannot use 'spaghetti' and 'surface' at the same time.")
        }
        spag <- out
        if (categorical) {
            spag <- do_call(cbind, array2list(spag))
        }
        sample <- rep(seq_rows(spag), each = ncol(spag))
        if (length(types) == 2L) {
            sample <- paste0(sample, "_", cond_data[[effects[2]]])
        }
        spag <- data.frame(as.numeric(t(spag)), factor(sample))
        colnames(spag) <- c("estimate__", "sample__")
        cond_data_spag <- repl(cond_data, nrow(spag)/nrow(cond_data))
        cond_data_spag <- Reduce(rbind, cond_data_spag)
        spag <- cbind(cond_data_spag, spag)
    }
    out <- posterior_summary(out, probs = probs, robust = robust)
    if (categorical || ordinal) {
        out <- do_call(rbind, array2list(out))
    }
    colnames(out) <- c("estimate__", "se__", "lower__", "upper__")
    out <- cbind(cond_data, out)
    if (!is.null(dpar)) {
        response <- dpar
    }
    else if (!is.null(nlpar)) {
        response <- nlpar
    }
    else {
        response <- as.character(x$formula[2])
    }
    attr(out, "effects") <- effects
    attr(out, "response") <- response
    attr(out, "surface") <- unname(both_numeric && surface)
    attr(out, "categorical") <- categorical
    attr(out, "catscale") <- catscale
    attr(out, "ordinal") <- ordinal
    attr(out, "spaghetti") <- spag
    attr(out, "points") <- make_point_frame(x, fit$data, effects, 
        ...)
    name <- paste0(usc(x$resp, "suffix"), paste0(effects, collapse = ":"))
    setNames(list(out), name)
}
