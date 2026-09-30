posterior_average.brmsfit <- 
function (x, ..., variable = NULL, pars = NULL, weights = "stacking", 
    ndraws = NULL, nsamples = NULL, missing = NULL, model_names = NULL, 
    control = list(), seed = NULL) 
{
    if (!is.null(seed)) {
        set.seed(seed)
    }
    variable <- use_alias(variable, pars)
    ndraws <- use_alias(ndraws, nsamples)
    models <- split_dots(x, ..., model_names = model_names, other = FALSE)
    vars_list <- lapply(models, variables)
    all_vars <- unique(unlist(vars_list))
    if (is.null(missing)) {
        common_vars <- lapply(vars_list, function(x) all_vars %in% 
            x)
        common_vars <- all_vars[Reduce("&", common_vars)]
        if (is.null(variable)) {
            variable <- setdiff(common_vars, "lp__")
        }
        variable <- as.character(variable)
        inv_vars <- setdiff(variable, common_vars)
        if (length(inv_vars)) {
            inv_vars <- collapse_comma(inv_vars)
            stop2("Parameters ", inv_vars, " cannot be found in all ", 
                "of the models. Consider using argument 'missing'.")
        }
    }
    else {
        if (is.null(variable)) {
            variable <- setdiff(all_vars, "lp__")
        }
        variable <- as.character(variable)
        inv_vars <- setdiff(variable, all_vars)
        if (length(inv_vars)) {
            inv_vars <- collapse_comma(inv_vars)
            stop2("Parameters ", inv_vars, " cannot be found in any of the models.")
        }
        if (is.list(missing)) {
            all_miss_vars <- unique(ulapply(models, function(m) setdiff(variable, 
                variables(m))))
            inv_vars <- setdiff(all_miss_vars, names(missing))
            if (length(inv_vars)) {
                stop2("Argument 'missing' has no value for parameters ", 
                  collapse_comma(inv_vars), ".")
            }
            missing <- lapply(missing, as_one_numeric, allow_na = TRUE)
        }
        else {
            missing <- as_one_numeric(missing, allow_na = TRUE)
            missing <- named_list(variable, missing)
        }
    }
    if (is.null(ndraws)) {
        ndraws <- ndraws(models[[1]])
    }
    ndraws <- as_one_integer(ndraws)
    weights <- validate_weights(weights, models, control)
    ndraws <- round_largest_remainder(weights * ndraws)
    names(weights) <- names(ndraws) <- names(models)
    out <- named_list(names(models))
    for (i in seq_along(out)) {
        if (ndraws[i] > 0) {
            draw <- sample(seq_len(ndraws(models[[i]])), ndraws[i])
            draw <- sort(draw)
            found_vars <- intersect(variable, variables(models[[i]]))
            if (length(found_vars)) {
                out[[i]] <- as.data.frame(models[[i]], variable = found_vars, 
                  draw = draw)
            }
            else {
                out[[i]] <- as.data.frame(matrix(numeric(0), 
                  nrow = ndraws[i], ncol = 0))
            }
            if (!is.null(missing)) {
                miss_vars <- setdiff(variable, names(out[[i]]))
                if (length(miss_vars)) {
                  out[[i]][miss_vars] <- missing[miss_vars]
                }
            }
        }
    }
    out <- do_call(rbind, out)
    rownames(out) <- NULL
    attr(out, "weights") <- weights
    attr(out, "ndraws") <- ndraws
    out
}
