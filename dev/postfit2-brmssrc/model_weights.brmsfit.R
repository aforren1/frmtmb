model_weights.brmsfit <- 
function (x, ..., weights = "stacking", model_names = NULL) 
{
    weights <- validate_weights_method(weights)
    args <- split_dots(x, ..., model_names = model_names)
    models <- args$models
    args$models <- NULL
    model_names <- names(models)
    if (weights %in% c("loo", "waic", "kfold")) {
        ics <- rep(NA, length(models))
        for (i in seq_along(ics)) {
            args$x <- models[[i]]
            args$model_names <- names(models)[i]
            ics[i] <- SW(do_call(weights, args))$estimates[3, 
                1]
        }
        ic_diffs <- ics - min(ics)
        out <- exp(-ic_diffs/2)
    }
    else if (weights %in% c("stacking", "pseudobma")) {
        args <- c(unname(models), args)
        args$method <- weights
        out <- do_call("loo_model_weights", args)
    }
    else if (weights %in% "bma") {
        args <- c(unname(models), args)
        out <- do_call("post_prob", args)
    }
    out <- as.numeric(out)
    out <- out/sum(out)
    names(out) <- model_names
    out
}
