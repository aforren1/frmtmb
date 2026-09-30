pp_average.brmsfit <- 
function (x, ..., weights = "stacking", method = "posterior_predict", 
    ndraws = NULL, nsamples = NULL, summary = TRUE, probs = c(0.025, 
        0.975), robust = FALSE, model_names = NULL, control = list(), 
    seed = NULL) 
{
    if (!is.null(seed)) {
        set.seed(seed)
    }
    method <- validate_pp_method(method)
    ndraws <- use_alias(ndraws, nsamples)
    if (any(c("draw_ids", "subset") %in% names(list(...)))) {
        stop2("Cannot use argument 'draw_ids' in pp_average.")
    }
    args <- split_dots(x, ..., model_names = model_names)
    args$summary <- FALSE
    models <- args$models
    args$models <- NULL
    if (!match_response(models)) {
        stop2("Can only average models predicting the same response.")
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
            args$object <- models[[i]]
            args$ndraws <- ndraws[i]
            out[[i]] <- do_call(method, args)
        }
    }
    out <- do_call(rbind, out)
    if (summary) {
        out <- posterior_summary(out, probs = probs, robust = robust)
    }
    attr(out, "weights") <- weights
    attr(out, "ndraws") <- ndraws
    out
}
