loo_model_weights.brmsfit <- 
function (x, ..., model_names = NULL) 
{
    args <- split_dots(x, ..., model_names = model_names)
    models <- args$models
    args$models <- NULL
    log_lik_list <- lapply(models, function(x) do_call(log_lik, 
        c(list(x), args)))
    args$x <- log_lik_list
    args$r_eff_list <- mapply(r_eff_log_lik, log_lik_list, fit = models, 
        SIMPLIFY = FALSE)
    out <- do_call(loo::loo_model_weights, args)
    names(out) <- names(models)
    out
}
