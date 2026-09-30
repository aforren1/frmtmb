prepare_predictions.mvbrmsframe <- 
function (x, draws, sdata, resp = NULL, ...) 
{
    resp <- validate_resp(resp, x$responses)
    if (length(resp) > 1) {
        if (has_subset(x)) {
            stop2("Argument 'resp' must be a single variable name ", 
                "for models using addition argument 'subset'.")
        }
        out <- list(ndraws = nrow(draws), nobs = sdata$N)
        out$resps <- named_list(resp)
        out$old_order <- attr(sdata, "old_order")
        for (r in resp) {
            out$resps[[r]] <- prepare_predictions(x$terms[[r]], 
                draws = draws, sdata = sdata, ...)
        }
        if (x$rescor) {
            out$family <- out$resps[[1]]$family
            out$family$fun <- paste0(out$family$family, "_mv")
            rescor <- get_cornames(resp, type = "rescor", brackets = FALSE)
            out$mvpars$rescor <- prepare_draws(draws, rescor)
            if (out$family$family == "student") {
                out$dpars$nu <- as.vector(prepare_draws(draws, 
                  "nu"))
            }
            out$data$N <- out$resps[[1]]$data$N
            out$data$weights <- out$resps[[1]]$data$weights
            Y <- lapply(out$resps, function(x) x$data$Y)
            out$data$Y <- do_call(cbind, Y)
        }
        out <- structure(out, class = "mvbrmsprep")
    }
    else {
        out <- prepare_predictions(x$terms[[resp]], draws = draws, 
            sdata = sdata, ...)
    }
    out
}
