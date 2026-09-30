validate_weights <- 
function (weights, models, control = list()) 
{
    if (!is.numeric(weights)) {
        weight_args <- c(unname(models), control)
        weight_args$weights <- weights
        weights <- do_call(model_weights, weight_args)
    }
    else {
        if (length(weights) != length(models)) {
            stop2("If numeric, 'weights' must have the same length ", 
                "as the number of models.")
        }
        if (any(weights < 0)) {
            stop2("If numeric, 'weights' must be positive.")
        }
    }
    weights/sum(weights)
}
