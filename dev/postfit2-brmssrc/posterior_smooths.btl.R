posterior_smooths.btl <- 
function (object, fit, smooth, newdata = NULL, ndraws = NULL, 
    draw_ids = NULL, nsamples = NULL, subset = NULL, ...) 
{
    smooth <- rm_wsp(as_one_character(smooth))
    ndraws <- use_alias(ndraws, nsamples)
    draw_ids <- use_alias(draw_ids, subset)
    object$frame$sm <- frame_sm(object, fit$data)
    class(object) <- c("bframel", class(object))
    smframe <- object$frame$sm
    smframe$term <- rm_wsp(smframe$term)
    smterms <- unique(smframe$term)
    if (!smooth %in% smterms) {
        stop2("Term '", smooth, "' cannot be found. Available ", 
            "smooth terms are: ", collapse_comma(smterms))
    }
    sub_smframe <- subset2(smframe, term = smooth)
    covars <- all_vars(sub_smframe$covars[[1]])
    byvars <- all_vars(sub_smframe$byvars[[1]])
    req_vars <- c(covars, byvars)
    sdata <- standata(fit, newdata, re_formula = NA, internal = TRUE, 
        check_response = FALSE, req_vars = req_vars)
    draw_ids <- validate_draw_ids(fit, draw_ids, ndraws)
    draws <- as_draws_matrix(fit)
    draws <- suppressMessages(subset_draws(draws, draw = draw_ids))
    prep_args <- nlist(x = object, draws, sdata, data = fit$data)
    prep <- do_call(prepare_predictions, prep_args)
    i <- which(smterms %in% smooth)[1]
    J <- which(smframe$termnum == i)
    scs <- unlist(attr(prep$sm$fe$Xs, "smcols")[J])
    prep$sm$fe$Xs <- prep$sm$fe$Xs[, scs, drop = FALSE]
    prep$sm$fe$bs <- prep$sm$fe$bs[, scs, drop = FALSE]
    prep$sm$re <- prep$sm$re[J]
    prep$family <- brmsfamily("gaussian")
    predictor(prep, i = NULL)
}
