prepare_predictions_thres <- 
function (bframe, draws, sdata, ...) 
{
    out <- list()
    if (!is_ordinal(bframe$family)) {
        return(out)
    }
    resp <- usc(bframe$resp)
    out$nthres <- sdata[[paste0("nthres", resp)]]
    out$Jthres <- sdata[[paste0("Jthres", resp)]]
    p <- usc(combine_prefix(bframe))
    thres_regex <- paste0("^b", p, "_Intercept\\[")
    out$thres <- prepare_draws(draws, thres_regex, regex = TRUE)
    out
}
