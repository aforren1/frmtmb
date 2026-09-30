prepare_predictions_cs <- 
function (bframe, draws, sdata, ...) 
{
    stopifnot(is.bframel(bframe))
    out <- list()
    if (!is_ordinal(bframe$family)) {
        return(out)
    }
    resp <- usc(bframe$resp)
    out$nthres <- sdata[[paste0("nthres", resp)]]
    csef <- bframe$frame$cs$vars
    if (length(csef)) {
        p <- usc(combine_prefix(bframe))
        cs_pars <- paste0("^bcs", p, "_", escape_all(csef), "\\[")
        out$bcs <- prepare_draws(draws, cs_pars, regex = TRUE)
        out$Xcs <- sdata[[paste0("Xcs", p)]]
    }
    out
}
