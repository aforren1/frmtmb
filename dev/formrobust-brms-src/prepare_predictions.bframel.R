prepare_predictions.bframel <- 
function (x, draws, sdata, ...) 
{
    ndraws <- nrow(draws)
    nobs <- sdata[[paste0("N", usc(x$resp))]]
    out <- nlist(family = x$family, ndraws, nobs)
    class(out) <- "bprepl"
    out$fe <- prepare_predictions_fe(x, draws, sdata, ...)
    out$sp <- prepare_predictions_sp(x, draws, sdata, ...)
    out$cs <- prepare_predictions_cs(x, draws, sdata, ...)
    out$sm <- prepare_predictions_sm(x, draws, sdata, ...)
    out$gp <- prepare_predictions_gp(x, draws, sdata, ...)
    out$re <- prepare_predictions_re(x, sdata, ...)
    out$ac <- prepare_predictions_ac(x, draws, sdata, nat_cov = FALSE, 
        ...)
    out$offset <- prepare_predictions_offset(x, sdata, ...)
    out
}
