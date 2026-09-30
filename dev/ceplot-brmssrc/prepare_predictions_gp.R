prepare_predictions_gp <- 
function (bframe, draws, sdata, new = FALSE, nug = NULL, ...) 
{
    stopifnot(is.bframel(bframe))
    gpframe <- bframe$frame$gp
    if (!has_rows(gpframe)) {
        return(list())
    }
    p <- usc(combine_prefix(bframe))
    if (is.null(nug)) {
        nug <- ifelse(new, 1e-08, 1e-12)
    }
    out <- named_list(gpframe$label)
    for (i in seq_along(out)) {
        cons <- gpframe$cons[[i]]
        if (length(cons)) {
            gp <- named_list(cons)
            for (j in seq_along(cons)) {
                gp[[j]] <- .prepare_predictions_gp(gpframe, draws = draws, 
                  sdata = sdata, nug = nug, new = new, byj = j, 
                  p = p, i = i)
            }
            attr(gp, "byfac") <- TRUE
        }
        else {
            gp <- .prepare_predictions_gp(gpframe, draws = draws, 
                sdata = sdata, nug = nug, new = new, p = p, i = i)
        }
        out[[i]] <- gp
    }
    out
}
