prepare_predictions_sm <- 
function (bframe, draws, sdata, ...) 
{
    stopifnot(is.bframel(bframe))
    out <- list()
    smframe <- bframe$frame$sm
    if (!has_rows(smframe)) {
        return(out)
    }
    p <- usc(combine_prefix(bframe))
    Xs_names <- attr(smframe, "Xs_names")
    if (length(Xs_names)) {
        out$fe$Xs <- sdata[[paste0("Xs", p)]]
        bspars <- paste0("^bs?", p, "_", escape_all(Xs_names), 
            "$")
        out$fe$bs <- prepare_draws(draws, bspars, regex = TRUE)
    }
    out$re <- named_list(smframe$label)
    for (i in seq_rows(smframe)) {
        sm <- list()
        for (j in seq_len(smframe$nbases[i])) {
            sm$Zs[[j]] <- sdata[[paste0("Zs", p, "_", i, "_", 
                j)]]
            spars <- paste0("^s", p, "_", smframe$label[i], "_", 
                j, "\\[")
            sm$s[[j]] <- prepare_draws(draws, spars, regex = TRUE)
        }
        out$re[[i]] <- sm
    }
    out
}
