prepare_predictions_bhaz <- 
function (bframe, draws, sdata, ...) 
{
    if (!is_cox(bframe$family)) {
        return(NULL)
    }
    out <- list()
    p <- usc(combine_prefix(bframe))
    Zbhaz <- sdata[[paste0("Zbhaz", p)]]
    Zcbhaz <- sdata[[paste0("Zcbhaz", p)]]
    if (has_bhaz_groups(bframe)) {
        groups <- get_bhaz_groups(bframe)
        Jgrbhaz <- sdata[[paste0("Jgrbhaz", p)]]
        out$bhaz <- out$cbhaz <- matrix(nrow = nrow(draws), ncol = nrow(Zbhaz))
        for (k in seq_along(groups)) {
            sbhaz_regex <- paste0("^sbhaz", p, "\\[", groups[k], 
                ",")
            sbhaz <- prepare_draws(draws, sbhaz_regex, regex = TRUE)
            take <- Jgrbhaz == k
            out$bhaz[, take] <- tcrossprod(sbhaz, Zbhaz[take, 
                ])
            out$cbhaz[, take] <- tcrossprod(sbhaz, Zcbhaz[take, 
                ])
        }
    }
    else {
        sbhaz_regex <- paste0("^sbhaz", p)
        sbhaz <- prepare_draws(draws, sbhaz_regex, regex = TRUE)
        out$bhaz <- tcrossprod(sbhaz, Zbhaz)
        out$cbhaz <- tcrossprod(sbhaz, Zcbhaz)
    }
    out
}
