validate_draw_ids <- 
function (x, draw_ids = NULL, ndraws = NULL) 
{
    ndraws_total <- ndraws(x)
    if (is.null(draw_ids) && !is.null(ndraws)) {
        ndraws <- as_one_integer(ndraws)
        if (ndraws < 1 || ndraws > ndraws_total) {
            stop2("Argument 'ndraws' should be between 1 and ", 
                "the maximum number of draws (", ndraws_total, 
                ").")
        }
        draw_ids <- sample(seq_len(ndraws_total), ndraws)
    }
    if (!is.null(draw_ids)) {
        draw_ids <- as.integer(draw_ids)
        if (any(draw_ids < 1) || any(draw_ids > ndraws_total)) {
            stop2("Some 'draw_ids' indices are out of range.")
        }
    }
    draw_ids
}
