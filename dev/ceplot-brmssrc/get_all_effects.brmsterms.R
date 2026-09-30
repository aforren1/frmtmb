get_all_effects.brmsterms <- 
function (x, rsv_vars = NULL, comb_all = FALSE, ...) 
{
    stopifnot(is_atomic_or_null(rsv_vars))
    out <- list()
    for (dp in names(x$dpars)) {
        out <- c(out, get_all_effects(x$dpars[[dp]]))
    }
    for (nlp in names(x$nlpars)) {
        out <- c(out, get_all_effects(x$nlpars[[nlp]]))
    }
    out <- rmNULL(lapply(out, setdiff, y = rsv_vars))
    if (comb_all) {
        out <- unique(unlist(out))
        out <- c(out, get_group_vars(x))
        if (length(out)) {
            int <- expand.grid(out, out, stringsAsFactors = FALSE)
            int <- int[int[, 1] != int[, 2], ]
            int <- as.list(as.data.frame(t(int), stringsAsFactors = FALSE))
            int <- unique(unname(lapply(int, sort)))
            out <- c(as.list(out), int)
        }
    }
    unique(out[lengths(out) <= 2])
}
