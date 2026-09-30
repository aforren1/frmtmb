function (..., flist = NULL, dpar = NULL, resp = NULL, center = NULL, 
    cmc = NULL, sparse = NULL, decomp = NULL) 
{
    out <- c(list(...), flist)
    warn_dpar(dpar)
    if (!is.null(resp)) {
        resp <- as_one_character(resp)
    }
    cmc <- if (!is.null(cmc)) 
        as_one_logical(cmc)
    center <- if (!is.null(center)) 
        as_one_logical(center)
    decomp <- if (!is.null(decomp)) 
        match.arg(decomp, decomp_opts())
    for (i in seq_along(out)) {
        if (!is.null(cmc)) {
            attr(out[[i]], "cmc") <- cmc
        }
        if (!is.null(center)) {
            attr(out[[i]], "center") <- center
        }
        if (!is.null(sparse)) {
            attr(out[[i]], "sparse") <- sparse
        }
        if (!is.null(decomp)) {
            attr(out[[i]], "decomp") <- decomp
        }
    }
    structure(out, resp = resp)
}
