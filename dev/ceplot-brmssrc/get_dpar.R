get_dpar <- 
function (prep, dpar, i = NULL, inv_link = NULL) 
{
    stopifnot(is.brmsprep(prep) || is.mvbrmsprep(prep))
    dpar <- as_one_character(dpar)
    x <- prep$dpars[[dpar]]
    stopifnot(!is.null(x))
    if (is.list(x)) {
        out <- predictor(x, i = i, fprep = prep)
        if (is.null(inv_link)) {
            inv_link <- apply_dpar_inv_link(dpar, family = prep$family)
        }
        else {
            inv_link <- as_one_logical(inv_link)
        }
        if (inv_link) {
            out <- inv_link(out, x$family$link)
        }
        if (length(i) == 1) {
            out <- slice_col(out, 1)
        }
    }
    else if (!is.null(i) && !is.null(dim(x))) {
        out <- slice_col(x, i)
    }
    else {
        out <- x
    }
    out
}
