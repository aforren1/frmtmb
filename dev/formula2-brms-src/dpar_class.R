function (dpar, family = NULL) 
{
    out <- sub("[[:digit:]]*$", "", dpar)
    if (!is.null(family)) {
        if (conv_cats_dpars(family)) {
            multi_dpars <- valid_dpars(family, type = "multi")
            for (dp in multi_dpars) {
                sel <- grepl(paste0("^", dp), out)
                out[sel] <- dp
            }
        }
    }
    out
}
