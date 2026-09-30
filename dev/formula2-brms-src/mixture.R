function (..., flist = NULL, nmix = 1, order = NULL) 
{
    dots <- c(list(...), flist)
    if (length(nmix) == 1L) {
        nmix <- rep(nmix, length(dots))
    }
    if (length(dots) != length(nmix)) {
        stop2("The length of 'nmix' should be the same ", "as the number of mixture components.")
    }
    dots <- dots[rep(seq_along(dots), nmix)]
    family <- list(family = "mixture", link = "identity", mix = lapply(dots, 
        validate_family))
    class(family) <- c("mixfamily", "brmsfamily", "family")
    if (length(family$mix) < 2L) {
        stop2("Expecting at least 2 mixture components.")
    }
    if (use_real(family) && use_int(family)) {
        stop2("Cannot mix families with real and integer support.")
    }
    is_ordinal <- ulapply(family$mix, is_ordinal)
    if (any(is_ordinal) && any(!is_ordinal)) {
        stop2("Cannot mix ordinal and non-ordinal families.")
    }
    no_mixture <- ulapply(family$mix, no_mixture)
    if (any(no_mixture)) {
        stop2("Some of the families are not allowed in mixture models.")
    }
    for (fam in family$mix) {
        if (is.customfamily(fam) && "theta" %in% fam$dpars) {
            stop2("Parameter name 'theta' is reserved in mixture models.")
        }
    }
    if (is.null(order)) {
        if (any(is_ordinal)) {
            family$order <- "none"
            message("Setting order = 'none' for mixtures of ordinal families.")
        }
        else if (length(unique(family_names(family))) == 1L) {
            family$order <- "mu"
            message("Setting order = 'mu' for mixtures of the same family.")
        }
        else {
            family$order <- "none"
            message("Setting order = 'none' for mixtures of different families.")
        }
    }
    else {
        if (length(order) != 1L) {
            stop2("Argument 'order' must be of length 1.")
        }
        if (is.character(order)) {
            valid_order <- c("none", "mu")
            if (!order %in% valid_order) {
                stop2("Argument 'order' is invalid. Valid options are: ", 
                  collapse_comma(valid_order))
            }
            family$order <- order
        }
        else {
            family$order <- ifelse(as.logical(order), "mu", "none")
        }
    }
    family
}
