function (..., flist = NULL, rescor = NULL) 
{
    dots <- c(list(...), flist)
    if (!length(dots)) {
        stop2("No objects passed to 'mvbrmsformula'.")
    }
    forms <- list()
    for (i in seq_along(dots)) {
        if (is.mvbrmsformula(dots[[i]])) {
            forms <- c(forms, dots[[i]]$forms)
            if (is.null(rescor)) {
                rescor <- dots[[i]]$rescor
            }
        }
        else {
            forms <- c(forms, list(bf(dots[[i]])))
        }
    }
    if (!is.null(rescor)) {
        rescor <- as_one_logical(rescor)
    }
    responses <- ufrom_list(forms, "resp")
    if (any(duplicated(responses))) {
        stop2("Cannot use the same response variable twice in the same model.")
    }
    names(forms) <- responses
    structure(nlist(forms, responses, rescor), class = c("mvbrmsformula", 
        "bform"))
}
