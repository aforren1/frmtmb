function (family, link = NULL, threshold = NULL) 
{
    if (is.function(family)) {
        family <- family()
    }
    if (!is(family, "brmsfamily")) {
        if (is.family(family)) {
            link <- family$link
            family <- family$family
        }
        if (is.character(family)) {
            if (is.null(link)) {
                link <- family[2]
            }
            family <- .brmsfamily(family[1], link = link)
        }
        else {
            stop2("Argument 'family' is invalid.")
        }
    }
    if (is_ordinal(family) && !is.null(threshold)) {
        threshold <- match.arg(threshold, c("flexible", "equidistant"))
        family$threshold <- threshold
    }
    family
}
