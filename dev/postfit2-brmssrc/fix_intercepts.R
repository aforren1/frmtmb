fix_intercepts <- 
function (bterms) 
{
    dpar <- dpar_class(bterms[["dpar"]])
    if (!length(dpar)) 
        dpar <- "mu"
    isTRUE(is_ordinal(bterms) && dpar %in% bterms$family[["order"]])
}
