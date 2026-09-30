has_ad_terms <- 
function (bterms, terms) 
{
    stopifnot(is.brmsterms(bterms), is.character(terms))
    any(ulapply(bterms$adforms[terms], is.formula))
}
