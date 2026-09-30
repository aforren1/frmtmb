+.bform <- 
function (e1, e2) 
{
    if (is.brmsformula(e1)) {
        out <- plus_brmsformula(e1, e2)
    }
    else if (is.mvbrmsformula(e1)) {
        out <- plus_mvbrmsformula(e1, e2)
    }
    else {
        stop2("Method '+.bform' not implemented for ", class(e1)[1], 
            " objects.")
    }
    out
}
