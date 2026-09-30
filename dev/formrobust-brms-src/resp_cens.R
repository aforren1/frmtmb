resp_cens <- 
function (x, y2 = NA) 
{
    cens <- deparse0(substitute(x))
    y2 <- deparse0(substitute(y2))
    class_resp_special("cens", call = match.call(), vars = nlist(cens, 
        y2))
}
